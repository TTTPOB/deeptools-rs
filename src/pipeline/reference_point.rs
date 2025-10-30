use std::collections::BTreeMap;
use std::sync::{Arc, mpsc};
use std::thread;

use anyhow::{Context, Result, anyhow, bail};
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

use crate::config::{GeneralOptions, IoOptions, ReferencePointOptions, SortRegions};
use crate::io::BedRecord;
use crate::io::writers::{self, StreamingMatrixWriter};
use crate::pipeline::core::{self, WorkerSamples};
use crate::pipeline::zones::ReferencePointPlan;

use super::RunOutcome;
use super::matrix::{MatrixData, MatrixHeader, MatrixRow};

struct RegionTask {
    index: usize,
    group_index: usize,
    record: BedRecord,
}

struct RegionResult {
    index: usize,
    group_index: usize,
    row: Option<MatrixRow>,
}

pub fn run(
    general: &GeneralOptions,
    io: &IoOptions,
    options: &ReferencePointOptions,
) -> Result<RunOutcome> {
    let bin_size = general.bin_size;
    ensure_positive(bin_size, "binSize")?;

    ensure_multiple(bin_size, options.upstream, "beforeRegionStartLength")?;
    ensure_multiple(bin_size, options.downstream, "afterRegionStartLength")?;

    let upstream_bins = (options.upstream / bin_size) as usize;
    let downstream_bins = (options.downstream / bin_size) as usize;
    let total_bins = upstream_bins + downstream_bins;

    if total_bins == 0 {
        bail!("Reference-point mode requires at least one upstream or downstream bin");
    }

    let sample_labels = core::derive_sample_labels(&io.scores, general)?;
    let sample_count = io.scores.len();

    let groups = core::load_groups(&io.regions)?;
    let group_labels: Vec<String> = groups.iter().map(|group| group.label.clone()).collect();

    let mut tasks = Vec::new();
    let mut group_capacity = Vec::with_capacity(groups.len());
    for (group_index, group) in groups.iter().enumerate() {
        group_capacity.push(group.records.len());
        for record in &group.records {
            let index = tasks.len();
            tasks.push(RegionTask {
                index,
                group_index,
                record: record.clone(),
            });
        }
    }

    let thread_count = std::cmp::max(1, general.number_of_processors.resolve() as usize);
    let row_count = tasks.len();
    let should_stream = writers::should_use_streaming_for_plan(
        row_count,
        sample_count,
        total_bins,
        matches!(general.sort_regions, SortRegions::Keep),
        io,
    );

    if should_stream {
        let header_estimate = build_reference_point_header(
            general,
            options,
            &sample_labels,
            &group_labels,
            &group_capacity,
            bin_size,
            total_bins,
            sample_count,
            thread_count,
        );
        writers::ensure_streaming_header_capacity(&header_estimate)?;

        return run_streaming_reference_point(
            general.clone(),
            options.clone(),
            io,
            sample_labels,
            group_labels,
            tasks,
            thread_count,
            bin_size,
            total_bins,
            upstream_bins,
            downstream_bins,
        );
    }

    let mut group_counts = vec![0usize; group_labels.len()];
    let mut rows = Vec::with_capacity(row_count);

    if row_count > 0 {
        let pool = ThreadPoolBuilder::new()
            .num_threads(thread_count)
            .build()
            .context("Failed to initialise rayon thread pool for reference-point scheduling")?;

        let sample_paths = Arc::new(io.scores.clone());
        let collected = pool.install(|| {
            tasks
                .into_par_iter()
                .map_init(
                    || WorkerSamples::new(sample_paths.clone()),
                    |state, task| {
                        let RegionTask {
                            index,
                            group_index,
                            record,
                        } = task;

                        let samples = state.samples()?;

                        let plan = ReferencePointPlan::reference_point(
                            &record,
                            options.reference_point,
                            bin_size,
                            upstream_bins,
                            downstream_bins,
                        );

                        let maybe_values = core::compute_row(
                            samples.as_mut_slice(),
                            &record,
                            &plan,
                            general,
                            options.nan_after_end,
                        )?;
                        let row = maybe_values.map(|values| MatrixRow { record, values });

                        Ok(RegionResult {
                            index,
                            group_index,
                            row,
                        })
                    },
                )
                .collect::<Vec<_>>()
        });

        let mut collected = collected.into_iter().collect::<Result<Vec<_>>>()?;
        collected.sort_by_key(|entry| entry.index);

        for entry in collected {
            if let Some(row) = entry.row {
                group_counts[entry.group_index] += 1;
                rows.push(row);
            }
        }
    }

    let header = build_reference_point_header(
        general,
        options,
        &sample_labels,
        &group_labels,
        &group_counts,
        bin_size,
        total_bins,
        sample_count,
        thread_count,
    );

    let sort_sample_indices =
        core::normalize_sort_sample_indices(general.sort_using_samples.as_ref(), sample_count)?;

    let mut matrix = MatrixData {
        header,
        rows,
        bin_count: total_bins,
        sample_count,
    };

    matrix.sort_groups(
        general.sort_regions,
        general.sort_using,
        sort_sample_indices.as_deref(),
    )?;

    matrix.prune_zero_rows();

    Ok(RunOutcome::Matrix(matrix))
}

fn run_streaming_reference_point(
    general: GeneralOptions,
    options: ReferencePointOptions,
    io: &IoOptions,
    sample_labels: Vec<String>,
    group_labels: Vec<String>,
    tasks: Vec<RegionTask>,
    thread_count: usize,
    bin_size: u32,
    total_bins: usize,
    upstream_bins: usize,
    downstream_bins: usize,
) -> Result<RunOutcome> {
    let sample_count = io.scores.len();
    let task_count = tasks.len();

    if task_count == 0 {
        let writer = StreamingMatrixWriter::start(&io.matrix_output)?;
        let empty_counts = vec![0usize; group_labels.len()];
        let header = build_reference_point_header(
            &general,
            &options,
            &sample_labels,
            &group_labels,
            &empty_counts,
            bin_size,
            total_bins,
            sample_count,
            thread_count,
        );
        writer.finish(&header)?;
        return Ok(RunOutcome::Streamed);
    }

    let writer = StreamingMatrixWriter::start(&io.matrix_output)?;
    let (tx, rx) = mpsc::channel();

    let aggregator_general = general.clone();
    let aggregator_options = options.clone();
    let aggregator_sample_labels = sample_labels;
    let aggregator_group_labels = group_labels;
    let aggregator_handle = thread::Builder::new()
        .name("matrix-streamer".into())
        .spawn(move || {
            stream_reference_point_rows(
                rx,
                writer,
                aggregator_general,
                aggregator_options,
                aggregator_sample_labels,
                aggregator_group_labels,
                bin_size,
                total_bins,
                sample_count,
                thread_count,
                task_count,
            )
        })
        .map_err(|err| anyhow!("Failed to spawn matrix streaming thread: {err}"))?;

    let pool = ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build()
        .context("Failed to initialise rayon thread pool for reference-point streaming")?;

    let sample_paths = Arc::new(io.scores.clone());
    let compute_result = pool.install(|| {
        tasks
            .into_par_iter()
            .map_init(
                || WorkerSamples::new(sample_paths.clone()),
                |state, task| {
                    let RegionTask {
                        index,
                        group_index,
                        record,
                    } = task;

                    let samples = state.samples()?;
                    let plan = ReferencePointPlan::reference_point(
                        &record,
                        options.reference_point,
                        bin_size,
                        upstream_bins,
                        downstream_bins,
                    );

                    let maybe_values = core::compute_row(
                        samples.as_mut_slice(),
                        &record,
                        &plan,
                        &general,
                        options.nan_after_end,
                    )?;
                    let row = maybe_values.map(|values| MatrixRow { record, values });

                    tx.send(RegionResult {
                        index,
                        group_index,
                        row,
                    })
                    .map_err(|err| anyhow!("Failed to stream computed row: {err}"))?;

                    Ok::<(), anyhow::Error>(())
                },
            )
            .try_reduce(|| (), |_, _| Ok::<(), anyhow::Error>(()))
    });

    drop(tx);

    let writer_result = aggregator_handle
        .join()
        .map_err(|_| anyhow!("Matrix streaming thread panicked"))?;

    compute_result?;
    writer_result?;

    Ok(RunOutcome::Streamed)
}

fn stream_reference_point_rows(
    rx: mpsc::Receiver<RegionResult>,
    mut writer: StreamingMatrixWriter,
    general: GeneralOptions,
    options: ReferencePointOptions,
    sample_labels: Vec<String>,
    group_labels: Vec<String>,
    bin_size: u32,
    total_bins: usize,
    sample_count: usize,
    thread_count: usize,
    task_count: usize,
) -> Result<()> {
    let mut group_counts = vec![0usize; group_labels.len()];
    let mut buffer = BTreeMap::new();
    let mut next_index = 0usize;

    while let Ok(result) = rx.recv() {
        buffer.insert(result.index, result);
        while let Some(entry) = buffer.remove(&next_index) {
            if let Some(row) = entry.row {
                writer.write_row(&row)?;
                group_counts[entry.group_index] += 1;
            }
            next_index += 1;
        }
    }

    while let Some(entry) = buffer.remove(&next_index) {
        if let Some(row) = entry.row {
            writer.write_row(&row)?;
            group_counts[entry.group_index] += 1;
        }
        next_index += 1;
    }

    if next_index != task_count {
        return Err(anyhow!(
            "Streamed matrix received {} of {} expected rows",
            next_index,
            task_count
        ));
    }

    let header = build_reference_point_header(
        &general,
        &options,
        &sample_labels,
        &group_labels,
        &group_counts,
        bin_size,
        total_bins,
        sample_count,
        thread_count,
    );
    writer.finish(&header)?;
    Ok(())
}

fn build_reference_point_header(
    general: &GeneralOptions,
    options: &ReferencePointOptions,
    sample_labels: &[String],
    group_labels: &[String],
    group_counts: &[usize],
    bin_size: u32,
    total_bins: usize,
    sample_count: usize,
    thread_count: usize,
) -> MatrixHeader {
    let group_boundaries = MatrixData::group_boundaries_from_counts(group_counts);
    let sample_boundaries = MatrixData::sample_boundaries_uniform(sample_count, total_bins);
    MatrixHeader {
        verbose: general.verbose,
        scale: general.scale_factor,
        skip_zeros: general.skip_zeros,
        nan_after_end: options.nan_after_end,
        sort_using: general.sort_using.to_string(),
        unscaled_5_prime: vec![0; sample_count],
        body: vec![0; sample_count],
        sample_labels: sample_labels.to_vec(),
        downstream: vec![options.downstream; sample_count],
        unscaled_3_prime: vec![0; sample_count],
        group_labels: group_labels.to_vec(),
        bin_size: vec![bin_size; sample_count],
        upstream: vec![options.upstream; sample_count],
        group_boundaries,
        sample_boundaries,
        missing_data_as_zero: general.missing_data_as_zero,
        ref_point: vec![Some(options.reference_point.to_string()); sample_count],
        min_threshold: general.min_threshold,
        sort_regions: general.sort_regions.to_string(),
        proc_number: thread_count as u32,
        bin_avg_type: general.average_type_bins.to_string(),
        max_threshold: general.max_threshold,
    }
}

fn ensure_positive(value: u32, flag: &str) -> Result<()> {
    if value == 0 {
        bail!("{flag} must be a positive integer");
    }
    Ok(())
}

fn ensure_multiple(bin_size: u32, distance: u32, flag: &str) -> Result<()> {
    if distance % bin_size != 0 {
        bail!(
            "{flag} ({distance}) must be a multiple of the bin size ({bin_size}) in reference-point mode"
        );
    }
    Ok(())
}
