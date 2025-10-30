use std::sync::{Arc, mpsc};

use anyhow::{Context, Result, anyhow, bail};
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

use crate::config::{GeneralOptions, IoOptions, ReferencePointOptions, SortRegions};
use crate::io::BedRecord;
use crate::io::writers;
use crate::pipeline::core::{
    self, AggregatedResult, InMemorySink, RegionResult, StreamingSink, WorkerSamples,
};
use crate::pipeline::zones::ReferencePointPlan;

use super::RunOutcome;
use super::matrix::{MatrixData, MatrixHeader, MatrixRow};

struct RegionTask {
    index: usize,
    group_index: usize,
    record: BedRecord,
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

    let sample_paths = Arc::new(io.scores.clone());
    let make_header_builder = |sample_labels: &[String], group_labels: &[String]| {
        let general = general.clone();
        let options = options.clone();
        let sample_labels = sample_labels.to_vec();
        let group_labels = group_labels.to_vec();
        move |group_counts: Vec<usize>| -> Result<MatrixHeader> {
            Ok(build_reference_point_header(
                &general,
                &options,
                &sample_labels,
                &group_labels,
                &group_counts,
                bin_size,
                total_bins,
                sample_count,
                thread_count,
            ))
        }
    };

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

        let writer = writers::StreamingMatrixWriter::start(&io.matrix_output)?;
        let sink = StreamingSink::new(writer);
        execute_reference_point(
            tasks,
            general,
            options,
            sample_paths,
            sink,
            thread_count,
            upstream_bins,
            downstream_bins,
            make_header_builder(&sample_labels, &group_labels),
            group_labels.len(),
        )?;
        return Ok(RunOutcome::Streamed);
    }

    let sink = InMemorySink::with_capacity(row_count);
    let aggregation = execute_reference_point(
        tasks,
        general,
        options,
        sample_paths,
        sink,
        thread_count,
        upstream_bins,
        downstream_bins,
        make_header_builder(&sample_labels, &group_labels),
        group_labels.len(),
    )?;

    let sort_sample_indices =
        core::normalize_sort_sample_indices(general.sort_using_samples.as_ref(), sample_count)?;

    let mut matrix = MatrixData {
        header: aggregation.header,
        rows: aggregation.output,
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

fn execute_reference_point<S, F>(
    tasks: Vec<RegionTask>,
    general: &GeneralOptions,
    options: &ReferencePointOptions,
    sample_paths: Arc<Vec<std::path::PathBuf>>,
    sink: S,
    thread_count: usize,
    upstream_bins: usize,
    downstream_bins: usize,
    header_builder: F,
    group_count: usize,
) -> Result<AggregatedResult<S::Output>>
where
    S: core::RowSink + Send + 'static,
    F: FnOnce(Vec<usize>) -> Result<MatrixHeader> + Send + 'static,
{
    let task_count = tasks.len();
    let (tx, rx) = mpsc::channel();

    let aggregator_handle =
        core::spawn_row_aggregator(rx, sink, group_count, task_count, header_builder)?;

    let pool = ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build()
        .context("Failed to initialise rayon thread pool for reference-point scheduling")?;

    let compute_result = pool.install(|| {
        tasks
            .into_par_iter()
            .map_init(
                || (WorkerSamples::new(sample_paths.clone()), tx.clone()),
                |state, task| {
                    let (worker_samples, sender) = state;
                    let RegionTask {
                        index,
                        group_index,
                        record,
                    } = task;

                    let samples = worker_samples.samples()?;
                    let plan = ReferencePointPlan::reference_point(
                        &record,
                        options.reference_point,
                        general.bin_size,
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

                    sender
                        .send(RegionResult {
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

    compute_result?;

    let aggregation = aggregator_handle
        .join()
        .map_err(|_| anyhow!("Matrix aggregation thread panicked"))??;

    Ok(aggregation)
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
