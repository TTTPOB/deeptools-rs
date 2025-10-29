use std::sync::Arc;

use anyhow::{Context, Result, bail};
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

use crate::config::{GeneralOptions, IoOptions, ReferencePointOptions};
use crate::io::BedRecord;
use crate::pipeline::core::{self, WorkerSamples};
use crate::pipeline::zones::ReferencePointPlan;

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
) -> Result<MatrixData> {
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
    let mut group_counts = vec![0usize; groups.len()];

    let mut tasks = Vec::new();
    for (group_index, group) in groups.iter().enumerate() {
        for record in &group.records {
            let index = tasks.len();
            tasks.push(RegionTask {
                index,
                group_index,
                record: record.clone(),
            });
        }
    }

    let task_count = tasks.len();
    let thread_count = std::cmp::max(1, general.number_of_processors.resolve() as usize);
    let mut rows = Vec::with_capacity(task_count);

    if task_count > 0 {
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

    let group_labels: Vec<String> = groups.iter().map(|group| group.label.clone()).collect();
    let group_boundaries = MatrixData::group_boundaries_from_counts(&group_counts);
    let sample_boundaries = MatrixData::sample_boundaries_uniform(sample_count, total_bins);

    let proc_number = thread_count as u32;

    let header = MatrixHeader {
        verbose: general.verbose,
        scale: general.scale_factor,
        skip_zeros: general.skip_zeros,
        nan_after_end: options.nan_after_end,
        sort_using: general.sort_using.to_string(),
        unscaled_5_prime: vec![0; sample_count],
        body: vec![0; sample_count],
        sample_labels: sample_labels.clone(),
        downstream: vec![options.downstream; sample_count],
        unscaled_3_prime: vec![0; sample_count],
        group_labels,
        bin_size: vec![bin_size; sample_count],
        upstream: vec![options.upstream; sample_count],
        group_boundaries,
        sample_boundaries,
        missing_data_as_zero: general.missing_data_as_zero,
        ref_point: vec![Some(options.reference_point.to_string()); sample_count],
        min_threshold: general.min_threshold,
        sort_regions: general.sort_regions.to_string(),
        proc_number,
        bin_avg_type: general.average_type_bins.to_string(),
        max_threshold: general.max_threshold,
    };

    let sort_sample_indices =
        normalize_sort_sample_indices(general.sort_using_samples.as_ref(), sample_count)?;

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

    Ok(matrix)
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

fn normalize_sort_sample_indices(
    raw: Option<&Vec<usize>>,
    sample_count: usize,
) -> Result<Option<Vec<usize>>> {
    let Some(raw_indices) = raw else {
        return Ok(None);
    };

    if raw_indices.is_empty() {
        return Ok(None);
    }

    let mut normalized = Vec::with_capacity(raw_indices.len());
    for &value in raw_indices {
        if value == 0 || value > sample_count {
            bail!(
                "The value {} for --sortUsingSamples is not valid. Only values from 1 to {} are allowed.",
                value,
                sample_count
            );
        }
        normalized.push(value - 1);
    }

    Ok(Some(normalized))
}
