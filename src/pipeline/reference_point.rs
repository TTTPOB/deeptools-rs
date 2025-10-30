use std::sync::Arc;

use anyhow::{Result, bail};

use crate::config::{GeneralOptions, IoOptions, ReferencePointOptions, SortRegions};
use crate::io::BedRecord;
use crate::io::writers;
use crate::pipeline::core::{self, FileCollector, InMemoryCollector, PipelineMode, RegionTask};
use crate::pipeline::zones::ReferencePointPlan;

use super::RunOutcome;
use super::matrix::{MatrixData, MatrixHeader, MatrixRow};

#[derive(Clone)]
struct ReferencePointMode {
    options: ReferencePointOptions,
}

impl ReferencePointMode {
    fn new(options: ReferencePointOptions) -> Self {
        Self { options }
    }
}

#[derive(Clone)]
struct ReferencePointMetadata {
    bin_size: u32,
    upstream_bins: usize,
    downstream_bins: usize,
    total_bins: usize,
}

impl PipelineMode for ReferencePointMode {
    type Plan = ReferencePointPlan;
    type Metadata = ReferencePointMetadata;

    fn validate(&self, general: &GeneralOptions) -> Result<Self::Metadata> {
        ensure_positive(general.bin_size, "binSize")?;
        ensure_multiple(
            general.bin_size,
            self.options.upstream,
            "beforeRegionStartLength",
        )?;
        ensure_multiple(
            general.bin_size,
            self.options.downstream,
            "afterRegionStartLength",
        )?;

        let upstream_bins = (self.options.upstream / general.bin_size) as usize;
        let downstream_bins = (self.options.downstream / general.bin_size) as usize;
        let total_bins = upstream_bins + downstream_bins;

        if total_bins == 0 {
            bail!("Reference-point mode requires at least one upstream or downstream bin");
        }

        Ok(ReferencePointMetadata {
            bin_size: general.bin_size,
            upstream_bins,
            downstream_bins,
            total_bins,
        })
    }

    fn total_bins(&self, metadata: &Self::Metadata) -> usize {
        metadata.total_bins
    }

    fn plan_for(&self, record: &BedRecord, metadata: &Self::Metadata) -> Self::Plan {
        ReferencePointPlan::reference_point(
            record,
            self.options.reference_point,
            metadata.bin_size,
            metadata.upstream_bins,
            metadata.downstream_bins,
        )
    }

    fn nan_after_end(&self, _metadata: &Self::Metadata) -> bool {
        self.options.nan_after_end
    }

    fn postprocess_row(
        &self,
        record: BedRecord,
        values: Vec<Vec<f32>>,
        _metadata: &Self::Metadata,
    ) -> MatrixRow {
        MatrixRow { record, values }
    }

    fn build_header(
        &self,
        general: &GeneralOptions,
        metadata: &Self::Metadata,
        sample_labels: &[String],
        group_labels: &[String],
        group_counts: &[usize],
        thread_count: usize,
        sample_count: usize,
    ) -> MatrixHeader {
        build_reference_point_header(
            general,
            &self.options,
            sample_labels,
            group_labels,
            group_counts,
            metadata.bin_size,
            metadata.total_bins,
            sample_count,
            thread_count,
        )
    }
}

pub fn run(
    general: &GeneralOptions,
    io: &IoOptions,
    options: &ReferencePointOptions,
) -> Result<RunOutcome> {
    let mode = ReferencePointMode::new(options.clone());
    let metadata = Arc::new(mode.validate(general)?);

    let sample_labels = core::derive_sample_labels(&io.scores, general)?;
    let sample_count = sample_labels.len();

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
    let total_bins = mode.total_bins(metadata.as_ref());
    let row_count = tasks.len();
    let should_stream = writers::should_use_streaming_for_plan(
        row_count,
        sample_count,
        total_bins,
        matches!(general.sort_regions, SortRegions::Keep),
        io,
    );

    let sample_paths = Arc::new(io.scores.clone());

    if should_stream {
        let header_estimate = mode.build_header(
            general,
            metadata.as_ref(),
            &sample_labels,
            &group_labels,
            &group_capacity,
            thread_count,
            sample_count,
        );
        writers::ensure_streaming_header_capacity(&header_estimate)?;

        let writer = writers::StreamingMatrixWriter::start(&io.matrix_output)?;
        let collector = FileCollector::new(writer);
        let header_builder = {
            let general = general.clone();
            let sample_labels = sample_labels.clone();
            let group_labels = group_labels.clone();
            let metadata = Arc::clone(&metadata);
            let mode = mode.clone();
            move |group_counts: Vec<usize>| -> Result<MatrixHeader> {
                Ok(mode.build_header(
                    &general,
                    metadata.as_ref(),
                    &sample_labels,
                    &group_labels,
                    &group_counts,
                    thread_count,
                    sample_count,
                ))
            }
        };

        core::execute_mode(
            tasks,
            general,
            Arc::clone(&sample_paths),
            collector,
            thread_count,
            &mode,
            Arc::clone(&metadata),
            header_builder,
            group_labels.len(),
        )?;
        return Ok(RunOutcome::Streamed);
    }

    let collector = InMemoryCollector::with_capacity(row_count, sample_count, total_bins);
    let header_builder = {
        let general = general.clone();
        let sample_labels = sample_labels.clone();
        let group_labels = group_labels.clone();
        let metadata = Arc::clone(&metadata);
        let mode = mode.clone();
        move |group_counts: Vec<usize>| -> Result<MatrixHeader> {
            Ok(mode.build_header(
                &general,
                metadata.as_ref(),
                &sample_labels,
                &group_labels,
                &group_counts,
                thread_count,
                sample_count,
            ))
        }
    };

    let mut matrix = core::execute_mode(
        tasks,
        general,
        sample_paths,
        collector,
        thread_count,
        &mode,
        metadata,
        header_builder,
        group_labels.len(),
    )?;

    let sort_sample_indices =
        core::normalize_sort_sample_indices(general.sort_using_samples.as_ref(), sample_count)?;

    matrix.sort_groups(
        general.sort_regions,
        general.sort_using,
        sort_sample_indices.as_deref(),
    )?;

    matrix.prune_zero_rows();

    Ok(RunOutcome::Matrix(matrix))
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
