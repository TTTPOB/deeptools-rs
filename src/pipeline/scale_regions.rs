use std::sync::Arc;

use anyhow::{Result, bail};

use crate::config::{GeneralOptions, IoOptions, ScaleRegionsOptions, SortRegions};
use crate::io::writers;
use crate::io::{BedRecord, Strand};
use crate::pipeline::core::{self, InMemorySink, PipelineMode, RegionTask, StreamingSink};
use crate::pipeline::matrix::{MatrixData, MatrixHeader, MatrixRow};
use crate::pipeline::zones::ScaleRegionsPlan;

use super::RunOutcome;

#[derive(Clone)]
struct ScaleRegionsMode {
    options: ScaleRegionsOptions,
}

impl ScaleRegionsMode {
    fn new(options: ScaleRegionsOptions) -> Self {
        Self { options }
    }
}

#[derive(Clone)]
struct ScaleRegionsMetadata {
    bin_size: u32,
    total_bins: usize,
}

impl PipelineMode for ScaleRegionsMode {
    type Plan = ScaleRegionsPlan;
    type Metadata = ScaleRegionsMetadata;

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
        ensure_multiple(
            general.bin_size,
            self.options.region_body_length,
            "regionBodyLength",
        )?;
        ensure_multiple(
            general.bin_size,
            self.options.unscaled_5_prime,
            "unscaled5prime",
        )?;
        ensure_multiple(
            general.bin_size,
            self.options.unscaled_3_prime,
            "unscaled3prime",
        )?;

        if self.options.region_body_length == 0
            && (self.options.unscaled_5_prime + self.options.unscaled_3_prime) > 0
        {
            bail!(
                "Unscaled 5- and 3-prime regions require a non-zero --regionBodyLength in scale-regions mode"
            );
        }

        let upstream_bins = (self.options.upstream / general.bin_size) as usize;
        let downstream_bins = (self.options.downstream / general.bin_size) as usize;
        let unscaled5_bins = (self.options.unscaled_5_prime / general.bin_size) as usize;
        let unscaled3_bins = (self.options.unscaled_3_prime / general.bin_size) as usize;
        let body_bins = (self.options.region_body_length / general.bin_size) as usize;

        let total_bins =
            upstream_bins + downstream_bins + unscaled5_bins + unscaled3_bins + body_bins;

        if total_bins == 0 {
            bail!("Scale-regions mode requires at least one bin to be generated");
        }

        Ok(ScaleRegionsMetadata {
            bin_size: general.bin_size,
            total_bins,
        })
    }

    fn total_bins(&self, metadata: &Self::Metadata) -> usize {
        metadata.total_bins
    }

    fn plan_for(&self, record: &BedRecord, metadata: &Self::Metadata) -> Self::Plan {
        ScaleRegionsPlan::scale_regions(record, &self.options, metadata.bin_size)
    }

    fn nan_after_end(&self, _metadata: &Self::Metadata) -> bool {
        false
    }

    fn postprocess_row(
        &self,
        record: BedRecord,
        mut values: Vec<Vec<f32>>,
        _metadata: &Self::Metadata,
    ) -> MatrixRow {
        if matches!(record.strand, Strand::Negative) {
            for sample_values in &mut values {
                sample_values.reverse();
            }
        }
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
        build_scale_regions_header(
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
    options: &ScaleRegionsOptions,
) -> Result<RunOutcome> {
    let mode = ScaleRegionsMode::new(options.clone());
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
        let sink = StreamingSink::new(writer);
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
            sink,
            thread_count,
            &mode,
            Arc::clone(&metadata),
            header_builder,
            group_labels.len(),
        )?;
        return Ok(RunOutcome::Streamed);
    }

    let sink = InMemorySink::with_capacity(row_count);
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

    let aggregation = core::execute_mode(
        tasks,
        general,
        sample_paths,
        sink,
        thread_count,
        &mode,
        metadata,
        header_builder,
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

fn build_scale_regions_header(
    general: &GeneralOptions,
    options: &ScaleRegionsOptions,
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
        nan_after_end: false,
        sort_using: general.sort_using.to_string(),
        unscaled_5_prime: vec![options.unscaled_5_prime; sample_count],
        body: vec![options.region_body_length; sample_count],
        sample_labels: sample_labels.to_vec(),
        downstream: vec![options.downstream; sample_count],
        unscaled_3_prime: vec![options.unscaled_3_prime; sample_count],
        group_labels: group_labels.to_vec(),
        bin_size: vec![bin_size; sample_count],
        upstream: vec![options.upstream; sample_count],
        group_boundaries,
        sample_boundaries,
        missing_data_as_zero: general.missing_data_as_zero,
        ref_point: vec![None; sample_count],
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

fn ensure_multiple(bin_size: u32, value: u32, flag: &str) -> Result<()> {
    if value % bin_size != 0 {
        bail!(
            "{flag} ({value}) must be a multiple of the bin size ({bin_size}) in scale-regions mode"
        );
    }
    Ok(())
}
