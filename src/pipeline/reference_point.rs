use std::sync::Arc;

use anyhow::{Result, bail};

use crate::config::{GeneralOptions, GtfOptions, IoOptions, ReferencePointOptions, SortRegions};
use crate::io::BedRecord;
use crate::io::writers;
use crate::pipeline::core::{
    self, FileCollector, InMemoryCollector, ModeTag, PipelineMode, RegionTask,
};
use crate::pipeline::matrix::{LayoutVectors, MatrixHeader, MatrixHeaderBuilder, MatrixRow};
use crate::pipeline::zones::ReferencePointPlan;

use super::RunOutcome;

#[derive(Clone)]
struct ReferencePointMode {
    options: ReferencePointOptions,
    keep_exons: bool,
}

impl ReferencePointMode {
    fn new(options: ReferencePointOptions, keep_exons: bool) -> Self {
        Self {
            options,
            keep_exons,
        }
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
        let mode = ModeTag::ReferencePoint;
        core::ensure_positive(general.bin_size, "binSize", mode)?;
        core::ensure_multiple(
            general.bin_size,
            self.options.upstream,
            "beforeRegionStartLength",
            mode,
        )?;
        core::ensure_multiple(
            general.bin_size,
            self.options.downstream,
            "afterRegionStartLength",
            mode,
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
            self.keep_exons,
            self.options.nan_after_end,
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
        // Extract exon coordinates for metagene mode output
        let exon_coords = if self.keep_exons {
            record.exons()
        } else {
            None
        };
        MatrixRow {
            record,
            values,
            exon_coords,
        }
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
        let layout = LayoutVectors::uniform(
            sample_count,
            metadata.bin_size,
            self.options.upstream,
            self.options.downstream,
            0,
            0,
            0,
            Some(self.options.reference_point.to_string()),
        );

        MatrixHeaderBuilder::new(
            general,
            sample_labels,
            group_labels,
            group_counts,
            thread_count,
            sample_count,
            self.options.nan_after_end,
        )
        .with_layout(layout)
        .with_uniform_sample_boundaries(metadata.total_bins)
        .build()
    }
}

pub fn run(
    general: &GeneralOptions,
    io: &IoOptions,
    gtf: &GtfOptions,
    options: &ReferencePointOptions,
) -> Result<RunOutcome> {
    let mode = ReferencePointMode::new(options.clone(), gtf.keep_exons);
    let metadata = Arc::new(mode.validate(general)?);

    let sample_labels = core::derive_sample_labels(&io.scores, general)?;
    let sample_count = sample_labels.len();

    let groups = core::load_groups(&io.regions, gtf)?;
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
