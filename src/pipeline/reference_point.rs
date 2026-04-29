use anyhow::{Result, bail};

use crate::config::{GeneralOptions, GtfOptions, IoOptions, ReferencePointOptions};
use crate::io::BedRecord;
use crate::pipeline::core::{self, ModeTag, PipelineMode};
use crate::pipeline::matrix::{LayoutVectors, MatrixHeader, MatrixHeaderBuilder, MatrixRow};
use crate::pipeline::zones::ReferencePointPlan;

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
        values: Vec<f64>,
        sample_count: usize,
        bin_count: usize,
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
            sample_count,
            bin_count,
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
) -> Result<()> {
    let mode = ReferencePointMode::new(options.clone(), gtf.keep_exons);
    super::run_pipeline(mode, general, io, gtf)
}
