use anyhow::{Result, bail};

use crate::config::{GeneralOptions, GtfOptions, IoOptions, ScaleRegionsOptions};
use crate::io::{BedRecord, Strand};
use crate::pipeline::core::{self, ModeTag, PipelineMode};
use crate::pipeline::matrix::{LayoutVectors, MatrixHeader, MatrixHeaderBuilder, MatrixRow};
use crate::pipeline::zones::ScaleRegionsPlan;

use super::RunOutcome;

#[derive(Clone)]
struct ScaleRegionsMode {
    options: ScaleRegionsOptions,
    keep_exons: bool,
}

impl ScaleRegionsMode {
    fn new(options: ScaleRegionsOptions, keep_exons: bool) -> Self {
        Self {
            options,
            keep_exons,
        }
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
        let mode = ModeTag::ScaleRegions;
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
        core::ensure_multiple(
            general.bin_size,
            self.options.region_body_length,
            "regionBodyLength",
            mode,
        )?;
        core::ensure_multiple(
            general.bin_size,
            self.options.unscaled_5_prime,
            "unscaled5prime",
            mode,
        )?;
        core::ensure_multiple(
            general.bin_size,
            self.options.unscaled_3_prime,
            "unscaled3prime",
            mode,
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
        ScaleRegionsPlan::scale_regions(record, &self.options, metadata.bin_size, self.keep_exons)
    }

    fn nan_after_end(&self, _metadata: &Self::Metadata) -> bool {
        false
    }

    fn postprocess_row(
        &self,
        record: BedRecord,
        mut values: Vec<f32>,
        sample_count: usize,
        bin_count: usize,
        _metadata: &Self::Metadata,
    ) -> MatrixRow {
        if matches!(record.strand, Strand::Negative) {
            for sample_idx in 0..sample_count {
                let start = sample_idx * bin_count;
                values[start..start + bin_count].reverse();
            }
        }
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
            self.options.region_body_length,
            self.options.unscaled_5_prime,
            self.options.unscaled_3_prime,
            None,
        );

        MatrixHeaderBuilder::new(
            general,
            sample_labels,
            group_labels,
            group_counts,
            thread_count,
            sample_count,
            false,
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
    options: &ScaleRegionsOptions,
) -> Result<RunOutcome> {
    let mode = ScaleRegionsMode::new(options.clone(), gtf.keep_exons);
    super::run_pipeline(mode, general, io, gtf)
}
