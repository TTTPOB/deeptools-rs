use std::{fmt, path::PathBuf};

use clap::ValueEnum;

#[derive(Debug, Clone)]
pub struct Config {
    pub mode: ModeConfig,
    pub io: IoOptions,
    pub general: GeneralOptions,
    pub gtf: GtfOptions,
}

#[derive(Debug, Clone)]
pub struct IoOptions {
    pub regions: Vec<PathBuf>,
    pub scores: Vec<PathBuf>,
    pub matrix_output: PathBuf,
    pub matrix_values_output: Option<PathBuf>,
    pub sorted_regions_output: Option<PathBuf>,
}

#[derive(Debug, Clone)]
pub struct GeneralOptions {
    pub bin_size: u32,
    pub sort_regions: SortRegions,
    pub sort_using: SortUsing,
    pub sort_using_samples: Option<Vec<usize>>,
    pub average_type_bins: AverageTypeBins,
    pub missing_data_as_zero: bool,
    pub skip_zeros: bool,
    pub min_threshold: Option<f64>,
    pub max_threshold: Option<f64>,
    pub blacklist: Option<PathBuf>,
    pub samples_label: Option<Vec<String>>,
    pub smart_labels: bool,
    pub quiet: bool,
    pub verbose: bool,
    pub scale_factor: f64,
    pub number_of_processors: ProcessorRequest,
}

#[derive(Debug, Clone)]
pub struct GtfOptions {
    pub keep_exons: bool,
    pub transcript_id: String,
    pub exon_id: String,
    pub transcript_id_designator: String,
}

impl Default for GtfOptions {
    fn default() -> Self {
        Self {
            keep_exons: false,
            transcript_id: "transcript".to_string(),
            exon_id: "exon".to_string(),
            transcript_id_designator: "transcript_id".to_string(),
        }
    }
}

#[derive(Debug, Clone)]
pub enum ModeConfig {
    ScaleRegions(ScaleRegionsOptions),
    ReferencePoint(ReferencePointOptions),
}

#[derive(Debug, Clone)]
pub struct ScaleRegionsOptions {
    pub region_body_length: u32,
    pub start_label: String,
    pub end_label: String,
    pub upstream: u32,
    pub downstream: u32,
    pub unscaled_5_prime: u32,
    pub unscaled_3_prime: u32,
}

#[derive(Debug, Clone)]
pub struct ReferencePointOptions {
    pub reference_point: ReferencePoint,
    pub upstream: u32,
    pub downstream: u32,
    pub nan_after_end: bool,
}

#[derive(Debug, Clone)]
pub enum ProcessorRequest {
    Max,
    MaxHalf,
    Fixed(u32),
}

#[derive(Debug, Clone, Copy, ValueEnum)]
#[clap(rename_all = "lower")]
pub enum SortRegions {
    Descend,
    Ascend,
    #[clap(name = "no")]
    No,
    Keep,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum SortUsing {
    Mean,
    Median,
    Max,
    Min,
    Sum,
    #[clap(name = "region_length")]
    RegionLength,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum AverageTypeBins {
    Mean,
    Median,
    Min,
    Max,
    Std,
    Sum,
}

/// Strategy for handling groups that become empty during the pipeline.
///
/// Two phases can produce empty groups:
///
/// **Pre-execution** (blacklist filtering): Occurs before output writers are
/// created.  The policy is determined by `SortRegions` — `Keep` uses `Error`,
/// all other modes use `Drop` (removing the empty group from `group_labels` and
/// `group_boundaries`, remapping subsequent indices).
///
/// **Runtime** (`--skipZeros`, `--minThreshold`, `--maxThreshold`): Occurs
/// while rows are being written.  Group structure is already fixed at output
/// time, so empty groups are implicitly preserved with a count of 0
/// (consecutive identical entries in `group_boundaries`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EmptyGroupPolicy {
    /// Remove empty groups from the header and remap subsequent group indices.
    /// Only valid for pre-execution filtering.
    Drop,
    /// Abort with an error if any group becomes empty.
    Error,
    /// Keep the group in the output, represented by consecutive identical
    /// boundaries (count = 0).  Implicit for runtime row filtering.
    PreserveWithZeroCount,
}

impl EmptyGroupPolicy {
    /// Policy for blacklist (pre-execution) filtering.
    pub fn for_sort_regions(sort: SortRegions) -> Self {
        match sort {
            SortRegions::Keep => EmptyGroupPolicy::Error,
            _ => EmptyGroupPolicy::Drop,
        }
    }

    /// Validate group counts against this policy.  Called when group counts
    /// are finalised, before building the output header.
    ///
    /// `Error` rejects any zero-count group.  `Drop` is invalid at runtime
    /// (indices are already fixed in output).  `PreserveWithZeroCount`
    /// explicitly allows zero-count groups.
    pub fn validate(&self, group_counts: &[usize], group_labels: &[String]) -> Result<(), String> {
        match self {
            EmptyGroupPolicy::Error => {
                if let Some(idx) = group_counts.iter().position(|c| *c == 0) {
                    return Err(format!("Group '{}' is empty", group_labels[idx]));
                }
            }
            EmptyGroupPolicy::Drop => {
                return Err(
                    "Drop policy cannot be applied at runtime — group indices are \
                     already fixed in output. Use PreserveWithZeroCount."
                        .into(),
                );
            }
            EmptyGroupPolicy::PreserveWithZeroCount => {
                // Zero-count groups are explicitly allowed.
            }
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy)]
pub enum ReferencePoint {
    Tss,
    Tes,
    Center,
}

impl ValueEnum for ReferencePoint {
    fn value_variants<'a>() -> &'a [Self] {
        static VARIANTS: [ReferencePoint; 3] = [
            ReferencePoint::Tss,
            ReferencePoint::Tes,
            ReferencePoint::Center,
        ];
        &VARIANTS
    }

    fn to_possible_value(&self) -> Option<clap::builder::PossibleValue> {
        match self {
            ReferencePoint::Tss => Some(clap::builder::PossibleValue::new("TSS")),
            ReferencePoint::Tes => Some(clap::builder::PossibleValue::new("TES")),
            ReferencePoint::Center => Some(clap::builder::PossibleValue::new("center")),
        }
    }
}

impl Default for SortRegions {
    fn default() -> Self {
        SortRegions::Keep
    }
}

impl Default for SortUsing {
    fn default() -> Self {
        SortUsing::Mean
    }
}

impl Default for AverageTypeBins {
    fn default() -> Self {
        AverageTypeBins::Mean
    }
}

impl Default for ReferencePoint {
    fn default() -> Self {
        ReferencePoint::Tss
    }
}

impl fmt::Display for SortRegions {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let value = match self {
            SortRegions::Descend => "descend",
            SortRegions::Ascend => "ascend",
            SortRegions::No => "no",
            SortRegions::Keep => "keep",
        };
        f.write_str(value)
    }
}

impl fmt::Display for SortUsing {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let value = match self {
            SortUsing::Mean => "mean",
            SortUsing::Median => "median",
            SortUsing::Max => "max",
            SortUsing::Min => "min",
            SortUsing::Sum => "sum",
            SortUsing::RegionLength => "region_length",
        };
        f.write_str(value)
    }
}

impl fmt::Display for AverageTypeBins {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let value = match self {
            AverageTypeBins::Mean => "mean",
            AverageTypeBins::Median => "median",
            AverageTypeBins::Min => "min",
            AverageTypeBins::Max => "max",
            AverageTypeBins::Std => "std",
            AverageTypeBins::Sum => "sum",
        };
        f.write_str(value)
    }
}

impl fmt::Display for ReferencePoint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let value = match self {
            ReferencePoint::Tss => "TSS",
            ReferencePoint::Tes => "TES",
            ReferencePoint::Center => "center",
        };
        f.write_str(value)
    }
}

impl fmt::Display for ProcessorRequest {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ProcessorRequest::Max => f.write_str("max"),
            ProcessorRequest::MaxHalf => f.write_str("max/2"),
            ProcessorRequest::Fixed(value) => write!(f, "{}", value),
        }
    }
}

impl ProcessorRequest {
    pub fn resolve(&self) -> u32 {
        match self {
            ProcessorRequest::Max => available_cpus(),
            ProcessorRequest::MaxHalf => {
                let cpus = available_cpus();
                std::cmp::max(1, cpus / 2)
            }
            ProcessorRequest::Fixed(value) => std::cmp::max(1, *value),
        }
    }
}

fn available_cpus() -> u32 {
    let count = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    u32::try_from(count).unwrap_or(u32::MAX)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resolve_max_returns_at_least_one() {
        assert!(ProcessorRequest::Max.resolve() >= 1);
    }

    #[test]
    fn resolve_max_half_returns_at_least_one() {
        assert!(ProcessorRequest::MaxHalf.resolve() >= 1);
    }

    #[test]
    fn resolve_fixed_zero_clamps_to_one() {
        assert_eq!(ProcessorRequest::Fixed(0).resolve(), 1);
    }

    #[test]
    fn resolve_fixed_value() {
        assert_eq!(ProcessorRequest::Fixed(4).resolve(), 4);
    }
}
