use std::path::PathBuf;

use anyhow::{Context, anyhow, bail};
use clap::{ArgAction, Args, Parser, Subcommand};

use crate::config::{
    AverageTypeBins, Config, GeneralOptions, GtfOptions, IoOptions, ModeConfig, ProcessorRequest,
    ReferencePoint, ReferencePointOptions, ScaleRegionsOptions, SortRegions, SortUsing,
};

#[derive(Debug, Parser)]
#[command(
    name = "computeMatrix",
    about = "Rust reimplementation of deeptools computeMatrix",
    version,
    disable_help_subcommand = true,
    arg_required_else_help = true
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Debug, Subcommand)]
pub enum Command {
    #[command(name = "scale-regions")]
    ScaleRegions(ScaleRegionsCli),
    #[command(name = "reference-point")]
    ReferencePoint(ReferencePointCli),
}

#[derive(Debug, Args)]
pub struct ScaleRegionsCli {
    #[command(flatten)]
    required: RequiredArgs,
    #[command(flatten)]
    io: OutputArgs,
    #[command(flatten)]
    general: GeneralArgs,
    #[command(flatten)]
    gtf: GtfArgs,
    #[command(flatten)]
    mode: ScaleRegionsModeArgs,
}

#[derive(Debug, Args)]
pub struct ReferencePointCli {
    #[command(flatten)]
    required: RequiredArgs,
    #[command(flatten)]
    io: OutputArgs,
    #[command(flatten)]
    general: GeneralArgs,
    #[command(flatten)]
    gtf: GtfArgs,
    #[command(flatten)]
    mode: ReferencePointModeArgs,
}

#[derive(Debug, Args)]
struct RequiredArgs {
    #[arg(
        short = 'R',
        long = "regionsFileName",
        value_name = "File",
        help = "File name or names, in BED or GTF format, containing the regions to plot. If multiple bed files are given, each one is considered a group that can be plotted separately. Also, adding a \"#\" symbol in the bed file causes all the regions until the previous \"#\" to be considered one group.",
        num_args = 1..,
        required = true
    )]
    regions: Vec<PathBuf>,

    #[arg(
        short = 'S',
        long = "scoreFileName",
        value_name = "File",
        help = "bigWig file(s) containing the scores to be plotted. Multiple files should be separated by spaces. BigWig files can be obtained by using the bamCoverage or bamCompare tools. More information about the bigWig file format can be found at http://genome.ucsc.edu/goldenPath/help/bigWig.html",
        num_args = 1..,
        required = true
    )]
    scores: Vec<PathBuf>,
}

#[derive(Debug, Args)]
struct OutputArgs {
    #[arg(
        short = 'o',
        long = "outFileName",
        value_name = "FILE",
        help = "File name to save the gzipped matrix file needed by the \"plotHeatmap\" and \"plotProfile\" tools."
    )]
    matrix_output: PathBuf,

    #[arg(
        long = "outFileNameMatrix",
        value_name = "FILE",
        help = "If this option is given, then the matrix of values underlying the heatmap will be saved using the indicated name, e.g. IndividualValues.tab. This matrix can easily be loaded into R or other programs."
    )]
    matrix_values: Option<PathBuf>,

    #[arg(
        long = "outFileSortedRegions",
        value_name = "BED file",
        help = "File name in which the regions are saved after skiping zeros or min/max threshold values. The order of the regions in the file follows the sorting order selected. This is useful, for example, to generate other heatmaps keeping the sorting of the first heatmap. Example: Heatmap1sortedRegions.bed"
    )]
    sorted_regions: Option<PathBuf>,
}

#[derive(Debug, Args)]
struct GeneralArgs {
    #[arg(
        long = "binSize",
        alias = "bs",
        value_name = "INT bp",
        default_value_t = 10,
        help = "Length, in bases, of the non-overlapping bins for averaging the score over the regions length. (Default: 10)"
    )]
    bin_size: u32,

    #[arg(
        long = "sortRegions",
        default_value_t = SortRegions::Keep,
        help = "Whether the output file should present the regions sorted. The default is to not sort the regions. Note that this is only useful if you plan to plot the results yourself and not, for example, with plotHeatmap, which will override this. Note also that unsorted output will be in whatever order the regions happen to be processed in and not match the order in the input files. If you require the output order to match that of the input regions, then either specify \"keep\" or use computeMatrixOperations to resort the results file. (Default: keep)"
    )]
    sort_regions: SortRegions,

    #[arg(
        long = "sortUsing",
        default_value_t = SortUsing::Mean,
        help = "Indicate which method should be used for sorting. The value is computed for each row. Note that the region_length option will lead to a dotted line within the heatmap that indicates the end of the regions. (Default: mean)"
    )]
    sort_using: SortUsing,

    #[arg(
        long = "sortUsingSamples",
        value_name = "INT",
        num_args = 1..,
        help = "List of sample numbers (order as in matrix), that are used for sorting by --sortUsing, no value uses all samples, example: --sortUsingSamples 1 3"
    )]
    sort_using_samples: Vec<usize>,

    #[arg(
        long = "averageTypeBins",
        default_value_t = AverageTypeBins::Mean,
        help = "Define the type of statistic that should be used over the bin size range. The options are: \"mean\", \"median\", \"min\", \"max\", \"sum\" and \"std\". The default is \"mean\". (Default: mean)"
    )]
    average_type_bins: AverageTypeBins,

    #[arg(
        long = "missingDataAsZero",
        action = ArgAction::SetTrue,
        help = "If set, missing data (NAs) will be treated as zeros. The default is to ignore such cases, which will be depicted as black areas in a heatmap. (see the --missingDataColor argument of the plotHeatmap command for additional options)."
    )]
    missing_data_as_zero: bool,

    #[arg(
        long = "skipZeros",
        action = ArgAction::SetTrue,
        help = "Whether regions with only scores of zero should be included or not. Default is to include them."
    )]
    skip_zeros: bool,

    #[arg(
        long = "minThreshold",
        value_name = "FLOAT",
        help = "Numeric value. Any region containing a value that is less than or equal to this will be skipped. This is useful to skip, for example, genes where the read count is zero for any of the bins. This could be the result of unmappable areas and can bias the overall results. (Default: None)"
    )]
    min_threshold: Option<f64>,

    #[arg(
        long = "maxThreshold",
        value_name = "FLOAT",
        help = "Numeric value. Any region containing a value greater than or equal to this will be skipped. The maxThreshold is useful to skip those few regions with very high read counts (e.g. micro satellites) that may bias the average values. (Default: None)"
    )]
    max_threshold: Option<f64>,

    #[arg(
        long = "blackListFileName",
        alias = "bl",
        value_name = "BED file",
        help = "A BED file containing regions that should be excluded from all analyses. Currently this works by rejecting genomic chunks that happen to overlap an entry. Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, then the read/fragment might still be considered."
    )]
    blacklist: Option<PathBuf>,

    #[arg(
        long = "samplesLabel",
        value_name = "LABEL",
        num_args = 1..,
        help = "Labels for the samples. This will then be passed to plotHeatmap and plotProfile. The default is to use the file name of the sample. The sample labels should be separated by spaces and quoted if a label itself contains a space E.g. --samplesLabel label-1 \"label 2\""
    )]
    samples_label: Vec<String>,

    #[arg(
        long = "smartLabels",
        action = ArgAction::SetTrue,
        help = "Instead of manually specifying labels for the input bigWig and BED/GTF files, this causes deepTools to use the file name after removing the path and extension."
    )]
    smart_labels: bool,

    #[arg(
        long = "quiet",
        short = 'q',
        action = ArgAction::SetTrue,
        help = "Set to remove any warning or processing messages."
    )]
    quiet: bool,

    #[arg(
        long = "verbose",
        action = ArgAction::SetTrue,
        help = "Being VERY verbose in the status messages. --quiet will disable this."
    )]
    verbose: bool,

    #[arg(
        long = "scale",
        value_name = "FLOAT",
        default_value_t = 1.0,
        help = "If set, all values are multiplied by this number. (Default: 1)"
    )]
    scale: f64,

    #[arg(
        long = "numberOfProcessors",
        short = 'p',
        value_name = "INT",
        default_value = "1",
        help = "Number of processors to use. Type \"max/2\" to use half the maximum number of processors or \"max\" to use all available processors. (Default: 1)"
    )]
    number_of_processors: String,
}

#[derive(Debug, Args)]
struct GtfArgs {
    #[arg(
        long = "metagene",
        action = ArgAction::SetTrue,
        help = "When either a BED12 or GTF file are used to provide regions, perform the computation on the merged exons, rather than using the genomic interval defined by the 5-prime and 3-prime most transcript bound (i.e., columns 2 and 3 of a BED file). If a BED3 or BED6 file is used as input, then columns 2 and 3 are used as an exon. (Default: False)"
    )]
    keep_exons: bool,

    #[arg(
        long = "transcriptID",
        default_value = "transcript",
        help = "When a GTF file is used to provide regions, only entries with this value as their feature (column 3) will be processed as transcripts. (Default: transcript)"
    )]
    transcript_id: String,

    #[arg(
        long = "exonID",
        default_value = "exon",
        help = "When a GTF file is used to provide regions, only entries with this value as their feature (column 3) will be processed as exons. CDS would be another common value for this. (Default: exon)"
    )]
    exon_id: String,

    #[arg(
        long = "transcript_id_designator",
        default_value = "transcript_id",
        help = "Each region has an ID (e.g., ACTB) assigned to it, which for BED files is either column 4 (if it exists) or the interval bounds. For GTF files this is instead stored in the last column as a key:value pair (e.g., as 'transcript_id \"ACTB\"', for a key of transcript_id and a value of ACTB). In some cases it can be convenient to use a different identifier. To do so, set this to the desired key. (Default: transcript_id)"
    )]
    transcript_id_designator: String,
}

#[derive(Debug, Args)]
struct ScaleRegionsModeArgs {
    #[arg(
        long = "regionBodyLength",
        short = 'm',
        value_name = "INT bp",
        default_value_t = 1000,
        help = "Distance in bases to which all regions will be fit. (Default: 1000)"
    )]
    region_body_length: u32,

    #[arg(
        long = "startLabel",
        default_value = "TSS",
        help = "Label shown in the plot for the start of the region. Default is TSS (transcription start site), but could be changed to anything, e.g. \"peak start\". Note that this is only useful if you plan to plot the results yourself and not, for example, with plotHeatmap, which will override this. (Default: TSS)"
    )]
    start_label: String,

    #[arg(
        long = "endLabel",
        default_value = "TES",
        help = "Label shown in the plot for the region end. Default is TES (transcription end site). See the --startLabel option for more information. (Default: TES)"
    )]
    end_label: String,

    #[arg(
        long = "beforeRegionStartLength",
        short = 'b',
        alias = "upstream",
        default_value_t = 0,
        value_name = "INT bp",
        help = "Distance upstream of the start site of the regions defined in the region file. If the regions are genes, this would be the distance upstream of the transcription start site. (Default: 0)"
    )]
    upstream: i64,

    #[arg(
        long = "afterRegionStartLength",
        short = 'a',
        alias = "downstream",
        default_value_t = 0,
        value_name = "INT bp",
        help = "Distance downstream of the end site of the given regions. If the regions are genes, this would be the distance downstream of the transcription end site. (Default: 0)"
    )]
    downstream: i64,

    #[arg(
        long = "unscaled5prime",
        default_value_t = 0,
        help = "Number of bases at the 5-prime end of the region to exclude from scaling. By default, each region is scaled to a given length (see the --regionBodyLength option). In some cases it is useful to look at unscaled signals around region boundaries, so this setting specifies the number of unscaled bases on the 5-prime end of each boundary. (Default: 0)"
    )]
    unscaled_5_prime: u32,

    #[arg(
        long = "unscaled3prime",
        default_value_t = 0,
        help = "Like --unscaled5prime, but for the 3-prime end. (Default: 0)"
    )]
    unscaled_3_prime: u32,
}

#[derive(Debug, Args)]
struct ReferencePointModeArgs {
    #[arg(
        long = "referencePoint",
        default_value = "TSS",
        value_enum,
        help = "The reference point for the plotting could be either the region start (TSS), the region end (TES) or the center of the region. Note that regardless of what you specify, plotHeatmap/plotProfile will default to using \"TSS\" as the label. (Default: TSS)"
    )]
    reference_point: ReferencePoint,

    #[arg(
        long = "beforeRegionStartLength",
        short = 'b',
        alias = "upstream",
        default_value_t = 500,
        value_name = "INT bp",
        help = "Distance upstream of the reference-point selected. (Default: 500)"
    )]
    upstream: i64,

    #[arg(
        long = "afterRegionStartLength",
        short = 'a',
        alias = "downstream",
        default_value_t = 1500,
        value_name = "INT bp",
        help = "Distance downstream of the reference-point selected. (Default: 1500)"
    )]
    downstream: i64,

    #[arg(
        long = "nanAfterEnd",
        action = ArgAction::SetTrue,
        help = "If set, any values after the region end are discarded. This is useful to visualize the region end when not using the scale-regions mode and when the reference-point is set to the TSS."
    )]
    nan_after_end: bool,
}

impl Cli {
    pub fn build_config(self) -> anyhow::Result<Config> {
        match self.command {
            Command::ScaleRegions(cli) => cli.try_into_config(),
            Command::ReferencePoint(cli) => cli.try_into_config(),
        }
    }
}

impl ScaleRegionsCli {
    fn try_into_config(self) -> anyhow::Result<Config> {
        let general = self.general.normalize()?;
        let io = self.io.into_io_options(self.required)?;
        let gtf = self.gtf.into_gtf();

        let upstream = positive_distance(self.mode.upstream, "beforeRegionStartLength")?;
        let downstream = positive_distance(self.mode.downstream, "afterRegionStartLength")?;

        let mode = ModeConfig::ScaleRegions(ScaleRegionsOptions {
            region_body_length: self.mode.region_body_length,
            start_label: self.mode.start_label,
            end_label: self.mode.end_label,
            upstream,
            downstream,
            unscaled_5_prime: self.mode.unscaled_5_prime,
            unscaled_3_prime: self.mode.unscaled_3_prime,
        });

        Ok(Config {
            mode,
            io,
            general,
            gtf,
        })
    }
}

impl ReferencePointCli {
    fn try_into_config(self) -> anyhow::Result<Config> {
        let general = self.general.normalize()?;
        let io = self.io.into_io_options(self.required)?;
        let gtf = self.gtf.into_gtf();

        let upstream = positive_distance(self.mode.upstream, "beforeRegionStartLength")?;
        let downstream = positive_distance(self.mode.downstream, "afterRegionStartLength")?;

        if upstream == 0 && downstream == 0 {
            bail!(
                "Upstream and downstream regions are both set to 0. Nothing to output. Maybe you want to use the scale-regions mode?"
            );
        }

        let mode = ModeConfig::ReferencePoint(ReferencePointOptions {
            reference_point: self.mode.reference_point,
            upstream,
            downstream,
            nan_after_end: self.mode.nan_after_end,
        });

        Ok(Config {
            mode,
            io,
            general,
            gtf,
        })
    }
}

impl GeneralArgs {
    fn normalize(mut self) -> anyhow::Result<GeneralOptions> {
        if self.quiet {
            self.verbose = false;
        }

        let processor_request =
            parse_processor_request(&self.number_of_processors).with_context(|| {
                format!(
                    "Invalid value for --numberOfProcessors: {}",
                    self.number_of_processors
                )
            })?;

        Ok(GeneralOptions {
            bin_size: self.bin_size,
            sort_regions: self.sort_regions,
            sort_using: self.sort_using,
            sort_using_samples: if self.sort_using_samples.is_empty() {
                None
            } else {
                Some(self.sort_using_samples)
            },
            average_type_bins: self.average_type_bins,
            missing_data_as_zero: self.missing_data_as_zero,
            skip_zeros: self.skip_zeros,
            min_threshold: self.min_threshold,
            max_threshold: self.max_threshold,
            blacklist: self.blacklist,
            samples_label: if self.samples_label.is_empty() {
                None
            } else {
                Some(self.samples_label)
            },
            smart_labels: self.smart_labels,
            quiet: self.quiet,
            verbose: self.verbose,
            scale_factor: self.scale,
            number_of_processors: processor_request,
        })
    }
}

impl OutputArgs {
    fn into_io_options(self, required: RequiredArgs) -> anyhow::Result<IoOptions> {
        Ok(IoOptions {
            regions: required.regions,
            scores: required.scores,
            matrix_output: self.matrix_output,
            matrix_values_output: self.matrix_values,
            sorted_regions_output: self.sorted_regions,
        })
    }
}

impl GtfArgs {
    fn into_gtf(self) -> GtfOptions {
        GtfOptions {
            keep_exons: self.keep_exons,
            transcript_id: self.transcript_id,
            exon_id: self.exon_id,
            transcript_id_designator: self.transcript_id_designator,
        }
    }
}

fn positive_distance(value: i64, flag: &str) -> anyhow::Result<u32> {
    if value < 0 {
        let adjusted = value.saturating_abs();
        eprintln!("{flag} changed from {value} into {adjusted}");
        u32::try_from(adjusted).map_err(|_| anyhow!("{flag} overflow when converting to unsigned"))
    } else {
        u32::try_from(value).map_err(|_| anyhow!("{flag} overflow when converting to unsigned"))
    }
}

fn parse_processor_request(raw: &str) -> anyhow::Result<ProcessorRequest> {
    match raw {
        "max" | "MAX" => Ok(ProcessorRequest::Max),
        "max/2" | "MAX/2" => Ok(ProcessorRequest::MaxHalf),
        _ => {
            let value: u32 = raw
                .parse()
                .context("expected integer, \"max\" or \"max/2\"")?;
            Ok(ProcessorRequest::Fixed(value))
        }
    }
}
