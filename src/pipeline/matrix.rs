use crate::io::BedRecord;
use serde::Serialize;

/// Serializable metadata header mirroring the JSON preamble written by the
/// Python implementation of `computeMatrix`.
#[derive(Debug, Clone, Serialize)]
pub struct MatrixHeader {
    #[serde(rename = "verbose")]
    pub verbose: bool,
    #[serde(rename = "scale")]
    pub scale: f64,
    #[serde(rename = "skip zeros")]
    pub skip_zeros: bool,
    #[serde(rename = "nan after end")]
    pub nan_after_end: bool,
    #[serde(rename = "sort using")]
    pub sort_using: String,
    #[serde(rename = "unscaled 5 prime")]
    pub unscaled_5_prime: Vec<u32>,
    #[serde(rename = "body")]
    pub body: Vec<u32>,
    #[serde(rename = "sample_labels")]
    pub sample_labels: Vec<String>,
    #[serde(rename = "downstream")]
    pub downstream: Vec<u32>,
    #[serde(rename = "unscaled 3 prime")]
    pub unscaled_3_prime: Vec<u32>,
    #[serde(rename = "group_labels")]
    pub group_labels: Vec<String>,
    #[serde(rename = "bin size")]
    pub bin_size: Vec<u32>,
    #[serde(rename = "upstream")]
    pub upstream: Vec<u32>,
    #[serde(rename = "group_boundaries")]
    pub group_boundaries: Vec<usize>,
    #[serde(rename = "sample_boundaries")]
    pub sample_boundaries: Vec<usize>,
    #[serde(rename = "missing data as zero")]
    pub missing_data_as_zero: bool,
    #[serde(rename = "ref point")]
    pub ref_point: Vec<Option<String>>,
    #[serde(rename = "min threshold")]
    pub min_threshold: Option<f64>,
    #[serde(rename = "sort regions")]
    pub sort_regions: String,
    #[serde(rename = "proc number")]
    pub proc_number: u32,
    #[serde(rename = "bin avg type")]
    pub bin_avg_type: String,
    #[serde(rename = "max threshold")]
    pub max_threshold: Option<f64>,
}

/// A single region row within the matrix output, tracking both the original
/// BED metadata and the per-sample binned signal values.
#[derive(Debug, Clone)]
pub struct MatrixRow {
    pub record: BedRecord,
    /// Matrix values organised as `sample -> bin`.
    pub values: Vec<Vec<f32>>,
}

impl MatrixRow {
    /// Returns a flattened view of the row values in sample-major order.
    pub fn flattened_values(&self) -> Vec<f32> {
        self.values
            .iter()
            .flat_map(|sample| sample.iter().copied())
            .collect()
    }
}

/// In-memory representation of the computeMatrix result required to serialise
/// the gzipped matrix as well as auxiliary artifacts such as the plain matrix
/// table or sorted BED output.
#[derive(Debug, Clone)]
pub struct MatrixData {
    pub header: MatrixHeader,
    pub rows: Vec<MatrixRow>,
    pub bin_count: usize,
    pub sample_count: usize,
}
