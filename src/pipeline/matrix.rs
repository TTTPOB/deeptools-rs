use std::cmp::Ordering;

use crate::config::{GeneralOptions, SortUsing};
use crate::io::BedRecord;
use serde::{Serialize, Serializer};

/// Serialize a f64 as an integer when its value is a whole number, otherwise
/// as a float.  This matches the Python JSON output where `scale=1` is written
/// as `1` rather than `1.0`.
fn serialize_scale<S: Serializer>(value: &f64, s: S) -> Result<S::Ok, S::Error> {
    if value.fract() == 0.0 && value.is_finite() {
        s.serialize_i64(*value as i64)
    } else {
        s.serialize_f64(*value)
    }
}

/// Serializable metadata header mirroring the JSON preamble written by the
/// Python implementation of `computeMatrix`.
///
/// Field order matches Python's dict insertion order so that serde produces
/// identical JSON key ordering.
#[derive(Debug, Clone, Serialize)]
pub struct MatrixHeader {
    #[serde(rename = "upstream")]
    pub upstream: Vec<u32>,
    #[serde(rename = "downstream")]
    pub downstream: Vec<u32>,
    #[serde(rename = "body")]
    pub body: Vec<u32>,
    #[serde(rename = "bin size")]
    pub bin_size: Vec<u32>,
    #[serde(rename = "ref point")]
    pub ref_point: Vec<Option<String>>,
    #[serde(rename = "verbose")]
    pub verbose: bool,
    #[serde(rename = "bin avg type")]
    pub bin_avg_type: String,
    #[serde(rename = "missing data as zero")]
    pub missing_data_as_zero: bool,
    #[serde(rename = "min threshold")]
    pub min_threshold: Option<f64>,
    #[serde(rename = "max threshold")]
    pub max_threshold: Option<f64>,
    #[serde(rename = "scale", serialize_with = "serialize_scale")]
    pub scale: f64,
    #[serde(rename = "skip zeros")]
    pub skip_zeros: bool,
    #[serde(rename = "nan after end")]
    pub nan_after_end: bool,
    #[serde(rename = "proc number")]
    pub proc_number: u32,
    #[serde(rename = "sort regions")]
    pub sort_regions: String,
    #[serde(rename = "sort using")]
    pub sort_using: String,
    #[serde(rename = "unscaled 5 prime")]
    pub unscaled_5_prime: Vec<u32>,
    #[serde(rename = "unscaled 3 prime")]
    pub unscaled_3_prime: Vec<u32>,
    #[serde(rename = "group_labels")]
    pub group_labels: Vec<String>,
    #[serde(rename = "group_boundaries")]
    pub group_boundaries: Vec<usize>,
    #[serde(rename = "sample_labels")]
    pub sample_labels: Vec<String>,
    #[serde(rename = "sample_boundaries")]
    pub sample_boundaries: Vec<usize>,
}

#[derive(Debug, Clone)]
pub struct LayoutVectors {
    pub upstream: Vec<u32>,
    pub downstream: Vec<u32>,
    pub body: Vec<u32>,
    pub unscaled_5_prime: Vec<u32>,
    pub unscaled_3_prime: Vec<u32>,
    pub bin_size: Vec<u32>,
    pub ref_point: Vec<Option<String>>,
}

impl LayoutVectors {
    pub fn uniform(
        sample_count: usize,
        bin_size: u32,
        upstream: u32,
        downstream: u32,
        body: u32,
        unscaled_5_prime: u32,
        unscaled_3_prime: u32,
        ref_point: Option<String>,
    ) -> Self {
        let mut ref_points = Vec::with_capacity(sample_count);
        for _ in 0..sample_count {
            ref_points.push(ref_point.as_ref().cloned());
        }

        Self {
            upstream: vec![upstream; sample_count],
            downstream: vec![downstream; sample_count],
            body: vec![body; sample_count],
            unscaled_5_prime: vec![unscaled_5_prime; sample_count],
            unscaled_3_prime: vec![unscaled_3_prime; sample_count],
            bin_size: vec![bin_size; sample_count],
            ref_point: ref_points,
        }
    }
}

pub struct MatrixHeaderBuilder<'a> {
    general: &'a GeneralOptions,
    sample_labels: &'a [String],
    group_labels: &'a [String],
    group_counts: &'a [usize],
    thread_count: usize,
    sample_count: usize,
    nan_after_end: bool,
    layout: Option<LayoutVectors>,
    sample_boundaries: Option<Vec<usize>>,
}

impl<'a> MatrixHeaderBuilder<'a> {
    pub fn new(
        general: &'a GeneralOptions,
        sample_labels: &'a [String],
        group_labels: &'a [String],
        group_counts: &'a [usize],
        thread_count: usize,
        sample_count: usize,
        nan_after_end: bool,
    ) -> Self {
        Self {
            general,
            sample_labels,
            group_labels,
            group_counts,
            thread_count,
            sample_count,
            nan_after_end,
            layout: None,
            sample_boundaries: None,
        }
    }

    pub fn with_layout(mut self, layout: LayoutVectors) -> Self {
        debug_assert_eq!(layout.upstream.len(), self.sample_count);
        debug_assert_eq!(layout.downstream.len(), self.sample_count);
        debug_assert_eq!(layout.body.len(), self.sample_count);
        debug_assert_eq!(layout.unscaled_5_prime.len(), self.sample_count);
        debug_assert_eq!(layout.unscaled_3_prime.len(), self.sample_count);
        debug_assert_eq!(layout.bin_size.len(), self.sample_count);
        debug_assert_eq!(layout.ref_point.len(), self.sample_count);
        self.layout = Some(layout);
        self
    }

    pub fn with_sample_boundaries(mut self, boundaries: Vec<usize>) -> Self {
        debug_assert!(
            boundaries.len() == self.sample_count + 1,
            "sample boundaries must include start and end markers"
        );
        self.sample_boundaries = Some(boundaries);
        self
    }

    pub fn with_uniform_sample_boundaries(self, bins_per_sample: usize) -> Self {
        let sample_count = self.sample_count;
        let boundaries = sample_boundaries_uniform(sample_count, bins_per_sample);
        self.with_sample_boundaries(boundaries)
    }

    pub fn build(self) -> MatrixHeader {
        let layout = self
            .layout
            .expect("MatrixHeaderBuilder requires layout before build");
        let sample_boundaries = self
            .sample_boundaries
            .expect("MatrixHeaderBuilder requires sample boundaries before build");
        let group_boundaries = group_boundaries_from_counts(self.group_counts);

        MatrixHeader {
            upstream: layout.upstream,
            downstream: layout.downstream,
            body: layout.body,
            bin_size: layout.bin_size,
            ref_point: layout.ref_point,
            verbose: self.general.verbose,
            bin_avg_type: self.general.average_type_bins.to_string(),
            missing_data_as_zero: self.general.missing_data_as_zero,
            min_threshold: self.general.min_threshold,
            max_threshold: self.general.max_threshold,
            scale: self.general.scale_factor,
            skip_zeros: self.general.skip_zeros,
            nan_after_end: self.nan_after_end,
            proc_number: self.thread_count as u32,
            sort_regions: self.general.sort_regions.to_string(),
            sort_using: self.general.sort_using.to_string(),
            unscaled_5_prime: layout.unscaled_5_prime,
            unscaled_3_prime: layout.unscaled_3_prime,
            group_labels: self.group_labels.to_vec(),
            group_boundaries,
            sample_labels: self.sample_labels.to_vec(),
            sample_boundaries,
        }
    }
}

/// A single region row within the matrix output, tracking both the original
/// BED metadata and the per-sample binned signal values.
#[derive(Debug, Clone)]
pub struct MatrixRow {
    pub record: BedRecord,
    /// Flattened values in sample-major order: sample 0 bins, sample 1 bins, ...
    pub values: Vec<f64>,
    pub sample_count: usize,
    pub bin_count: usize,
    /// When metagene mode is used, stores the exon coordinates as (start, end) pairs
    /// for writing comma-separated coordinates in the output.
    pub exon_coords: Option<Vec<(u32, u32)>>,
}

/// Compute cumulative group boundaries from per-group row counts.
pub fn group_boundaries_from_counts(counts: &[usize]) -> Vec<usize> {
    let mut boundaries = Vec::with_capacity(counts.len() + 1);
    let mut running = 0usize;
    boundaries.push(0);
    for count in counts {
        running += *count;
        boundaries.push(running);
    }
    boundaries
}

/// Compute cumulative sample boundaries for a uniform number of bins per sample.
pub fn sample_boundaries_uniform(sample_count: usize, bins_per_sample: usize) -> Vec<usize> {
    let mut boundaries = Vec::with_capacity(sample_count + 1);
    for index in 0..=sample_count {
        boundaries.push(index * bins_per_sample);
    }
    boundaries
}

pub(crate) fn compute_sort_metric(
    row: &MatrixRow,
    sort_using: SortUsing,
    sample_list: Option<&[usize]>,
) -> f64 {
    match sort_using {
        SortUsing::RegionLength => {
            // Python sums the exon interval lengths for metagene regions;
            // for plain regions the single (start, end) span is identical to
            // `record.length()`.
            if let Some(ref exons) = row.exon_coords {
                exons.iter().map(|(s, e)| (e - s) as f64).sum()
            } else {
                row.record.length() as f64
            }
        }
        SortUsing::Sum => {
            let values = collect_values(row, sample_list);
            values.into_iter().fold(0.0f64, |acc, value| acc + value)
        }
        SortUsing::Mean => {
            let values = collect_values(row, sample_list);
            if values.is_empty() {
                f64::NAN
            } else {
                values.iter().copied().sum::<f64>() / values.len() as f64
            }
        }
        SortUsing::Median => {
            let mut values = collect_values(row, sample_list);
            if values.is_empty() {
                f64::NAN
            } else {
                values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
                let mid = values.len() / 2;
                if values.len() % 2 == 1 {
                    values[mid]
                } else {
                    (values[mid - 1] + values[mid]) / 2.0
                }
            }
        }
        SortUsing::Max => {
            let values = collect_values(row, sample_list);
            if values.is_empty() {
                f64::NAN
            } else {
                values
                    .into_iter()
                    .fold(f64::NEG_INFINITY, |acc, value| acc.max(value))
            }
        }
        SortUsing::Min => {
            let values = collect_values(row, sample_list);
            if values.is_empty() {
                f64::NAN
            } else {
                values
                    .into_iter()
                    .fold(f64::INFINITY, |acc, value| acc.min(value))
            }
        }
    }
}

fn collect_values(row: &MatrixRow, sample_list: Option<&[usize]>) -> Vec<f64> {
    let bin_count = row.bin_count;
    let mut values = Vec::new();
    match sample_list {
        Some(indices) => {
            for &sample_index in indices {
                let start = sample_index * bin_count;
                if start < row.values.len() {
                    let end = (start + bin_count).min(row.values.len());
                    values.extend(
                        row.values[start..end]
                            .iter()
                            .copied()
                            .filter(|v| !v.is_nan()),
                    );
                }
            }
        }
        None => {
            values.extend(row.values.iter().copied().filter(|v| !v.is_nan()));
        }
    }
    values
}

pub(crate) fn compare_ascending(left: f64, right: f64) -> Ordering {
    match (left.is_nan(), right.is_nan()) {
        (true, true) => Ordering::Equal,
        (true, false) => Ordering::Greater,
        (false, true) => Ordering::Less,
        (false, false) => left.partial_cmp(&right).unwrap_or_else(|| Ordering::Equal),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::{BedRecord, Strand};
    use std::sync::Arc;

    /// Helper to create a minimal MatrixRow for testing sort metrics.
    fn make_row(
        start: u32,
        end: u32,
        exon_coords: Option<Vec<(u32, u32)>>,
        values: Vec<f64>,
    ) -> MatrixRow {
        MatrixRow {
            record: BedRecord {
                chrom: Arc::from("chr1"),
                start,
                end,
                name: None,
                score: None,
                score_raw: None,
                strand: Strand::Unstranded,
                strand_raw: None,
                extra_fields: Vec::new(),
            },
            sample_count: 1,
            bin_count: values.len(),
            values,
            exon_coords,
        }
    }

    // ---------------------------------------------------------------
    // Issue 3a: sortUsing=region_length must sum exon lengths in metagene mode
    // ---------------------------------------------------------------

    #[test]
    fn region_length_without_exons_uses_genomic_span() {
        // Plain region: 1000-5000 => length 4000
        let row = make_row(1000, 5000, None, vec![1.0; 10]);
        let metric = compute_sort_metric(&row, SortUsing::RegionLength, None);
        assert_eq!(metric, 4000.0);
    }

    #[test]
    fn region_length_with_exons_sums_exon_lengths() {
        // Genomic span 1000-5000 (4000 bp), but only two exons:
        //   exon1: 1000-2000 (1000 bp)
        //   exon2: 3000-4000 (1000 bp)
        // Total exon length = 2000, NOT 4000.
        let exons = vec![(1000, 2000), (3000, 4000)];
        let row = make_row(1000, 5000, Some(exons), vec![1.0; 10]);
        let metric = compute_sort_metric(&row, SortUsing::RegionLength, None);
        assert_eq!(metric, 2000.0);
    }

    #[test]
    fn region_length_single_exon_matches_span() {
        // When there's a single exon spanning the whole region, the sum
        // equals the genomic span.
        let exons = vec![(1000, 5000)];
        let row = make_row(1000, 5000, Some(exons), vec![1.0; 10]);
        let metric = compute_sort_metric(&row, SortUsing::RegionLength, None);
        assert_eq!(metric, 4000.0);
    }

    #[test]
    fn region_length_empty_exon_list_returns_zero() {
        // Edge case: exon_coords is Some but empty.
        let row = make_row(1000, 5000, Some(vec![]), vec![1.0; 10]);
        let metric = compute_sort_metric(&row, SortUsing::RegionLength, None);
        assert_eq!(metric, 0.0);
    }

    // ---------------------------------------------------------------
    // Issue 3b: group_boundaries_from_counts with zero-count groups
    // ---------------------------------------------------------------

    #[test]
    fn group_boundaries_basic() {
        // Two groups with 3 and 5 rows
        let boundaries = group_boundaries_from_counts(&[3, 5]);
        assert_eq!(boundaries, vec![0, 3, 8]);
    }

    #[test]
    fn group_boundaries_with_empty_first_group() {
        // First group has 0 rows (e.g. all filtered by skipZeros),
        // second group has 5 rows.
        let boundaries = group_boundaries_from_counts(&[0, 5]);
        assert_eq!(boundaries, vec![0, 0, 5]);
    }

    #[test]
    fn group_boundaries_with_empty_middle_group() {
        let boundaries = group_boundaries_from_counts(&[3, 0, 5]);
        assert_eq!(boundaries, vec![0, 3, 3, 8]);
    }

    #[test]
    fn group_boundaries_all_empty() {
        let boundaries = group_boundaries_from_counts(&[0, 0, 0]);
        assert_eq!(boundaries, vec![0, 0, 0, 0]);
    }

    #[test]
    fn group_boundaries_single_group() {
        let boundaries = group_boundaries_from_counts(&[10]);
        assert_eq!(boundaries, vec![0, 10]);
    }

    #[test]
    fn group_boundaries_no_groups() {
        let boundaries = group_boundaries_from_counts(&[]);
        assert_eq!(boundaries, vec![0]);
    }
}

#[cfg(test)]
impl MatrixHeader {
    pub fn default_for_test(group_counts: Vec<usize>) -> Self {
        let mut boundaries = vec![0usize];
        for &count in &group_counts {
            boundaries.push(boundaries.last().unwrap() + count);
        }
        Self {
            upstream: vec![0],
            downstream: vec![0],
            body: vec![0],
            bin_size: vec![10],
            ref_point: vec![None],
            verbose: false,
            bin_avg_type: "mean".into(),
            missing_data_as_zero: false,
            min_threshold: None,
            max_threshold: None,
            scale: 1.0,
            skip_zeros: false,
            nan_after_end: false,
            proc_number: 1,
            sort_regions: "keep".into(),
            sort_using: "mean".into(),
            unscaled_5_prime: vec![0],
            unscaled_3_prime: vec![0],
            group_labels: (0..group_counts.len())
                .map(|i| format!("group{i}"))
                .collect(),
            group_boundaries: boundaries,
            sample_labels: vec!["test".into()],
            sample_boundaries: vec![0, 1],
        }
    }
}
