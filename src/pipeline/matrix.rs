use std::cmp::Ordering;

use anyhow::{Result, bail};

use crate::config::{GeneralOptions, SortRegions, SortUsing};
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
        let boundaries = MatrixData::sample_boundaries_uniform(sample_count, bins_per_sample);
        self.with_sample_boundaries(boundaries)
    }

    pub fn build(self) -> MatrixHeader {
        let layout = self
            .layout
            .expect("MatrixHeaderBuilder requires layout before build");
        let sample_boundaries = self
            .sample_boundaries
            .expect("MatrixHeaderBuilder requires sample boundaries before build");
        let group_boundaries = MatrixData::group_boundaries_from_counts(self.group_counts);

        MatrixHeader {
            verbose: self.general.verbose,
            scale: self.general.scale_factor,
            skip_zeros: self.general.skip_zeros,
            nan_after_end: self.nan_after_end,
            sort_using: self.general.sort_using.to_string(),
            unscaled_5_prime: layout.unscaled_5_prime,
            body: layout.body,
            sample_labels: self.sample_labels.to_vec(),
            downstream: layout.downstream,
            unscaled_3_prime: layout.unscaled_3_prime,
            group_labels: self.group_labels.to_vec(),
            bin_size: layout.bin_size,
            upstream: layout.upstream,
            group_boundaries,
            sample_boundaries,
            missing_data_as_zero: self.general.missing_data_as_zero,
            ref_point: layout.ref_point,
            min_threshold: self.general.min_threshold,
            sort_regions: self.general.sort_regions.to_string(),
            proc_number: self.thread_count as u32,
            bin_avg_type: self.general.average_type_bins.to_string(),
            max_threshold: self.general.max_threshold,
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

impl MatrixRow {
    /// Returns a clone of the flat values in sample-major order.
    pub fn flattened_values(&self) -> Vec<f64> {
        self.values.clone()
    }
}

/// Summary information for each group to support diagnostics and future features
/// like clustering or reporting.
#[derive(Debug, Clone)]
pub struct GroupStats {
    pub label: String,
    pub row_count: usize,
    pub sample_count: usize,
    pub start_row: usize,
    pub end_row: usize,
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

impl MatrixData {
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

    /// Compute cumulative sample boundaries from per-sample bin counts.
    pub fn sample_boundaries_from_counts(bin_counts: &[usize]) -> Vec<usize> {
        let mut boundaries = Vec::with_capacity(bin_counts.len() + 1);
        let mut running = 0usize;
        boundaries.push(0);
        for count in bin_counts {
            running += *count;
            boundaries.push(running);
        }
        boundaries
    }

    /// Sort rows within each group according to the configured method, mirroring
    /// the behaviour of DeepTools' `_matrix.sort_groups`.
    pub fn sort_groups(
        &mut self,
        sort_method: SortRegions,
        sort_using: SortUsing,
        sample_list: Option<&[usize]>,
    ) -> Result<()> {
        if matches!(sort_method, SortRegions::No | SortRegions::Keep) {
            return Ok(());
        }

        if self.rows.is_empty() {
            return Ok(());
        }

        let sample_list = sample_list.filter(|indices| !indices.is_empty());
        if let Some(indices) = sample_list {
            for &index in indices {
                if index >= self.sample_count {
                    bail!(
                        "The value {} for --sortUsingSamples is not valid. Only values from 1 to {} are allowed.",
                        index + 1,
                        self.sample_count
                    );
                }
            }
        }

        let metrics: Vec<f64> = self
            .rows
            .iter()
            .map(|row| compute_sort_metric(row, sort_using, sample_list))
            .collect();

        // Move rows out of self and wrap in Option so we can take() individual
        // entries without cloning any MatrixRow (avoids duplicating ~400 MB for
        // realistic inputs).
        let old_rows = std::mem::take(&mut self.rows);
        let mut takeable: Vec<Option<MatrixRow>> = old_rows.into_iter().map(Some).collect();
        let old_len = takeable.len();
        let mut reordered = Vec::with_capacity(old_len);

        for window in self.header.group_boundaries.windows(2) {
            let start = window[0];
            let end = window[1];
            if start >= end {
                continue;
            }

            let mut indices: Vec<usize> = (start..end).collect();
            indices.sort_by(|&left, &right| compare_ascending(metrics[left], metrics[right]));
            if matches!(sort_method, SortRegions::Descend) {
                indices.reverse();
            }

            for index in indices {
                if let Some(row) = takeable[index].take() {
                    reordered.push(row);
                }
            }
        }

        self.rows = reordered;
        Ok(())
    }

    /// Remove rows containing only zeros (ignoring NaNs) when skip-zero behaviour
    /// is requested. Updates group boundaries to reflect the filtered matrix.
    pub fn prune_zero_rows(&mut self) {
        if !self.header.skip_zeros || self.rows.is_empty() {
            return;
        }

        let mut filtered_rows = Vec::with_capacity(self.rows.len());
        let mut group_counts = Vec::with_capacity(self.header.group_labels.len());

        for window in self.header.group_boundaries.windows(2) {
            let start = window[0];
            let end = window[1];
            if start >= end {
                group_counts.push(0);
                continue;
            }

            let mut retained = 0usize;
            for row in &self.rows[start..end] {
                if row_is_all_zero(row) {
                    continue;
                }
                filtered_rows.push(row.clone());
                retained += 1;
            }
            group_counts.push(retained);
        }

        self.rows = filtered_rows;
        self.header.group_boundaries = MatrixData::group_boundaries_from_counts(&group_counts);
    }

    /// Provide high-level statistics per group for diagnostics or logging.
    pub fn group_stats(&self) -> Vec<GroupStats> {
        let mut stats = Vec::with_capacity(self.header.group_labels.len());
        for (idx, window) in self.header.group_boundaries.windows(2).enumerate() {
            let start = window[0];
            let end = window[1];
            let row_count = end.saturating_sub(start);
            stats.push(GroupStats {
                label: self
                    .header
                    .group_labels
                    .get(idx)
                    .cloned()
                    .unwrap_or_else(|| format!("group {}", idx)),
                row_count,
                sample_count: self.sample_count,
                start_row: start,
                end_row: end,
            });
        }
        stats
    }
}

fn compute_sort_metric(
    row: &MatrixRow,
    sort_using: SortUsing,
    sample_list: Option<&[usize]>,
) -> f64 {
    match sort_using {
        SortUsing::RegionLength => row.record.length() as f64,
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
                    values.extend(row.values[start..end].iter().copied().filter(|v| !v.is_nan()));
                }
            }
        }
        None => {
            values.extend(row.values.iter().copied().filter(|v| !v.is_nan()));
        }
    }
    values
}

fn compare_ascending(left: f64, right: f64) -> Ordering {
    match (left.is_nan(), right.is_nan()) {
        (true, true) => Ordering::Equal,
        (true, false) => Ordering::Greater,
        (false, true) => Ordering::Less,
        (false, false) => left.partial_cmp(&right).unwrap_or_else(|| Ordering::Equal),
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
            verbose: false,
            scale: 1.0,
            skip_zeros: false,
            nan_after_end: false,
            sort_using: "mean".into(),
            unscaled_5_prime: vec![0],
            body: vec![0],
            sample_labels: vec!["test".into()],
            downstream: vec![0],
            unscaled_3_prime: vec![0],
            group_labels: (0..group_counts.len()).map(|i| format!("group{i}")).collect(),
            bin_size: vec![10],
            upstream: vec![0],
            group_boundaries: boundaries,
            sample_boundaries: vec![0, 1],
            missing_data_as_zero: false,
            ref_point: vec![None],
            min_threshold: None,
            sort_regions: "keep".into(),
            proc_number: 1,
            bin_avg_type: "mean".into(),
            max_threshold: None,
        }
    }
}

fn row_is_all_zero(row: &MatrixRow) -> bool {
    for &value in &row.values {
        if value.is_nan() {
            continue;
        }
        if value != 0.0 {
            return false;
        }
    }
    true
}
