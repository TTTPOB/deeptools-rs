use std::collections::{HashMap, HashSet};
use std::convert::TryFrom;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context, Result, anyhow, bail};
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

use crate::config::{AverageTypeBins, GeneralOptions, IoOptions, ReferencePointOptions};
use crate::io::{BedReadError, BedRecord, BigWigReader};
use crate::pipeline::zones::ReferencePointPlan;

use super::matrix::{MatrixData, MatrixHeader, MatrixRow};

struct RegionTask {
    index: usize,
    group_index: usize,
    record: BedRecord,
}

struct RegionResult {
    index: usize,
    group_index: usize,
    row: Option<MatrixRow>,
}

pub fn run(
    general: &GeneralOptions,
    io: &IoOptions,
    options: &ReferencePointOptions,
) -> Result<MatrixData> {
    let bin_size = general.bin_size;
    ensure_positive(bin_size, "binSize")?;

    ensure_multiple(bin_size, options.upstream, "beforeRegionStartLength")?;
    ensure_multiple(bin_size, options.downstream, "afterRegionStartLength")?;

    let upstream_bins = (options.upstream / bin_size) as usize;
    let downstream_bins = (options.downstream / bin_size) as usize;
    let total_bins = upstream_bins + downstream_bins;

    if total_bins == 0 {
        bail!("Reference-point mode requires at least one upstream or downstream bin");
    }

    let sample_labels = derive_sample_labels(&io.scores, general)?;
    let sample_count = io.scores.len();

    let groups = load_groups(&io.regions)?;
    let mut group_counts = vec![0usize; groups.len()];

    let mut tasks = Vec::new();
    for (group_index, group) in groups.iter().enumerate() {
        for record in &group.records {
            let index = tasks.len();
            tasks.push(RegionTask {
                index,
                group_index,
                record: record.clone(),
            });
        }
    }

    let task_count = tasks.len();
    let thread_count = std::cmp::max(1, general.number_of_processors.resolve() as usize);
    let mut rows = Vec::with_capacity(task_count);

    if task_count > 0 {
        let pool = ThreadPoolBuilder::new()
            .num_threads(thread_count)
            .build()
            .context("Failed to initialise rayon thread pool for reference-point scheduling")?;

        let sample_paths = Arc::new(io.scores.clone());
        let collected = pool.install(|| {
            tasks
                .into_par_iter()
                .map_init(
                    || WorkerSamples::new(sample_paths.clone()),
                    |state, task| {
                        let RegionTask {
                            index,
                            group_index,
                            record,
                        } = task;

                        let samples = state.samples()?;

                        let plan = ReferencePointPlan::reference_point(
                            &record,
                            options.reference_point,
                            bin_size,
                            upstream_bins,
                            downstream_bins,
                        );

                        let maybe_values =
                            compute_row(samples.as_mut_slice(), &record, &plan, options, general)?;
                        let row = maybe_values.map(|values| MatrixRow { record, values });

                        Ok(RegionResult {
                            index,
                            group_index,
                            row,
                        })
                    },
                )
                .collect::<Vec<_>>()
        });

        let mut collected = collected.into_iter().collect::<Result<Vec<_>>>()?;
        collected.sort_by_key(|entry| entry.index);

        for entry in collected {
            if let Some(row) = entry.row {
                group_counts[entry.group_index] += 1;
                rows.push(row);
            }
        }
    }

    let group_labels: Vec<String> = groups.iter().map(|group| group.label.clone()).collect();
    let group_boundaries = MatrixData::group_boundaries_from_counts(&group_counts);
    let sample_boundaries = MatrixData::sample_boundaries_uniform(sample_count, total_bins);

    let proc_number = thread_count as u32;

    let header = MatrixHeader {
        verbose: general.verbose,
        scale: general.scale_factor,
        skip_zeros: general.skip_zeros,
        nan_after_end: options.nan_after_end,
        sort_using: general.sort_using.to_string(),
        unscaled_5_prime: vec![0; sample_count],
        body: vec![0; sample_count],
        sample_labels: sample_labels.clone(),
        downstream: vec![options.downstream; sample_count],
        unscaled_3_prime: vec![0; sample_count],
        group_labels,
        bin_size: vec![bin_size; sample_count],
        upstream: vec![options.upstream; sample_count],
        group_boundaries,
        sample_boundaries,
        missing_data_as_zero: general.missing_data_as_zero,
        ref_point: vec![Some(options.reference_point.to_string()); sample_count],
        min_threshold: general.min_threshold,
        sort_regions: general.sort_regions.to_string(),
        proc_number,
        bin_avg_type: general.average_type_bins.to_string(),
        max_threshold: general.max_threshold,
    };

    let sort_sample_indices =
        normalize_sort_sample_indices(general.sort_using_samples.as_ref(), sample_count)?;

    let mut matrix = MatrixData {
        header,
        rows,
        bin_count: total_bins,
        sample_count,
    };

    matrix.sort_groups(
        general.sort_regions,
        general.sort_using,
        sort_sample_indices.as_deref(),
    )?;

    matrix.prune_zero_rows();

    Ok(matrix)
}

struct Group {
    label: String,
    records: Vec<BedRecord>,
}

struct Sample {
    path: PathBuf,
    reader: BigWigReader,
    chrom_lengths: HashMap<String, u32>,
}

impl Sample {
    fn open(path: &Path) -> Result<Self> {
        let reader = BigWigReader::open(path)
            .with_context(|| format!("Failed to open bigWig file '{}'", path.display()))?;
        let chrom_lengths = reader
            .chroms()
            .iter()
            .map(|chrom| (chrom.name.clone(), chrom.length))
            .collect::<HashMap<_, _>>();

        Ok(Self {
            path: path.to_path_buf(),
            reader,
            chrom_lengths,
        })
    }

    fn chrom_length(&self, chrom: &str) -> Option<u32> {
        self.chrom_lengths.get(chrom).copied()
    }
}

struct WorkerSamples {
    samples: Result<Vec<Sample>, String>,
}

impl WorkerSamples {
    fn new(paths: Arc<Vec<PathBuf>>) -> Self {
        let samples = open_samples(paths.as_ref()).map_err(|err| err.to_string());
        Self { samples }
    }

    fn samples(&mut self) -> Result<&mut Vec<Sample>> {
        match &mut self.samples {
            Ok(samples) => Ok(samples),
            Err(message) => Err(anyhow!(message.clone())),
        }
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

fn load_groups(paths: &[PathBuf]) -> Result<Vec<Group>> {
    let mut groups = Vec::new();
    let mut seen_labels = HashSet::new();
    for path in paths {
        let mut file_groups = parse_grouped_bed(path, &mut seen_labels)
            .map_err(anyhow::Error::new)
            .with_context(|| format!("Failed to parse regions file '{}'", path.display()))?;
        groups.append(&mut file_groups);
    }

    Ok(groups)
}

fn parse_grouped_bed(
    path: &Path,
    seen_labels: &mut HashSet<String>,
) -> Result<Vec<Group>, BedReadError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let default_label = label_from_path(path);
    let mut groups = Vec::new();
    let mut current_records = Vec::new();

    for (line_number, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with('#') {
            finalize_group(
                trimmed.strip_prefix('#').unwrap_or("").trim(),
                &default_label,
                &mut current_records,
                &mut groups,
                seen_labels,
            );
            continue;
        }

        match BedRecord::parse(trimmed) {
            Ok(record) => current_records.push(record),
            Err(message) => {
                return Err(BedReadError::Parse {
                    line_number: line_number + 1,
                    message,
                    line,
                });
            }
        }
    }

    if !current_records.is_empty() {
        finalize_group(
            "",
            &default_label,
            &mut current_records,
            &mut groups,
            seen_labels,
        );
    } else if groups.is_empty() {
        let label = next_unique_label("", &default_label, seen_labels);
        groups.push(Group {
            label,
            records: Vec::new(),
        });
    }

    Ok(groups)
}

fn finalize_group(
    raw_label: &str,
    default_label: &str,
    current_records: &mut Vec<BedRecord>,
    groups: &mut Vec<Group>,
    seen_labels: &mut HashSet<String>,
) {
    if current_records.is_empty() {
        return;
    }

    let label = next_unique_label(raw_label, default_label, seen_labels);
    let records = std::mem::take(current_records);
    groups.push(Group { label, records });
}

fn next_unique_label(
    raw_label: &str,
    default_label: &str,
    seen_labels: &mut HashSet<String>,
) -> String {
    let candidate = if raw_label.trim().is_empty() {
        default_label.to_string()
    } else {
        raw_label.trim().to_string()
    };

    if seen_labels.insert(candidate.clone()) {
        return candidate;
    }

    let mut suffix = 1;
    loop {
        let proposal = format!("{}_{}", candidate, suffix);
        if seen_labels.insert(proposal.clone()) {
            return proposal;
        }
        suffix += 1;
    }
}

fn open_samples(paths: &[PathBuf]) -> Result<Vec<Sample>> {
    let mut samples = Vec::with_capacity(paths.len());
    for path in paths {
        samples.push(Sample::open(path)?);
    }
    Ok(samples)
}

fn normalize_sort_sample_indices(
    raw: Option<&Vec<usize>>,
    sample_count: usize,
) -> Result<Option<Vec<usize>>> {
    let Some(raw_indices) = raw else {
        return Ok(None);
    };

    if raw_indices.is_empty() {
        return Ok(None);
    }

    let mut normalized = Vec::with_capacity(raw_indices.len());
    for &value in raw_indices {
        if value == 0 || value > sample_count {
            bail!(
                "The value {} for --sortUsingSamples is not valid. Only values from 1 to {} are allowed.",
                value,
                sample_count
            );
        }
        normalized.push(value - 1);
    }

    Ok(Some(normalized))
}

fn derive_sample_labels(paths: &[PathBuf], general: &GeneralOptions) -> Result<Vec<String>> {
    if let Some(labels) = &general.samples_label {
        if labels.len() != paths.len() {
            bail!(
                "--samplesLabel expects {} entries but {} were provided",
                paths.len(),
                labels.len()
            );
        }
        return Ok(labels.clone());
    }

    if general.smart_labels {
        return Ok(paths.iter().map(|path| smart_label(path)).collect());
    }

    Ok(paths.iter().map(|path| default_label(path)).collect())
}

fn default_label(path: &Path) -> String {
    path.file_name()
        .and_then(|stem| stem.to_str())
        .map(|value| value.to_string())
        .unwrap_or_else(|| path.display().to_string())
}

fn smart_label(path: &Path) -> String {
    let mut label = default_label(path);
    if let Some(stripped) = label.strip_suffix(".bw") {
        label = stripped.to_string();
    }
    label
}

fn label_from_path(path: &Path) -> String {
    default_label(path)
}

fn compute_row(
    samples: &mut [Sample],
    record: &BedRecord,
    plan: &ReferencePointPlan,
    options: &ReferencePointOptions,
    general: &GeneralOptions,
) -> Result<Option<Vec<Vec<f32>>>> {
    let mut per_sample = Vec::with_capacity(samples.len());
    for sample in samples.iter_mut() {
        let values = compute_sample_bins(sample, record, plan, options, general)?;
        per_sample.push(values);
    }

    if should_skip_row(&per_sample, general) {
        return Ok(None);
    }

    Ok(Some(per_sample))
}

fn should_skip_row(values: &[Vec<f32>], general: &GeneralOptions) -> bool {
    if general.skip_zeros {
        let mut all_zero = true;
        for value in values.iter().flat_map(|sample| sample.iter()) {
            if value.is_nan() {
                continue;
            }
            if *value != 0.0 {
                all_zero = false;
                break;
            }
        }
        if all_zero {
            return true;
        }
    }

    if let Some(min_threshold) = general.min_threshold {
        if values
            .iter()
            .flat_map(|sample| sample.iter())
            .filter(|value| !value.is_nan())
            .any(|value| (*value as f64) <= min_threshold)
        {
            return true;
        }
    }

    if let Some(max_threshold) = general.max_threshold {
        if values
            .iter()
            .flat_map(|sample| sample.iter())
            .filter(|value| !value.is_nan())
            .any(|value| (*value as f64) >= max_threshold)
        {
            return true;
        }
    }

    false
}

fn compute_sample_bins(
    sample: &mut Sample,
    record: &BedRecord,
    plan: &ReferencePointPlan,
    options: &ReferencePointOptions,
    general: &GeneralOptions,
) -> Result<Vec<f32>> {
    let bin_count = plan.bins.len();
    let chrom_length = match sample.chrom_length(&record.chrom) {
        Some(length) => length,
        None => {
            return Ok(vec![f32::NAN; bin_count]);
        }
    };

    let window_span = plan.window_end - plan.window_start;
    if window_span <= 0 {
        return Ok(vec![f32::NAN; bin_count]);
    }

    let window_len =
        usize::try_from(window_span).expect("reference-point window span exceeds usize");
    let default_fill = if general.missing_data_as_zero {
        0.0f32
    } else {
        f32::NAN
    };

    let mut coverage = vec![default_fill; window_len];
    let fetch_start = clamp_coordinate(plan.window_start, chrom_length);
    let fetch_end = clamp_coordinate(plan.window_end, chrom_length);

    if fetch_start < fetch_end {
        let intervals = sample
            .reader
            .values(&record.chrom, fetch_start, fetch_end)
            .map_err(anyhow::Error::new)
            .with_context(|| {
                format!(
                    "Failed to read bigWig intervals for '{}' in '{}'",
                    record.chrom,
                    sample.path.display()
                )
            })?;

        let base_offset = plan.window_start;
        for interval in intervals {
            let overlap_start = interval.start.max(fetch_start);
            let overlap_end = interval.end.min(fetch_end);
            if overlap_start >= overlap_end {
                continue;
            }
            let rel_start = (overlap_start as i64 - base_offset) as usize;
            let rel_end = (overlap_end as i64 - base_offset) as usize;
            coverage[rel_start..rel_end].fill(interval.value);
        }
    }

    let mut bins = Vec::with_capacity(bin_count);
    for bin in &plan.bins {
        let start_idx = index_from_coordinate(bin.start, plan.window_start, window_len);
        let end_idx = index_from_coordinate(bin.end, plan.window_start, window_len);

        let mut value = if start_idx < end_idx {
            aggregate_slice(&coverage[start_idx..end_idx], general.average_type_bins)
        } else {
            None
        };

        if value.is_none() && general.missing_data_as_zero {
            value = Some(0.0);
        }

        let mut value = value.unwrap_or(f32::NAN);

        if options.nan_after_end && bin.beyond_region {
            value = f32::NAN;
        }

        if value.is_finite() {
            value *= general.scale_factor as f32;
        }

        bins.push(value);
    }

    Ok(bins)
}

fn index_from_coordinate(value: i64, base: i64, window_len: usize) -> usize {
    if value <= base {
        return 0;
    }
    let diff = value - base;
    let idx = usize::try_from(diff).unwrap_or(window_len);
    idx.min(window_len)
}

fn aggregate_slice(slice: &[f32], average_type: AverageTypeBins) -> Option<f32> {
    let mut values: Vec<f64> = slice
        .iter()
        .copied()
        .filter(|value| !value.is_nan())
        .map(|value| value as f64)
        .collect();

    if values.is_empty() {
        return None;
    }

    match average_type {
        AverageTypeBins::Mean => {
            let mean = values.iter().sum::<f64>() / values.len() as f64;
            Some(mean as f32)
        }
        AverageTypeBins::Median => {
            values.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let mid = values.len() / 2;
            if values.len() % 2 == 0 {
                Some(((values[mid - 1] + values[mid]) / 2.0) as f32)
            } else {
                Some(values[mid] as f32)
            }
        }
        AverageTypeBins::Min => values
            .into_iter()
            .reduce(f64::min)
            .map(|value| value as f32),
        AverageTypeBins::Max => values
            .into_iter()
            .reduce(f64::max)
            .map(|value| value as f32),
        AverageTypeBins::Sum => Some(values.into_iter().sum::<f64>() as f32),
        AverageTypeBins::Std => {
            let mean = values.iter().sum::<f64>() / values.len() as f64;
            let variance = values
                .iter()
                .map(|value| {
                    let delta = value - mean;
                    delta * delta
                })
                .sum::<f64>()
                / values.len() as f64;
            Some(variance.sqrt() as f32)
        }
    }
}

fn clamp_coordinate(value: i64, chrom_length: u32) -> u32 {
    value
        .max(0)
        .min(chrom_length as i64)
        .try_into()
        .unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn index_from_coordinate_bounds_checks() {
        let base = 100;
        let window_len = 50;
        assert_eq!(index_from_coordinate(90, base, window_len), 0);
        assert_eq!(index_from_coordinate(100, base, window_len), 0);
        assert_eq!(index_from_coordinate(125, base, window_len), 25);
        assert_eq!(index_from_coordinate(200, base, window_len), window_len);
    }

    #[test]
    fn aggregate_slice_ignores_nans() {
        let data = [1.0, f32::NAN, 3.0, 5.0];
        let mean = aggregate_slice(&data, AverageTypeBins::Mean).unwrap();
        assert!((mean - 3.0).abs() < 1e-6);

        let max = aggregate_slice(&data, AverageTypeBins::Max).unwrap();
        assert_eq!(max, 5.0);

        let median = aggregate_slice(&data, AverageTypeBins::Median).unwrap();
        assert_eq!(median, 3.0);
    }
}
