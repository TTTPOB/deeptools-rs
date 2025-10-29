use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};

use crate::config::{GeneralOptions, IoOptions, ReferencePoint, ReferencePointOptions};
use crate::io::{BedReadError, BedRecord, BigWigReader, BigWigValue, Strand};

use super::matrix::{MatrixData, MatrixHeader, MatrixRow};

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
    let mut samples = open_samples(&io.scores)?;
    let sample_count = samples.len();

    let groups = load_groups(&io.regions)?;
    let mut rows = Vec::new();
    let mut group_counts = vec![0usize; groups.len()];

    for (group_index, group) in groups.iter().enumerate() {
        for record in &group.records {
            match compute_row(
                &mut samples,
                record,
                options,
                general,
                bin_size,
                upstream_bins,
                downstream_bins,
            )? {
                Some(values) => {
                    group_counts[group_index] += 1;
                    rows.push(MatrixRow {
                        record: record.clone(),
                        values,
                    });
                }
                None => {
                    // Region skipped due to thresholds or zero masking.
                }
            }
        }
    }

    let group_labels: Vec<String> = groups.iter().map(|group| group.label.clone()).collect();
    let group_boundaries = compute_group_boundaries(&group_counts);
    let sample_boundaries = compute_sample_boundaries(sample_count, total_bins);

    let proc_number = general.number_of_processors.resolve();

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

    Ok(MatrixData {
        header,
        rows,
        bin_count: total_bins,
        sample_count,
    })
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

fn compute_group_boundaries(counts: &[usize]) -> Vec<usize> {
    let mut boundaries = Vec::with_capacity(counts.len() + 1);
    let mut running = 0usize;
    boundaries.push(0);
    for count in counts {
        running += *count;
        boundaries.push(running);
    }
    boundaries
}

fn compute_sample_boundaries(sample_count: usize, bins_per_sample: usize) -> Vec<usize> {
    let mut boundaries = Vec::with_capacity(sample_count + 1);
    for index in 0..=sample_count {
        boundaries.push(index * bins_per_sample);
    }
    boundaries
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
    path.file_stem()
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
    options: &ReferencePointOptions,
    general: &GeneralOptions,
    bin_size: u32,
    upstream_bins: usize,
    downstream_bins: usize,
) -> Result<Option<Vec<Vec<f32>>>> {
    let mut per_sample = Vec::with_capacity(samples.len());
    for sample in samples.iter_mut() {
        let values = compute_sample_bins(
            sample,
            record,
            options,
            general,
            bin_size,
            upstream_bins,
            downstream_bins,
        )?;
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
    options: &ReferencePointOptions,
    general: &GeneralOptions,
    bin_size: u32,
    upstream_bins: usize,
    downstream_bins: usize,
) -> Result<Vec<f32>> {
    let bin_count = upstream_bins + downstream_bins;
    let reference = reference_coordinate(record, &options.reference_point);
    let (window_start, window_end) = window_bounds(reference, record.strand, options);

    let chrom_length = match sample.chrom_length(&record.chrom) {
        Some(length) => length,
        None => {
            return Ok(vec![f32::NAN; bin_count]);
        }
    };

    let fetch_start = clamp_coordinate(window_start, chrom_length);
    let fetch_end = clamp_coordinate(window_end, chrom_length);

    let intervals = if fetch_start < fetch_end {
        sample
            .reader
            .values(&record.chrom, fetch_start, fetch_end)
            .map_err(anyhow::Error::new)
            .with_context(|| {
                format!(
                    "Failed to read bigWig intervals for '{}' in '{}'",
                    record.chrom,
                    sample.path.display()
                )
            })?
    } else {
        Vec::new()
    };

    let mut bins = Vec::with_capacity(bin_count);
    for bin_index in 0..bin_count {
        let (bin_start, bin_end) =
            bin_boundaries(bin_index, upstream_bins, bin_size, reference, record.strand);

        let clamped_start = clamp_coordinate(bin_start, chrom_length);
        let clamped_end = clamp_coordinate(bin_end, chrom_length);

        let (sum, covered) = if clamped_start < clamped_end {
            summarise_bin(&intervals, clamped_start, clamped_end)
        } else {
            (0.0, 0)
        };

        let mut value = if covered == 0 {
            if general.missing_data_as_zero {
                0.0
            } else {
                f32::NAN
            }
        } else {
            (sum / covered as f64) as f32
        };

        if options.nan_after_end && bin_beyond_region(record, bin_start, bin_end) {
            value = f32::NAN;
        }

        if value.is_finite() {
            value *= general.scale_factor as f32;
        }

        bins.push(value);
    }

    Ok(bins)
}

fn reference_coordinate(record: &BedRecord, reference_point: &ReferencePoint) -> i64 {
    let start = record.start as i64;
    let end = record.end as i64;

    match reference_point {
        ReferencePoint::Tss => match record.strand {
            Strand::Negative => end,
            _ => start,
        },
        ReferencePoint::Tes => match record.strand {
            Strand::Negative => start,
            _ => end,
        },
        ReferencePoint::Center => (start + end) / 2,
    }
}

fn window_bounds(reference: i64, strand: Strand, options: &ReferencePointOptions) -> (i64, i64) {
    let upstream = options.upstream as i64;
    let downstream = options.downstream as i64;
    match strand {
        Strand::Negative => (reference - downstream, reference + upstream),
        _ => (reference - upstream, reference + downstream),
    }
}

fn bin_boundaries(
    bin_index: usize,
    upstream_bins: usize,
    bin_size: u32,
    reference: i64,
    strand: Strand,
) -> (i64, i64) {
    let bin_size = bin_size as i64;
    let offset_start = (bin_index as i64 - upstream_bins as i64) * bin_size;
    let offset_end = offset_start + bin_size;

    match strand {
        Strand::Negative => {
            let start = reference - offset_end;
            let end = reference - offset_start;
            (start.min(end), start.max(end))
        }
        _ => {
            let start = reference + offset_start;
            let end = reference + offset_end;
            (start.min(end), start.max(end))
        }
    }
}

fn bin_beyond_region(record: &BedRecord, bin_start: i64, bin_end: i64) -> bool {
    let region_start = record.start as i64;
    let region_end = record.end as i64;

    match record.strand {
        Strand::Negative => bin_end <= region_start,
        _ => bin_start >= region_end,
    }
}

fn summarise_bin(intervals: &[BigWigValue], start: u32, end: u32) -> (f64, u32) {
    let mut sum = 0.0f64;
    let mut covered = 0u32;
    for value in intervals {
        if value.end <= start {
            continue;
        }
        if value.start >= end {
            break;
        }
        let overlap_start = value.start.max(start);
        let overlap_end = value.end.min(end);
        if overlap_start >= overlap_end {
            continue;
        }
        let length = overlap_end - overlap_start;
        sum += value.value as f64 * length as f64;
        covered += length;
    }
    (sum, covered)
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

    fn build_record(strand: Strand, start: u32, end: u32) -> BedRecord {
        BedRecord {
            chrom: "chr1".to_string(),
            start,
            end,
            name: None,
            score: None,
            strand,
            extra_fields: Vec::new(),
        }
    }

    #[test]
    fn reference_coordinate_respects_strand_for_tss() {
        let positive = build_record(Strand::Positive, 100, 200);
        let negative = build_record(Strand::Negative, 100, 200);

        assert_eq!(reference_coordinate(&positive, &ReferencePoint::Tss), 100);
        assert_eq!(reference_coordinate(&negative, &ReferencePoint::Tss), 200);
    }

    #[test]
    fn reference_coordinate_respects_strand_for_tes() {
        let positive = build_record(Strand::Positive, 100, 200);
        let negative = build_record(Strand::Negative, 100, 200);

        assert_eq!(reference_coordinate(&positive, &ReferencePoint::Tes), 200);
        assert_eq!(reference_coordinate(&negative, &ReferencePoint::Tes), 100);
    }

    #[test]
    fn reference_coordinate_center_midpoint() {
        let record = build_record(Strand::Positive, 5, 15);
        assert_eq!(reference_coordinate(&record, &ReferencePoint::Center), 10);
    }

    #[test]
    fn bin_boundaries_handle_positive_strand() {
        let reference = 100;
        let (start_upstream, end_upstream) = bin_boundaries(0, 2, 10, reference, Strand::Positive);
        assert_eq!((start_upstream, end_upstream), (80, 90));

        let (start_downstream, end_downstream) =
            bin_boundaries(2, 2, 10, reference, Strand::Positive);
        assert_eq!((start_downstream, end_downstream), (100, 110));
    }

    #[test]
    fn bin_boundaries_handle_negative_strand() {
        let reference = 200;
        let (start_upstream, end_upstream) = bin_boundaries(0, 2, 10, reference, Strand::Negative);
        assert_eq!((start_upstream, end_upstream), (210, 220));

        let (start_downstream, end_downstream) =
            bin_boundaries(2, 2, 10, reference, Strand::Negative);
        assert_eq!((start_downstream, end_downstream), (190, 200));
    }

    #[test]
    fn bin_beyond_region_detects_downstream_tail() {
        let record = build_record(Strand::Positive, 100, 200);
        assert!(bin_beyond_region(&record, 205, 215));
        assert!(!bin_beyond_region(&record, 195, 205));

        let record_neg = build_record(Strand::Negative, 100, 200);
        assert!(bin_beyond_region(&record_neg, 80, 90));
        assert!(!bin_beyond_region(&record_neg, 95, 105));
    }
}
