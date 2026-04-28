use std::collections::{HashMap, HashSet};
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::{Arc, mpsc};

use anyhow::{Context, Result, anyhow, bail};
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

use crate::config::{AverageTypeBins, GeneralOptions, GtfOptions};
use crate::io::writers::StreamingMatrixWriter;
use crate::io::{BedReadError, BedRecord, BigWigReader, SharedBigWigReader, load_gtf_records};
use crate::pipeline::matrix::{MatrixData, MatrixHeader, MatrixRow};

pub trait SignalBin {
    fn start(&self) -> i64;
    fn end(&self) -> i64;
    fn beyond_region(&self) -> bool;
}

#[derive(Debug, Clone, Copy)]
pub enum ModeTag {
    ReferencePoint,
    ScaleRegions,
}

impl fmt::Display for ModeTag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let label = match self {
            ModeTag::ReferencePoint => "reference-point",
            ModeTag::ScaleRegions => "scale-regions",
        };
        f.write_str(label)
    }
}

pub fn ensure_positive(value: u32, flag: &str, mode: ModeTag) -> Result<()> {
    if value == 0 {
        bail!("[{mode}] {flag} must be a positive integer");
    }
    Ok(())
}

pub fn ensure_multiple(bin_size: u32, distance: u32, flag: &str, mode: ModeTag) -> Result<()> {
    if distance % bin_size != 0 {
        bail!("[{mode}] {flag} ({distance}) must be a multiple of the bin size ({bin_size})");
    }
    Ok(())
}

pub trait RegionPlan {
    type Bin: SignalBin;

    fn window_start(&self) -> i64;
    fn window_end(&self) -> i64;
    fn bins(&self) -> &[Self::Bin];

    /// Optional list of intervals that should contribute signal to this plan.
    /// When present (metagene mode), coverage outside these intervals is
    /// treated as missing data to avoid counting intronic signal.
    fn included_intervals(&self) -> Option<&[(i64, i64)]> {
        None
    }
}

pub struct Sample {
    path: PathBuf,
    reader: BigWigReader,
    chrom_lengths: HashMap<String, u32>,
}

impl Sample {
    pub fn open(path: &Path) -> Result<Self> {
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

    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn reader(&self) -> &BigWigReader {
        &self.reader
    }

    pub fn reader_mut(&mut self) -> &mut BigWigReader {
        &mut self.reader
    }

    /// Create a Sample from an already-opened shared reader.  The underlying
    /// mmap and metadata are shared via Arc; only the per-worker caches are
    /// fresh.
    pub fn from_shared(path: PathBuf, shared: Arc<SharedBigWigReader>) -> Self {
        let chrom_lengths = shared
            .chroms()
            .iter()
            .map(|chrom| (chrom.name.clone(), chrom.length))
            .collect::<HashMap<_, _>>();
        Self {
            path,
            reader: BigWigReader::from_shared(shared),
            chrom_lengths,
        }
    }

    pub fn chrom_length(&self, chrom: &str) -> Option<u32> {
        self.chrom_lengths.get(chrom).copied()
    }
}

pub struct WorkerSamples {
    samples: Result<Vec<Sample>, String>,
}

impl WorkerSamples {
    pub fn new(paths: Arc<Vec<PathBuf>>) -> Self {
        let samples = open_samples(paths.as_ref()).map_err(|err| err.to_string());
        Self { samples }
    }

    /// Create per-worker Sample instances from pre-opened shared readers.
    /// Each worker gets its own caches but shares the mmap-backed immutable
    /// state, avoiding redundant mmap entries per thread.
    pub fn from_shared(
        paths: Arc<Vec<PathBuf>>,
        shared_readers: Arc<Vec<Arc<SharedBigWigReader>>>,
    ) -> Self {
        let samples = paths
            .iter()
            .zip(shared_readers.iter())
            .map(|(path, shared)| Sample::from_shared(path.clone(), Arc::clone(shared)))
            .collect();
        Self { samples: Ok(samples) }
    }

    pub fn samples(&mut self) -> Result<&mut Vec<Sample>> {
        match &mut self.samples {
            Ok(samples) => Ok(samples),
            Err(message) => Err(anyhow!(message.clone())),
        }
    }
}

pub struct Group {
    pub label: String,
    pub records: Vec<BedRecord>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RegionFormat {
    Bed,
    Gtf,
}

#[derive(Clone)]
pub struct RegionTask {
    pub index: usize,
    pub group_index: usize,
    pub record: Arc<BedRecord>,
}

pub fn load_groups(paths: &[PathBuf], gtf: &GtfOptions) -> Result<Vec<Group>> {
    let mut groups = Vec::new();
    let mut seen_labels = HashSet::new();
    // When there's only one file, Python uses "genes" as the default label
    let use_default_genes_label = paths.len() == 1;
    for path in paths {
        match infer_region_format(path) {
            RegionFormat::Bed => {
                let mut file_groups =
                    parse_grouped_bed(path, use_default_genes_label, &mut seen_labels)
                        .map_err(anyhow::Error::new)
                        .with_context(|| {
                            format!("Failed to parse regions file '{}'", path.display())
                        })?;
                groups.append(&mut file_groups);
            }
            RegionFormat::Gtf => {
                let mut file_groups =
                    parse_grouped_gtf(path, gtf, use_default_genes_label, &mut seen_labels)
                        .with_context(|| {
                            format!("Failed to parse regions file '{}'", path.display())
                        })?;
                groups.append(&mut file_groups);
            }
        }
    }

    Ok(groups)
}

pub fn derive_sample_labels(paths: &[PathBuf], general: &GeneralOptions) -> Result<Vec<String>> {
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

    Ok(paths
        .iter()
        .map(|path| label_from_path(path, general.smart_labels))
        .collect())
}

pub fn compute_row<P: RegionPlan>(
    samples: &mut [Sample],
    record: &BedRecord,
    plan: &P,
    general: &GeneralOptions,
    nan_after_end: bool,
) -> Result<Option<(Vec<f32>, usize, usize)>> {
    let sample_count = samples.len();
    let bin_count = plan.bins().len();
    let mut all_values = Vec::with_capacity(sample_count * bin_count);
    for sample in samples.iter_mut() {
        let values = compute_sample_bins(sample, record, plan, general, nan_after_end)?;
        all_values.extend(values);
    }

    if should_skip_row_flat(&all_values, general) {
        return Ok(None);
    }

    Ok(Some((all_values, sample_count, bin_count)))
}

fn should_skip_row_flat(values: &[f32], general: &GeneralOptions) -> bool {
    if general.skip_zeros {
        let mut all_zero = true;
        for &value in values {
            if value.is_nan() {
                continue;
            }
            if value != 0.0 {
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
            .filter(|value| !value.is_nan())
            .any(|value| (*value as f64) <= min_threshold)
        {
            return true;
        }
    }

    if let Some(max_threshold) = general.max_threshold {
        if values
            .iter()
            .filter(|value| !value.is_nan())
            .any(|value| (*value as f64) >= max_threshold)
        {
            return true;
        }
    }

    false
}

pub fn normalize_sort_sample_indices(
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

fn compute_sample_bins<P: RegionPlan>(
    sample: &mut Sample,
    record: &BedRecord,
    plan: &P,
    general: &GeneralOptions,
    nan_after_end: bool,
) -> Result<Vec<f32>> {
    let bin_count = plan.bins().len();
    let chrom_length = match sample.chrom_length(&record.chrom) {
        Some(length) => length,
        None => {
            return Ok(vec![f32::NAN; bin_count]);
        }
    };

    let window_span = plan.window_end() - plan.window_start();
    if window_span <= 0 {
        return Ok(vec![f32::NAN; bin_count]);
    }

    let window_len = usize::try_from(window_span).expect("region plan window span exceeds usize");
    let default_fill = if general.missing_data_as_zero {
        0.0f32
    } else {
        f32::NAN
    };

    thread_local! {
        static COVERAGE_BUF: std::cell::RefCell<Vec<f32>> = std::cell::RefCell::new(Vec::new());
    }

    let bins = COVERAGE_BUF.with(|cell| {
        let mut buf = cell.borrow_mut();
        buf.clear();
        buf.resize(window_len, default_fill);

        if let Some(allowed) = plan.included_intervals() {
            let base_offset = plan.window_start();
            for (seg_start, seg_end) in allowed {
                let fetch_start = clamp_coordinate(*seg_start, chrom_length);
                let fetch_end = clamp_coordinate(*seg_end, chrom_length);
                if fetch_start >= fetch_end {
                    continue;
                }

                let intervals = sample
                    .reader_mut()
                    .values(&record.chrom, fetch_start, fetch_end)
                    .map_err(anyhow::Error::new)
                    .with_context(|| {
                        format!(
                            "Failed to read bigWig intervals for '{}' in '{}'",
                            record.chrom,
                            sample.path().display()
                        )
                    })?;

                for interval in intervals {
                    let overlap_start = i64::from(interval.start).max(i64::from(fetch_start));
                    let overlap_end = i64::from(interval.end).min(i64::from(fetch_end));
                    if overlap_start >= overlap_end {
                        continue;
                    }
                    let rel_start = usize::try_from(overlap_start - base_offset)
                        .expect("relative start offset exceeded usize");
                    let rel_end = usize::try_from(overlap_end - base_offset)
                        .expect("relative end offset exceeded usize");
                    buf[rel_start..rel_end].fill(interval.value);
                }
            }
        } else {
            let fetch_start = clamp_coordinate(plan.window_start(), chrom_length);
            let fetch_end = clamp_coordinate(plan.window_end(), chrom_length);

            if fetch_start < fetch_end {
                let intervals = sample
                    .reader_mut()
                    .values(&record.chrom, fetch_start, fetch_end)
                    .map_err(anyhow::Error::new)
                    .with_context(|| {
                        format!(
                            "Failed to read bigWig intervals for '{}' in '{}'",
                            record.chrom,
                            sample.path().display()
                        )
                    })?;

                let base_offset = plan.window_start();
                for interval in intervals {
                    let overlap_start = i64::from(interval.start).max(i64::from(fetch_start));
                    let overlap_end = i64::from(interval.end).min(i64::from(fetch_end));
                    if overlap_start >= overlap_end {
                        continue;
                    }
                    let rel_start = usize::try_from(overlap_start - base_offset)
                        .expect("relative start offset exceeded usize");
                    let rel_end = usize::try_from(overlap_end - base_offset)
                        .expect("relative end offset exceeded usize");
                    buf[rel_start..rel_end].fill(interval.value);
                }
            }
        }

        let mut bins = Vec::with_capacity(bin_count);
        for bin in plan.bins() {
            let start_idx = index_from_coordinate(bin.start(), plan.window_start(), window_len);
            let end_idx = index_from_coordinate(bin.end(), plan.window_start(), window_len);

            let mut value = if start_idx < end_idx {
                aggregate_slice(&buf[start_idx..end_idx], general.average_type_bins)
            } else {
                None
            };

            if value.is_none() && general.missing_data_as_zero {
                value = Some(0.0);
            }

            let mut value = value.unwrap_or(f32::NAN);

            if nan_after_end && bin.beyond_region() {
                value = f32::NAN;
            }

            if value.is_finite() {
                value *= general.scale_factor as f32;
            }

            bins.push(value);
        }

        Ok::<Vec<f32>, anyhow::Error>(bins)
    })?;

    Ok(bins)
}

fn aggregate_slice(slice: &[f32], average_type: AverageTypeBins) -> Option<f32> {
    let len = slice.len();
    if len == 0 {
        return None;
    }

    match average_type {
        AverageTypeBins::Mean => {
            let mut sum = 0.0f32;
            let mut count = 0u32;
            for &value in slice {
                if !value.is_nan() {
                    sum += value;
                    count += 1;
                }
            }
            if count == 0 { None } else { Some(sum / count as f32) }
        }
        AverageTypeBins::Sum => {
            let mut sum = 0.0f32;
            let mut found = false;
            for &value in slice {
                if !value.is_nan() {
                    sum += value;
                    found = true;
                }
            }
            if found { Some(sum) } else { None }
        }
        AverageTypeBins::Min => {
            let mut min = f32::INFINITY;
            let mut found = false;
            for &value in slice {
                if !value.is_nan() {
                    min = min.min(value);
                    found = true;
                }
            }
            if found { Some(min) } else { None }
        }
        AverageTypeBins::Max => {
            let mut max = f32::NEG_INFINITY;
            let mut found = false;
            for &value in slice {
                if !value.is_nan() {
                    max = max.max(value);
                    found = true;
                }
            }
            if found { Some(max) } else { None }
        }
        AverageTypeBins::Std => {
            let mut sum = 0.0f32;
            let mut count = 0u32;
            for &value in slice {
                if !value.is_nan() {
                    sum += value;
                    count += 1;
                }
            }
            if count == 0 {
                return None;
            }
            let mean = sum / count as f32;
            let mut variance_sum = 0.0f64;
            for &value in slice {
                if !value.is_nan() {
                    let delta = value as f64 - mean as f64;
                    variance_sum += delta * delta;
                }
            }
            Some((variance_sum / count as f64).sqrt() as f32)
        }
        AverageTypeBins::Median => {
            let mut values: Vec<f32> = slice
                .iter()
                .copied()
                .filter(|v| !v.is_nan())
                .collect();
            if values.is_empty() {
                return None;
            }
            values.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let mid = values.len() / 2;
            if values.len() % 2 == 0 {
                Some((values[mid - 1] + values[mid]) / 2.0)
            } else {
                Some(values[mid])
            }
        }
    }
}

fn index_from_coordinate(value: i64, base: i64, window_len: usize) -> usize {
    if value <= base {
        return 0;
    }
    let diff = value - base;
    let idx = usize::try_from(diff).unwrap_or(window_len);
    idx.min(window_len)
}

fn clamp_coordinate(value: i64, chrom_length: u32) -> u32 {
    value
        .max(0)
        .min(chrom_length as i64)
        .try_into()
        .unwrap_or(0)
}

fn open_samples(paths: &[PathBuf]) -> Result<Vec<Sample>> {
    let mut samples = Vec::with_capacity(paths.len());
    for path in paths {
        samples.push(Sample::open(path)?);
    }
    Ok(samples)
}

fn parse_grouped_bed(
    path: &Path,
    use_default_genes_label: bool,
    seen_labels: &mut HashSet<String>,
) -> Result<Vec<Group>, BedReadError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let default_label = if use_default_genes_label {
        "genes".to_string()
    } else {
        bed_file_label(path)
    };
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

fn parse_grouped_gtf(
    path: &Path,
    options: &GtfOptions,
    use_default_genes_label: bool,
    seen_labels: &mut HashSet<String>,
) -> Result<Vec<Group>> {
    let default_label = if use_default_genes_label {
        "genes".to_string()
    } else {
        bed_file_label(path)
    };
    let mut groups = Vec::new();

    let records = load_gtf_records(path, options)?;
    let label = next_unique_label("", &default_label, seen_labels);
    groups.push(Group { label, records });

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

fn label_from_path(path: &Path, use_stem: bool) -> String {
    if use_stem {
        if let Some(stem) = path.file_stem().and_then(|s| s.to_str()) {
            return stem.to_string();
        }
    }

    path.file_name()
        .and_then(|name| name.to_str())
        .map(|s| s.to_string())
        .unwrap_or_else(|| path.display().to_string())
}

fn bed_file_label(path: &Path) -> String {
    label_from_path(path, false)
}

fn infer_region_format(path: &Path) -> RegionFormat {
    let lower = path.to_string_lossy().to_ascii_lowercase();
    if lower.ends_with(".gtf")
        || lower.ends_with(".gtf.gz")
        || lower.ends_with(".gff")
        || lower.ends_with(".gff.gz")
        || lower.ends_with(".gff3")
        || lower.ends_with(".gff3.gz")
    {
        RegionFormat::Gtf
    } else {
        RegionFormat::Bed
    }
}

pub trait PipelineMode: Sync {
    type Plan: RegionPlan;
    type Metadata: Send + Sync;

    fn validate(&self, general: &GeneralOptions) -> Result<Self::Metadata>;
    fn total_bins(&self, metadata: &Self::Metadata) -> usize;
    fn plan_for(&self, record: &BedRecord, metadata: &Self::Metadata) -> Self::Plan;
    fn nan_after_end(&self, metadata: &Self::Metadata) -> bool;
    fn postprocess_row(
        &self,
        record: BedRecord,
        values: Vec<f32>,
        sample_count: usize,
        bin_count: usize,
        metadata: &Self::Metadata,
    ) -> MatrixRow;
    fn build_header(
        &self,
        general: &GeneralOptions,
        metadata: &Self::Metadata,
        sample_labels: &[String],
        group_labels: &[String],
        group_counts: &[usize],
        thread_count: usize,
        sample_count: usize,
    ) -> MatrixHeader;
}

pub trait RowCollector: Send {
    type Output: Send;

    fn on_row(&mut self, row: MatrixRow) -> Result<()>;
    fn finalize(self, header: MatrixHeader) -> Result<Self::Output>;
}

pub struct InMemoryCollector {
    rows: Vec<MatrixRow>,
    sample_count: usize,
    bin_count: usize,
}

impl InMemoryCollector {
    pub fn with_capacity(capacity: usize, sample_count: usize, bin_count: usize) -> Self {
        Self {
            rows: Vec::with_capacity(capacity),
            sample_count,
            bin_count,
        }
    }
}

impl RowCollector for InMemoryCollector {
    type Output = MatrixData;

    fn on_row(&mut self, row: MatrixRow) -> Result<()> {
        self.rows.push(row);
        Ok(())
    }

    fn finalize(self, header: MatrixHeader) -> Result<Self::Output> {
        Ok(MatrixData {
            header,
            rows: self.rows,
            bin_count: self.bin_count,
            sample_count: self.sample_count,
        })
    }
}

pub struct FileCollector {
    writer: StreamingMatrixWriter,
}

impl FileCollector {
    pub fn new(writer: StreamingMatrixWriter) -> Self {
        Self { writer }
    }
}

impl RowCollector for FileCollector {
    type Output = ();

    fn on_row(&mut self, row: MatrixRow) -> Result<()> {
        self.writer.write_row(&row)
    }

    fn finalize(self, header: MatrixHeader) -> Result<Self::Output> {
        self.writer.finish(&header)
    }
}

type BatchResult = (usize, usize, Option<MatrixRow>);

/// Internal work item carrying the original index and I/O sort key.
struct WorkItem {
    orig_idx: usize,
    group_index: usize,
    record: Arc<BedRecord>,
    query_start: i64,
    query_end: i64,
}

// ── Query coalescing ──────────────────────────────────────────────────────

const COALESCE_CLAMP_MAX: i64 = 2000;

enum CoalesceStrategy {
    Coalesce(i64),
    NoCoalesce,
}

/// Estimate a coalescing gap threshold from the actual distribution of gaps
/// between consecutive same-chromosome items in the sorted work list.
///
/// Returns a threshold in `[100, 2000]` based on the 75th percentile of
/// observed gaps.  Falls back to 500 bp when there are too few gaps (< 10).
fn estimate_coalesce_gap(work_items: &[WorkItem]) -> i64 {
    if work_items.len() < 2 {
        return 500;
    }

    let mut gaps: Vec<i64> = Vec::new();
    for w in work_items.windows(2) {
        if w[0].record.chrom == w[1].record.chrom {
            let gap = w[1].query_start - w[0].query_end;
            if gap > 0 {
                gaps.push(gap);
            }
        }
    }

    if gaps.len() < 10 {
        eprintln!(
            "[coalesce-gap] {} gaps (< 10), using default threshold: 500",
            gaps.len()
        );
        return 500;
    }

    gaps.sort_unstable();

    let n = gaps.len();
    let p50 = gaps[n / 2];
    let p75 = gaps[(n * 3) / 4];
    let threshold = p75.clamp(100, 2000);

    eprintln!(
        "[coalesce-gap] n_gaps={} p50={} p75={} threshold={}",
        n, p50, p75, threshold
    );
    threshold
}

/// A batch of consecutive work items on the same chromosome whose query
/// windows overlap or are separated by at most the caller-supplied
/// `coalesce_gap` threshold.  Records are **moved** (not cloned) from
/// `WorkItem`s, and `work_items` is consumed.
struct CoalescedBatch {
    /// Items in original sorted order: (orig_idx, group_index, record).
    items: Vec<(usize, usize, Arc<BedRecord>)>,
    /// Start of the merged query window (minimum of all item windows).
    query_start: i64,
    /// End of the merged query window (maximum of all item windows).
    query_end: i64,
}

/// Create batches according to the chosen strategy.
fn create_batches(work_items: Vec<WorkItem>, strategy: &CoalesceStrategy) -> Vec<CoalescedBatch> {
    match strategy {
        CoalesceStrategy::Coalesce(coalesce_gap) => {
            create_coalesced_batches(work_items, *coalesce_gap)
        }
        CoalesceStrategy::NoCoalesce => create_per_item_batches(work_items),
    }
}

fn create_per_item_batches(work_items: Vec<WorkItem>) -> Vec<CoalescedBatch> {
    work_items
        .into_iter()
        .map(|item| CoalescedBatch {
            query_start: item.query_start,
            query_end: item.query_end,
            items: vec![(item.orig_idx, item.group_index, item.record)],
        })
        .collect()
}

/// Scan the sorted `work_items`, group consecutive same-chromosome items
/// whose query windows overlap or are gapped by at most `coalesce_gap`,
/// and move records into [`CoalescedBatch`]es.  `work_items` is consumed.
fn create_coalesced_batches(work_items: Vec<WorkItem>, coalesce_gap: i64) -> Vec<CoalescedBatch> {
    let mut batches = Vec::new();
    let mut current_chrom: Arc<str> = Arc::from("");
    let mut current_items: Vec<(usize, usize, Arc<BedRecord>)> = Vec::new();
    let mut batch_start: i64 = 0;
    let mut batch_end: i64 = 0;

    for item in work_items {
        if current_items.is_empty() {
            current_chrom = item.record.chrom.clone();
            batch_start = item.query_start;
            batch_end = item.query_end;
            current_items.push((item.orig_idx, item.group_index, item.record));
        } else if item.record.chrom != current_chrom
            || item.query_start > batch_end.saturating_add(coalesce_gap)
        {
            batches.push(CoalescedBatch {
                items: std::mem::take(&mut current_items),
                query_start: batch_start,
                query_end: batch_end,
            });
            current_chrom = item.record.chrom.clone();
            batch_start = item.query_start;
            batch_end = item.query_end;
            current_items.push((item.orig_idx, item.group_index, item.record));
        } else {
            batch_end = batch_end.max(item.query_end);
            current_items.push((item.orig_idx, item.group_index, item.record));
        }
    }

    if !current_items.is_empty() {
        batches.push(CoalescedBatch {
            items: current_items,
            query_start: batch_start,
            query_end: batch_end,
        });
    }

    batches
}

thread_local! {
    static COVERAGE_POOL: std::cell::RefCell<Vec<Vec<f32>>> =
        std::cell::RefCell::new(Vec::new());
}

fn take_coverage_buffers(sample_count: usize, window_len: usize, default_fill: f32) -> Vec<Vec<f32>> {
    COVERAGE_POOL.with(|pool| {
        let mut bufs = pool.borrow_mut();
        bufs.resize_with(sample_count, Vec::new);
        for buf in bufs.iter_mut() {
            buf.clear();
            buf.resize(window_len, default_fill);
        }
        std::mem::take(&mut *bufs)
    })
}

fn return_coverage_buffers(bufs: Vec<Vec<f32>>) {
    COVERAGE_POOL.with(|pool| {
        *pool.borrow_mut() = bufs;
    });
}

/// Process a single coalesced batch.
///
/// Performs one bigWig read per sample for the batch's merged query window,
/// then extracts per-region bins from the pre-read coverage buffers.
/// Items in metagene mode (where `included_intervals()` returns `Some`) fall
/// back to the original per-item `compute_row` path for correctness.
///
/// Records are **moved** (not cloned) out of the batch items, so once a
/// batch is processed its records are transferred into the result rows.
fn process_batch<M: PipelineMode>(
    samples: &mut [Sample],
    batch: CoalescedBatch,
    mode: &M,
    general: &GeneralOptions,
    metadata: &M::Metadata,
) -> Result<Vec<(usize, usize, Option<MatrixRow>)>> {
    if batch.items.is_empty() {
        return Ok(Vec::new());
    }

    let item_count = batch.items.len();
    let window_span = batch.query_end - batch.query_start;
    let sample_count = samples.len();

    // Zero or negative window span — delegate to per-item path.
    if window_span <= 0 {
        let nan_after_end = mode.nan_after_end(metadata);
        let mut results = Vec::with_capacity(item_count);
        for (orig_idx, group_index, record) in batch.items {
            let plan = mode.plan_for(&record, metadata);
            let maybe_values =
                compute_row(samples, &record, &plan, general, nan_after_end)?;
            let row = maybe_values
                .map(|(flat, sc, bc)| mode.postprocess_row(Arc::unwrap_or_clone(record), flat, sc, bc, metadata));
            results.push((orig_idx, group_index, row));
        }
        return Ok(results);
    }

    let window_len =
        usize::try_from(window_span).context("batch window span exceeds usize")?;

    let default_fill = if general.missing_data_as_zero {
        0.0f32
    } else {
        f32::NAN
    };

    let chrom = &batch.items[0].2.chrom;

    // ── ONE bigWig read per sample for the entire merged window ────────
    let mut sample_coverages = take_coverage_buffers(sample_count, window_len, default_fill);
    for (si, sample) in samples.iter_mut().enumerate() {
        let chrom_length = match sample.chrom_length(chrom) {
            Some(l) => l,
            None => {
                continue;
            }
        };

        let fetch_start = clamp_coordinate(batch.query_start, chrom_length);
        let fetch_end = clamp_coordinate(batch.query_end, chrom_length);

        if fetch_start >= fetch_end {
            continue;
        }

        let intervals = sample
            .reader_mut()
            .values(chrom, fetch_start, fetch_end)
            .map_err(anyhow::Error::new)
            .with_context(|| {
                format!(
                    "Failed to read bigWig intervals for '{}' in '{}'",
                    chrom,
                    sample.path().display()
                )
            })?;

        let cov = &mut sample_coverages[si];
        for v in intervals {
            let rs = i64::from(v.start)
                .saturating_sub(batch.query_start)
                .max(0);
            let re = i64::from(v.end)
                .saturating_sub(batch.query_start)
                .min(window_span)
                .max(0);
            if rs < re {
                cov[rs as usize..re as usize].fill(v.value);
            }
        }
    }

    // ── Extract per-region bins from the pre-read coverage buffers ─────
    let nan_after_end = mode.nan_after_end(metadata);
    let mut results = Vec::with_capacity(item_count);

    for (orig_idx, group_index, record) in batch.items {
        let plan = mode.plan_for(&record, metadata);

        // Metagene fallback: items with explicit included_intervals
        // (intron-skipping) must read individual exon intervals; we
        // delegate to the original per-item compute_row path.
        if plan.included_intervals().is_some() {
            let maybe_values =
                compute_row(samples, &record, &plan, general, nan_after_end)?;
            let row = maybe_values
                .map(|(flat, sc, bc)| mode.postprocess_row(Arc::unwrap_or_clone(record), flat, sc, bc, metadata));
            results.push((orig_idx, group_index, row));
            continue;
        }

        let bins = plan.bins();
        let bin_count = bins.len();
        let mut all_values = Vec::with_capacity(sample_count * bin_count);

        for si in 0..sample_count {
            let cov = &sample_coverages[si];
            for bin in bins {
                let bs =
                    ((bin.start() - batch.query_start).max(0) as usize).min(window_len);
                let be =
                    ((bin.end() - batch.query_start).max(0) as usize).min(window_len);

                let mut value = if bs < be {
                    aggregate_slice(&cov[bs..be], general.average_type_bins)
                } else {
                    None
                };

                if value.is_none() && general.missing_data_as_zero {
                    value = Some(0.0);
                }

                let mut value = value.unwrap_or(f32::NAN);

                if nan_after_end && bin.beyond_region() {
                    value = f32::NAN;
                }

                if value.is_finite() {
                    value *= general.scale_factor as f32;
                }

                all_values.push(value);
            }
        }

        let row = if should_skip_row_flat(&all_values, general) {
            None
        } else {
            Some(mode.postprocess_row(
                Arc::unwrap_or_clone(record),
                all_values,
                sample_count,
                bin_count,
                metadata,
            ))
        };
        results.push((orig_idx, group_index, row));
    }

    return_coverage_buffers(sample_coverages);
    Ok(results)
}

pub fn execute_mode<M, C, F>(
    tasks: Vec<RegionTask>,
    general: &GeneralOptions,
    sample_paths: Arc<Vec<PathBuf>>,
    collector: C,
    thread_count: usize,
    mode: &M,
    metadata: Arc<M::Metadata>,
    header_builder: F,
    group_count: usize,
) -> Result<C::Output>
where
    M: PipelineMode,
    C: RowCollector + Send + 'static,
    F: FnOnce(Vec<usize>) -> Result<MatrixHeader> + Send + 'static,
{
    let task_count = tasks.len();

    // ── Phase 1: Pre-compute sort keys for I/O locality ──────────────────
    let mut work_items: Vec<WorkItem> = tasks
        .into_iter()
        .map(|task| {
            let plan = mode.plan_for(&task.record, metadata.as_ref());
            WorkItem {
                orig_idx: task.index,
                group_index: task.group_index,
                record: task.record,
                query_start: plan.window_start(),
                query_end: plan.window_end(),
            }
        })
        .collect();

    // ── Phase 2: Sort by (chrom, window_start, window_end) ────────────────
    work_items.sort_by(|a, b| {
        a.record
            .chrom
            .cmp(&b.record.chrom)
            .then(a.query_start.cmp(&b.query_start))
            .then(a.query_end.cmp(&b.query_end))
    });

    // Empty input: build an empty header and return.
    if task_count == 0 {
        let header = header_builder(vec![0; group_count])?;
        return collector.finalize(header);
    }

    // ── Phase 3: Open shared bigWig readers once ────────────────────────
    // Open one set of readers per file and share them across rayon workers
    // via Arc.  The shared readers use pread-based I/O (not mmap) so RSS
    // stays low — file pages live in the kernel page cache, not in the
    // process address space.  We still drop the Arc before the thread-pool
    // exits for clean resource management.
    let shared_readers = Arc::new(
        sample_paths
            .iter()
            .map(|path| {
                SharedBigWigReader::open(path)
                    .map(Arc::new)
                    .with_context(|| {
                        format!("Failed to open bigWig file '{}'", path.display())
                    })
            })
            .collect::<Result<Vec<_>>>()?,
    );

    // ── Phase 3.5: Create coalesced batches ─────────────────────────────
    // Estimate a coalescing gap from the actual gap distribution, then
    // decide whether to coalesce or skip it for sparse datasets.  When
    // the estimated gap exceeds COALESCE_CLAMP_MAX the data is sparse
    // enough that coalescing would not merge many items, so we skip it.
    // Records are moved (not cloned) from work_items, so work_items is
    // consumed here.
    let coalesce_gap = estimate_coalesce_gap(&work_items);
    let strategy = if coalesce_gap >= COALESCE_CLAMP_MAX {
        CoalesceStrategy::NoCoalesce
    } else {
        CoalesceStrategy::Coalesce(coalesce_gap)
    };
    let batches = create_batches(work_items, &strategy);
    eprintln!(
        "[coalesce-gap] strategy={:?} batches={} items={} ratio={:.2}",
        match &strategy {
            CoalesceStrategy::Coalesce(g) => format!("coalesce({g})"),
            CoalesceStrategy::NoCoalesce => "no-coalesce".into(),
        },
        batches.len(),
        task_count,
        batches.len() as f64 / task_count as f64
    );

    // ── Phase 4: Spawn writer thread ────────────────────────────────────
    // Channel buffers up to 256 results between compute workers and the
    // writer.  The writer reorders via a BTreeMap so rows are emitted in
    // the same order as the original input (orig_idx).
    let (tx, rx) = mpsc::sync_channel::<BatchResult>(256);

    let writer_handle = std::thread::Builder::new()
        .name("matrix-writer".into())
        .spawn(move || {
            let mut next_idx: usize = 0;
            let mut pending: std::collections::BTreeMap<usize, (usize, Option<MatrixRow>)> =
                std::collections::BTreeMap::new();
            let mut collector = collector;
            let mut group_counts = vec![0usize; group_count];

            for (orig_idx, group_index, row) in rx {
                pending.insert(orig_idx, (group_index, row));
                while let Some(entry) = pending.remove(&next_idx) {
                    let (grp, row_opt) = entry;
                    if let Some(row) = row_opt {
                        collector.on_row(row)?;
                        group_counts[grp] += 1;
                    }
                    next_idx += 1;
                }
            }
            let header = header_builder(group_counts)?;
            collector.finalize(header)
        })
        .context("Failed to spawn writer thread")?;

    // ── Phase 5: Parallel processing ────────────────────────────────────
    // Compute workers process batches in parallel and stream results to
    // the writer thread via the sync channel.  This avoids buffering all
    // results in a result_slots Vec, reducing peak RSS and letting the
    // writer start I/O while compute is still in progress.
    let pool = ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build()
        .context("Failed to initialise rayon thread pool for pipeline scheduling")?;

    let sample_paths_for_workers = Arc::clone(&sample_paths);
    let shared_for_workers = Arc::clone(&shared_readers);
    let metadata_ref = metadata.as_ref();

    let batch_errors: Vec<Result<()>> = pool.install(|| {
        batches
            .into_par_iter()
            .map_init(
                move || {
                    WorkerSamples::from_shared(
                        Arc::clone(&sample_paths_for_workers),
                        Arc::clone(&shared_for_workers),
                    )
                },
                |worker_samples, batch| {
                    let samples = worker_samples.samples()?;
                    let results = process_batch(
                        samples.as_mut_slice(),
                        batch,
                        mode,
                        general,
                        metadata_ref,
                    )?;
                    for result in results {
                        if tx.send(result).is_err() {
                            break;
                        }
                    }
                    Ok(())
                },
            )
            .collect()
    });

    drop(tx);
    drop(shared_readers);

    // Propagate any batch processing errors
    for result in batch_errors {
        result?;
    }

    // ── Phase 6: Join writer thread ─────────────────────────────────────
    match writer_handle.join() {
        Ok(Ok(output)) => Ok(output),
        Ok(Err(e)) => Err(e),
        Err(panic) => std::panic::resume_unwind(panic),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Clone)]
    struct TestBin {
        start: i64,
        end: i64,
        beyond_region: bool,
    }

    impl SignalBin for TestBin {
        fn start(&self) -> i64 {
            self.start
        }

        fn end(&self) -> i64 {
            self.end
        }

        fn beyond_region(&self) -> bool {
            self.beyond_region
        }
    }

    struct TestPlan {
        start: i64,
        end: i64,
        bins: Vec<TestBin>,
    }

    impl RegionPlan for TestPlan {
        type Bin = TestBin;

        fn window_start(&self) -> i64 {
            self.start
        }

        fn window_end(&self) -> i64 {
            self.end
        }

        fn bins(&self) -> &[Self::Bin] {
            &self.bins
        }
    }

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
