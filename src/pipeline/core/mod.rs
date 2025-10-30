use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::{Arc, mpsc};
use std::thread;

use anyhow::{Context, Result, anyhow, bail};
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

use crate::config::{AverageTypeBins, GeneralOptions};
use crate::io::writers::StreamingMatrixWriter;
use crate::io::{BedReadError, BedRecord, BigWigReader};
use crate::pipeline::matrix::{MatrixData, MatrixHeader, MatrixRow};

pub trait SignalBin {
    fn start(&self) -> i64;
    fn end(&self) -> i64;
    fn beyond_region(&self) -> bool;
}

pub trait RegionPlan {
    type Bin: SignalBin;

    fn window_start(&self) -> i64;
    fn window_end(&self) -> i64;
    fn bins(&self) -> &[Self::Bin];
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

#[derive(Clone)]
pub struct RegionTask {
    pub index: usize,
    pub group_index: usize,
    pub record: BedRecord,
}

pub fn load_groups(paths: &[PathBuf]) -> Result<Vec<Group>> {
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
) -> Result<Option<Vec<Vec<f32>>>> {
    let mut per_sample = Vec::with_capacity(samples.len());
    for sample in samples.iter_mut() {
        let values = compute_sample_bins(sample, record, plan, general, nan_after_end)?;
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

    let mut coverage = vec![default_fill; window_len];
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
            coverage[rel_start..rel_end].fill(interval.value);
        }
    }

    let mut bins = Vec::with_capacity(bin_count);
    for bin in plan.bins() {
        let start_idx = index_from_coordinate(bin.start(), plan.window_start(), window_len);
        let end_idx = index_from_coordinate(bin.end(), plan.window_start(), window_len);

        let mut value = if start_idx < end_idx {
            aggregate_slice(&coverage[start_idx..end_idx], general.average_type_bins)
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

    Ok(bins)
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
    seen_labels: &mut HashSet<String>,
) -> Result<Vec<Group>, BedReadError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let default_label = bed_file_label(path);
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
        values: Vec<Vec<f32>>,
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

#[derive(Debug)]
pub struct RegionResult {
    pub index: usize,
    pub group_index: usize,
    pub row: Option<MatrixRow>,
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

pub fn spawn_row_aggregator<C, F>(
    rx: mpsc::Receiver<RegionResult>,
    collector: C,
    group_count: usize,
    task_count: usize,
    header_builder: F,
) -> Result<thread::JoinHandle<Result<C::Output>>>
where
    C: RowCollector + 'static,
    F: FnOnce(Vec<usize>) -> Result<MatrixHeader> + Send + 'static,
{
    thread::Builder::new()
        .name("matrix-aggregator".into())
        .spawn(move || consume_row_results(rx, collector, group_count, task_count, header_builder))
        .map_err(|err| anyhow!("Failed to spawn matrix streaming thread: {err}"))
}

fn consume_row_results<C, F>(
    rx: mpsc::Receiver<RegionResult>,
    mut collector: C,
    group_count: usize,
    task_count: usize,
    header_builder: F,
) -> Result<C::Output>
where
    C: RowCollector,
    F: FnOnce(Vec<usize>) -> Result<MatrixHeader>,
{
    let mut group_counts = vec![0usize; group_count];
    let mut buffer = BTreeMap::new();
    let mut next_index = 0usize;

    while let Ok(result) = rx.recv() {
        buffer.insert(result.index, result);
        flush_ready_entries(
            &mut buffer,
            &mut collector,
            &mut group_counts,
            &mut next_index,
        )?;
    }

    flush_ready_entries(
        &mut buffer,
        &mut collector,
        &mut group_counts,
        &mut next_index,
    )?;

    if next_index != task_count {
        return Err(anyhow!(
            "Streamed matrix received {} of {} expected rows",
            next_index,
            task_count
        ));
    }

    let header = header_builder(group_counts)?;
    collector.finalize(header)
}

fn flush_ready_entries<C: RowCollector>(
    buffer: &mut BTreeMap<usize, RegionResult>,
    collector: &mut C,
    group_counts: &mut [usize],
    next_index: &mut usize,
) -> Result<()> {
    loop {
        let key = *next_index;
        let Some(entry) = buffer.remove(&key) else {
            break;
        };

        if let Some(row) = entry.row {
            collector.on_row(row)?;
            if let Some(count) = group_counts.get_mut(entry.group_index) {
                *count += 1;
            }
        }

        *next_index += 1;
    }

    Ok(())
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
    let (tx, rx) = mpsc::channel();
    let aggregator_handle =
        spawn_row_aggregator(rx, collector, group_count, task_count, header_builder)?;

    if task_count == 0 {
        drop(tx);
        return aggregator_handle
            .join()
            .map_err(|_| anyhow!("Matrix aggregation thread panicked"))?;
    }

    let pool = ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build()
        .context("Failed to initialise rayon thread pool for pipeline scheduling")?;

    let sample_paths_for_workers = Arc::clone(&sample_paths);
    let metadata_for_workers = Arc::clone(&metadata);
    let tx_template = tx.clone();

    let compute_result = pool.install(|| {
        tasks
            .into_par_iter()
            .map_init(
                move || {
                    (
                        WorkerSamples::new(Arc::clone(&sample_paths_for_workers)),
                        tx_template.clone(),
                        Arc::clone(&metadata_for_workers),
                    )
                },
                |state, task| {
                    let (worker_samples, sender, metadata) = state;
                    let metadata_ref = metadata.as_ref();

                    let RegionTask {
                        index,
                        group_index,
                        record,
                    } = task;

                    let samples = worker_samples.samples()?;
                    let plan = mode.plan_for(&record, metadata_ref);
                    let maybe_values = compute_row(
                        samples.as_mut_slice(),
                        &record,
                        &plan,
                        general,
                        mode.nan_after_end(metadata_ref),
                    )?;
                    let row = maybe_values
                        .map(|values| mode.postprocess_row(record, values, metadata_ref));

                    sender
                        .send(RegionResult {
                            index,
                            group_index,
                            row,
                        })
                        .map_err(|err| anyhow!("Failed to stream computed row: {err}"))?;

                    Ok::<(), anyhow::Error>(())
                },
            )
            .try_reduce(|| (), |_, _| Ok::<(), anyhow::Error>(()))
    });

    drop(tx);

    let aggregation = aggregator_handle
        .join()
        .map_err(|_| anyhow!("Matrix aggregation thread panicked"))??;

    compute_result?;

    Ok(aggregation)
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
