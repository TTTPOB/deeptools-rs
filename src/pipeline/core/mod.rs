use std::collections::{HashMap, HashSet};
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context, Result, anyhow, bail};
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

use crate::config::{GeneralOptions, GtfOptions, SortRegions};
use crate::io::{BedReadError, BedRecord, BigWigReader, SharedBigWigReader, load_gtf_records};
use crate::pipeline::matrix::{MatrixHeader, MatrixRow};

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

mod collector;
pub use collector::{RowCollector, InMemoryCollector, FileCollector, GroupBucketCollector};

mod coalesce;
use coalesce::{
    CoalesceStrategy, COALESCE_CLAMP_MAX, create_batches, estimate_coalesce_gap,
};

mod worker;
use worker::process_batch;
pub use worker::compute_row;

/// Output dispatch strategy for the chunk-collect pipeline.
enum OutputStrategy {
    /// Input order already matches compute-sorted order (or SortRegions::No).
    /// Results can be streamed directly to the collector in chunk order.
    StreamOrdered,
    /// SortRegions::Keep but input order differs from compute-sorted order.
    /// Collect all results in memory, sort by orig_idx, then emit.
    InMemoryKeep,
    /// SortRegions::Ascend or Descend — bucket rows by group, then let the
    /// downstream sort handle ordering within each group.
    InMemoryGroupBucket,
}

fn into_chunks<T>(items: Vec<T>, chunk_size: usize) -> Vec<Vec<T>> {
    let mut chunks = Vec::new();
    let mut iter = items.into_iter();
    loop {
        let chunk: Vec<T> = iter.by_ref().take(chunk_size).collect();
        if chunk.is_empty() {
            break;
        }
        chunks.push(chunk);
    }
    chunks
}

type BatchResult = (usize, usize, Option<MatrixRow>);

/// Internal work item carrying the original index and I/O sort key.
pub(crate) struct WorkItem {
    orig_idx: usize,
    group_index: usize,
    record: Arc<BedRecord>,
    query_start: i64,
    query_end: i64,
}

fn input_order_is_compute_sorted(items: &[WorkItem]) -> bool {
    items.windows(2).all(|pair| {
        let a = &pair[0];
        let b = &pair[1];
        (a.record.chrom.as_ref(), a.query_start, a.query_end)
            <= (b.record.chrom.as_ref(), b.query_start, b.query_end)
    })
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

    // ── Phase 2: Determine output strategy, then sort for I/O locality ───
    let already_sorted = input_order_is_compute_sorted(&work_items);
    let output_strategy = match general.sort_regions {
        SortRegions::Keep if already_sorted => OutputStrategy::StreamOrdered,
        SortRegions::No => OutputStrategy::StreamOrdered,
        SortRegions::Keep => OutputStrategy::InMemoryKeep,
        SortRegions::Ascend | SortRegions::Descend => OutputStrategy::InMemoryGroupBucket,
    };

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

    // ── Phase 4: Create coalesced batches ───────────────────────────────
    let coalesce_gap = estimate_coalesce_gap(&work_items);
    let coalesce_strategy = if coalesce_gap >= COALESCE_CLAMP_MAX {
        CoalesceStrategy::NoCoalesce
    } else {
        CoalesceStrategy::Coalesce(coalesce_gap)
    };
    let batches = create_batches(work_items, &coalesce_strategy);
    eprintln!(
        "[coalesce-gap] strategy={:?} batches={} items={} ratio={:.2}",
        match &coalesce_strategy {
            CoalesceStrategy::Coalesce(g) => format!("coalesce({g})"),
            CoalesceStrategy::NoCoalesce => "no-coalesce".into(),
        },
        batches.len(),
        task_count,
        batches.len() as f64 / task_count as f64
    );

    // ── Phase 5: Build thread pool ──────────────────────────────────────
    let pool = ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build()
        .context("Failed to initialise rayon thread pool for pipeline scheduling")?;

    // ── Phase 6: Split batches into chunks ──────────────────────────────
    let chunk_size = std::cmp::max(256, batches.len() / (thread_count * 4));
    let chunks = into_chunks(batches, chunk_size);

    let metadata_ref = metadata.as_ref();

    // ── Phase 7: Dispatch based on output strategy ──────────────────────
    match output_strategy {
        OutputStrategy::StreamOrdered => {
            let mut collector = collector;
            let mut group_counts = vec![0usize; group_count];

            for chunk in chunks {
                let sample_paths_c = Arc::clone(&sample_paths);
                let shared_c = Arc::clone(&shared_readers);

                let chunk_results: Vec<Result<Vec<BatchResult>>> = pool.install(|| {
                    chunk
                        .into_par_iter()
                        .map_init(
                            move || {
                                WorkerSamples::from_shared(
                                    Arc::clone(&sample_paths_c),
                                    Arc::clone(&shared_c),
                                )
                            },
                            |worker_samples, batch| {
                                let samples = worker_samples.samples()?;
                                process_batch(
                                    samples.as_mut_slice(),
                                    batch,
                                    mode,
                                    general,
                                    metadata_ref,
                                )
                            },
                        )
                        .collect()
                });

                // Emit results in order — par_iter preserves input order
                for batch_result in chunk_results {
                    let rows = batch_result?;
                    for (_orig_idx, group_index, row) in rows {
                        if let Some(row) = row {
                            collector.on_row(row)?;
                            group_counts[group_index] += 1;
                        }
                    }
                }
            }

            drop(shared_readers);
            let header = header_builder(group_counts)?;
            collector.finalize(header)
        }

        OutputStrategy::InMemoryKeep => {
            let mut all_results: Vec<BatchResult> = Vec::with_capacity(task_count);

            for chunk in chunks {
                let sample_paths_c = Arc::clone(&sample_paths);
                let shared_c = Arc::clone(&shared_readers);

                let chunk_results: Vec<Result<Vec<BatchResult>>> = pool.install(|| {
                    chunk
                        .into_par_iter()
                        .map_init(
                            move || {
                                WorkerSamples::from_shared(
                                    Arc::clone(&sample_paths_c),
                                    Arc::clone(&shared_c),
                                )
                            },
                            |worker_samples, batch| {
                                let samples = worker_samples.samples()?;
                                process_batch(
                                    samples.as_mut_slice(),
                                    batch,
                                    mode,
                                    general,
                                    metadata_ref,
                                )
                            },
                        )
                        .collect()
                });

                for batch_result in chunk_results {
                    let rows = batch_result?;
                    all_results.extend(rows);
                }
            }

            drop(shared_readers);

            // Restore original input order
            all_results.sort_by_key(|r| r.0);

            let mut collector = collector;
            let mut group_counts = vec![0usize; group_count];
            for (_orig_idx, group_index, row) in all_results {
                if let Some(row) = row {
                    collector.on_row(row)?;
                    group_counts[group_index] += 1;
                }
            }

            let header = header_builder(group_counts)?;
            collector.finalize(header)
        }

        OutputStrategy::InMemoryGroupBucket => {
            // Use GroupBucketCollector to bucket rows by group; downstream
            // sort (Ascend/Descend) is applied by the caller on the
            // returned MatrixData.
            let mut all_results: Vec<BatchResult> = Vec::with_capacity(task_count);

            for chunk in chunks {
                let sample_paths_c = Arc::clone(&sample_paths);
                let shared_c = Arc::clone(&shared_readers);

                let chunk_results: Vec<Result<Vec<BatchResult>>> = pool.install(|| {
                    chunk
                        .into_par_iter()
                        .map_init(
                            move || {
                                WorkerSamples::from_shared(
                                    Arc::clone(&sample_paths_c),
                                    Arc::clone(&shared_c),
                                )
                            },
                            |worker_samples, batch| {
                                let samples = worker_samples.samples()?;
                                process_batch(
                                    samples.as_mut_slice(),
                                    batch,
                                    mode,
                                    general,
                                    metadata_ref,
                                )
                            },
                        )
                        .collect()
                });

                for batch_result in chunk_results {
                    let rows = batch_result?;
                    all_results.extend(rows);
                }
            }

            drop(shared_readers);

            // Feed results through the generic collector interface.
            // For Ascend/Descend, the caller passes an InMemoryCollector
            // and handles sorting on the returned MatrixData, so we just
            // need to emit rows in a reasonable order.  Sort by orig_idx
            // so the per-group order is deterministic.
            all_results.sort_by_key(|r| r.0);

            let mut collector = collector;
            let mut group_counts = vec![0usize; group_count];
            for (_orig_idx, group_index, row) in all_results {
                if let Some(row) = row {
                    collector.on_row(row)?;
                    group_counts[group_index] += 1;
                }
            }

            let header = header_builder(group_counts)?;
            collector.finalize(header)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::worker::{aggregate_slice, index_from_coordinate};
    use crate::config::AverageTypeBins;
    use crate::io::readers::bed::Strand;

    fn work_item(idx: usize, chrom: &str, start: i64, end: i64) -> WorkItem {
        WorkItem {
            orig_idx: idx,
            group_index: 0,
            record: Arc::new(BedRecord {
                chrom: Arc::from(chrom),
                start: start as u32,
                end: end as u32,
                name: None,
                score: None,
                score_raw: None,
                strand: Strand::Unstranded,
                strand_raw: None,
                extra_fields: vec![],
            }),
            query_start: start,
            query_end: end,
        }
    }

    #[test]
    fn test_input_is_compute_sorted_true() {
        let items = vec![
            work_item(0, "chr1", 100, 200),
            work_item(1, "chr1", 300, 400),
            work_item(2, "chr2", 50, 150),
        ];
        assert!(input_order_is_compute_sorted(&items));
    }

    #[test]
    fn test_input_is_compute_sorted_false() {
        let items = vec![
            work_item(0, "chr2", 100, 200),
            work_item(1, "chr1", 300, 400),
        ];
        assert!(!input_order_is_compute_sorted(&items));
    }

    #[test]
    fn test_input_is_compute_sorted_empty() {
        let items: Vec<WorkItem> = vec![];
        assert!(input_order_is_compute_sorted(&items));
    }

    #[test]
    fn test_input_is_compute_sorted_single() {
        let items = vec![work_item(0, "chr1", 100, 200)];
        assert!(input_order_is_compute_sorted(&items));
    }

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
