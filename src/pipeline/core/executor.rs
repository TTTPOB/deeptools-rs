use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{Context, Result};
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

use crate::config::{GeneralOptions, SortRegions};
use crate::io::readers::block_cache::compute_per_file_block_cache_capacity;
use crate::io::{BedRecord, BigWigFile};
use crate::pipeline::matrix::{MatrixHeader, compute_sort_metric};

use super::coalesce::{
    COALESCE_CLAMP_MAX, CoalesceStrategy, create_batches, estimate_coalesce_gap,
};
use super::collector::FileCollector;
use super::regions::{RegionTask, normalize_sort_sample_indices};
use super::samples::WorkerSamples;
use super::spill::HybridBucketCollector;
use super::traits::{PipelineMode, RegionPlan};
use super::worker::process_batch;

const MIN_CHUNK_SIZE: usize = 256;
const CHUNKS_PER_THREAD: usize = 4;

/// Output dispatch strategy for the chunk-collect pipeline.
enum OutputStrategy {
    /// Input order already matches compute-sorted order (or SortRegions::No
    /// with per-group I/O sort). Results can be streamed directly to the
    /// collector in chunk order.
    StreamOrdered,
    /// Everything else: keep+coalesced, ascend, descend. Rows are bucketed
    /// by group and emitted via HybridBucketCollector.
    HybridBucket,
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

type BatchResult = (usize, usize, Option<crate::pipeline::matrix::MatrixRow>);

/// Internal work item carrying the original index and I/O sort key.
pub(crate) struct WorkItem {
    pub(crate) orig_idx: usize,
    pub(crate) group_index: usize,
    pub(crate) record: Arc<BedRecord>,
    pub(crate) query_start: i64,
    pub(crate) query_end: i64,
}

pub(super) fn input_order_is_compute_sorted(items: &[WorkItem]) -> bool {
    items.windows(2).all(|pair| {
        let a = &pair[0];
        let b = &pair[1];
        (a.record.chrom.as_ref(), a.query_start, a.query_end)
            <= (b.record.chrom.as_ref(), b.query_start, b.query_end)
    })
}

pub fn execute_mode<M, F>(
    tasks: Vec<RegionTask>,
    general: &GeneralOptions,
    sample_paths: Arc<Vec<PathBuf>>,
    collector: FileCollector,
    thread_count: usize,
    mode: &M,
    metadata: Arc<M::Metadata>,
    header_builder: F,
    group_count: usize,
    task_count: usize,
    sample_count: usize,
) -> Result<()>
where
    M: PipelineMode,
    F: FnOnce(Vec<usize>) -> Result<MatrixHeader> + Send + 'static,
{
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
        SortRegions::Keep => OutputStrategy::HybridBucket,
        SortRegions::Ascend | SortRegions::Descend => OutputStrategy::HybridBucket,
    };

    // Compute per-group item counts before sorting (needed for per-group I/O sort).
    let group_item_counts: Vec<usize> = {
        let mut counts = vec![0usize; group_count];
        for item in &work_items {
            counts[item.group_index] += 1;
        }
        counts
    };

    // Sort for I/O locality.
    match general.sort_regions {
        SortRegions::No if !already_sorted => {
            // Per-group I/O sort: sort within each group's span.
            // Work items are laid out contiguously by group (group 0 items, then group 1, etc.).
            // We need to figure out group spans. Since work_items are in original order
            // (group0 items, group1 items, ...), we can use group_item_counts.
            let mut start = 0;
            for &count in &group_item_counts {
                let end = start + count;
                work_items[start..end].sort_by(|a, b| {
                    a.record
                        .chrom
                        .cmp(&b.record.chrom)
                        .then(a.query_start.cmp(&b.query_start))
                        .then(a.query_end.cmp(&b.query_end))
                });
                start = end;
            }
        }
        _ if !already_sorted => {
            // Global sort for keep/ascend/descend.
            work_items.sort_by(|a, b| {
                a.record
                    .chrom
                    .cmp(&b.record.chrom)
                    .then(a.query_start.cmp(&b.query_start))
                    .then(a.query_end.cmp(&b.query_end))
            });
        }
        _ => {}
    }

    // Empty input: build an empty header and return.
    if task_count == 0 {
        let header = header_builder(vec![0; group_count])?;
        return collector.finalize(header);
    }

    // ── Phase 3: Open shared bigWig readers once ────────────────────────
    // Each file gets its own block cache whose capacity is a fair share of
    // the global budget, preventing cross-sample block collisions.
    let sample_count_for_cache = sample_paths.len();
    let shared_readers = Arc::new(
        sample_paths
            .iter()
            .enumerate()
            .map(|(sample_index, path)| {
                let cache_capacity =
                    compute_per_file_block_cache_capacity(sample_count_for_cache, sample_index);
                BigWigFile::open_with_block_cache_capacity(path, cache_capacity)
                    .map(Arc::new)
                    .with_context(|| format!("Failed to open bigWig file '{}'", path.display()))
            })
            .collect::<Result<Vec<_>>>()?,
    );

    // ── Phase 4: Create coalesced batches ───────────────────────────────
    let coalesce_gap = estimate_coalesce_gap(&work_items, general.verbose);
    let coalesce_strategy = if coalesce_gap >= COALESCE_CLAMP_MAX {
        CoalesceStrategy::NoCoalesce
    } else {
        CoalesceStrategy::Coalesce(coalesce_gap)
    };
    let batches = create_batches(work_items, &coalesce_strategy);
    if general.verbose {
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
    }

    // ── Phase 5: Build thread pool ──────────────────────────────────────
    let pool = ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build()
        .context("Failed to initialise rayon thread pool for pipeline scheduling")?;

    // ── Phase 6: Split batches into chunks ──────────────────────────────
    let chunk_size = std::cmp::max(
        MIN_CHUNK_SIZE,
        batches.len() / (thread_count * CHUNKS_PER_THREAD),
    );
    let chunks = into_chunks(batches, chunk_size);

    let metadata_ref = metadata.as_ref();

    // Helper closure: dispatch a single chunk to the thread pool and collect
    // the per-batch results.  Extracted to eliminate identical code in both
    // StreamOrdered and HybridBucket arms.
    let dispatch_chunk = |chunk: Vec<_>| -> Vec<Result<Vec<BatchResult>>> {
        let sample_paths_c = Arc::clone(&sample_paths);
        let shared_c = Arc::clone(&shared_readers);
        pool.install(|| {
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
                        let (samples, work_buf, decode_buf) = worker_samples.samples_and_bufs()?;
                        process_batch(
                            samples.as_mut_slice(),
                            batch,
                            mode,
                            general,
                            metadata_ref,
                            work_buf,
                            decode_buf,
                        )
                    },
                )
                .collect()
        })
    };

    // ── Phase 7: Dispatch based on output strategy ──────────────────────
    match output_strategy {
        OutputStrategy::StreamOrdered => {
            let mut collector = collector;
            let mut group_counts = vec![0usize; group_count];

            let compute_result: Result<()> = (|| {
                for chunk in chunks {
                    let chunk_results = dispatch_chunk(chunk);

                    // Emit results in order — par_iter preserves input order
                    for batch_result in chunk_results {
                        let rows = batch_result?;
                        for (_orig_idx, group_index, row) in rows {
                            if let Some(row) = row {
                                collector.on_row(group_index, row)?;
                                group_counts[group_index] += 1;
                            }
                        }
                    }
                }
                Ok(())
            })();

            match compute_result {
                Ok(()) => {
                    let header = header_builder(group_counts)?;
                    collector.finalize(header)
                }
                Err(e) => {
                    collector.abort();
                    Err(e)
                }
            }
        }

        OutputStrategy::HybridBucket => {
            let total_bins = mode.total_bins(metadata_ref);

            // Compute sort sample indices for sort key computation.
            let sort_sample_indices =
                normalize_sort_sample_indices(general.sort_using_samples.as_ref(), sample_count)?;

            let mut bucket_collector =
                HybridBucketCollector::new(group_count, sample_count, total_bins);

            let compute_result: Result<()> = (|| {
                for chunk in chunks {
                    let chunk_results = dispatch_chunk(chunk);

                    for batch_result in chunk_results {
                        for (orig_idx, group_index, row) in batch_result? {
                            if let Some(row) = row {
                                let sort_key = compute_sort_metric(
                                    &row,
                                    general.sort_using,
                                    sort_sample_indices.as_deref(),
                                );
                                bucket_collector.push(row, orig_idx, group_index, sort_key)?;
                            }
                        }
                    }
                }
                Ok(())
            })();

            match compute_result {
                Ok(()) => {
                    let mut collector = collector;
                    let header = match general.sort_regions {
                        SortRegions::Ascend => {
                            bucket_collector.finalize_sorted(true, header_builder, |gi, row| {
                                collector.on_row(gi, row)
                            })?
                        }
                        SortRegions::Descend => {
                            bucket_collector.finalize_sorted(false, header_builder, |gi, row| {
                                collector.on_row(gi, row)
                            })?
                        }
                        SortRegions::Keep => bucket_collector.finalize_keep_order(
                            task_count,
                            header_builder,
                            |gi, row| collector.on_row(gi, row),
                        )?,
                        SortRegions::No => {
                            // sort=No should use StreamOrdered, not HybridBucket.
                            unreachable!("sort=No should use StreamOrdered strategy")
                        }
                    };
                    collector.finalize(header)
                }
                Err(e) => {
                    collector.abort();
                    Err(e)
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use crate::io::BedRecord;
    use crate::io::readers::bed::Strand;

    use super::*;

    /// Build a minimal WorkItem for testing sort/order checks.
    fn make_item(chrom: &str, query_start: i64, query_end: i64) -> WorkItem {
        WorkItem {
            orig_idx: 0,
            group_index: 0,
            record: Arc::new(BedRecord {
                chrom: Arc::from(chrom),
                start: query_start.max(0) as u32,
                end: query_end.max(0) as u32,
                name: None,
                score: None,
                score_raw: None,
                strand: Strand::Unstranded,
                strand_raw: None,
                extra_fields: Vec::new(),
            }),
            query_start,
            query_end,
        }
    }

    // ── into_chunks ────────────────────────────────────────────────────────

    #[test]
    fn into_chunks_evenly_divisible() {
        // 6 items, chunk_size 2 → 3 chunks of 2
        let items = vec![1, 2, 3, 4, 5, 6];
        let chunks = into_chunks(items, 2);
        assert_eq!(chunks, vec![vec![1, 2], vec![3, 4], vec![5, 6]]);
    }

    #[test]
    fn into_chunks_last_chunk_smaller() {
        // 5 items, chunk_size 2 → 2 full chunks + 1 remainder chunk
        let items = vec![1, 2, 3, 4, 5];
        let chunks = into_chunks(items, 2);
        assert_eq!(chunks, vec![vec![1, 2], vec![3, 4], vec![5]]);
    }

    #[test]
    fn into_chunks_size_larger_than_items() {
        // chunk_size exceeds item count → single chunk containing all items
        let items = vec![1, 2, 3];
        let chunks = into_chunks(items, 10);
        assert_eq!(chunks, vec![vec![1, 2, 3]]);
    }

    #[test]
    fn into_chunks_empty_input() {
        // empty input → no chunks
        let items: Vec<i32> = vec![];
        let chunks = into_chunks(items, 5);
        assert!(chunks.is_empty());
    }

    #[test]
    fn into_chunks_size_one() {
        // chunk_size 1 → each item is its own chunk
        let items = vec![10, 20, 30];
        let chunks = into_chunks(items, 1);
        assert_eq!(chunks, vec![vec![10], vec![20], vec![30]]);
    }

    // ── input_order_is_compute_sorted ──────────────────────────────────────

    #[test]
    fn sorted_already_sorted() {
        // Items in correct (chrom, start, end) order → true
        let items = vec![
            make_item("chr1", 100, 200),
            make_item("chr1", 200, 300),
            make_item("chr2", 50, 150),
        ];
        assert!(input_order_is_compute_sorted(&items));
    }

    #[test]
    fn sorted_empty_slice() {
        // Empty slice is trivially sorted
        assert!(input_order_is_compute_sorted(&[]));
    }

    #[test]
    fn sorted_single_item() {
        // A single item is always sorted
        let items = vec![make_item("chr1", 100, 200)];
        assert!(input_order_is_compute_sorted(&items));
    }

    #[test]
    fn sorted_unsorted_by_chrom() {
        // chr2 comes before chr1 → not sorted
        let items = vec![make_item("chr2", 100, 200), make_item("chr1", 100, 200)];
        assert!(!input_order_is_compute_sorted(&items));
    }

    #[test]
    fn sorted_same_chrom_unsorted_by_start() {
        // Same chrom, but start decreases → not sorted
        let items = vec![make_item("chr1", 300, 400), make_item("chr1", 100, 200)];
        assert!(!input_order_is_compute_sorted(&items));
    }

    #[test]
    fn sorted_same_chrom_and_start_unsorted_by_end() {
        // Same chrom and start, but end decreases → not sorted
        let items = vec![make_item("chr1", 100, 400), make_item("chr1", 100, 200)];
        assert!(!input_order_is_compute_sorted(&items));
    }

    #[test]
    fn sorted_equal_items() {
        // Identical items satisfy the <= ordering → true
        let items = vec![make_item("chr1", 100, 200), make_item("chr1", 100, 200)];
        assert!(input_order_is_compute_sorted(&items));
    }
}
