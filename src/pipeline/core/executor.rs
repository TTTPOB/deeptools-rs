use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{Context, Result};
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

use crate::config::{GeneralOptions, SortRegions};
use crate::io::{BedRecord, SharedBigWigReader, SharedBlockCache};
use crate::pipeline::matrix::{MatrixHeader, compute_sort_metric};

use super::collector::{FileCollector, RowCollector};
use super::coalesce::{
    CoalesceStrategy, COALESCE_CLAMP_MAX, create_batches, estimate_coalesce_gap,
};
use super::regions::{RegionTask, normalize_sort_sample_indices};
use super::samples::WorkerSamples;
use super::spill::HybridBucketCollector;
use super::traits::{PipelineMode, RegionPlan};
use super::worker::process_batch;

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
    let block_cache = Arc::new(SharedBlockCache::new());
    let shared_readers = Arc::new(
        sample_paths
            .iter()
            .map(|path| {
                SharedBigWigReader::open_with_cache(path, Arc::clone(&block_cache))
                    .map(Arc::new)
                    .with_context(|| {
                        format!("Failed to open bigWig file '{}'", path.display())
                    })
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
    let chunk_size = std::cmp::max(256, batches.len() / (thread_count * 4));
    let chunks = into_chunks(batches, chunk_size);

    let metadata_ref = metadata.as_ref();

    // ── Phase 7: Dispatch based on output strategy ──────────────────────
    match output_strategy {
        OutputStrategy::StreamOrdered => {
            let mut collector = collector;
            let mut group_counts = vec![0usize; group_count];

            let compute_result: Result<()> = (|| {
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
                        SortRegions::Ascend => bucket_collector.finalize_sorted(
                            true,
                            header_builder,
                            |gi, row| collector.on_row(gi, row),
                        )?,
                        SortRegions::Descend => bucket_collector.finalize_sorted(
                            false,
                            header_builder,
                            |gi, row| collector.on_row(gi, row),
                        )?,
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
