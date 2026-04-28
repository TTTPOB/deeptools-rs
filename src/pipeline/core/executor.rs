use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{Context, Result};
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

use crate::config::{GeneralOptions, SortRegions};
use crate::io::{BedRecord, SharedBigWigReader, SharedBlockCache};
use crate::pipeline::matrix::MatrixHeader;

use super::collector::{GroupBucketCollector, RowCollector};
use super::coalesce::{
    CoalesceStrategy, COALESCE_CLAMP_MAX, create_batches, estimate_coalesce_gap,
};
use super::regions::RegionTask;
use super::samples::WorkerSamples;
use super::traits::{PipelineMode, RegionPlan};
use super::worker::process_batch;

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

    if !already_sorted {
        work_items.sort_by(|a, b| {
            a.record
                .chrom
                .cmp(&b.record.chrom)
                .then(a.query_start.cmp(&b.query_start))
                .then(a.query_end.cmp(&b.query_end))
        });
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
                                collector.on_row(row)?;
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
            let total_bins = mode.total_bins(metadata_ref);
            let sample_count = sample_paths.len();
            let mut bucket_collector =
                GroupBucketCollector::new(group_count, sample_count, total_bins);

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
                    for (_orig_idx, group_index, row) in batch_result? {
                        if let Some(row) = row {
                            bucket_collector.on_row_with_group(group_index, row)?;
                        }
                    }
                }
            }

            let matrix_data = bucket_collector.finalize_grouped(header_builder)?;
            let mut collector = collector;
            let header = matrix_data.header;
            for row in matrix_data.rows {
                collector.on_row(row)?;
            }
            collector.finalize(header)
        }
    }
}
