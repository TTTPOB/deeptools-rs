use std::sync::Arc;

use crate::io::BedRecord;

use super::executor::WorkItem;

// ── Query coalescing ──────────────────────────────────────────────────────

pub(super) const COALESCE_CLAMP_MAX: i64 = 2000;

pub(super) enum CoalesceStrategy {
    Coalesce(i64),
    NoCoalesce,
}

/// Estimate a coalescing gap threshold from the actual distribution of gaps
/// between consecutive same-chromosome items in the sorted work list.
///
/// Returns a threshold in `[100, 2000]` based on the 75th percentile of
/// observed gaps.  Falls back to 500 bp when there are too few gaps (< 10).
pub(super) fn estimate_coalesce_gap(work_items: &[WorkItem], verbose: bool) -> i64 {
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
        if verbose {
            eprintln!(
                "[coalesce-gap] {} gaps (< 10), using default threshold: 500",
                gaps.len()
            );
        }
        return 500;
    }

    gaps.sort_unstable();

    let n = gaps.len();
    let p50 = gaps[n / 2];
    let p75 = gaps[(n * 3) / 4];
    let threshold = p75.clamp(100, 2000);

    if verbose {
        eprintln!(
            "[coalesce-gap] n_gaps={} p50={} p75={} threshold={}",
            n, p50, p75, threshold
        );
    }
    threshold
}

/// A batch of consecutive work items on the same chromosome whose query
/// windows overlap or are separated by at most the caller-supplied
/// `coalesce_gap` threshold.  Records are **moved** (not cloned) from
/// `WorkItem`s, and `work_items` is consumed.
pub(super) struct CoalescedBatch {
    /// Items in original sorted order: (orig_idx, group_index, record).
    pub(super) items: Vec<(usize, usize, Arc<BedRecord>)>,
    /// Start of the merged query window (minimum of all item windows).
    pub(super) query_start: i64,
    /// End of the merged query window (maximum of all item windows).
    pub(super) query_end: i64,
}

/// Create batches according to the chosen strategy.
pub(super) fn create_batches(
    work_items: Vec<WorkItem>,
    strategy: &CoalesceStrategy,
) -> Vec<CoalescedBatch> {
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
