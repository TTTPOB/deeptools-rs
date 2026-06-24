use std::sync::Arc;

use crate::io::BedRecord;

use super::executor::WorkItem;

// ── Query coalescing ──────────────────────────────────────────────────────

pub(super) const COALESCE_CLAMP_MAX: i64 = 2000;
pub(super) const COALESCE_CLAMP_MIN: i64 = 100;
pub(super) const COALESCE_DEFAULT_GAP: i64 = 500;
pub(super) const COALESCE_MIN_GAP_SAMPLES: usize = 10;

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
        return COALESCE_DEFAULT_GAP;
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

    if gaps.len() < COALESCE_MIN_GAP_SAMPLES {
        if verbose {
            eprintln!(
                "[coalesce-gap] {} gaps (< {COALESCE_MIN_GAP_SAMPLES}), using default threshold: {COALESCE_DEFAULT_GAP}",
                gaps.len()
            );
        }
        return COALESCE_DEFAULT_GAP;
    }

    gaps.sort_unstable();

    let n = gaps.len();
    let p50 = gaps[n / 2];
    let p75 = gaps[(n * 3) / 4];
    let threshold = p75.clamp(COALESCE_CLAMP_MIN, COALESCE_CLAMP_MAX);

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::BedRecord;
    use crate::io::readers::bed::Strand;
    use std::sync::Arc;

    use super::super::executor::WorkItem;

    // Helper to build a WorkItem with minimal fields.
    fn make_item(chrom: &str, query_start: i64, query_end: i64, orig_idx: usize) -> WorkItem {
        WorkItem {
            orig_idx,
            group_index: 0,
            record: Arc::new(BedRecord {
                chrom: Arc::from(chrom),
                start: query_start.max(0) as u32,
                end: query_end.max(0) as u32,
                name: None,
                bed_field_count: None,
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

    // ── estimate_coalesce_gap ─────────────────────────────────────────────

    #[test]
    fn test_estimate_gap_empty() {
        let items: Vec<WorkItem> = vec![];
        assert_eq!(estimate_coalesce_gap(&items, false), 500);
    }

    #[test]
    fn test_estimate_gap_single_item() {
        let items = vec![make_item("chr1", 0, 100, 0)];
        assert_eq!(estimate_coalesce_gap(&items, false), 500);
    }

    #[test]
    fn test_estimate_gap_fewer_than_10_gaps() {
        // 5 items on the same chromosome → 4 positive gaps, < 10 → fallback 500
        let items: Vec<WorkItem> = (0..5)
            .map(|i| make_item("chr1", i * 200, i * 200 + 100, i as usize))
            .collect();
        assert_eq!(estimate_coalesce_gap(&items, false), 500);
    }

    #[test]
    fn test_estimate_gap_normal_case() {
        // 11 items spaced 300 bp apart → 10 gaps of 200 each (end=100, next start=300).
        // p75 of [200]*10 = 200, clamped to [100,2000] → 200.
        let items: Vec<WorkItem> = (0..11)
            .map(|i| make_item("chr1", i * 300, i * 300 + 100, i as usize))
            .collect();
        let gap = estimate_coalesce_gap(&items, false);
        assert!(
            (100..=2000).contains(&gap),
            "expected gap in [100,2000], got {}",
            gap
        );
        assert_eq!(gap, 200);
    }

    #[test]
    fn test_estimate_gap_all_small_clamped_to_100() {
        // 11 items spaced 1 bp apart → gaps of 1, p75 = 1, clamped to 100.
        let items: Vec<WorkItem> = (0..11)
            .map(|i| make_item("chr1", i * 2, i * 2 + 1, i as usize))
            .collect();
        assert_eq!(estimate_coalesce_gap(&items, false), 100);
    }

    #[test]
    fn test_estimate_gap_all_large_clamped_to_2000() {
        // 11 items spaced 10_000 bp apart → gaps of 9_900, p75 > 2000, clamped to 2000.
        let items: Vec<WorkItem> = (0..11)
            .map(|i| make_item("chr1", i * 10_000, i * 10_000 + 100, i as usize))
            .collect();
        assert_eq!(estimate_coalesce_gap(&items, false), 2000);
    }

    // ── create_batches – NoCoalesce ───────────────────────────────────────

    #[test]
    fn test_no_coalesce_each_item_own_batch() {
        let items = vec![
            make_item("chr1", 0, 100, 0),
            make_item("chr1", 50, 150, 1),
            make_item("chr2", 0, 100, 2),
        ];
        let batches = create_batches(items, &CoalesceStrategy::NoCoalesce);
        assert_eq!(batches.len(), 3);
        for (i, batch) in batches.iter().enumerate() {
            assert_eq!(batch.items.len(), 1);
            assert_eq!(batch.items[0].0, i); // orig_idx
        }
    }

    #[test]
    fn test_no_coalesce_empty() {
        let batches = create_batches(vec![], &CoalesceStrategy::NoCoalesce);
        assert!(batches.is_empty());
    }

    // ── create_batches – Coalesce (exercises create_coalesced_batches) ────

    #[test]
    fn test_coalesce_adjacent_overlapping_merged() {
        // Items overlap: 0-200 and 100-300 → one batch.
        let items = vec![make_item("chr1", 0, 200, 0), make_item("chr1", 100, 300, 1)];
        let batches = create_batches(items, &CoalesceStrategy::Coalesce(500));
        assert_eq!(batches.len(), 1);
        assert_eq!(batches[0].items.len(), 2);
        assert_eq!(batches[0].query_start, 0);
        assert_eq!(batches[0].query_end, 300);
    }

    #[test]
    fn test_coalesce_separated_by_exactly_gap_merged() {
        // Gap between items == coalesce_gap → merged.
        // item 0: 0-100, item 1: 600-700, gap = 500 == coalesce_gap 500.
        let items = vec![make_item("chr1", 0, 100, 0), make_item("chr1", 600, 700, 1)];
        let batches = create_batches(items, &CoalesceStrategy::Coalesce(500));
        assert_eq!(batches.len(), 1);
        assert_eq!(batches[0].query_start, 0);
        assert_eq!(batches[0].query_end, 700);
    }

    #[test]
    fn test_coalesce_separated_more_than_gap_split() {
        // Gap = 501 > coalesce_gap 500 → two separate batches.
        let items = vec![make_item("chr1", 0, 100, 0), make_item("chr1", 601, 700, 1)];
        let batches = create_batches(items, &CoalesceStrategy::Coalesce(500));
        assert_eq!(batches.len(), 2);
        assert_eq!(batches[0].query_start, 0);
        assert_eq!(batches[0].query_end, 100);
        assert_eq!(batches[1].query_start, 601);
        assert_eq!(batches[1].query_end, 700);
    }

    #[test]
    fn test_coalesce_different_chroms_always_split() {
        let items = vec![make_item("chr1", 0, 100, 0), make_item("chr2", 0, 100, 1)];
        let batches = create_batches(items, &CoalesceStrategy::Coalesce(10_000));
        assert_eq!(batches.len(), 2);
        assert_eq!(*batches[0].items[0].2.chrom, *"chr1");
        assert_eq!(*batches[1].items[0].2.chrom, *"chr2");
    }

    #[test]
    fn test_coalesce_single_item() {
        let items = vec![make_item("chr1", 50, 150, 0)];
        let batches = create_batches(items, &CoalesceStrategy::Coalesce(500));
        assert_eq!(batches.len(), 1);
        assert_eq!(batches[0].query_start, 50);
        assert_eq!(batches[0].query_end, 150);
    }

    #[test]
    fn test_coalesce_empty_input() {
        let batches = create_batches(vec![], &CoalesceStrategy::Coalesce(500));
        assert!(batches.is_empty());
    }

    #[test]
    fn test_coalesce_batch_query_bounds_are_min_max() {
        // Three items; middle one extends the window furthest right.
        let items = vec![
            make_item("chr1", 100, 200, 0),
            make_item("chr1", 150, 500, 1),
            make_item("chr1", 300, 400, 2),
        ];
        let batches = create_batches(items, &CoalesceStrategy::Coalesce(500));
        assert_eq!(batches.len(), 1);
        assert_eq!(batches[0].query_start, 100);
        assert_eq!(batches[0].query_end, 500);
        assert_eq!(batches[0].items.len(), 3);
    }
}
