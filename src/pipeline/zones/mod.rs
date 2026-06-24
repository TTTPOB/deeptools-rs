use crate::config::{ReferencePoint, ScaleRegionsOptions};
use crate::io::BedRecord;
use crate::io::Strand;
use crate::pipeline::core::{RegionPlan, SignalBin};

pub(crate) mod metagene;

/// Represents a single bin span within the reference-point window.
#[derive(Debug, Clone, Copy)]
pub struct ReferenceBin {
    pub start: i64,
    pub end: i64,
    pub beyond_region: bool,
}

/// Precomputed layout for reference-point binning.
#[derive(Debug, Clone)]
pub struct ReferencePointPlan {
    pub reference: i64,
    pub window_start: i64,
    pub window_end: i64,
    pub bins: Vec<ReferenceBin>,
    pub included_intervals: Option<Vec<(i64, i64)>>,
}

impl SignalBin for ReferenceBin {
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

impl RegionPlan for ReferencePointPlan {
    type Bin = ReferenceBin;

    fn window_start(&self) -> i64 {
        self.window_start
    }

    fn window_end(&self) -> i64 {
        self.window_end
    }

    fn bins(&self) -> &[Self::Bin] {
        &self.bins
    }

    fn included_intervals(&self) -> Option<&[(i64, i64)]> {
        self.included_intervals.as_deref()
    }
}

impl ReferencePointPlan {
    pub fn reference_point(
        record: &BedRecord,
        reference_point: ReferencePoint,
        bin_size: u32,
        upstream_bins: usize,
        downstream_bins: usize,
        keep_exons: bool,
        nan_after_end: bool,
    ) -> Self {
        let reference = reference_coordinate(record, reference_point);
        let (bins, included_intervals) = if keep_exons {
            metagene::reference_bins(
                record,
                reference_point,
                bin_size,
                upstream_bins,
                downstream_bins,
                nan_after_end,
            )
            .map(|(b, intervals)| (b, Some(intervals)))
            .unwrap_or_else(|| {
                (
                    build_bins(record, reference, bin_size, upstream_bins, downstream_bins),
                    None,
                )
            })
        } else {
            (
                build_bins(record, reference, bin_size, upstream_bins, downstream_bins),
                None,
            )
        };

        let window_start = bins.iter().map(|bin| bin.start).min().unwrap_or(reference);
        let window_end = bins.iter().map(|bin| bin.end).max().unwrap_or(reference);

        Self {
            reference,
            window_start,
            window_end,
            bins,
            included_intervals,
        }
    }
}

/// Represents a single bin span within the scale-regions window.
#[derive(Debug, Clone, Copy)]
pub struct ScaleBin {
    pub start: i64,
    pub end: i64,
    pub beyond_region: bool,
}

/// Precomputed layout for scale-regions binning.
#[derive(Debug, Clone)]
pub struct ScaleRegionsPlan {
    pub window_start: i64,
    pub window_end: i64,
    pub bins: Vec<ScaleBin>,
    pub included_intervals: Option<Vec<(i64, i64)>>,
    /// When true the region body (after subtracting unscaled regions) is
    /// shorter than a single bin.  Python short-circuits the entire row
    /// to zeros/NaN in this case (heatmapper.py:402-411).
    pub body_too_short: bool,
}

impl SignalBin for ScaleBin {
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

impl RegionPlan for ScaleRegionsPlan {
    type Bin = ScaleBin;

    fn window_start(&self) -> i64 {
        self.window_start
    }

    fn window_end(&self) -> i64 {
        self.window_end
    }

    fn bins(&self) -> &[Self::Bin] {
        &self.bins
    }

    fn included_intervals(&self) -> Option<&[(i64, i64)]> {
        self.included_intervals.as_deref()
    }

    fn body_too_short(&self) -> bool {
        self.body_too_short
    }
}

impl ScaleRegionsPlan {
    pub fn scale_regions(
        record: &BedRecord,
        options: &ScaleRegionsOptions,
        bin_size: u32,
        keep_exons: bool,
    ) -> Self {
        if keep_exons {
            if let Some(plan) = metagene::scale_bins(record, options, bin_size) {
                return plan;
            }
        }

        let start = record.start as i64;
        let end = record.end as i64;
        let region_len = (end - start).max(0);

        let upstream_bins = (options.upstream / bin_size) as usize;
        let downstream_bins = (options.downstream / bin_size) as usize;
        let unscaled5_bins = (options.unscaled_5_prime / bin_size) as usize;
        let unscaled3_bins = (options.unscaled_3_prime / bin_size) as usize;
        let body_bins = (options.region_body_length / bin_size) as usize;

        let mut bins = Vec::with_capacity(
            upstream_bins + downstream_bins + unscaled5_bins + unscaled3_bins + body_bins,
        );

        let unscaled5_len = (options.unscaled_5_prime as i64).min(region_len);
        let unscaled3_len =
            (options.unscaled_3_prime as i64).min(region_len.saturating_sub(unscaled5_len));
        let body_len = region_len.saturating_sub(unscaled5_len + unscaled3_len);

        match record.strand {
            Strand::Negative => {
                let upstream_len = options.downstream as i64;
                let downstream_len = options.upstream as i64;

                append_bins(
                    &mut bins,
                    start - upstream_len,
                    upstream_len,
                    downstream_bins,
                );
                append_bins(&mut bins, start, unscaled3_len, unscaled3_bins);
                append_bins(&mut bins, start + unscaled3_len, body_len, body_bins);
                append_bins(
                    &mut bins,
                    end - unscaled5_len,
                    unscaled5_len,
                    unscaled5_bins,
                );
                append_bins(&mut bins, end, downstream_len, upstream_bins);
            }
            _ => {
                let upstream_len = options.upstream as i64;
                let downstream_len = options.downstream as i64;

                append_bins(&mut bins, start - upstream_len, upstream_len, upstream_bins);
                append_bins(&mut bins, start, unscaled5_len, unscaled5_bins);
                append_bins(&mut bins, start + unscaled5_len, body_len, body_bins);
                append_bins(
                    &mut bins,
                    end - unscaled3_len,
                    unscaled3_len,
                    unscaled3_bins,
                );
                append_bins(&mut bins, end, downstream_len, downstream_bins);
            }
        }

        let window_start = bins
            .iter()
            .map(|bin| bin.start.min(bin.end))
            .min()
            .unwrap_or(start);
        let window_end = bins
            .iter()
            .map(|bin| bin.start.max(bin.end))
            .max()
            .unwrap_or(end);

        // Python: if body > 0 and body_length < bin_size, skip entire row
        let scalable_body = region_len - unscaled5_len - unscaled3_len;
        let body_too_short = options.region_body_length > 0 && scalable_body < bin_size as i64;

        Self {
            window_start,
            window_end,
            bins,
            included_intervals: None,
            body_too_short,
        }
    }
}

fn build_bins(
    record: &BedRecord,
    reference: i64,
    bin_size: u32,
    upstream_bins: usize,
    downstream_bins: usize,
) -> Vec<ReferenceBin> {
    let total_bins = upstream_bins + downstream_bins;
    let mut bins = Vec::with_capacity(total_bins);
    for bin_index in 0..total_bins {
        let (start, end) =
            bin_boundaries(bin_index, upstream_bins, bin_size, reference, record.strand);
        let beyond_region = bin_beyond_region(record, start, end);
        bins.push(ReferenceBin {
            start,
            end,
            beyond_region,
        });
    }
    bins
}

fn reference_coordinate(record: &BedRecord, reference_point: ReferencePoint) -> i64 {
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
            if start <= end {
                (start, end)
            } else {
                (end, start)
            }
        }
        _ => {
            let start = reference + offset_start;
            let end = reference + offset_end;
            if start <= end {
                (start, end)
            } else {
                (end, start)
            }
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

fn append_bins(target: &mut Vec<ScaleBin>, start: i64, length: i64, count: usize) {
    if count == 0 {
        return;
    }

    if length <= 0 {
        for _ in 0..count {
            target.push(ScaleBin {
                start,
                end: start,
                beyond_region: false,
            });
        }
        return;
    }

    let mut positions = Vec::with_capacity(count + 1);
    for idx in 0..count {
        let pos = (length * idx as i64) / count as i64;
        positions.push(pos);
    }
    positions.push(length);

    for idx in 0..count {
        let rel_start = positions[idx];
        let mut rel_end = positions[idx + 1];
        if rel_end <= rel_start {
            rel_end = rel_start + 1;
        }

        target.push(ScaleBin {
            start: start + rel_start,
            end: start + rel_end,
            beyond_region: false,
        });
    }
}

fn coordinate_from_offset(
    intervals: &[(i64, i64)],
    mut offset: i64,
    prefer_next_interval: bool,
) -> i64 {
    for (start, end) in intervals {
        let length = (end - start).max(0);
        if offset < length {
            return start + offset;
        }
        if offset == length {
            if prefer_next_interval {
                offset = 0;
                continue;
            } else {
                return *end;
            }
        }
        offset -= length;
    }
    intervals.last().map(|(_, end)| *end).unwrap_or_default()
}

pub(crate) fn intervals_total_length(intervals: &[(i64, i64)]) -> i64 {
    intervals
        .iter()
        .map(|(start, end)| (end - start).max(0))
        .sum()
}

pub(crate) fn intervals_to_bins(
    intervals: &[(i64, i64)],
    bin_count: usize,
) -> Vec<(i64, i64, bool)> {
    let mut bins = Vec::with_capacity(bin_count);
    if bin_count == 0 {
        return bins;
    }

    let total_len = intervals_total_length(intervals);
    if total_len <= 0 {
        let anchor = intervals.first().map(|(start, _)| *start).unwrap_or(0);
        for _ in 0..bin_count {
            bins.push((anchor, anchor, true));
        }
        return bins;
    }

    for idx in 0..bin_count {
        let start_offset = total_len * idx as i64 / bin_count as i64;
        let mut end_offset = total_len * (idx + 1) as i64 / bin_count as i64;
        if end_offset <= start_offset {
            end_offset = start_offset + 1;
        }

        let start = coordinate_from_offset(intervals, start_offset, true);
        let mut end = coordinate_from_offset(intervals, end_offset, false);
        if end <= start {
            end = start + 1;
        }
        bins.push((start, end, false));
    }

    bins
}

pub(crate) fn collect_window_bounds(
    intervals: &[(i64, i64)],
    mut window: (i64, i64),
) -> (i64, i64) {
    if intervals.is_empty() {
        return window;
    }

    let min_coord = intervals
        .iter()
        .map(|(start, end)| start.min(end))
        .copied()
        .min()
        .unwrap();
    let max_coord = intervals
        .iter()
        .map(|(start, end)| start.max(end))
        .copied()
        .max()
        .unwrap();

    window.0 = window.0.min(min_coord);
    window.1 = window.1.max(max_coord);
    window
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    fn build_record(strand: Strand, start: u32, end: u32) -> BedRecord {
        BedRecord {
            chrom: Arc::from("chr1"),
            start,
            end,
            name: None,
            score: None,
            score_raw: None,
            strand,
            strand_raw: None,
            extra_fields: Vec::new(),
        }
    }

    fn scale_options() -> ScaleRegionsOptions {
        ScaleRegionsOptions {
            region_body_length: 100,
            start_label: "TSS".to_string(),
            end_label: "TES".to_string(),
            upstream: 40,
            downstream: 30,
            unscaled_5_prime: 20,
            unscaled_3_prime: 30,
        }
    }

    #[test]
    fn reference_plan_positive_strand() {
        let record = build_record(Strand::Positive, 100, 200);
        let plan = ReferencePointPlan::reference_point(
            &record,
            ReferencePoint::Tss,
            10,
            2,
            2,
            false,
            false,
        );
        assert_eq!(plan.reference, 100);
        assert_eq!(plan.bins.len(), 4);
        assert_eq!(plan.bins[0].start, 80);
        assert_eq!(plan.bins[0].end, 90);
        assert_eq!(plan.bins[3].start, 110);
        assert_eq!(plan.bins[3].end, 120);
        assert!(plan.included_intervals().is_none());
    }

    #[test]
    fn reference_plan_negative_strand() {
        let record = build_record(Strand::Negative, 100, 200);
        let plan = ReferencePointPlan::reference_point(
            &record,
            ReferencePoint::Tss,
            10,
            2,
            2,
            false,
            false,
        );
        assert_eq!(plan.reference, 200);
        assert_eq!(plan.bins[0].start, 210);
        assert_eq!(plan.bins[0].end, 220);
        assert_eq!(plan.bins[2].start, 190);
        assert_eq!(plan.bins[2].end, 200);
    }

    #[test]
    fn scale_plan_positive_strand() {
        let record = build_record(Strand::Positive, 100, 200);
        let options = scale_options();
        let plan = ScaleRegionsPlan::scale_regions(&record, &options, 10, false);

        assert_eq!(plan.bins.len(), 22);
        assert_eq!(plan.window_start, 60);
        assert_eq!(plan.window_end, 230);

        // First four bins correspond to upstream region.
        for (idx, bin) in plan.bins.iter().take(4).enumerate() {
            assert_eq!(bin.start, 60 + (idx as i64 * 10));
            assert_eq!(bin.end, 60 + ((idx as i64 + 1) * 10));
        }

        // Body bins start at 120 and end before the 3' unscaled block.
        let body_start = 4 + 2; // upstream + unscaled 5 prime bins
        let first_body_bin = &plan.bins[body_start];
        assert_eq!(first_body_bin.start, 120);
        let last_body_bin = &plan.bins[body_start + 9];
        assert!(last_body_bin.end <= 170);
    }

    #[test]
    fn scale_plan_negative_strand() {
        let record = build_record(Strand::Negative, 100, 200);
        let options = scale_options();
        let plan = ScaleRegionsPlan::scale_regions(&record, &options, 10, false);

        assert_eq!(plan.bins.len(), 22);
        assert_eq!(plan.window_start, 70);
        assert_eq!(plan.window_end, 240);

        // Upstream bins (derived from downstream distance) occupy the lowest coordinates.
        for (idx, bin) in plan.bins.iter().take(3).enumerate() {
            assert_eq!(bin.start, 70 + (idx as i64 * 10));
            assert_eq!(bin.end, 70 + ((idx as i64 + 1) * 10));
        }

        // The last downstream bins extend beyond the region end using the upstream distance.
        let downstream_bins = (options.upstream / 10) as usize;
        let tail = &plan.bins[plan.bins.len() - downstream_bins..];
        for (idx, bin) in tail.iter().enumerate() {
            assert_eq!(bin.start, 200 + (idx as i64 * 10));
            assert_eq!(bin.end, 200 + ((idx as i64 + 1) * 10));
        }
    }

    fn bed12_record(strand: Strand, exons: &[(u32, u32)]) -> BedRecord {
        let start = exons.iter().map(|(s, _)| *s).min().unwrap_or(0);
        let end = exons.iter().map(|(_, e)| *e).max().unwrap_or(start);
        let block_count = exons.len().to_string();
        let block_sizes = {
            let mut values: Vec<String> = exons.iter().map(|(s, e)| (e - s).to_string()).collect();
            values.push(String::new());
            values.join(",")
        };
        let block_starts = {
            let mut values: Vec<String> =
                exons.iter().map(|(s, _)| (s - start).to_string()).collect();
            values.push(String::new());
            values.join(",")
        };

        BedRecord {
            chrom: Arc::from("chr1"),
            start,
            end,
            name: None,
            score: None,
            score_raw: None,
            strand,
            strand_raw: None,
            extra_fields: vec![
                start.to_string(),
                end.to_string(),
                "0".to_string(),
                block_count,
                block_sizes,
                block_starts,
            ],
        }
    }

    #[test]
    fn metagene_scale_regions_collapses_introns() {
        let record = bed12_record(Strand::Positive, &[(100, 150), (250, 300)]);
        let mut options = scale_options();
        options.upstream = 0;
        options.downstream = 0;
        options.unscaled_5_prime = 0;
        options.unscaled_3_prime = 0;
        options.region_body_length = 100;

        let plan = ScaleRegionsPlan::scale_regions(&record, &options, 50, true);
        assert_eq!(plan.bins.len(), 2);
        assert_eq!(plan.bins[0].start, 100);
        assert_eq!(plan.bins[0].end, 150);
        assert!(!plan.bins[0].beyond_region);
        assert_eq!(plan.bins[1].start, 250);
        assert_eq!(plan.bins[1].end, 300);
        assert!(!plan.bins[1].beyond_region);
    }

    #[test]
    fn metagene_reference_point_skips_introns() {
        let record = bed12_record(Strand::Positive, &[(100, 150), (250, 300)]);
        let plan =
            ReferencePointPlan::reference_point(&record, ReferencePoint::Tss, 50, 0, 2, true, true);
        assert_eq!(plan.bins.len(), 2);
        assert_eq!(plan.bins[0].start, 100);
        assert_eq!(plan.bins[0].end, 150);
        assert_eq!(plan.bins[1].start, 250);
        assert_eq!(plan.bins[1].end, 300);

        let included = plan.included_intervals().unwrap();
        assert_eq!(included.len(), 2);
        assert_eq!(included[0], (100, 150));
        assert_eq!(included[1], (250, 300));
    }

    #[test]
    fn metagene_center_pads_empty_left_from_right_anchor() {
        let record = bed12_record(Strand::Positive, &[(100, 101)]);
        let plan = ReferencePointPlan::reference_point(
            &record,
            ReferencePoint::Center,
            10,
            10,
            10,
            true,
            false,
        );

        assert_eq!(plan.bins.len(), 20);
        assert_eq!(plan.bins[0].start, 0);
        assert_eq!(plan.bins[0].end, 10);
        assert_eq!(plan.bins[9].start, 90);
        assert_eq!(plan.bins[9].end, 100);
        assert_eq!(plan.bins[10].start, 100);
        assert_eq!(plan.bins[10].end, 110);
        assert_eq!(plan.bins[19].start, 190);
        assert_eq!(plan.bins[19].end, 200);

        let included = plan.included_intervals().unwrap();
        assert_eq!(included[0], (0, 100));
        assert_eq!(included[1], (100, 101));
        assert_eq!(included[2], (101, 200));
    }

    // ── Issue 2a: body_too_short flag ────────────────────────────────────

    #[test]
    fn body_too_short_when_scalable_body_less_than_bin_size() {
        // region: [100, 115), length = 15
        // unscaled5 = 5, unscaled3 = 5
        // scalable body = 15 - 5 - 5 = 5, bin_size = 10 → too short
        let record = build_record(Strand::Positive, 100, 115);
        let options = ScaleRegionsOptions {
            region_body_length: 50,
            start_label: "TSS".to_string(),
            end_label: "TES".to_string(),
            upstream: 0,
            downstream: 0,
            unscaled_5_prime: 5,
            unscaled_3_prime: 5,
        };
        let plan = ScaleRegionsPlan::scale_regions(&record, &options, 10, false);
        assert!(
            plan.body_too_short,
            "expected body_too_short when scalable body (5) < bin_size (10)"
        );
    }

    #[test]
    fn body_not_too_short_when_scalable_body_equals_bin_size() {
        // region: [100, 120), length = 20
        // unscaled5 = 5, unscaled3 = 5
        // scalable body = 20 - 5 - 5 = 10, bin_size = 10 → NOT too short
        let record = build_record(Strand::Positive, 100, 120);
        let options = ScaleRegionsOptions {
            region_body_length: 50,
            start_label: "TSS".to_string(),
            end_label: "TES".to_string(),
            upstream: 0,
            downstream: 0,
            unscaled_5_prime: 5,
            unscaled_3_prime: 5,
        };
        let plan = ScaleRegionsPlan::scale_regions(&record, &options, 10, false);
        assert!(
            !plan.body_too_short,
            "body should NOT be too short when scalable body (10) == bin_size (10)"
        );
    }

    #[test]
    fn body_not_too_short_when_region_body_length_is_zero() {
        // region_body_length = 0 → body parameter not set, should NOT
        // trigger the short-circuit even if region is tiny
        let record = build_record(Strand::Positive, 100, 101);
        let options = ScaleRegionsOptions {
            region_body_length: 0,
            start_label: "TSS".to_string(),
            end_label: "TES".to_string(),
            upstream: 0,
            downstream: 0,
            unscaled_5_prime: 0,
            unscaled_3_prime: 0,
        };
        let plan = ScaleRegionsPlan::scale_regions(&record, &options, 10, false);
        assert!(
            !plan.body_too_short,
            "body_too_short should be false when region_body_length == 0"
        );
    }

    #[test]
    fn body_too_short_negative_strand() {
        // Same logic should apply to negative strand
        let record = build_record(Strand::Negative, 100, 115);
        let options = ScaleRegionsOptions {
            region_body_length: 50,
            start_label: "TSS".to_string(),
            end_label: "TES".to_string(),
            upstream: 0,
            downstream: 0,
            unscaled_5_prime: 5,
            unscaled_3_prime: 5,
        };
        let plan = ScaleRegionsPlan::scale_regions(&record, &options, 10, false);
        assert!(
            plan.body_too_short,
            "body_too_short should apply to negative strand as well"
        );
    }

    #[test]
    fn body_too_short_no_unscaled_regions() {
        // region: [100, 105), length = 5, no unscaled regions
        // scalable body = 5, bin_size = 10 → too short
        let record = build_record(Strand::Positive, 100, 105);
        let options = ScaleRegionsOptions {
            region_body_length: 50,
            start_label: "TSS".to_string(),
            end_label: "TES".to_string(),
            upstream: 20,
            downstream: 20,
            unscaled_5_prime: 0,
            unscaled_3_prime: 0,
        };
        let plan = ScaleRegionsPlan::scale_regions(&record, &options, 10, false);
        assert!(
            plan.body_too_short,
            "body_too_short should be true when region (5bp) < bin_size (10)"
        );
    }

    #[test]
    fn body_too_short_metagene_exon_total_less_than_bin_size() {
        // Two exons: [100,103) and [200,202) → total = 5
        // unscaled5 = 0, unscaled3 = 0
        // scalable body = 5, bin_size = 10 → too short
        let record = bed12_record(Strand::Positive, &[(100, 103), (200, 202)]);
        let options = ScaleRegionsOptions {
            region_body_length: 50,
            start_label: "TSS".to_string(),
            end_label: "TES".to_string(),
            upstream: 0,
            downstream: 0,
            unscaled_5_prime: 0,
            unscaled_3_prime: 0,
        };
        let plan = ScaleRegionsPlan::scale_regions(&record, &options, 10, true);
        assert!(
            plan.body_too_short,
            "metagene: body_too_short when exon total (5) < bin_size (10)"
        );
    }

    #[test]
    fn body_not_too_short_metagene_exon_total_sufficient() {
        // Two exons: [100,150) and [200,250) → total = 100
        // unscaled5 = 10, unscaled3 = 10
        // scalable body = 100 - 10 - 10 = 80, bin_size = 10 → NOT too short
        let record = bed12_record(Strand::Positive, &[(100, 150), (200, 250)]);
        let options = ScaleRegionsOptions {
            region_body_length: 50,
            start_label: "TSS".to_string(),
            end_label: "TES".to_string(),
            upstream: 0,
            downstream: 0,
            unscaled_5_prime: 10,
            unscaled_3_prime: 10,
        };
        let plan = ScaleRegionsPlan::scale_regions(&record, &options, 10, true);
        assert!(
            !plan.body_too_short,
            "metagene: body should NOT be too short when exon total (80) >= bin_size (10)"
        );
    }

    #[test]
    fn reference_point_plan_body_too_short_default_false() {
        // ReferencePointPlan uses the default trait impl which returns false
        use crate::pipeline::core::RegionPlan;
        let record = build_record(Strand::Positive, 100, 200);
        let plan = ReferencePointPlan::reference_point(
            &record,
            ReferencePoint::Tss,
            10,
            2,
            2,
            false,
            false,
        );
        assert!(
            !plan.body_too_short(),
            "ReferencePointPlan should always return false for body_too_short"
        );
    }
}
