use crate::config::{ReferencePoint, ScaleRegionsOptions};
use crate::io::BedRecord;
use crate::io::Strand;
use crate::pipeline::core::{RegionPlan, SignalBin};

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
        let bins = if keep_exons {
            metagene::reference_bins(
                record,
                reference_point,
                bin_size,
                upstream_bins,
                downstream_bins,
                nan_after_end,
            )
            .unwrap_or_else(|| {
                build_bins(record, reference, bin_size, upstream_bins, downstream_bins)
            })
        } else {
            build_bins(record, reference, bin_size, upstream_bins, downstream_bins)
        };

        let window_start = bins.iter().map(|bin| bin.start).min().unwrap_or(reference);
        let window_end = bins.iter().map(|bin| bin.end).max().unwrap_or(reference);

        Self {
            reference,
            window_start,
            window_end,
            bins,
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

        Self {
            window_start,
            window_end,
            bins,
            included_intervals: None,
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

fn intervals_total_length(intervals: &[(i64, i64)]) -> i64 {
    intervals
        .iter()
        .map(|(start, end)| (end - start).max(0))
        .sum()
}

fn intervals_to_bins(intervals: &[(i64, i64)], bin_count: usize) -> Vec<(i64, i64, bool)> {
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

fn collect_window_bounds(intervals: &[(i64, i64)], mut window: (i64, i64)) -> (i64, i64) {
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

mod metagene {
    use super::{
        ReferenceBin, ReferencePoint, ScaleBin, ScaleRegionsPlan, collect_window_bounds,
        intervals_to_bins, intervals_total_length,
    };
    use crate::config::ScaleRegionsOptions;
    use crate::io::{BedRecord, Strand};

    pub fn reference_bins(
        record: &BedRecord,
        reference_point: ReferencePoint,
        bin_size: u32,
        upstream_bins: usize,
        downstream_bins: usize,
        nan_after_end: bool,
    ) -> Option<Vec<ReferenceBin>> {
        let mut exons = record
            .exons()?
            .into_iter()
            .map(|(start, end)| (start as i64, end as i64))
            .collect::<Vec<_>>();
        if exons.is_empty() {
            return None;
        }
        exons.sort_by_key(|(start, _)| *start);

        let mut bins = Vec::with_capacity(upstream_bins + downstream_bins);

        match reference_point {
            ReferencePoint::Tss => build_tss(
                record,
                &exons,
                bin_size,
                upstream_bins,
                downstream_bins,
                nan_after_end,
                &mut bins,
            ),
            ReferencePoint::Tes => build_tes(
                record,
                &exons,
                bin_size,
                upstream_bins,
                downstream_bins,
                nan_after_end,
                &mut bins,
            ),
            ReferencePoint::Center => build_center(
                record,
                &exons,
                bin_size,
                upstream_bins,
                downstream_bins,
                nan_after_end,
                &mut bins,
            ),
        }

        if bins.len() == upstream_bins + downstream_bins {
            Some(bins)
        } else {
            None
        }
    }

    pub fn scale_bins(
        record: &BedRecord,
        options: &ScaleRegionsOptions,
        bin_size: u32,
    ) -> Option<ScaleRegionsPlan> {
        let mut exons = record
            .exons()?
            .into_iter()
            .map(|(start, end)| (start as i64, end as i64))
            .collect::<Vec<_>>();
        if exons.is_empty() {
            return None;
        }
        exons.sort_by_key(|(start, _)| *start);

        let upstream_bins = (options.upstream / bin_size) as usize;
        let downstream_bins = (options.downstream / bin_size) as usize;
        let unscaled5_bins = (options.unscaled_5_prime / bin_size) as usize;
        let unscaled3_bins = (options.unscaled_3_prime / bin_size) as usize;
        let body_bins = (options.region_body_length / bin_size) as usize;

        let mut bins = Vec::with_capacity(
            upstream_bins + downstream_bins + unscaled5_bins + unscaled3_bins + body_bins,
        );
        let mut included_intervals = Vec::new();
        let mut window = (i64::MAX, i64::MIN);

        match record.strand {
            Strand::Negative => build_scale_negative(
                &exons,
                options,
                upstream_bins,
                downstream_bins,
                unscaled5_bins,
                unscaled3_bins,
                body_bins,
                &mut bins,
                &mut window,
                &mut included_intervals,
            ),
            _ => build_scale_positive(
                &exons,
                options,
                upstream_bins,
                downstream_bins,
                unscaled5_bins,
                unscaled3_bins,
                body_bins,
                &mut bins,
                &mut window,
                &mut included_intervals,
            ),
        }

        if window.0 == i64::MAX {
            window.0 = record.start as i64;
        }
        if window.1 == i64::MIN {
            window.1 = record.end as i64;
        }

        Some(ScaleRegionsPlan {
            window_start: window.0,
            window_end: window.1,
            bins,
            included_intervals: Some(included_intervals),
        })
    }

    fn build_tss(
        record: &BedRecord,
        exons: &[(i64, i64)],
        bin_size: u32,
        upstream_bins: usize,
        downstream_bins: usize,
        nan_after_end: bool,
        bins: &mut Vec<ReferenceBin>,
    ) {
        let feature_start = exons.first().map(|(start, _)| *start).unwrap_or(0);
        let feature_end = exons.last().map(|(_, end)| *end).unwrap_or(0);

        match record.strand {
            Strand::Negative => {
                let upstream_intervals = if upstream_bins > 0 {
                    vec![(
                        feature_end,
                        feature_end + (bin_size as usize * upstream_bins) as i64,
                    )]
                } else {
                    Vec::new()
                };
                append_reference_bins(
                    bins,
                    &upstream_intervals,
                    upstream_bins,
                    bin_size,
                    nan_after_end,
                );

                let (mut downstream_intervals, pad) =
                    take_from_end(exons, downstream_bins as u32 * bin_size);
                if pad > 0 && !downstream_intervals.is_empty() && !nan_after_end {
                    let start = downstream_intervals
                        .first()
                        .map(|(start, _)| *start)
                        .unwrap_or(feature_start);
                    downstream_intervals.insert(0, (start - pad as i64, start));
                }
                append_reference_bins(
                    bins,
                    &downstream_intervals,
                    downstream_bins,
                    bin_size,
                    nan_after_end,
                );
            }
            _ => {
                let upstream_intervals = if upstream_bins > 0 {
                    vec![(
                        feature_start - (bin_size as usize * upstream_bins) as i64,
                        feature_start,
                    )]
                } else {
                    Vec::new()
                };
                append_reference_bins(
                    bins,
                    &upstream_intervals,
                    upstream_bins,
                    bin_size,
                    nan_after_end,
                );

                let (mut downstream_intervals, pad) =
                    take_from_start(exons, downstream_bins as u32 * bin_size);
                if pad > 0 && !downstream_intervals.is_empty() && !nan_after_end {
                    let end = downstream_intervals
                        .last()
                        .map(|(_, end)| *end)
                        .unwrap_or(feature_end);
                    downstream_intervals.push((end, end + pad as i64));
                }
                append_reference_bins(
                    bins,
                    &downstream_intervals,
                    downstream_bins,
                    bin_size,
                    nan_after_end,
                );
            }
        }
    }

    fn build_tes(
        record: &BedRecord,
        exons: &[(i64, i64)],
        bin_size: u32,
        upstream_bins: usize,
        downstream_bins: usize,
        nan_after_end: bool,
        bins: &mut Vec<ReferenceBin>,
    ) {
        let feature_start = exons.first().map(|(start, _)| *start).unwrap_or(0);
        let feature_end = exons.last().map(|(_, end)| *end).unwrap_or(0);

        match record.strand {
            Strand::Negative => {
                let upstream_intervals = if upstream_bins > 0 {
                    vec![(
                        feature_end,
                        feature_end + (bin_size as usize * upstream_bins) as i64,
                    )]
                } else {
                    Vec::new()
                };
                append_reference_bins(
                    bins,
                    &upstream_intervals,
                    upstream_bins,
                    bin_size,
                    nan_after_end,
                );

                let (mut downstream_intervals, pad) =
                    take_from_end(exons, upstream_bins as u32 * bin_size);
                if pad > 0 && !downstream_intervals.is_empty() && !nan_after_end {
                    let start = downstream_intervals
                        .first()
                        .map(|(start, _)| *start)
                        .unwrap_or(feature_start);
                    downstream_intervals.insert(0, (start - pad as i64, start));
                }
                append_reference_bins(
                    bins,
                    &downstream_intervals,
                    downstream_bins,
                    bin_size,
                    nan_after_end,
                );
            }
            _ => {
                let upstream_intervals = if upstream_bins > 0 {
                    vec![(
                        feature_start - (bin_size as usize * upstream_bins) as i64,
                        feature_start,
                    )]
                } else {
                    Vec::new()
                };
                append_reference_bins(
                    bins,
                    &upstream_intervals,
                    upstream_bins,
                    bin_size,
                    nan_after_end,
                );

                let (mut downstream_intervals, pad) =
                    take_from_start(exons, upstream_bins as u32 * bin_size);
                if pad > 0 && !downstream_intervals.is_empty() && !nan_after_end {
                    let end = downstream_intervals
                        .last()
                        .map(|(_, end)| *end)
                        .unwrap_or(feature_end);
                    downstream_intervals.push((end, end + pad as i64));
                }
                append_reference_bins(
                    bins,
                    &downstream_intervals,
                    downstream_bins,
                    bin_size,
                    nan_after_end,
                );
            }
        }
    }

    fn build_center(
        record: &BedRecord,
        exons: &[(i64, i64)],
        bin_size: u32,
        upstream_bins: usize,
        downstream_bins: usize,
        nan_after_end: bool,
        bins: &mut Vec<ReferenceBin>,
    ) {
        let (mut left, mut right, pad_left, pad_right) = chop_regions_from_middle(
            exons,
            (upstream_bins as u32) * bin_size,
            (downstream_bins as u32) * bin_size,
        );

        if record.strand == Strand::Negative {
            std::mem::swap(&mut left, &mut right);
        }

        if pad_left > 0 && !left.is_empty() && !nan_after_end {
            let start = left
                .first()
                .map(|(start, _)| *start)
                .unwrap_or(right.first().map(|(start, _)| *start).unwrap_or(0));
            left.insert(0, (start - pad_left as i64, start));
        }
        if pad_right > 0 && !right.is_empty() && !nan_after_end {
            let end = right
                .last()
                .map(|(_, end)| *end)
                .unwrap_or(left.last().map(|(_, end)| *end).unwrap_or(0));
            right.push((end, end + pad_right as i64));
        }

        append_reference_bins(bins, &left, upstream_bins, bin_size, nan_after_end);
        append_reference_bins(bins, &right, downstream_bins, bin_size, nan_after_end);
    }

    fn append_reference_bins(
        target: &mut Vec<ReferenceBin>,
        intervals: &[(i64, i64)],
        expected_bins: usize,
        bin_size: u32,
        nan_after_end: bool,
    ) {
        if expected_bins == 0 {
            return;
        }

        let available_bins = (intervals_total_length(intervals) / bin_size as i64).max(0) as usize;
        let bins_to_use = available_bins.min(expected_bins);
        let bins = intervals_to_bins(intervals, bins_to_use);
        for (start, end, beyond_region) in bins {
            target.push(ReferenceBin {
                start,
                end,
                beyond_region,
            });
        }

        let missing = expected_bins.saturating_sub(bins_to_use);
        if missing > 0 {
            let mut anchor = target
                .last()
                .map(|bin| bin.end)
                .or_else(|| intervals.first().map(|(start, _)| *start))
                .unwrap_or(0);
            for _ in 0..missing {
                if nan_after_end {
                    target.push(ReferenceBin {
                        start: anchor,
                        end: anchor,
                        beyond_region: true,
                    });
                } else {
                    target.push(ReferenceBin {
                        start: anchor,
                        end: anchor + bin_size as i64,
                        beyond_region: false,
                    });
                    anchor += bin_size as i64;
                }
            }
        }
    }

    fn build_scale_positive(
        exons: &[(i64, i64)],
        options: &ScaleRegionsOptions,
        upstream_bins: usize,
        downstream_bins: usize,
        unscaled5_bins: usize,
        unscaled3_bins: usize,
        body_bins: usize,
        bins: &mut Vec<ScaleBin>,
        window: &mut (i64, i64),
        included: &mut Vec<(i64, i64)>,
    ) {
        let feature_start = exons.first().map(|(start, _)| *start).unwrap_or(0);
        let feature_end = exons.last().map(|(_, end)| *end).unwrap_or(0);

        let upstream_intervals = if upstream_bins > 0 {
            vec![(feature_start - options.upstream as i64, feature_start)]
        } else {
            Vec::new()
        };
        included.extend_from_slice(&upstream_intervals);
        *window = collect_window_bounds(&upstream_intervals, *window);
        append_scale_bins(bins, &upstream_intervals, upstream_bins);

        let (unscaled5, body, unscaled3, _, _) =
            chop_regions(exons, options.unscaled_5_prime, options.unscaled_3_prime);

        included.extend_from_slice(&unscaled5);
        *window = collect_window_bounds(&unscaled5, *window);
        append_scale_bins(bins, &unscaled5, unscaled5_bins);

        included.extend_from_slice(&body);
        *window = collect_window_bounds(&body, *window);
        append_scale_bins(bins, &body, body_bins);

        included.extend_from_slice(&unscaled3);
        *window = collect_window_bounds(&unscaled3, *window);
        append_scale_bins(bins, &unscaled3, unscaled3_bins);

        let downstream_intervals = if downstream_bins > 0 {
            vec![(feature_end, feature_end + options.downstream as i64)]
        } else {
            Vec::new()
        };
        included.extend_from_slice(&downstream_intervals);
        *window = collect_window_bounds(&downstream_intervals, *window);
        append_scale_bins(bins, &downstream_intervals, downstream_bins);
    }

    fn build_scale_negative(
        exons: &[(i64, i64)],
        options: &ScaleRegionsOptions,
        upstream_bins: usize,
        downstream_bins: usize,
        unscaled5_bins: usize,
        unscaled3_bins: usize,
        body_bins: usize,
        bins: &mut Vec<ScaleBin>,
        window: &mut (i64, i64),
        included: &mut Vec<(i64, i64)>,
    ) {
        let feature_start = exons.first().map(|(start, _)| *start).unwrap_or(0);
        let feature_end = exons.last().map(|(_, end)| *end).unwrap_or(0);

        let upstream_intervals = if downstream_bins > 0 {
            vec![(feature_start - options.downstream as i64, feature_start)]
        } else {
            Vec::new()
        };
        included.extend_from_slice(&upstream_intervals);
        *window = collect_window_bounds(&upstream_intervals, *window);
        append_scale_bins(bins, &upstream_intervals, downstream_bins);

        let (unscaled5, body, unscaled3, _, _) =
            chop_regions(exons, options.unscaled_3_prime, options.unscaled_5_prime);

        // For - strand, the "left" chop (unscaled5) contains biological 3' unscaled region
        // and should use unscaled3_bins. The "right" chop (unscaled3) contains biological
        // 5' unscaled region and should use unscaled5_bins.
        // Python uses: zones = [(upstream, a), (unscaled5prime, b), (body, c), (unscaled3prime, d), (downstream, e)]
        // where b = unscaled_3_prime // bin_size and d = unscaled_5_prime // bin_size
        included.extend_from_slice(&unscaled5);
        *window = collect_window_bounds(&unscaled5, *window);
        append_scale_bins(bins, &unscaled5, unscaled3_bins);

        included.extend_from_slice(&body);
        *window = collect_window_bounds(&body, *window);
        append_scale_bins(bins, &body, body_bins);

        included.extend_from_slice(&unscaled3);
        *window = collect_window_bounds(&unscaled3, *window);
        append_scale_bins(bins, &unscaled3, unscaled5_bins);

        let downstream_intervals = if upstream_bins > 0 {
            vec![(feature_end, feature_end + options.upstream as i64)]
        } else {
            Vec::new()
        };
        included.extend_from_slice(&downstream_intervals);
        *window = collect_window_bounds(&downstream_intervals, *window);
        append_scale_bins(bins, &downstream_intervals, upstream_bins);
    }

    fn append_scale_bins(target: &mut Vec<ScaleBin>, intervals: &[(i64, i64)], count: usize) {
        if count == 0 {
            return;
        }
        let bins = intervals_to_bins(intervals, count);
        for (start, end, beyond_region) in bins {
            target.push(ScaleBin {
                start,
                end,
                beyond_region,
            });
        }
    }

    fn take_from_start(intervals: &[(i64, i64)], length: u32) -> (Vec<(i64, i64)>, u32) {
        let mut remaining = length as i64;
        let mut output = Vec::new();

        for (start, end) in intervals {
            if remaining <= 0 {
                break;
            }
            let width = (end - start).max(0);
            if width <= remaining {
                output.push((*start, *end));
                remaining -= width;
            } else {
                output.push((*start, start + remaining));
                remaining = 0;
            }
        }

        let pad = remaining.max(0) as u32;
        (output, pad)
    }

    fn take_from_end(intervals: &[(i64, i64)], length: u32) -> (Vec<(i64, i64)>, u32) {
        let mut remaining = length as i64;
        let mut output = Vec::new();

        for (start, end) in intervals.iter().rev() {
            if remaining <= 0 {
                break;
            }
            let width = (end - start).max(0);
            if width <= remaining {
                output.push((*start, *end));
                remaining -= width;
            } else {
                output.push((end - remaining, *end));
                remaining = 0;
            }
        }
        output.reverse();

        let pad = remaining.max(0) as u32;
        (output, pad)
    }

    fn chop_regions(
        exons: &[(i64, i64)],
        left: u32,
        right: u32,
    ) -> (Vec<(i64, i64)>, Vec<(i64, i64)>, Vec<(i64, i64)>, u32, u32) {
        let mut remaining_left = left as i64;
        let mut left_bins = Vec::new();
        let mut body = exons.to_vec();

        while !body.is_empty() && remaining_left > 0 {
            let (start, end) = body[0];
            let width = (end - start).max(0);
            if width <= remaining_left {
                left_bins.push((start, end));
                body.remove(0);
                remaining_left -= width;
            } else {
                left_bins.push((start, start + remaining_left));
                body[0] = (start + remaining_left, end);
                remaining_left = 0;
            }
        }
        let pad_left = remaining_left.max(0) as u32;

        let mut remaining_right = right as i64;
        let mut right_bins = Vec::new();
        while !body.is_empty() && remaining_right > 0 {
            let idx = body.len() - 1;
            let (start, end) = body[idx];
            let width = (end - start).max(0);
            if width <= remaining_right {
                right_bins.push((start, end));
                body.remove(idx);
                remaining_right -= width;
            } else {
                right_bins.push((end - remaining_right, end));
                body[idx] = (start, end - remaining_right);
                remaining_right = 0;
            }
        }
        right_bins.reverse();
        let pad_right = remaining_right.max(0) as u32;

        (left_bins, body, right_bins, pad_left, pad_right)
    }

    fn chop_regions_from_middle(
        exons: &[(i64, i64)],
        left: u32,
        right: u32,
    ) -> (Vec<(i64, i64)>, Vec<(i64, i64)>, u32, u32) {
        let total = intervals_total_length(exons);
        if total <= 0 {
            return (Vec::new(), Vec::new(), left, right);
        }

        let middle = total / 2;
        let mut cumulative = 0;
        let mut left_bins = Vec::new();
        let mut right_bins = Vec::new();

        for (start, end) in exons {
            let width = (end - start).max(0);
            if cumulative >= middle {
                right_bins.push((*start, *end));
            } else if cumulative + width < middle {
                left_bins.push((*start, *end));
            } else {
                let split = start + (middle - cumulative);
                if *start < split {
                    left_bins.push((*start, split));
                }
                if split < *end {
                    right_bins.push((split, *end));
                }
            }
            cumulative += width;
        }

        let mut pad_left = 0;
        let mut pad_right = 0;

        let left_sum = intervals_total_length(&left_bins);
        if left_sum > left as i64 {
            let mut acc = 0;
            for idx in (0..left_bins.len()).rev() {
                let (start, end) = left_bins[idx];
                let width = (end - start).max(0);
                if acc + width > left as i64 {
                    left_bins[idx].0 = end - (left as i64 - acc);
                    left_bins = left_bins[idx..].to_vec();
                    break;
                }
                acc += width;
                if acc == left as i64 {
                    left_bins = left_bins[idx..].to_vec();
                    break;
                }
            }
        } else if left_sum < left as i64 {
            pad_left = left - left_sum as u32;
            left_bins.clear();
        }

        let right_sum = intervals_total_length(&right_bins);
        if right_sum > right as i64 {
            let mut acc = 0;
            for idx in 0..right_bins.len() {
                let (start, end) = right_bins[idx];
                let width = (end - start).max(0);
                if acc + width > right as i64 {
                    right_bins[idx].1 = start + (right as i64 - acc);
                    right_bins.truncate(idx + 1);
                    break;
                }
                acc += width;
                if acc == right as i64 {
                    right_bins.truncate(idx + 1);
                    break;
                }
            }
        } else if right_sum < right as i64 {
            pad_right = right - right_sum as u32;
            right_bins.clear();
        }

        (left_bins, right_bins, pad_left, pad_right)
    }
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
            chrom: "chr1".to_string(),
            start,
            end,
            name: None,
            score: None,
            strand,
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
    }
}
