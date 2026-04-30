use super::{
    ReferenceBin, ScaleBin, ScaleRegionsPlan, collect_window_bounds, intervals_to_bins,
    intervals_total_length,
};
use crate::config::{ReferencePoint, ScaleRegionsOptions};
use crate::io::{BedRecord, Strand};

pub fn reference_bins(
    record: &BedRecord,
    reference_point: ReferencePoint,
    bin_size: u32,
    upstream_bins: usize,
    downstream_bins: usize,
    nan_after_end: bool,
) -> Option<(Vec<ReferenceBin>, Vec<(i64, i64)>)> {
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
    let mut included = Vec::new();

    match reference_point {
        ReferencePoint::Tss => build_tss(
            record,
            &exons,
            bin_size,
            upstream_bins,
            downstream_bins,
            nan_after_end,
            &mut bins,
            &mut included,
        ),
        ReferencePoint::Tes => build_tes(
            record,
            &exons,
            bin_size,
            upstream_bins,
            downstream_bins,
            nan_after_end,
            &mut bins,
            &mut included,
        ),
        ReferencePoint::Center => build_center(
            record,
            &exons,
            bin_size,
            upstream_bins,
            downstream_bins,
            nan_after_end,
            &mut bins,
            &mut included,
        ),
    }

    if bins.len() == upstream_bins + downstream_bins {
        Some((bins, included))
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

    // Python: if body > 0 and body_length < bin_size, skip entire row.
    // For metagene, body_length = total exon length - unscaled regions.
    let exon_total = intervals_total_length(&exons);
    let scalable_body =
        exon_total - options.unscaled_5_prime as i64 - options.unscaled_3_prime as i64;
    let body_too_short = options.region_body_length > 0 && scalable_body < bin_size as i64;

    Some(ScaleRegionsPlan {
        window_start: window.0,
        window_end: window.1,
        bins,
        included_intervals: Some(included_intervals),
        body_too_short,
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
    included: &mut Vec<(i64, i64)>,
) {
    let feature_start = exons.first().map(|(start, _)| *start).unwrap_or(0);
    let feature_end = exons.last().map(|(_, end)| *end).unwrap_or(0);

    match record.strand {
        Strand::Negative => {
            // Python builds zones in ascending genomic order for negative-
            // strand TSS, then reverses the whole coverage array.  Match
            // that layout so postprocess_row can do the same reversal.
            //
            // Python zones: [(exon_right_chopped, downstream_bins),
            //                (feature_end_extension, upstream_bins)]
            let (mut downstream_intervals, pad) =
                take_from_end(exons, downstream_bins as u32 * bin_size);
            if pad > 0 && !downstream_intervals.is_empty() && !nan_after_end {
                let start = downstream_intervals
                    .first()
                    .map(|(start, _)| *start)
                    .unwrap_or(feature_start);
                downstream_intervals.insert(0, (start - pad as i64, start));
            }
            included.extend_from_slice(&downstream_intervals);
            append_reference_bins(
                bins,
                &downstream_intervals,
                downstream_bins,
                bin_size,
                nan_after_end,
            );

            let upstream_intervals = if upstream_bins > 0 {
                vec![(
                    feature_end,
                    feature_end + (bin_size as usize * upstream_bins) as i64,
                )]
            } else {
                Vec::new()
            };
            included.extend_from_slice(&upstream_intervals);
            append_reference_bins(
                bins,
                &upstream_intervals,
                upstream_bins,
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
            included.extend_from_slice(&upstream_intervals);
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
            included.extend_from_slice(&downstream_intervals);
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
    included: &mut Vec<(i64, i64)>,
) {
    let feature_start = exons.first().map(|(start, _)| *start).unwrap_or(0);
    let feature_end = exons.last().map(|(_, end)| *end).unwrap_or(0);

    match record.strand {
        Strand::Negative => {
            // Python builds zones in ascending genomic order for negative-
            // strand TES, then reverses the whole coverage array.
            // TES for negative strand = biological start = genomic left (feature_start).
            //
            // Python zones: [(feature_start - downstream, feature_start),
            //                (left_exon_portion by upstream)]
            let upstream_intervals = if downstream_bins > 0 {
                vec![(
                    feature_start - (bin_size as usize * downstream_bins) as i64,
                    feature_start,
                )]
            } else {
                Vec::new()
            };
            included.extend_from_slice(&upstream_intervals);
            append_reference_bins(
                bins,
                &upstream_intervals,
                downstream_bins,
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
            included.extend_from_slice(&downstream_intervals);
            append_reference_bins(
                bins,
                &downstream_intervals,
                upstream_bins,
                bin_size,
                nan_after_end,
            );
        }
        _ => {
            // Python TES positive: upstream = chopRegions(exons, right=upstream)
            // (right portion of exons), downstream = extension beyond feature_end.
            // zones = [(upstream, a), (downstream, e)] in ascending genomic order.
            let (mut upstream_intervals, pad) =
                take_from_end(exons, upstream_bins as u32 * bin_size);
            if pad > 0 && !upstream_intervals.is_empty() && !nan_after_end {
                let start = upstream_intervals
                    .first()
                    .map(|(start, _)| *start)
                    .unwrap_or(feature_start);
                upstream_intervals.insert(0, (start - pad as i64, start));
            }
            included.extend_from_slice(&upstream_intervals);
            append_reference_bins(
                bins,
                &upstream_intervals,
                upstream_bins,
                bin_size,
                nan_after_end,
            );

            let downstream_intervals = if downstream_bins > 0 {
                vec![(
                    feature_end,
                    feature_end + (bin_size as usize * downstream_bins) as i64,
                )]
            } else {
                Vec::new()
            };
            included.extend_from_slice(&downstream_intervals);
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
    included: &mut Vec<(i64, i64)>,
) {
    // Python for negative strand passes (left=downstream, right=upstream) to
    // chopRegionsFromMiddle, whereas positive uses (left=upstream, right=downstream).
    // Both produce zones in ascending genomic order; negative is later reversed.
    let (left_param, right_param) = if record.strand == Strand::Negative {
        (
            (downstream_bins as u32) * bin_size,
            (upstream_bins as u32) * bin_size,
        )
    } else {
        (
            (upstream_bins as u32) * bin_size,
            (downstream_bins as u32) * bin_size,
        )
    };
    let (mut left, mut right, pad_left, pad_right) =
        chop_regions_from_middle(exons, left_param, right_param);

    // Determine bin counts for left and right halves.
    // For negative strand, left holds downstream_bins and right holds upstream_bins.
    let (left_bins, right_bins) = if record.strand == Strand::Negative {
        (downstream_bins, upstream_bins)
    } else {
        (upstream_bins, downstream_bins)
    };

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

    included.extend_from_slice(&left);
    append_reference_bins(bins, &left, left_bins, bin_size, nan_after_end);
    included.extend_from_slice(&right);
    append_reference_bins(bins, &right, right_bins, bin_size, nan_after_end);
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
        // Do NOT clear left_bins — keep the exon intervals for coverage
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
        // Do NOT clear right_bins — keep the exon intervals for coverage
    }

    (left_bins, right_bins, pad_left, pad_right)
}
