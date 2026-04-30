use std::sync::Arc;

use anyhow::{Context, Result};

use crate::config::{AverageTypeBins, GeneralOptions};
use crate::io::BedRecord;
use crate::pipeline::matrix::MatrixRow;

use super::aggregation::{
    aggregate_slice, direct_mean_bins, direct_sum_bins, index_from_coordinate,
};
use super::coalesce::CoalescedBatch;
use super::samples::Sample;
use super::traits::{PipelineMode, RegionPlan, SignalBin};

thread_local! {
    static COVERAGE_POOL: std::cell::RefCell<Vec<Vec<f64>>> =
        std::cell::RefCell::new(Vec::new());
}

fn take_coverage_buffers(
    sample_count: usize,
    window_len: usize,
    default_fill: f64,
) -> Vec<Vec<f64>> {
    COVERAGE_POOL.with(|pool| {
        let mut bufs = pool.borrow_mut();
        bufs.resize_with(sample_count, Vec::new);
        for buf in bufs.iter_mut() {
            buf.clear();
            buf.resize(window_len, default_fill);
        }
        std::mem::take(&mut *bufs)
    })
}

fn return_coverage_buffers(bufs: Vec<Vec<f64>>) {
    COVERAGE_POOL.with(|pool| {
        *pool.borrow_mut() = bufs;
    });
}

fn clamp_coordinate(value: i64, chrom_length: u32) -> u32 {
    value
        .max(0)
        .min(chrom_length as i64)
        .try_into()
        .unwrap_or(0)
}

fn should_skip_row_flat(values: &[f64], general: &GeneralOptions) -> bool {
    if general.skip_zeros {
        // Python computes np.mean(row) across all values; if the mean
        // is 0 or there are no non-NaN values the row is skipped.
        let mut sum = 0.0f64;
        let mut count = 0u64;
        for &value in values {
            if value.is_nan() {
                continue;
            }
            sum += value;
            count += 1;
        }
        if count == 0 || sum / (count as f64) == 0.0 {
            return true;
        }
    }

    if let Some(min_threshold) = general.min_threshold {
        if values
            .iter()
            .filter(|value| !value.is_nan())
            .any(|value| *value <= min_threshold)
        {
            return true;
        }
    }

    if let Some(max_threshold) = general.max_threshold {
        if values
            .iter()
            .filter(|value| !value.is_nan())
            .any(|value| *value >= max_threshold)
        {
            return true;
        }
    }

    false
}

fn compute_sample_bins<P: RegionPlan>(
    sample: &mut Sample,
    record: &BedRecord,
    plan: &P,
    general: &GeneralOptions,
    nan_after_end: bool,
) -> Result<Vec<f64>> {
    let sample_path = sample.path().to_path_buf();
    let bin_count = plan.bins().len();
    let chrom_length = match sample.chrom_length(&record.chrom) {
        Some(length) => length,
        None => {
            return Ok(vec![f64::NAN; bin_count]);
        }
    };

    let window_span = plan.window_end() - plan.window_start();
    if window_span <= 0 {
        return Ok(vec![f64::NAN; bin_count]);
    }

    let window_len = usize::try_from(window_span).expect("region plan window span exceeds usize");
    let default_fill = if general.missing_data_as_zero {
        0.0f64
    } else {
        f64::NAN
    };

    thread_local! {
        static COVERAGE_BUF: std::cell::RefCell<Vec<f64>> = std::cell::RefCell::new(Vec::new());
    }

    let bins = COVERAGE_BUF.with(|cell| {
        let mut buf = cell.borrow_mut();
        buf.clear();
        // In metagene mode, intron positions must stay NaN even with
        // missingDataAsZero — only exonic gaps get 0.0.
        let fill = if plan.included_intervals().is_some() {
            f64::NAN
        } else {
            default_fill
        };
        buf.resize(window_len, fill);

        if let Some(allowed) = plan.included_intervals() {
            let base_offset = plan.window_start();
            for (seg_start, seg_end) in allowed {
                let fetch_start = clamp_coordinate(*seg_start, chrom_length);
                let fetch_end = clamp_coordinate(*seg_end, chrom_length);
                if fetch_start >= fetch_end {
                    continue;
                }

                // Pre-fill the full included-interval range with 0.0
                // so uncovered positions (including those off the
                // chromosome) get 0.0 while introns stay NaN.
                if general.missing_data_as_zero {
                    let rs = (*seg_start - base_offset).max(0) as usize;
                    let re = ((*seg_end - base_offset).max(0) as usize).min(window_len);
                    if rs < re {
                        buf[rs..re].fill(0.0);
                    }
                }

                let intervals = sample
                    .reader_mut()
                    .values(&record.chrom, fetch_start, fetch_end)
                    .map_err(anyhow::Error::new)
                    .with_context(|| {
                        format!(
                            "Failed to read bigWig intervals for '{}' in '{}'",
                            record.chrom,
                            sample_path.display()
                        )
                    })?;

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
                    buf[rel_start..rel_end].fill(f64::from(interval.value));
                }
            }
        } else {
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
                            sample_path.display()
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
                    buf[rel_start..rel_end].fill(f64::from(interval.value));
                }
            }
        }

        let mut bins = Vec::with_capacity(bin_count);
        for bin in plan.bins() {
            let start_idx = index_from_coordinate(bin.start(), plan.window_start(), window_len);
            let end_idx = index_from_coordinate(bin.end(), plan.window_start(), window_len);

            let mut value = if start_idx < end_idx {
                aggregate_slice(&buf[start_idx..end_idx], general.average_type_bins)
            } else {
                None
            };

            if value.is_none() && general.missing_data_as_zero {
                value = Some(0.0);
            }

            let mut value = value.unwrap_or(f64::NAN);

            if nan_after_end && bin.beyond_region() {
                value = f64::NAN;
            }

            if value.is_finite() {
                value *= general.scale_factor;
            }

            bins.push(value);
        }

        Ok::<Vec<f64>, anyhow::Error>(bins)
    })?;

    Ok(bins)
}

pub fn compute_row<P: RegionPlan>(
    samples: &mut [Sample],
    record: &BedRecord,
    plan: &P,
    general: &GeneralOptions,
    nan_after_end: bool,
) -> Result<Option<(Vec<f64>, usize, usize)>> {
    let sample_count = samples.len();
    let bin_count = plan.bins().len();
    let mut all_values = Vec::with_capacity(sample_count * bin_count);
    for sample in samples.iter_mut() {
        let values = compute_sample_bins(sample, record, plan, general, nan_after_end)?;
        all_values.extend(values);
    }

    if should_skip_row_flat(&all_values, general) {
        return Ok(None);
    }

    Ok(Some((all_values, sample_count, bin_count)))
}

/// Process a single coalesced batch.
///
/// Performs one bigWig read per sample for the batch's merged query window,
/// then extracts per-region bins from the pre-read coverage buffers.
/// Items in metagene mode (where `included_intervals()` returns `Some`) fall
/// back to the original per-item `compute_row` path for correctness.
///
/// Records are **moved** (not cloned) out of the batch items, so once a
/// batch is processed its records are transferred into the result rows.
pub(super) fn process_batch<M: PipelineMode>(
    samples: &mut [Sample],
    batch: CoalescedBatch,
    mode: &M,
    general: &GeneralOptions,
    metadata: &M::Metadata,
    work_buf: &mut Vec<u8>,
    decode_buf: &mut Vec<u8>,
) -> Result<Vec<(usize, usize, Option<MatrixRow>)>> {
    if batch.items.is_empty() {
        return Ok(Vec::new());
    }

    let item_count = batch.items.len();
    let window_span = batch.query_end - batch.query_start;
    let sample_count = samples.len();

    // Zero or negative window span — delegate to per-item path.
    if window_span <= 0 {
        let nan_after_end = mode.nan_after_end(metadata);
        let mut results = Vec::with_capacity(item_count);
        for (orig_idx, group_index, record) in batch.items {
            let plan = mode.plan_for(&record, metadata);
            let maybe_values = compute_row(samples, &record, &plan, general, nan_after_end)?;
            let row = maybe_values.map(|(flat, sc, bc)| {
                mode.postprocess_row(Arc::unwrap_or_clone(record), flat, sc, bc, metadata)
            });
            results.push((orig_idx, group_index, row));
        }
        return Ok(results);
    }

    let window_len = usize::try_from(window_span).context("batch window span exceeds usize")?;

    let default_fill = if general.missing_data_as_zero {
        0.0f64
    } else {
        f64::NAN
    };

    let chrom = &batch.items[0].2.chrom;

    // Determine if direct aggregation is possible (Mean or Sum only)
    let use_direct = matches!(
        general.average_type_bins,
        AverageTypeBins::Mean | AverageTypeBins::Sum
    );

    // Store raw intervals per sample when using the direct path
    let mut sample_raw_intervals: Vec<Vec<(i64, i64, f64)>> = if use_direct {
        vec![Vec::new(); sample_count]
    } else {
        Vec::new()
    };

    // ── ONE bigWig read per sample for the entire merged window ────────
    let sample_paths: Vec<_> = samples.iter().map(|s| s.path().to_path_buf()).collect();

    // Only allocate coverage buffers when NOT using the direct path;
    // the direct path works from raw intervals and never needs per-base expansion.
    let mut sample_coverages = if use_direct {
        Vec::new()
    } else {
        take_coverage_buffers(sample_count, window_len, default_fill)
    };

    for (si, sample) in samples.iter_mut().enumerate() {
        let chrom_length = match sample.chrom_length(chrom) {
            Some(l) => l,
            None => {
                continue;
            }
        };

        let fetch_start = clamp_coordinate(batch.query_start, chrom_length);
        let fetch_end = clamp_coordinate(batch.query_end, chrom_length);

        if fetch_start >= fetch_end {
            continue;
        }

        let intervals = sample
            .values_with_bufs(chrom, fetch_start, fetch_end, work_buf, decode_buf)
            .with_context(|| {
                format!(
                    "Failed to read bigWig intervals for '{}' in '{}'",
                    chrom,
                    sample_paths[si].display()
                )
            })?;

        if use_direct {
            // Collect raw intervals for direct aggregation path;
            // no per-base coverage buffer fill needed.
            sample_raw_intervals[si] = intervals
                .iter()
                .map(|v| (i64::from(v.start), i64::from(v.end), f64::from(v.value)))
                .collect();
        } else {
            // Fill per-base coverage buffer for non-direct aggregation types.
            let cov = &mut sample_coverages[si];
            for v in intervals {
                let rs = i64::from(v.start).saturating_sub(batch.query_start).max(0);
                let re = i64::from(v.end)
                    .saturating_sub(batch.query_start)
                    .min(window_span)
                    .max(0);
                if rs < re {
                    cov[rs as usize..re as usize].fill(f64::from(v.value));
                }
            }
        }
    }

    // ── Extract per-region bins from the pre-read coverage buffers ─────
    let nan_after_end = mode.nan_after_end(metadata);
    let mut results = Vec::with_capacity(item_count);

    // Pre-allocate bin_coords outside the loop to avoid per-item allocation
    let mut bin_coords: Vec<(i64, i64)> = Vec::new();

    for (orig_idx, group_index, record) in batch.items {
        let plan = mode.plan_for(&record, metadata);

        // Metagene fallback: items with explicit included_intervals
        // (intron-skipping) must read individual exon intervals; we
        // delegate to the original per-item compute_row path.
        if plan.included_intervals().is_some() {
            let maybe_values = compute_row(samples, &record, &plan, general, nan_after_end)?;
            let row = maybe_values.map(|(flat, sc, bc)| {
                mode.postprocess_row(Arc::unwrap_or_clone(record), flat, sc, bc, metadata)
            });
            results.push((orig_idx, group_index, row));
            continue;
        }

        let bins = plan.bins();
        let bin_count = bins.len();
        let mut all_values = Vec::with_capacity(sample_count * bin_count);

        if use_direct {
            // Direct aggregation path: compute mean/sum directly from
            // raw intervals without expanding per-base coverage buffer.
            bin_coords.clear();
            bin_coords.extend(bins.iter().map(|bin| (bin.start(), bin.end())));

            for si in 0..sample_count {
                let direct_values = match general.average_type_bins {
                    AverageTypeBins::Mean => direct_mean_bins(
                        &bin_coords,
                        &sample_raw_intervals[si],
                        general.missing_data_as_zero,
                    ),
                    AverageTypeBins::Sum => direct_sum_bins(
                        &bin_coords,
                        &sample_raw_intervals[si],
                        general.missing_data_as_zero,
                    ),
                    _ => unreachable!(),
                };

                for (bi, value_option) in direct_values.into_iter().enumerate() {
                    // missing_data_as_zero is already handled inside direct_*_bins
                    let mut value = value_option.unwrap_or(f64::NAN);

                    if nan_after_end && bins[bi].beyond_region() {
                        value = f64::NAN;
                    }

                    if value.is_finite() {
                        value *= general.scale_factor;
                    }

                    all_values.push(value);
                }
            }
        } else {
            // Coverage buffer path: used for Median, Std, Min, Max
            for si in 0..sample_count {
                let cov = &sample_coverages[si];
                for bin in bins {
                    let bs = ((bin.start() - batch.query_start).max(0) as usize).min(window_len);
                    let be = ((bin.end() - batch.query_start).max(0) as usize).min(window_len);

                    let mut value = if bs < be {
                        aggregate_slice(&cov[bs..be], general.average_type_bins)
                    } else {
                        None
                    };

                    if value.is_none() && general.missing_data_as_zero {
                        value = Some(0.0);
                    }

                    let mut value = value.unwrap_or(f64::NAN);

                    if nan_after_end && bin.beyond_region() {
                        value = f64::NAN;
                    }

                    if value.is_finite() {
                        value *= general.scale_factor;
                    }

                    all_values.push(value);
                }
            }
        }

        let row = if should_skip_row_flat(&all_values, general) {
            None
        } else {
            Some(mode.postprocess_row(
                Arc::unwrap_or_clone(record),
                all_values,
                sample_count,
                bin_count,
                metadata,
            ))
        };
        results.push((orig_idx, group_index, row));
    }

    if !use_direct {
        return_coverage_buffers(sample_coverages);
    }
    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::AverageTypeBins;

    fn default_general() -> GeneralOptions {
        GeneralOptions {
            bin_size: 10,
            sort_regions: crate::config::SortRegions::Keep,
            sort_using: crate::config::SortUsing::Mean,
            sort_using_samples: None,
            average_type_bins: AverageTypeBins::Mean,
            missing_data_as_zero: false,
            skip_zeros: false,
            min_threshold: None,
            max_threshold: None,
            blacklist: None,
            samples_label: None,
            smart_labels: false,
            quiet: true,
            verbose: false,
            scale_factor: 1.0,
            number_of_processors: crate::config::ProcessorRequest::Fixed(1),
        }
    }

    // ── clamp_coordinate ──────────────────────────────────────────────────

    #[test]
    fn clamp_negative_to_zero() {
        assert_eq!(clamp_coordinate(-5, 1000), 0);
    }

    #[test]
    fn clamp_exceeding_chrom_length() {
        assert_eq!(clamp_coordinate(2000, 1000), 1000);
    }

    #[test]
    fn clamp_normal_value() {
        assert_eq!(clamp_coordinate(500, 1000), 500);
    }

    #[test]
    fn clamp_exactly_zero() {
        assert_eq!(clamp_coordinate(0, 1000), 0);
    }

    #[test]
    fn clamp_exactly_at_chrom_length() {
        assert_eq!(clamp_coordinate(1000, 1000), 1000);
    }

    // ── should_skip_row_flat ──────────────────────────────────────────────

    #[test]
    fn skip_zeros_all_zeros_returns_true() {
        let general = GeneralOptions {
            skip_zeros: true,
            ..default_general()
        };
        assert!(should_skip_row_flat(&[0.0, 0.0, 0.0], &general));
    }

    #[test]
    fn skip_zeros_has_nonzero_returns_false() {
        let general = GeneralOptions {
            skip_zeros: true,
            ..default_general()
        };
        assert!(!should_skip_row_flat(&[0.0, 1.0, 0.0], &general));
    }

    #[test]
    fn skip_zeros_all_nan_returns_true() {
        // NaN values are skipped; no non-zero found, so all_zero stays true
        let general = GeneralOptions {
            skip_zeros: true,
            ..default_general()
        };
        assert!(should_skip_row_flat(&[f64::NAN, f64::NAN], &general));
    }

    #[test]
    fn min_threshold_below_threshold_returns_true() {
        let general = GeneralOptions {
            min_threshold: Some(5.0),
            ..default_general()
        };
        assert!(should_skip_row_flat(&[3.0, 7.0], &general));
    }

    #[test]
    fn min_threshold_at_threshold_returns_true() {
        let general = GeneralOptions {
            min_threshold: Some(5.0),
            ..default_general()
        };
        assert!(should_skip_row_flat(&[5.0, 7.0], &general));
    }

    #[test]
    fn max_threshold_above_threshold_returns_true() {
        let general = GeneralOptions {
            max_threshold: Some(10.0),
            ..default_general()
        };
        assert!(should_skip_row_flat(&[8.0, 12.0], &general));
    }

    #[test]
    fn max_threshold_at_threshold_returns_true() {
        let general = GeneralOptions {
            max_threshold: Some(10.0),
            ..default_general()
        };
        assert!(should_skip_row_flat(&[8.0, 10.0], &general));
    }

    #[test]
    fn skip_row_both_thresholds_min_triggers() {
        let general = GeneralOptions {
            min_threshold: Some(5.0),
            max_threshold: Some(100.0),
            ..default_general()
        };
        assert!(should_skip_row_flat(&[3.0, 50.0], &general));
    }

    #[test]
    fn no_thresholds_returns_false() {
        let general = default_general();
        assert!(!should_skip_row_flat(&[1.0, 2.0, 3.0], &general));
    }

    // ── Issue 1a: skipZeros uses mean semantics ──────────────────────────

    #[test]
    fn skip_zeros_mean_zero_positive_negative_cancel() {
        // Python: mean([1.0, -1.0]) == 0.0 → skip
        let general = GeneralOptions {
            skip_zeros: true,
            ..default_general()
        };
        assert!(should_skip_row_flat(&[1.0, -1.0], &general));
    }

    #[test]
    fn skip_zeros_mean_nonzero_not_skipped() {
        // mean([2.0, -1.0]) == 0.5 → don't skip
        let general = GeneralOptions {
            skip_zeros: true,
            ..default_general()
        };
        assert!(!should_skip_row_flat(&[2.0, -1.0], &general));
    }

    #[test]
    fn skip_zeros_nan_and_zeros_skipped() {
        // mean of non-NaN values [0.0, 0.0] == 0.0 → skip
        let general = GeneralOptions {
            skip_zeros: true,
            ..default_general()
        };
        assert!(should_skip_row_flat(&[f64::NAN, 0.0, 0.0], &general));
    }
}
