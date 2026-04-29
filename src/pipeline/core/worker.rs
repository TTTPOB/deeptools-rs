use std::sync::Arc;

use anyhow::{Context, Result};

use crate::config::{AverageTypeBins, GeneralOptions};
use crate::io::BedRecord;
use crate::pipeline::matrix::MatrixRow;

use super::traits::{PipelineMode, RegionPlan, SignalBin};
use super::samples::Sample;
use super::coalesce::CoalescedBatch;

thread_local! {
    static COVERAGE_POOL: std::cell::RefCell<Vec<Vec<f64>>> =
        std::cell::RefCell::new(Vec::new());
}

fn take_coverage_buffers(sample_count: usize, window_len: usize, default_fill: f64) -> Vec<Vec<f64>> {
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

pub(super) fn index_from_coordinate(value: i64, base: i64, window_len: usize) -> usize {
    if value <= base {
        return 0;
    }
    let diff = value - base;
    let idx = usize::try_from(diff).unwrap_or(window_len);
    idx.min(window_len)
}

pub(super) fn aggregate_slice(slice: &[f64], average_type: AverageTypeBins) -> Option<f64> {
    let len = slice.len();
    if len == 0 {
        return None;
    }

    match average_type {
        AverageTypeBins::Mean => {
            let mut sum = 0.0f64;
            let mut count = 0u32;
            for &value in slice {
                if !value.is_nan() {
                    sum += value;
                    count += 1;
                }
            }
            if count == 0 { None } else { Some(sum / count as f64) }
        }
        AverageTypeBins::Sum => {
            let mut sum = 0.0f64;
            let mut found = false;
            for &value in slice {
                if !value.is_nan() {
                    sum += value;
                    found = true;
                }
            }
            if found { Some(sum) } else { None }
        }
        AverageTypeBins::Min => {
            let mut min = f64::INFINITY;
            let mut found = false;
            for &value in slice {
                if !value.is_nan() {
                    min = min.min(value);
                    found = true;
                }
            }
            if found { Some(min) } else { None }
        }
        AverageTypeBins::Max => {
            let mut max = f64::NEG_INFINITY;
            let mut found = false;
            for &value in slice {
                if !value.is_nan() {
                    max = max.max(value);
                    found = true;
                }
            }
            if found { Some(max) } else { None }
        }
        AverageTypeBins::Std => {
            let mut sum = 0.0f64;
            let mut count = 0u32;
            for &value in slice {
                if !value.is_nan() {
                    sum += value;
                    count += 1;
                }
            }
            if count == 0 {
                return None;
            }
            let mean = sum / count as f64;
            let mut variance_sum = 0.0f64;
            for &value in slice {
                if !value.is_nan() {
                    let delta = value - mean;
                    variance_sum += delta * delta;
                }
            }
            Some((variance_sum / count as f64).sqrt())
        }
        AverageTypeBins::Median => {
            let mut values: Vec<f64> = slice
                .iter()
                .copied()
                .filter(|v| !v.is_nan())
                .collect();
            if values.is_empty() {
                return None;
            }
            values.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let mid = values.len() / 2;
            if values.len() % 2 == 0 {
                Some((values[mid - 1] + values[mid]) / 2.0)
            } else {
                Some(values[mid])
            }
        }
    }
}

fn should_skip_row_flat(values: &[f64], general: &GeneralOptions) -> bool {
    if general.skip_zeros {
        let mut all_zero = true;
        for &value in values {
            if value.is_nan() {
                continue;
            }
            if value != 0.0 {
                all_zero = false;
                break;
            }
        }
        if all_zero {
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
        buf.resize(window_len, default_fill);

        if let Some(allowed) = plan.included_intervals() {
            let base_offset = plan.window_start();
            for (seg_start, seg_end) in allowed {
                let fetch_start = clamp_coordinate(*seg_start, chrom_length);
                let fetch_end = clamp_coordinate(*seg_end, chrom_length);
                if fetch_start >= fetch_end {
                    continue;
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
            let maybe_values =
                compute_row(samples, &record, &plan, general, nan_after_end)?;
            let row = maybe_values
                .map(|(flat, sc, bc)| mode.postprocess_row(Arc::unwrap_or_clone(record), flat, sc, bc, metadata));
            results.push((orig_idx, group_index, row));
        }
        return Ok(results);
    }

    let window_len =
        usize::try_from(window_span).context("batch window span exceeds usize")?;

    let default_fill = if general.missing_data_as_zero {
        0.0f64
    } else {
        f64::NAN
    };

    let chrom = &batch.items[0].2.chrom;

    // ── ONE bigWig read per sample for the entire merged window ────────
    let sample_paths: Vec<_> = samples.iter().map(|s| s.path().to_path_buf()).collect();
    let mut sample_coverages = take_coverage_buffers(sample_count, window_len, default_fill);
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
            .reader_mut()
            .values(chrom, fetch_start, fetch_end)
            .map_err(anyhow::Error::new)
            .with_context(|| {
                format!(
                    "Failed to read bigWig intervals for '{}' in '{}'",
                    chrom,
                    sample_paths[si].display()
                )
            })?;

        let cov = &mut sample_coverages[si];
        for v in intervals {
            let rs = i64::from(v.start)
                .saturating_sub(batch.query_start)
                .max(0);
            let re = i64::from(v.end)
                .saturating_sub(batch.query_start)
                .min(window_span)
                .max(0);
            if rs < re {
                cov[rs as usize..re as usize].fill(f64::from(v.value));
            }
        }
    }

    // ── Extract per-region bins from the pre-read coverage buffers ─────
    let nan_after_end = mode.nan_after_end(metadata);
    let mut results = Vec::with_capacity(item_count);

    for (orig_idx, group_index, record) in batch.items {
        let plan = mode.plan_for(&record, metadata);

        // Metagene fallback: items with explicit included_intervals
        // (intron-skipping) must read individual exon intervals; we
        // delegate to the original per-item compute_row path.
        if plan.included_intervals().is_some() {
            let maybe_values =
                compute_row(samples, &record, &plan, general, nan_after_end)?;
            let row = maybe_values
                .map(|(flat, sc, bc)| mode.postprocess_row(Arc::unwrap_or_clone(record), flat, sc, bc, metadata));
            results.push((orig_idx, group_index, row));
            continue;
        }

        let bins = plan.bins();
        let bin_count = bins.len();
        let mut all_values = Vec::with_capacity(sample_count * bin_count);

        for si in 0..sample_count {
            let cov = &sample_coverages[si];
            for bin in bins {
                let bs =
                    ((bin.start() - batch.query_start).max(0) as usize).min(window_len);
                let be =
                    ((bin.end() - batch.query_start).max(0) as usize).min(window_len);

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

    return_coverage_buffers(sample_coverages);
    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

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

    // ── aggregate_slice: Mean ──────────────────────────────────────────────

    #[test]
    fn aggregate_mean_normal() {
        let result = aggregate_slice(&[1.0, 2.0, 3.0], AverageTypeBins::Mean);
        assert_eq!(result, Some(2.0));
    }

    #[test]
    fn aggregate_mean_all_nan() {
        let result = aggregate_slice(&[f64::NAN, f64::NAN], AverageTypeBins::Mean);
        assert_eq!(result, None);
    }

    #[test]
    fn aggregate_mean_empty() {
        let result = aggregate_slice(&[], AverageTypeBins::Mean);
        assert_eq!(result, None);
    }

    #[test]
    fn aggregate_mean_mixed_nan() {
        let result = aggregate_slice(&[f64::NAN, 4.0, 6.0], AverageTypeBins::Mean);
        assert_eq!(result, Some(5.0));
    }

    #[test]
    fn aggregate_mean_single() {
        let result = aggregate_slice(&[7.0], AverageTypeBins::Mean);
        assert_eq!(result, Some(7.0));
    }

    // ── aggregate_slice: Sum ───────────────────────────────────────────────

    #[test]
    fn aggregate_sum_normal() {
        let result = aggregate_slice(&[1.0, 2.0, 3.0], AverageTypeBins::Sum);
        assert_eq!(result, Some(6.0));
    }

    #[test]
    fn aggregate_sum_all_nan() {
        let result = aggregate_slice(&[f64::NAN], AverageTypeBins::Sum);
        assert_eq!(result, None);
    }

    #[test]
    fn aggregate_sum_empty() {
        let result = aggregate_slice(&[], AverageTypeBins::Sum);
        assert_eq!(result, None);
    }

    #[test]
    fn aggregate_sum_mixed_nan() {
        let result = aggregate_slice(&[f64::NAN, 3.0, 5.0], AverageTypeBins::Sum);
        assert_eq!(result, Some(8.0));
    }

    #[test]
    fn aggregate_sum_single() {
        let result = aggregate_slice(&[9.0], AverageTypeBins::Sum);
        assert_eq!(result, Some(9.0));
    }

    // ── aggregate_slice: Min ───────────────────────────────────────────────

    #[test]
    fn aggregate_min_normal() {
        let result = aggregate_slice(&[3.0, 1.0, 2.0], AverageTypeBins::Min);
        assert_eq!(result, Some(1.0));
    }

    #[test]
    fn aggregate_min_all_nan() {
        let result = aggregate_slice(&[f64::NAN], AverageTypeBins::Min);
        assert_eq!(result, None);
    }

    #[test]
    fn aggregate_min_empty() {
        let result = aggregate_slice(&[], AverageTypeBins::Min);
        assert_eq!(result, None);
    }

    #[test]
    fn aggregate_min_mixed_nan() {
        let result = aggregate_slice(&[f64::NAN, 5.0, 2.0], AverageTypeBins::Min);
        assert_eq!(result, Some(2.0));
    }

    #[test]
    fn aggregate_min_single() {
        assert_eq!(
            aggregate_slice(&[42.0], AverageTypeBins::Min),
            Some(42.0)
        );
    }

    // ── aggregate_slice: Max ───────────────────────────────────────────────

    #[test]
    fn aggregate_max_normal() {
        let result = aggregate_slice(&[3.0, 1.0, 5.0], AverageTypeBins::Max);
        assert_eq!(result, Some(5.0));
    }

    #[test]
    fn aggregate_max_all_nan() {
        let result = aggregate_slice(&[f64::NAN], AverageTypeBins::Max);
        assert_eq!(result, None);
    }

    #[test]
    fn aggregate_max_empty() {
        let result = aggregate_slice(&[], AverageTypeBins::Max);
        assert_eq!(result, None);
    }

    #[test]
    fn aggregate_max_mixed_nan() {
        let result = aggregate_slice(&[f64::NAN, 4.0, 9.0], AverageTypeBins::Max);
        assert_eq!(result, Some(9.0));
    }

    #[test]
    fn aggregate_max_single() {
        assert_eq!(
            aggregate_slice(&[42.0], AverageTypeBins::Max),
            Some(42.0)
        );
    }

    // ── aggregate_slice: Std ───────────────────────────────────────────────

    #[test]
    fn aggregate_std_normal() {
        // values [2, 4, 4, 4, 5, 5, 7, 9], mean=5, population std=2
        let values = [2.0f64, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let result = aggregate_slice(&values, AverageTypeBins::Std).unwrap();
        assert!((result - 2.0).abs() < 1e-10);
    }

    #[test]
    fn aggregate_std_all_nan() {
        let result = aggregate_slice(&[f64::NAN], AverageTypeBins::Std);
        assert_eq!(result, None);
    }

    #[test]
    fn aggregate_std_empty() {
        let result = aggregate_slice(&[], AverageTypeBins::Std);
        assert_eq!(result, None);
    }

    #[test]
    fn aggregate_std_single() {
        // std of a single value is 0
        let result = aggregate_slice(&[5.0], AverageTypeBins::Std);
        assert_eq!(result, Some(0.0));
    }

    #[test]
    fn aggregate_std_mixed_nan() {
        // NaN values are ignored; std of [3, 3] is 0
        let result = aggregate_slice(&[f64::NAN, 3.0, 3.0], AverageTypeBins::Std);
        assert_eq!(result, Some(0.0));
    }

    // ── aggregate_slice: Median ────────────────────────────────────────────

    #[test]
    fn aggregate_median_odd_count() {
        let result = aggregate_slice(&[3.0, 1.0, 2.0], AverageTypeBins::Median);
        assert_eq!(result, Some(2.0));
    }

    #[test]
    fn aggregate_median_even_count() {
        let result = aggregate_slice(&[1.0, 2.0, 3.0, 4.0], AverageTypeBins::Median);
        assert_eq!(result, Some(2.5));
    }

    #[test]
    fn aggregate_median_all_nan() {
        let result = aggregate_slice(&[f64::NAN], AverageTypeBins::Median);
        assert_eq!(result, None);
    }

    #[test]
    fn aggregate_median_empty() {
        let result = aggregate_slice(&[], AverageTypeBins::Median);
        assert_eq!(result, None);
    }

    #[test]
    fn aggregate_median_mixed_nan() {
        // NaN filtered out; median of [1, 3] = 2.0
        let result = aggregate_slice(&[f64::NAN, 1.0, 3.0], AverageTypeBins::Median);
        assert_eq!(result, Some(2.0));
    }

    #[test]
    fn aggregate_median_single() {
        let result = aggregate_slice(&[42.0], AverageTypeBins::Median);
        assert_eq!(result, Some(42.0));
    }

    // ── index_from_coordinate ─────────────────────────────────────────────

    #[test]
    fn index_at_base_is_zero() {
        assert_eq!(index_from_coordinate(100, 100, 50), 0);
    }

    #[test]
    fn index_below_base_is_zero() {
        assert_eq!(index_from_coordinate(50, 100, 50), 0);
    }

    #[test]
    fn index_above_base_correct_offset() {
        assert_eq!(index_from_coordinate(110, 100, 50), 10);
    }

    #[test]
    fn index_beyond_window_len_clamped() {
        assert_eq!(index_from_coordinate(200, 100, 50), 50);
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
    fn no_thresholds_returns_false() {
        let general = default_general();
        assert!(!should_skip_row_flat(&[1.0, 2.0, 3.0], &general));
    }
}
