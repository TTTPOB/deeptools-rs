use crate::config::AverageTypeBins;

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::AverageTypeBins;

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
    fn aggregate_min_negative_values() {
        let result = aggregate_slice(&[-3.0, -1.0, -2.0], AverageTypeBins::Min);
        assert_eq!(result, Some(-3.0));
    }

    #[test]
    fn aggregate_min_single() {
        let result = aggregate_slice(&[42.0], AverageTypeBins::Min);
        assert_eq!(result, Some(42.0));
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
        let result = aggregate_slice(&[42.0], AverageTypeBins::Max);
        assert_eq!(result, Some(42.0));
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

    #[test]
    fn aggregate_std_mixed_nan_nontrivial() {
        let result = aggregate_slice(&[f64::NAN, 1.0, 3.0], AverageTypeBins::Std);
        assert!((result.unwrap() - 1.0).abs() < 1e-10);
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
    fn index_one_past_base_is_one() {
        assert_eq!(index_from_coordinate(101, 100, 50), 1);
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
}
