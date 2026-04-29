use crate::parse::MatrixFileRow;

/// Summary of a single value-level difference.
#[derive(Debug, Clone)]
pub struct ValueDiff {
    pub row: usize,
    pub col: usize,
    pub left: f64,
    pub right: f64,
    pub abs_diff: f64,
}

/// Overall result of a value comparison.
#[derive(Debug)]
pub struct ComparisonResult {
    /// True when row counts, column counts, and all values match within tolerance.
    pub matches: bool,
    pub row_count_left: usize,
    pub row_count_right: usize,
    pub row_count_match: bool,
    /// True only when every row has the same column count in both files.
    pub col_count_match: bool,
    /// Column count mismatches: (row_index, left_col_count, right_col_count)
    pub col_count_diffs: Vec<(usize, usize, usize)>,
    /// Value diffs up to `max_diffs_reported`.
    pub value_diffs: Vec<ValueDiff>,
    /// Total number of value diffs (may exceed the reported list).
    pub total_value_diffs: usize,
    pub tolerance: f64,
}

/// Compare numeric values of two sets of matrix rows.
///
/// * Checks row count equality.
/// * Checks column count equality for **all** rows (not just the first).
/// * Treats (NaN, NaN) pairs as equal.
/// * Reports diffs up to `max_diffs_reported`.
pub fn compare_values(
    left: &[MatrixFileRow],
    right: &[MatrixFileRow],
    tolerance: f64,
    max_diffs_reported: usize,
) -> ComparisonResult {
    let row_count_left = left.len();
    let row_count_right = right.len();
    let row_count_match = row_count_left == row_count_right;

    let mut col_count_diffs: Vec<(usize, usize, usize)> = Vec::new();
    let mut value_diffs: Vec<ValueDiff> = Vec::new();
    let mut total_value_diffs: usize = 0;

    let common_rows = row_count_left.min(row_count_right);

    for row_idx in 0..common_rows {
        let lrow = &left[row_idx];
        let rrow = &right[row_idx];

        let lcols = lrow.values.len();
        let rcols = rrow.values.len();

        if lcols != rcols {
            col_count_diffs.push((row_idx, lcols, rcols));
            // Still compare the overlapping columns so the caller gets value diffs too
        }

        let common_cols = lcols.min(rcols);
        for col_idx in 0..common_cols {
            let lv = lrow.values[col_idx];
            let rv = rrow.values[col_idx];

            if values_equal(lv, rv, tolerance) {
                continue;
            }

            total_value_diffs += 1;
            if value_diffs.len() < max_diffs_reported {
                value_diffs.push(ValueDiff {
                    row: row_idx,
                    col: col_idx,
                    left: lv,
                    right: rv,
                    abs_diff: (lv - rv).abs(),
                });
            }
        }

        // Count extra columns (beyond the shorter row) as diffs
        let extra = if lcols > rcols {
            lcols - rcols
        } else {
            rcols - lcols
        };
        if extra > 0 {
            // Each extra column in the longer row is effectively a mismatch vs a missing value.
            // We don't emit individual ValueDiff entries for these since they're structural,
            // but we count them toward the total to signal something is wrong.
            total_value_diffs += extra;
        }
    }

    let col_count_match = col_count_diffs.is_empty();

    let matches = row_count_match && col_count_match && total_value_diffs == 0;

    ComparisonResult {
        matches,
        row_count_left,
        row_count_right,
        row_count_match,
        col_count_match,
        col_count_diffs,
        value_diffs,
        total_value_diffs,
        tolerance,
    }
}

/// Two values are equal when:
/// - Both are NaN (NaN == NaN by convention here).
/// - Or their absolute difference is within tolerance.
#[inline]
fn values_equal(left: f64, right: f64, tolerance: f64) -> bool {
    if left.is_nan() && right.is_nan() {
        return true;
    }
    (left - right).abs() <= tolerance
}
