use serde_json::Value;

use crate::matrix_compare::header::{HeaderDiff, compare_headers};
use crate::matrix_compare::parse::MatrixFileRow;
use crate::matrix_compare::values::{ComparisonResult, compare_values};

/// Describes whether rows appear to be reordered between left and right.
#[derive(Debug)]
pub struct BedFieldDiff {
    pub reordered: bool,
    /// Rows present in left but missing in right (by key).
    pub left_only: Vec<String>,
    /// Rows present in right but missing in left (by key).
    pub right_only: Vec<String>,
    /// Row count difference when the totals differ.
    pub row_count_left: usize,
    pub row_count_right: usize,
}

/// Full diff result combining header, BED fields, and numeric values.
#[derive(Debug)]
pub struct FullDiff {
    pub header_diffs: Vec<HeaderDiff>,
    pub bed_diff: BedFieldDiff,
    pub value_result: ComparisonResult,
    pub matches: bool,
}

/// Run a full comparison: header + BED field ordering + numeric values.
pub fn full_diff(
    left_header: &Value,
    right_header: &Value,
    left_rows: &[MatrixFileRow],
    right_rows: &[MatrixFileRow],
    tolerance: f64,
    ignore_keys: &[String],
    max_diffs_reported: usize,
) -> FullDiff {
    let header_diffs = compare_headers(left_header, right_header, ignore_keys);
    let bed_diff = compare_bed_fields(left_rows, right_rows);
    let value_result = compare_values(left_rows, right_rows, tolerance, max_diffs_reported);

    let matches = header_diffs.is_empty()
        && !bed_diff.reordered
        && bed_diff.left_only.is_empty()
        && bed_diff.right_only.is_empty()
        && value_result.matches;

    FullDiff {
        header_diffs,
        bed_diff,
        value_result,
        matches,
    }
}

/// Compare BED fields (chrom:start:end:name tuples) to detect row reordering or
/// missing rows between the two files.
fn compare_bed_fields(left: &[MatrixFileRow], right: &[MatrixFileRow]) -> BedFieldDiff {
    use std::collections::HashMap;

    let mut left_counts: HashMap<String, usize> = HashMap::new();
    for row in left {
        *left_counts.entry(row.key()).or_insert(0) += 1;
    }

    let mut right_counts: HashMap<String, usize> = HashMap::new();
    for row in right {
        *right_counts.entry(row.key()).or_insert(0) += 1;
    }

    // Rows in left but not right
    let mut left_only: Vec<String> = left_counts
        .iter()
        .filter_map(|(key, &lcount)| {
            let rcount = right_counts.get(key).copied().unwrap_or(0);
            if lcount > rcount {
                Some(key.clone())
            } else {
                None
            }
        })
        .collect();
    left_only.sort();

    // Rows in right but not left
    let mut right_only: Vec<String> = right_counts
        .iter()
        .filter_map(|(key, &rcount)| {
            let lcount = left_counts.get(key).copied().unwrap_or(0);
            if rcount > lcount {
                Some(key.clone())
            } else {
                None
            }
        })
        .collect();
    right_only.sort();

    // Detect reordering: same set of rows but in different order
    let reordered = left_only.is_empty() && right_only.is_empty() && left.len() == right.len() && {
        // Check positional order differs
        left.iter()
            .zip(right.iter())
            .any(|(l, r)| l.key() != r.key())
    };

    BedFieldDiff {
        reordered,
        left_only,
        right_only,
        row_count_left: left.len(),
        row_count_right: right.len(),
    }
}

/// Print a human-readable summary of a FullDiff to stdout.
pub fn print_full_diff(diff: &FullDiff, left_path: &str, right_path: &str) {
    println!("=== Full diff: {} vs {} ===", left_path, right_path);

    // Header section
    if diff.header_diffs.is_empty() {
        println!("[header] OK — no differences");
    } else {
        println!("[header] {} difference(s):", diff.header_diffs.len());
        for d in &diff.header_diffs {
            println!("  key {:?}:", d.key);
            println!("    left:  {}", d.left);
            println!("    right: {}", d.right);
        }
    }

    // BED field section
    println!();
    print_bed_diff(&diff.bed_diff);

    // Values section
    println!();
    print_value_result(&diff.value_result);

    // Overall
    println!();
    if diff.matches {
        println!("RESULT: MATCH");
    } else {
        println!("RESULT: MISMATCH");
    }
}

/// Print a human-readable summary of BED field differences.
pub fn print_bed_diff(diff: &BedFieldDiff) {
    if diff.row_count_left != diff.row_count_right {
        println!(
            "[bed] Row count mismatch: left={}, right={}",
            diff.row_count_left, diff.row_count_right
        );
    }
    if diff.reordered {
        println!("[bed] Rows are reordered (same set, different order)");
    }
    if !diff.left_only.is_empty() {
        println!("[bed] {} row(s) in left only:", diff.left_only.len());
        for key in diff.left_only.iter().take(10) {
            println!("  {}", key);
        }
        if diff.left_only.len() > 10 {
            println!("  ... ({} more)", diff.left_only.len() - 10);
        }
    }
    if !diff.right_only.is_empty() {
        println!("[bed] {} row(s) in right only:", diff.right_only.len());
        for key in diff.right_only.iter().take(10) {
            println!("  {}", key);
        }
        if diff.right_only.len() > 10 {
            println!("  ... ({} more)", diff.right_only.len() - 10);
        }
    }
    if diff.reordered || !diff.left_only.is_empty() || !diff.right_only.is_empty() {
        // already printed something
    } else if diff.row_count_left == diff.row_count_right {
        println!("[bed] OK — {} rows, same order", diff.row_count_left);
    }
}

/// Print a human-readable summary of a value ComparisonResult.
pub fn print_value_result(result: &ComparisonResult) {
    if !result.row_count_match {
        println!(
            "[values] Row count mismatch: left={}, right={}",
            result.row_count_left, result.row_count_right
        );
    }
    if !result.col_count_match {
        println!(
            "[values] {} row(s) have column count mismatches:",
            result.col_count_diffs.len()
        );
        for &(row, lcols, rcols) in result.col_count_diffs.iter().take(5) {
            println!("  row {}: left={} cols, right={} cols", row, lcols, rcols);
        }
        if result.col_count_diffs.len() > 5 {
            println!("  ... ({} more)", result.col_count_diffs.len() - 5);
        }
    }
    if result.total_value_diffs == 0 && result.row_count_match && result.col_count_match {
        println!(
            "[values] OK — {} rows, tolerance={}",
            result.row_count_left, result.tolerance
        );
    } else {
        println!(
            "[values] {} value diff(s) found (tolerance={}):",
            result.total_value_diffs, result.tolerance
        );
        for d in &result.value_diffs {
            println!(
                "  row={} col={}: left={:.6} right={:.6} |diff|={:.2e}",
                d.row, d.col, d.left, d.right, d.abs_diff
            );
        }
        if result.total_value_diffs > result.value_diffs.len() {
            println!(
                "  ... ({} more not shown)",
                result.total_value_diffs - result.value_diffs.len()
            );
        }
    }
}
