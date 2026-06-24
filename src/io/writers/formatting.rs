use std::cell::RefCell;
use std::io::{self, Write};

use anyhow::{Context, Result};
use itoa::Buffer;

use crate::pipeline::matrix::MatrixRow;

thread_local! {
    static ROW_BUFFER: RefCell<Vec<u8>> = RefCell::new(Vec::with_capacity(32768));
}

pub fn write_matrix_row<W: Write>(writer: &mut W, row: &MatrixRow) -> Result<()> {
    ROW_BUFFER.with(|cell| -> Result<()> {
        let mut buffer = cell.borrow_mut();
        buffer.clear();

        buffer.extend_from_slice(row.record.chrom.as_bytes());
        buffer.push(b'\t');

        let mut int_buffer = Buffer::new();

        if let Some(ref exon_coords) = row.exon_coords {
            for (i, (start, _)) in exon_coords.iter().enumerate() {
                if i > 0 {
                    buffer.push(b',');
                }
                buffer.extend_from_slice(int_buffer.format(*start).as_bytes());
            }
            buffer.push(b'\t');
            for (i, (_, end)) in exon_coords.iter().enumerate() {
                if i > 0 {
                    buffer.push(b',');
                }
                buffer.extend_from_slice(int_buffer.format(*end).as_bytes());
            }
        } else {
            buffer.extend_from_slice(int_buffer.format(row.record.start).as_bytes());
            buffer.push(b'\t');
            buffer.extend_from_slice(int_buffer.format(row.record.end).as_bytes());
        }
        buffer.push(b'\t');

        if matches!(row.record.bed_field_count, Some(3..=5)) {
            write_bed_coordinate_name(
                &mut *buffer,
                row.record.chrom.as_ref(),
                row.record.start,
                row.record.end,
            )?;
        } else {
            let name = row.record.name.as_deref().unwrap_or(".");
            buffer.extend_from_slice(name.as_bytes());
        }
        buffer.push(b'\t');

        if matches!(row.record.bed_field_count, Some(5)) {
            buffer.push(b'.');
        } else if row.record.score_raw.is_some() {
            write_score_value(&mut *buffer, 0.0)?;
        } else if let Some(score) = row.record.score {
            write_score_value(&mut *buffer, f64::from(score))?;
        } else {
            buffer.push(b'.');
        }
        buffer.push(b'\t');

        if row.record.strand_raw.is_some() {
            buffer.push(b'.');
        } else {
            buffer.push(row.record.strand.as_char() as u8);
        }

        for value in &row.values {
            buffer.push(b'\t');
            write_matrix_value(&mut *buffer, *value)?;
        }

        buffer.push(b'\n');

        writer
            .write_all(&buffer)
            .context("Failed to write matrix row")?;
        Ok(())
    })
}

pub fn write_plain_row<W: Write>(writer: &mut W, row: &MatrixRow) -> Result<()> {
    let mut first = true;
    for value in &row.values {
        if !first {
            writer.write_all(b"\t")?;
        } else {
            first = false;
        }
        writer.write_all(format_plain_value(*value).as_bytes())?;
    }
    writer.write_all(b"\n")?;
    Ok(())
}

fn write_bed_coordinate_name<W: Write>(
    writer: &mut W,
    chrom: &str,
    start: u32,
    end: u32,
) -> io::Result<()> {
    let mut int_buffer = Buffer::new();
    writer.write_all(chrom.as_bytes())?;
    writer.write_all(b":")?;
    writer.write_all(int_buffer.format(start).as_bytes())?;
    writer.write_all(b"-")?;
    writer.write_all(int_buffer.format(end).as_bytes())
}

/// Format an f64 to match Python's `str(float(x))` output.
///
/// Python 3 uses shortest-representation formatting, always including a
/// decimal point: `0` → `"0.0"`, `1.5` → `"1.5"`, `0.123456789` →
/// `"0.123456789"`.  Rust's `Display` for `f64` omits the fractional part
/// for whole numbers, so we append `.0` when needed.
pub fn write_score_value<W: Write>(writer: &mut W, value: f64) -> io::Result<()> {
    let s = format!("{value}");
    writer.write_all(s.as_bytes())?;
    if !s.contains('.') {
        writer.write_all(b".0")?;
    }
    Ok(())
}

pub fn write_matrix_value<W: Write>(writer: &mut W, value: f64) -> io::Result<()> {
    if value.is_nan() {
        return writer.write_all(b"nan");
    }

    if value.is_infinite() {
        return if value.is_sign_negative() {
            writer.write_all(b"-inf")
        } else {
            writer.write_all(b"inf")
        };
    }

    if value > -1e7 && value < 1e7 {
        let scaled = (value * 1_000_000.0).round_ties_even();
        return write_scaled_i64(writer, scaled as i64);
    }

    let scaled = (value * 1_000_000.0).round_ties_even();
    if !scaled.is_finite() || scaled.abs() > i128::MAX as f64 {
        let fallback = format!("{value:.6}");
        return writer.write_all(fallback.as_bytes());
    }

    write_scaled_i128(writer, scaled as i128)
}

#[inline]
pub fn write_scaled_i64<W: Write>(writer: &mut W, scaled: i64) -> io::Result<()> {
    let mut buffer = itoa::Buffer::new();

    if scaled == 0 {
        return writer.write_all(b"0.000000");
    }

    if scaled < 0 {
        writer.write_all(b"-")?;
    }
    let abs = scaled.unsigned_abs();
    let integer_part = abs / 1_000_000;
    let fractional_part = (abs % 1_000_000) as u32;

    writer.write_all(buffer.format(integer_part).as_bytes())?;
    writer.write_all(b".")?;

    let mut frac_digits = [b'0'; 6];
    let mut remainder = fractional_part;
    for slot in frac_digits.iter_mut().rev() {
        *slot = b'0' + (remainder % 10) as u8;
        remainder /= 10;
    }
    writer.write_all(&frac_digits)
}

#[inline]
pub fn write_scaled_i128<W: Write>(writer: &mut W, scaled: i128) -> io::Result<()> {
    let mut buffer = itoa::Buffer::new();
    let sign_negative = scaled < 0;

    if scaled == 0 {
        if sign_negative {
            writer.write_all(b"-")?;
        }
        return writer.write_all(b"0.000000");
    }

    let abs = if sign_negative { -scaled } else { scaled };
    let integer_part = (abs / 1_000_000) as u128;
    let fractional_part = (abs % 1_000_000) as u32;

    if sign_negative {
        writer.write_all(b"-")?;
    }
    writer.write_all(buffer.format(integer_part).as_bytes())?;
    writer.write_all(b".")?;

    let mut frac_digits = [b'0'; 6];
    let mut remainder = fractional_part;
    for slot in frac_digits.iter_mut().rev() {
        *slot = b'0' + (remainder % 10) as u8;
        remainder /= 10;
    }
    writer.write_all(&frac_digits)
}

/// Format an f64 to match Python's `%.4g` format (4 significant digits,
/// g-format: switches between fixed and scientific, strips trailing zeros).
///
/// Rules (matching C/Python `%.4g`):
/// - NaN  → `"nan"`
/// - ±Inf → `"inf"` / `"-inf"`
/// - 0    → `"0"`
/// - If exponent `e` satisfies `-4 <= e < 4`: fixed notation with
///   `max(0, 3 - e)` decimal places, trailing zeros and trailing dot stripped.
/// - Otherwise: scientific notation `X.XXXe±NN`, trailing zeros in the
///   mantissa stripped, exponent always two-digit with explicit sign.
pub fn format_plain_value(value: f64) -> String {
    if value.is_nan() {
        return "nan".to_string();
    }
    if value.is_infinite() {
        return if value.is_sign_negative() {
            "-inf".to_string()
        } else {
            "inf".to_string()
        };
    }
    if value == 0.0 {
        return "0".to_string();
    }

    let abs = value.abs();

    // Use Rust's built-in formatting engine to round to 4 significant digits.
    // This matches C printf rounding semantics, avoiding manual multiply/round/divide
    // which can introduce floating-point errors on rounding boundaries.
    let sci = format!("{abs:.3e}");

    // Parse the exponent from the formatted string to determine format bucket.
    let exp: i32 = sci
        .split_once('e')
        .expect("expected 'e' in scientific format")
        .1
        .parse()
        .expect("expected integer exponent");

    if exp >= -4 && exp < 4 {
        // Fixed notation: number of decimal places = max(0, precision - 1 - exp)
        let decimals = (3 - exp).max(0) as usize;
        // Re-round value to the correct number of decimal places for fixed output.
        // Parse the rounded mantissa from the scientific string to avoid double-rounding.
        let rounded_abs: f64 =
            sci.split_once('e').unwrap().0.parse::<f64>().unwrap() * 10f64.powi(exp);
        let rounded_value = if value.is_sign_negative() {
            -rounded_abs
        } else {
            rounded_abs
        };
        let s = format!("{rounded_value:.prec$}", prec = decimals);
        strip_trailing_zeros_fixed(&s)
    } else {
        // Scientific notation — reuse the already-formatted string
        let s = if value.is_sign_negative() {
            format!("-{sci}")
        } else {
            sci
        };
        normalize_scientific_notation(&s)
    }
}

/// Strip trailing zeros after the decimal point in a fixed-notation string.
/// Also strips the decimal point if it becomes the last character.
fn strip_trailing_zeros_fixed(s: &str) -> String {
    if !s.contains('.') {
        return s.to_string();
    }
    let trimmed = s.trim_end_matches('0');
    let trimmed = trimmed.trim_end_matches('.');
    trimmed.to_string()
}

/// Normalize Rust's scientific notation output to match Python's `%.4g`.
///
/// Rust produces e.g. `1.235e4` or `1.235e-5`.
/// Python produces `1.235e+04` or `1.235e-05`.
///
/// This function:
/// 1. Strips trailing zeros from the mantissa (before 'e').
/// 2. Strips a trailing decimal point from the mantissa.
/// 3. Normalizes the exponent to always have a sign and at least 2 digits.
fn normalize_scientific_notation(s: &str) -> String {
    let (mantissa, exp_str) = s
        .split_once('e')
        .expect("expected 'e' in scientific format");

    // Strip trailing zeros from mantissa
    let mantissa = if mantissa.contains('.') {
        let m = mantissa.trim_end_matches('0');
        m.trim_end_matches('.')
    } else {
        mantissa
    };

    // Parse exponent and format with sign and at least 2 digits
    let exp_val: i32 = exp_str.parse().expect("expected integer exponent");
    let exp_formatted = if exp_val >= 0 {
        format!("e+{:02}", exp_val)
    } else {
        format!("e-{:02}", exp_val.unsigned_abs())
    };

    format!("{mantissa}{exp_formatted}")
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use crate::io::{BedRecord, Strand};
    use crate::pipeline::matrix::MatrixRow;

    use super::*;

    fn fmt_row(row: &MatrixRow) -> String {
        let mut buf = Vec::new();
        write_matrix_row(&mut buf, row).unwrap();
        String::from_utf8(buf).unwrap()
    }

    fn matrix_row_with_raw_bed_fields() -> MatrixRow {
        MatrixRow {
            record: BedRecord {
                chrom: Arc::from("chr1"),
                start: 10,
                end: 20,
                name: Some("raw".to_string()),
                bed_field_count: Some(6),
                score: None,
                score_raw: Some("abc".to_string()),
                strand: Strand::Unstranded,
                strand_raw: Some("strandx".to_string()),
                extra_fields: Vec::new(),
            },
            sample_count: 1,
            bin_count: 1,
            values: vec![1.0],
            exon_coords: None,
        }
    }

    fn minimal_matrix_row(field_count: usize) -> MatrixRow {
        MatrixRow {
            record: BedRecord {
                chrom: Arc::from("chr1"),
                start: 10,
                end: 20,
                name: (field_count >= 4).then(|| "named_only".to_string()),
                bed_field_count: Some(field_count),
                score: (field_count >= 5).then_some(5.0),
                score_raw: None,
                strand: Strand::Unstranded,
                strand_raw: None,
                extra_fields: Vec::new(),
            },
            sample_count: 1,
            bin_count: 1,
            values: vec![1.0],
            exon_coords: None,
        }
    }

    fn bed6_matrix_row_with_empty_score_and_strand() -> MatrixRow {
        MatrixRow {
            record: BedRecord {
                chrom: Arc::from("chr1"),
                start: 10,
                end: 20,
                name: Some("foo".to_string()),
                bed_field_count: Some(6),
                score: None,
                score_raw: None,
                strand: Strand::Unstranded,
                strand_raw: None,
                extra_fields: Vec::new(),
            },
            sample_count: 1,
            bin_count: 1,
            values: vec![1.0],
            exon_coords: None,
        }
    }

    fn fmt_score(value: f64) -> String {
        let mut buf = Vec::new();
        write_score_value(&mut buf, value).unwrap();
        String::from_utf8(buf).unwrap()
    }

    fn fmt_value(value: f64) -> String {
        let mut buf = Vec::new();
        write_matrix_value(&mut buf, value).unwrap();
        String::from_utf8(buf).unwrap()
    }

    #[test]
    fn score_integers() {
        assert_eq!(fmt_score(0.0), "0.0");
        assert_eq!(fmt_score(1.0), "1.0");
        assert_eq!(fmt_score(5.0), "5.0");
        assert_eq!(fmt_score(-3.0), "-3.0");
    }

    #[test]
    fn score_fractional() {
        assert_eq!(fmt_score(1.5), "1.5");
        assert_eq!(fmt_score(0.123456789), "0.123456789");
        assert_eq!(fmt_score(-0.5), "-0.5");
    }

    #[test]
    fn score_special() {
        // NaN and inf don't contain '.', so they get ".0" appended.
        // These are not realistic BED scores, but verify the output is stable.
        assert_eq!(fmt_score(f64::NAN), "NaN.0");
        assert_eq!(fmt_score(f64::INFINITY), "inf.0");
        assert_eq!(fmt_score(f64::NEG_INFINITY), "-inf.0");
    }

    #[test]
    fn matrix_row_normalizes_invalid_bed_score_and_strand() {
        let line = fmt_row(&matrix_row_with_raw_bed_fields());
        let cols: Vec<&str> = line.trim_end().split('\t').collect();

        assert_eq!(cols[4], "0.0");
        assert_eq!(cols[5], ".");
        assert_eq!(cols[6], "1.000000");
    }

    #[test]
    fn matrix_row_uses_coordinate_name_for_minimal_bed_fields() {
        for field_count in 3..=5 {
            let line = fmt_row(&minimal_matrix_row(field_count));
            let cols: Vec<&str> = line.trim_end().split('\t').collect();

            assert_eq!(cols[3], "chr1:10-20");
            assert_eq!(cols[4], ".");
            assert_eq!(cols[5], ".");
        }
    }

    #[test]
    fn matrix_row_preserves_bed6_name_with_empty_score_and_strand() {
        let line = fmt_row(&bed6_matrix_row_with_empty_score_and_strand());
        let cols: Vec<&str> = line.trim_end().split('\t').collect();

        assert_eq!(cols[3], "foo");
        assert_eq!(cols[4], ".");
        assert_eq!(cols[5], ".");
    }

    #[test]
    fn nan_and_infinity() {
        assert_eq!(fmt_value(f64::NAN), "nan");
        assert_eq!(fmt_value(f64::INFINITY), "inf");
        assert_eq!(fmt_value(f64::NEG_INFINITY), "-inf");
    }

    #[test]
    fn zero() {
        assert_eq!(fmt_value(0.0), "0.000000");
    }

    #[test]
    fn small_positive() {
        assert_eq!(fmt_value(0.5), "0.500000");
        assert_eq!(fmt_value(1.0), "1.000000");
        assert_eq!(fmt_value(1.5), "1.500000");
        assert_eq!(fmt_value(0.001), "0.001000");
    }

    #[test]
    fn small_negative() {
        assert_eq!(fmt_value(-0.5), "-0.500000");
        assert_eq!(fmt_value(-1.0), "-1.000000");
        assert_eq!(fmt_value(-1.5), "-1.500000");
    }

    #[test]
    fn large_still_in_i64_range() {
        let v = 9_999_999.0f64;
        let result = fmt_value(v);
        assert!(result.starts_with("9999999"));
        assert!(result.ends_with("000000"));
    }

    #[test]
    fn very_large_falls_back_to_i128() {
        let v = 1e8f64;
        let result = fmt_value(v);
        assert!(!result.contains("nan") && !result.contains("inf"));
    }

    #[test]
    fn rounding_ties_to_even() {
        assert_eq!(fmt_value(0.0000005), "0.000000");
        assert_eq!(fmt_value(0.0000015), "0.000002");
    }

    #[test]
    fn max_precision_values() {
        assert_eq!(fmt_value(0.123456), "0.123456");
        assert_eq!(fmt_value(0.000001), "0.000001");
        assert_eq!(fmt_value(0.999999), "0.999999");
    }

    // ---------------------------------------------------------------
    // format_plain_value: Python %.4g parity tests
    // ---------------------------------------------------------------

    #[test]
    fn plain_value_nan() {
        assert_eq!(format_plain_value(f64::NAN), "nan");
    }

    #[test]
    fn plain_value_inf() {
        assert_eq!(format_plain_value(f64::INFINITY), "inf");
        assert_eq!(format_plain_value(f64::NEG_INFINITY), "-inf");
    }

    #[test]
    fn plain_value_zero() {
        assert_eq!(format_plain_value(0.0), "0");
    }

    #[test]
    fn plain_value_whole_numbers() {
        assert_eq!(format_plain_value(1.0), "1");
        assert_eq!(format_plain_value(2.0), "2");
        assert_eq!(format_plain_value(9999.0), "9999");
    }

    #[test]
    fn plain_value_simple_fractions() {
        assert_eq!(format_plain_value(0.5), "0.5");
        assert_eq!(format_plain_value(1.5), "1.5");
        assert_eq!(format_plain_value(-0.5), "-0.5");
    }

    #[test]
    fn plain_value_small_fixed() {
        assert_eq!(format_plain_value(0.001), "0.001");
        assert_eq!(format_plain_value(0.0001234), "0.0001234");
    }

    #[test]
    fn plain_value_scientific_small() {
        // exp = -5, outside [-4, 4) range → scientific
        assert_eq!(format_plain_value(0.00009999), "9.999e-05");
    }

    #[test]
    fn plain_value_scientific_large() {
        // exp = 4, outside [-4, 4) range → scientific
        assert_eq!(format_plain_value(10000.0), "1e+04");
        assert_eq!(format_plain_value(12345.0), "1.234e+04");
    }

    #[test]
    fn plain_value_rounding() {
        // 123.456 → 4 sig digits → 123.5
        assert_eq!(format_plain_value(123.456), "123.5");
        // 12345.678 → 4 sig digits → 1.235e+04
        assert_eq!(format_plain_value(12345.678), "1.235e+04");
    }

    #[test]
    fn plain_value_negative() {
        assert_eq!(format_plain_value(-1.0), "-1");
        assert_eq!(format_plain_value(-123.456), "-123.5");
        assert_eq!(format_plain_value(-0.00009999), "-9.999e-05");
    }

    #[test]
    fn plain_value_strip_trailing_zeros() {
        // 1200.0 → 4 sig digits → "1200" (no decimal)
        assert_eq!(format_plain_value(1200.0), "1200");
        // 1.500 → "1.5" (trailing zero stripped)
        assert_eq!(format_plain_value(1.500), "1.5");
    }

    // ---------------------------------------------------------------
    // Rounding across exponent boundary tests
    // ---------------------------------------------------------------

    #[test]
    fn plain_value_rounding_crosses_fixed_to_scientific() {
        // 9999.5 rounds to 10000 → exp crosses from 3 to 4 → scientific
        assert_eq!(format_plain_value(9999.5), "1e+04");
        // Negative version
        assert_eq!(format_plain_value(-9999.5), "-1e+04");
    }

    #[test]
    fn plain_value_rounding_crosses_scientific_to_fixed() {
        // 0.000099999 rounds to 0.0001 → exp crosses from -5 to -4 → fixed
        assert_eq!(format_plain_value(0.000099999), "0.0001");
    }

    #[test]
    fn plain_value_rounding_stays_fixed_no_boundary_cross() {
        // 99.95 stays at exp=1 after rounding → fixed
        assert_eq!(format_plain_value(99.95), "99.95");
    }

    #[test]
    fn plain_value_rounding_small_crosses_into_fixed() {
        // 0.00099995 → rounds to 0.001, exp from -4 to -3 → stays in fixed range
        assert_eq!(format_plain_value(0.00099995), "0.001");
    }

    #[test]
    fn plain_value_small_scientific_no_boundary_cross() {
        // 0.00009999 stays at exp=-5 → scientific
        assert_eq!(format_plain_value(0.00009999), "9.999e-05");
        // 0.000099995 IEEE 754 representation is 9.99949999...e-05, rounds down → scientific
        assert_eq!(format_plain_value(0.000099995), "9.999e-05");
    }
}
