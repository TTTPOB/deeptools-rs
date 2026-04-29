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

        let name = row.record.name.as_deref().unwrap_or(".");
        buffer.extend_from_slice(name.as_bytes());
        buffer.push(b'\t');

        if let Some(raw) = row.record.score_raw.as_deref() {
            buffer.extend_from_slice(raw.as_bytes());
        } else if let Some(score) = row.record.score {
            write_matrix_value(&mut *buffer, f64::from(score))?;
        } else {
            buffer.push(b'.');
        }
        buffer.push(b'\t');

        if let Some(raw) = row.record.strand_raw.as_deref() {
            buffer.extend_from_slice(raw.as_bytes());
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

pub fn format_plain_value(value: f64) -> String {
    if value.is_nan() {
        "nan".to_string()
    } else {
        format!("{value:.4}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fmt_value(value: f64) -> String {
        let mut buf = Vec::new();
        write_matrix_value(&mut buf, value).unwrap();
        String::from_utf8(buf).unwrap()
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
}
