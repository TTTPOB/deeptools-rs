use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::Compression;
use flate2::write::GzEncoder;

use crate::config::IoOptions;
use crate::pipeline::matrix::{MatrixData, MatrixRow};

/// Persist all requested outputs derived from the matrix computation.
pub fn write_outputs(matrix: &MatrixData, io: &IoOptions) -> Result<()> {
    write_matrix_gz(&io.matrix_output, matrix)?;

    if let Some(path) = &io.matrix_values_output {
        write_matrix_values(path, matrix)?;
    }

    if let Some(path) = &io.sorted_regions_output {
        write_sorted_regions(path, matrix)?;
    }

    Ok(())
}

fn write_matrix_gz(path: &Path, matrix: &MatrixData) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Failed to create matrix file '{}'", path.display()))?;
    let mut encoder = GzEncoder::new(file, Compression::default());

    let header = serde_json::to_string(&matrix.header)?;
    writeln!(encoder, "{}", header)?;

    for row in &matrix.rows {
        write_matrix_row(&mut encoder, row)?;
    }

    encoder.finish()?;
    Ok(())
}

fn write_matrix_row<W: Write>(writer: &mut W, row: &MatrixRow) -> Result<()> {
    let name = row.record.name.as_deref().unwrap_or(".");
    let score = row
        .record
        .score
        .map(|score| format!("{score:.6}"))
        .unwrap_or_else(|| ".".to_string());

    write!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}",
        row.record.chrom,
        row.record.start,
        row.record.end,
        name,
        score,
        row.record.strand.as_char()
    )?;

    for field in &row.record.extra_fields {
        write!(writer, "\t{}", field)?;
    }

    for sample_values in &row.values {
        for value in sample_values {
            write!(writer, "\t{}", format_matrix_value(*value))?;
        }
    }

    writer.write_all(b"\n")?;
    Ok(())
}

fn write_matrix_values(path: &Path, matrix: &MatrixData) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Failed to create matrix values file '{}'", path.display()))?;
    let mut writer = BufWriter::new(file);

    write_matrix_values_header(&mut writer, matrix)?;

    for row in &matrix.rows {
        write_plain_row(&mut writer, row)?;
    }

    writer.flush()?;
    Ok(())
}

fn write_matrix_values_header<W: Write>(writer: &mut W, matrix: &MatrixData) -> Result<()> {
    let header = &matrix.header;
    let group_lengths = diff(&header.group_boundaries);
    let mut group_info = Vec::new();
    for (label, length) in header.group_labels.iter().zip(group_lengths.iter()) {
        group_info.push(format!("{}:{}", label, length));
    }
    writeln!(writer, "#{}", group_info.join("\t"))?;

    let downstream = header.downstream.get(0).copied().unwrap_or_default();
    let upstream = header.upstream.get(0).copied().unwrap_or_default();
    let body = header.body.get(0).copied().unwrap_or_default();
    let bin_size = header.bin_size.get(0).copied().unwrap_or_default();
    let unscaled5 = header.unscaled_5_prime.get(0).copied().unwrap_or_default();
    let unscaled3 = header.unscaled_3_prime.get(0).copied().unwrap_or_default();

    writeln!(
        writer,
        "#downstream:{}\tupstream:{}\tbody:{}\tbin size:{}\tunscaled 5 prime:{}\tunscaled 3 prime:{}",
        downstream, upstream, body, bin_size, unscaled5, unscaled3
    )?;

    let sample_lengths = diff(&header.sample_boundaries);
    let mut labels_expanded = Vec::new();
    for (label, length) in header.sample_labels.iter().zip(sample_lengths.iter()) {
        for _ in 0..*length {
            labels_expanded.push(label.clone());
        }
    }
    writeln!(writer, "{}", labels_expanded.join("\t"))?;

    Ok(())
}

fn write_plain_row<W: Write>(writer: &mut W, row: &MatrixRow) -> Result<()> {
    let mut first = true;
    for sample_values in &row.values {
        for value in sample_values {
            if !first {
                writer.write_all(b"\t")?;
            } else {
                first = false;
            }
            writer.write_all(format_plain_value(*value).as_bytes())?;
        }
    }
    writer.write_all(b"\n")?;
    Ok(())
}

fn write_sorted_regions(path: &Path, matrix: &MatrixData) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Failed to create sorted regions file '{}'", path.display()))?;
    let mut writer = BufWriter::new(file);
    writer.write_all(b"#chrom\tstart\tend\tname\tscore\tstrand\tgroup\n")?;

    for (index, row) in matrix.rows.iter().enumerate() {
        let name = row.record.name.as_deref().unwrap_or(".");
        let score = row
            .record
            .score
            .map(|score| format!("{score:.6}"))
            .unwrap_or_else(|| ".".to_string());
        let group = group_label_for_index(matrix, index).unwrap_or(".");
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.record.chrom,
            row.record.start,
            row.record.end,
            name,
            score,
            row.record.strand.as_char(),
            group
        )?;
    }

    writer.flush()?;
    Ok(())
}

fn group_label_for_index<'a>(matrix: &'a MatrixData, index: usize) -> Option<&'a str> {
    let boundaries = &matrix.header.group_boundaries;
    let labels = &matrix.header.group_labels;
    for (idx, window) in boundaries.windows(2).enumerate() {
        if index >= window[0] && index < window[1] {
            return labels.get(idx).map(|label| label.as_str());
        }
    }
    None
}

fn diff(values: &[usize]) -> Vec<usize> {
    if values.len() < 2 {
        return Vec::new();
    }
    values.windows(2).map(|pair| pair[1] - pair[0]).collect()
}

fn format_matrix_value(value: f32) -> String {
    if value.is_nan() {
        "nan".to_string()
    } else {
        format!("{value:.6}")
    }
}

fn format_plain_value(value: f32) -> String {
    if value.is_nan() {
        "nan".to_string()
    } else {
        format!("{value:.4}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diff_computes_segment_lengths() {
        assert_eq!(diff(&[0, 2, 5, 9]), vec![2, 3, 4]);
        assert!(diff(&[0]).is_empty());
    }

    #[test]
    fn formatters_handle_nan() {
        assert_eq!(format_matrix_value(f32::NAN), "nan");
        assert_eq!(format_plain_value(f32::NAN), "nan");
    }
}
