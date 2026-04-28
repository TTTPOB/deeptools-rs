use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};

use crate::pipeline::matrix::MatrixData;
use super::formatting::write_plain_row;

pub fn write_matrix_values(path: &Path, matrix: &MatrixData) -> Result<()> {
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

pub fn write_sorted_regions(path: &Path, matrix: &MatrixData) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Failed to create sorted regions file '{}'", path.display()))?;
    let mut writer = BufWriter::new(file);
    writer.write_all(b"#chrom\tstart\tend\tname\tscore\tstrand\tgroup\n")?;

    for (index, row) in matrix.rows.iter().enumerate() {
        let name = row.record.name.as_deref().unwrap_or(".");
        let score = row
            .record
            .score_raw
            .clone()
            .or_else(|| row.record.score.map(|score| format!("{score:.6}")))
            .unwrap_or_else(|| ".".to_string());
        let strand = row
            .record
            .strand_raw
            .clone()
            .unwrap_or_else(|| row.record.strand.as_char().to_string());
        let group = group_label_for_index(matrix, index).unwrap_or(".");
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.record.chrom,
            row.record.start,
            row.record.end,
            name,
            score,
            strand,
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
