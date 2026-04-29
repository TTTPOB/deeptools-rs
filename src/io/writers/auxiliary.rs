use std::io::Write;

use anyhow::Result;

use super::formatting::write_plain_row;
use crate::pipeline::matrix::MatrixRow;

/// Write a single sorted-region BED row for streaming output.
pub fn write_sorted_region_row<W: Write>(
    writer: &mut W,
    row: &MatrixRow,
    group_label: &str,
) -> Result<()> {
    let name = row.record.name.as_deref().unwrap_or(".");
    let score = row
        .record
        .score_raw
        .as_deref()
        .map(|s| s.to_string())
        .or_else(|| row.record.score.map(|score| format!("{score:.6}")))
        .unwrap_or_else(|| ".".to_string());
    let strand = row
        .record
        .strand_raw
        .as_deref()
        .map(|s| s.to_string())
        .unwrap_or_else(|| row.record.strand.as_char().to_string());
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
        row.record.chrom, row.record.start, row.record.end, name, score, strand, group_label
    )?;
    Ok(())
}

/// Write a single plain-values row for streaming output.
pub fn write_plain_values_row<W: Write>(writer: &mut W, row: &MatrixRow) -> Result<()> {
    write_plain_row(writer, row)
}
