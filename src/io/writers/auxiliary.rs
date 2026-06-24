use std::io::Write;

use anyhow::Result;

use super::formatting::{write_plain_row, write_score_value};
use crate::pipeline::matrix::MatrixRow;

/// Header for the 13-column sorted-regions BED12+group output.
pub const SORTED_REGIONS_HEADER: &[u8] =
    b"#chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRGB\tblockCount\tblockSizes\tblockStart\tdeepTools_group\n";

/// Write a single sorted-region BED12+group row for streaming output.
///
/// Produces 13 columns matching Python deepTools:
///   chrom, start, end, name, score, strand,
///   thickStart, thickEnd, itemRGB, blockCount, blockSizes, blockStarts,
///   deepTools_group
pub fn write_sorted_region_row<W: Write>(
    writer: &mut W,
    row: &MatrixRow,
    group_label: &str,
) -> Result<()> {
    let name = row.record.name.as_deref().unwrap_or(".");
    let strand = row
        .record
        .strand_raw
        .as_deref()
        .map(|s| s.to_string())
        .unwrap_or_else(|| row.record.strand.as_char().to_string());

    let start = row.record.start;
    let end = row.record.end;

    // BED6 columns: chrom, start, end, name
    write!(
        writer,
        "{}\t{}\t{}\t{}\t",
        row.record.chrom, start, end, name
    )?;

    // Score: use write_score_value for Python-matching float formatting
    if let Some(raw) = row.record.score_raw.as_deref() {
        writer.write_all(raw.as_bytes())?;
    } else if let Some(score) = row.record.score {
        write_score_value(writer, f64::from(score))?;
    } else {
        writer.write_all(b".")?;
    }

    // strand + BED12 synthetic fields + deepTools_group
    if let Some(ref exon_coords) = row.exon_coords {
        // Metagene (BED12/GTF) input: use actual exon blocks
        let block_count = exon_coords.len();
        let block_sizes: Vec<String> = exon_coords
            .iter()
            .map(|(s, e)| (e - s).to_string())
            .collect();
        // blockStarts: standard BED12 uses offsets relative to chromStart (0 for
        // single-block regions).  Python deepTools (via deeptoolsintervals) outputs
        // absolute genomic coordinates here instead, which is non-standard.  We
        // follow the BED12 specification; downstream tools that parse this column
        // will get correct results.
        let block_starts: Vec<String> = exon_coords
            .iter()
            .map(|(s, _)| (s - start).to_string())
            .collect();
        writeln!(
            writer,
            "\t{}\t{}\t{}\t0\t{}\t{}\t{}\t{}",
            strand,
            start,
            end,
            block_count,
            block_sizes.join(","),
            block_starts.join(","),
            group_label
        )?;
    } else {
        // Non-metagene (BED6) input: single synthetic block
        let block_size = end - start;
        writeln!(
            writer,
            "\t{}\t{}\t{}\t0\t1\t{}\t0\t{}",
            strand, start, end, block_size, group_label
        )?;
    }

    Ok(())
}

/// Write a single plain-values row for streaming output.
pub fn write_plain_values_row<W: Write>(writer: &mut W, row: &MatrixRow) -> Result<()> {
    write_plain_row(writer, row)
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::io::readers::bed::{BedRecord, Strand};

    fn make_row(
        start: u32,
        end: u32,
        score: Option<f32>,
        exon_coords: Option<Vec<(u32, u32)>>,
    ) -> MatrixRow {
        MatrixRow {
            record: BedRecord {
                chrom: Arc::from("ch1"),
                start,
                end,
                name: Some("gene1".to_string()),
                bed_field_count: None,
                score,
                score_raw: None,
                strand: Strand::Positive,
                strand_raw: None,
                extra_fields: vec![],
            },
            values: vec![],
            sample_count: 1,
            bin_count: 0,
            exon_coords,
        }
    }

    #[test]
    fn sorted_region_header_has_13_columns() {
        let header = std::str::from_utf8(SORTED_REGIONS_HEADER).unwrap();
        let cols: Vec<&str> = header.trim().split('\t').collect();
        assert_eq!(cols.len(), 13);
        assert_eq!(cols[0], "#chrom");
        assert_eq!(cols[12], "deepTools_group");
    }

    #[test]
    fn sorted_region_bed6_single_block() {
        let row = make_row(100, 150, Some(0.0), None);
        let mut buf = Vec::new();
        write_sorted_region_row(&mut buf, &row, "Group 1").unwrap();
        let line = String::from_utf8(buf).unwrap();
        let cols: Vec<&str> = line.trim().split('\t').collect();
        assert_eq!(cols.len(), 13, "expected 13 columns, got: {cols:?}");
        assert_eq!(cols[0], "ch1");
        assert_eq!(cols[1], "100");
        assert_eq!(cols[2], "150");
        assert_eq!(cols[3], "gene1");
        assert_eq!(cols[4], "0.0"); // Python-style float
        assert_eq!(cols[5], "+");
        assert_eq!(cols[6], "100"); // thickStart = start
        assert_eq!(cols[7], "150"); // thickEnd = end
        assert_eq!(cols[8], "0"); // itemRGB
        assert_eq!(cols[9], "1"); // blockCount
        assert_eq!(cols[10], "50"); // blockSizes = end - start
        // blockStarts follows BED12 spec (relative to chromStart), not the
        // absolute-coordinate convention used by Python's deeptoolsintervals.
        assert_eq!(cols[11], "0"); // blockStarts
        assert_eq!(cols[12], "Group 1");
    }

    #[test]
    fn sorted_region_metagene_exon_blocks() {
        // Two exon blocks: (100, 120) and (140, 150) within region 100..150
        let row = make_row(100, 150, Some(5.0), Some(vec![(100, 120), (140, 150)]));
        let mut buf = Vec::new();
        write_sorted_region_row(&mut buf, &row, "Group 2").unwrap();
        let line = String::from_utf8(buf).unwrap();
        let cols: Vec<&str> = line.trim().split('\t').collect();
        assert_eq!(cols.len(), 13);
        assert_eq!(cols[4], "5.0");
        assert_eq!(cols[9], "2"); // blockCount = 2 exons
        assert_eq!(cols[10], "20,10"); // blockSizes: 120-100, 150-140
        // blockStarts follows BED12 spec (relative to chromStart), not the
        // absolute-coordinate convention used by Python's deeptoolsintervals.
        assert_eq!(cols[11], "0,40"); // blockStarts: 100-100, 140-100
        assert_eq!(cols[12], "Group 2");
    }

    #[test]
    fn sorted_region_score_raw_passthrough() {
        let mut row = make_row(0, 100, None, None);
        row.record.score_raw = Some("3.14".to_string());
        let mut buf = Vec::new();
        write_sorted_region_row(&mut buf, &row, "G").unwrap();
        let line = String::from_utf8(buf).unwrap();
        let cols: Vec<&str> = line.trim().split('\t').collect();
        assert_eq!(cols[4], "3.14");
    }

    #[test]
    fn sorted_region_no_score_outputs_dot() {
        let row = make_row(0, 100, None, None);
        let mut buf = Vec::new();
        write_sorted_region_row(&mut buf, &row, "G").unwrap();
        let line = String::from_utf8(buf).unwrap();
        let cols: Vec<&str> = line.trim().split('\t').collect();
        assert_eq!(cols[4], ".");
    }
}
