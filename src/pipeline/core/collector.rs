use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::path::Path;

use anyhow::{Context, Result};

use crate::io::writers::StreamingMatrixWriter;
use crate::io::writers::auxiliary::{
    SORTED_REGIONS_HEADER, write_plain_values_row, write_sorted_region_row,
};
use crate::pipeline::matrix::{MatrixHeader, MatrixRow};

/// Information needed to rewrite values header lines at finalize time.
struct ValuesHeaderInfo {
    /// Byte length of the original header line 1 (including newline).
    line1_byte_len: usize,
    /// Byte length of header line 2 (including newline).
    /// Used to compute the seek offset for line 3 rewriting.
    line2_byte_len: usize,
    /// Byte length of the original header line 3 (including newline).
    line3_byte_len: usize,
}

pub struct FileCollector {
    writer: StreamingMatrixWriter,
    values_writer: Option<BufWriter<File>>,
    regions_writer: Option<BufWriter<File>>,
    group_labels: Vec<String>,
    values_header_info: Option<ValuesHeaderInfo>,
}

impl FileCollector {
    /// Create a new FileCollector with optional auxiliary writers.
    pub fn new(
        writer: StreamingMatrixWriter,
        group_labels: &[String],
        header_estimate: &MatrixHeader,
        sorted_regions_path: Option<&Path>,
        matrix_values_path: Option<&Path>,
    ) -> Result<Self> {
        // Optionally open sorted-regions writer.
        let regions_writer = if let Some(path) = sorted_regions_path {
            let file = File::create(path).with_context(|| {
                format!("Failed to create sorted regions file '{}'", path.display())
            })?;
            let mut w = BufWriter::new(file);
            w.write_all(SORTED_REGIONS_HEADER)?;
            Some(w)
        } else {
            None
        };

        // Optionally open values writer.
        let (values_writer, values_header_info) = if let Some(path) = matrix_values_path {
            let file = File::create(path).with_context(|| {
                format!("Failed to create matrix values file '{}'", path.display())
            })?;
            let mut w = BufWriter::new(file);

            // Write header line 1 using group_capacity counts (placeholder only —
            // will be rewritten at finalize time with actual counts).
            let line1 = build_values_header_line1(header_estimate);
            let line1_len = line1.len();
            w.write_all(line1.as_bytes())?;

            // Write header lines 2 and 3.
            // Line 3 contains group counts and will be rewritten at finalize time.
            let (line2_len, line3_len) = write_values_header_lines_2_3(&mut w, header_estimate)?;

            w.flush()?;

            let info = ValuesHeaderInfo {
                line1_byte_len: line1_len,
                line2_byte_len: line2_len,
                line3_byte_len: line3_len,
            };
            (Some(w), Some(info))
        } else {
            (None, None)
        };

        Ok(Self {
            writer,
            values_writer,
            regions_writer,
            group_labels: group_labels.to_vec(),
            values_header_info,
        })
    }

    /// Discard the underlying writer without finalising the gzip stream.
    pub fn abort(self) {
        self.writer.abort();
    }

    pub fn on_row(&mut self, group_index: usize, row: MatrixRow) -> Result<()> {
        // Write main gzip row.
        self.writer.write_row(&row)?;

        // Write optional plain values row.
        if let Some(ref mut w) = self.values_writer {
            write_plain_values_row(w, &row)?;
        }

        // Write optional sorted region row.
        if let Some(ref mut w) = self.regions_writer {
            let label = self
                .group_labels
                .get(group_index)
                .map(|s| s.as_str())
                .unwrap_or(".");
            write_sorted_region_row(w, &row, label)?;
        }

        Ok(())
    }

    pub fn finalize(mut self, header: MatrixHeader) -> Result<()> {
        // Rewrite values header lines 1 and 3 if needed (group counts may have
        // changed due to --skipZeros / threshold filtering).
        if let (Some(w), Some(info)) = (&mut self.values_writer, &self.values_header_info) {
            w.flush()?;
            let inner = w.get_mut();

            // Rewrite line 1.
            inner.seek(SeekFrom::Start(0))?;
            let actual_line1 = build_values_header_line1(&header);
            let padded = pad_line_to_length(&actual_line1, info.line1_byte_len);
            inner.write_all(padded.as_bytes())?;

            // Rewrite line 3.
            let line3_offset = (info.line1_byte_len + info.line2_byte_len) as u64;
            inner.seek(SeekFrom::Start(line3_offset))?;
            let actual_line3 = build_values_header_line3(&header);
            let padded_line3 = pad_line_to_length(&actual_line3, info.line3_byte_len);
            inner.write_all(padded_line3.as_bytes())?;

            inner.flush()?;
        }

        // Flush auxiliary writers.
        if let Some(mut w) = self.values_writer {
            w.flush()?;
        }
        if let Some(mut w) = self.regions_writer {
            w.flush()?;
        }

        // Finalize main gzip.
        self.writer.finish(&header)
    }
}

/// Build values header line 1: `#Group1:N\tGroup2:M\n`
fn build_values_header_line1(header: &MatrixHeader) -> String {
    let counts = header.group_counts();
    let mut parts = Vec::with_capacity(header.group_labels.len());
    for (label, count) in header.group_labels.iter().zip(counts.iter()) {
        parts.push(format!("{}:{}", label, count));
    }
    format!("#{}\n", parts.join("\t"))
}

/// Pad a line (which must end with `\n`) to exactly `target_len` bytes by
/// inserting spaces before the trailing newline.
fn pad_line_to_length(line: &str, target_len: usize) -> String {
    debug_assert!(line.ends_with('\n'));
    if line.len() >= target_len {
        return line.to_string();
    }
    let padding = target_len - line.len();
    let without_newline = &line[..line.len() - 1];
    let mut result = String::with_capacity(target_len);
    result.push_str(without_newline);
    for _ in 0..padding {
        result.push(' ');
    }
    result.push('\n');
    result
}

/// Write values header lines 2 and 3.  Returns `(line2_byte_len, line3_byte_len)`.
fn write_values_header_lines_2_3<W: Write>(
    writer: &mut W,
    header: &MatrixHeader,
) -> Result<(usize, usize)> {
    let downstream = header.downstream.first().copied().unwrap_or_default();
    let upstream = header.upstream.first().copied().unwrap_or_default();
    let body = header.body.first().copied().unwrap_or_default();
    let bin_size = header.bin_size.first().copied().unwrap_or_default();
    let unscaled5 = header.unscaled_5_prime.first().copied().unwrap_or_default();
    let unscaled3 = header.unscaled_3_prime.first().copied().unwrap_or_default();

    let line2 = format!(
        "#downstream:{}\tupstream:{}\tbody:{}\tbin size:{}\tunscaled 5 prime:{}\tunscaled 3 prime:{}",
        downstream, upstream, body, bin_size, unscaled5, unscaled3
    );
    writeln!(writer, "{}", line2)?;
    let line2_byte_len = line2.len() + 1; // +1 for the newline from writeln!

    // Line 3: group labels with counts, then sample labels repeated per bin.
    let line3 = build_values_header_line3(header);
    writer.write_all(line3.as_bytes())?;
    let line3_byte_len = line3.len();

    Ok((line2_byte_len, line3_byte_len))
}

/// Build header line 3 for the matrix values file.
///
/// Format: `Group1:N\tGroup2:M\tsample1\tsample1\t...\tsample2\tsample2\t...\n`
fn build_values_header_line3(header: &MatrixHeader) -> String {
    let group_counts = header.group_counts();
    let mut parts = Vec::new();
    for (label, count) in header.group_labels.iter().zip(group_counts.iter()) {
        parts.push(format!("{}:{}", label, count));
    }

    let sample_lengths = diff(&header.sample_boundaries);
    for (label, length) in header.sample_labels.iter().zip(sample_lengths.iter()) {
        for _ in 0..*length {
            parts.push(label.clone());
        }
    }
    format!("{}\n", parts.join("\t"))
}

fn diff(values: &[usize]) -> Vec<usize> {
    if values.len() < 2 {
        return Vec::new();
    }
    values.windows(2).map(|pair| pair[1] - pair[0]).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_values_header_line1_basic() {
        let header = MatrixHeader::default_for_test(vec![100, 50]);
        let line = build_values_header_line1(&header);
        assert_eq!(line, "#group0:100\tgroup1:50\n");
    }

    #[test]
    fn pad_line_to_length_exact() {
        let line = "#A:10\n";
        let padded = pad_line_to_length(line, line.len());
        assert_eq!(padded, line);
    }

    #[test]
    fn pad_line_to_length_needs_padding() {
        let line = "#A:5\n";
        let padded = pad_line_to_length(line, 10);
        assert_eq!(padded.len(), 10);
        assert!(padded.ends_with('\n'));
        assert_eq!(&padded[..4], "#A:5");
        // 5 padding spaces + newline
        assert_eq!(&padded[4..9], "     ");
    }

    #[test]
    fn pad_line_to_length_shorter_is_noop() {
        let line = "#A:10000\n";
        let padded = pad_line_to_length(line, 5);
        assert_eq!(padded, line);
    }
}
