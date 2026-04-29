use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::path::Path;

use anyhow::{Context, Result};

use crate::io::writers::StreamingMatrixWriter;
use crate::io::writers::auxiliary::{write_plain_values_row, write_sorted_region_row};
use crate::pipeline::matrix::{MatrixHeader, MatrixRow};

pub trait RowCollector: Send {
    type Output: Send;

    fn on_row(&mut self, group_index: usize, row: MatrixRow) -> Result<()>;
    fn finalize(self, header: MatrixHeader) -> Result<Self::Output>;

    /// Discard this collector without finalising.  The default implementation
    /// simply drops `self`, which is appropriate for in-memory collectors.
    /// File-backed collectors should override this to release I/O resources.
    fn abort(self)
    where
        Self: Sized,
    {
    }
}

/// Information needed to rewrite the values header line 1 at finalize time.
struct ValuesHeaderInfo {
    /// Byte length of the original header line 1 (including newline).
    reserved_byte_len: usize,
}

pub struct FileCollector {
    writer: StreamingMatrixWriter,
    values_writer: Option<BufWriter<File>>,
    regions_writer: Option<BufWriter<File>>,
    group_labels: Vec<String>,
    values_header_info: Option<ValuesHeaderInfo>,
    /// Running per-group row counts (for rewriting values header).
    group_counts: Vec<usize>,
}

impl FileCollector {
    /// Create a new FileCollector with optional auxiliary writers.
    pub fn new(
        writer: StreamingMatrixWriter,
        group_labels: &[String],
        group_capacity: &[usize],
        header_estimate: &MatrixHeader,
        sorted_regions_path: Option<&Path>,
        matrix_values_path: Option<&Path>,
    ) -> Result<Self> {
        let group_count = group_labels.len();

        // Optionally open sorted-regions writer.
        let regions_writer = if let Some(path) = sorted_regions_path {
            let file = File::create(path).with_context(|| {
                format!("Failed to create sorted regions file '{}'", path.display())
            })?;
            let mut w = BufWriter::new(file);
            w.write_all(b"#chrom\tstart\tend\tname\tscore\tstrand\tgroup\n")?;
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

            // Write header line 1 using group_capacity counts (placeholder).
            let line1 = build_values_header_line1(group_labels, group_capacity);
            let reserved_byte_len = line1.len();
            w.write_all(line1.as_bytes())?;

            // Write header lines 2 and 3 (these are final, not rewritten).
            write_values_header_lines_2_3(&mut w, header_estimate)?;

            w.flush()?;

            let info = ValuesHeaderInfo { reserved_byte_len };
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
            group_counts: vec![0usize; group_count],
        })
    }

    /// Discard the underlying writer without finalising the gzip stream.
    pub fn abort(self) {
        self.writer.abort();
    }
}

impl RowCollector for FileCollector {
    type Output = ();

    fn on_row(&mut self, group_index: usize, row: MatrixRow) -> Result<()> {
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

        // Track per-group counts.
        if group_index < self.group_counts.len() {
            self.group_counts[group_index] += 1;
        }

        Ok(())
    }

    fn finalize(mut self, header: MatrixHeader) -> Result<Self::Output> {
        // Rewrite values header line 1 if needed.
        if let (Some(w), Some(info)) = (&mut self.values_writer, &self.values_header_info) {
            w.flush()?;
            let inner = w.get_mut();
            inner.seek(SeekFrom::Start(0))?;

            let actual_line1 = build_values_header_line1(&self.group_labels, &self.group_counts);
            // Pad to reserved length with spaces (before the newline).
            let padded = pad_line_to_length(&actual_line1, info.reserved_byte_len);
            inner.write_all(padded.as_bytes())?;
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

    fn abort(self) {
        self.writer.abort();
    }
}

/// Build values header line 1: `#Group1:N\tGroup2:M\n`
fn build_values_header_line1(group_labels: &[String], group_counts: &[usize]) -> String {
    let mut parts = Vec::with_capacity(group_labels.len());
    for (label, &count) in group_labels.iter().zip(group_counts.iter()) {
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

/// Write values header lines 2 and 3 (these depend only on the header estimate
/// and are not rewritten later).
fn write_values_header_lines_2_3<W: Write>(writer: &mut W, header: &MatrixHeader) -> Result<()> {
    let downstream = header.downstream.first().copied().unwrap_or_default();
    let upstream = header.upstream.first().copied().unwrap_or_default();
    let body = header.body.first().copied().unwrap_or_default();
    let bin_size = header.bin_size.first().copied().unwrap_or_default();
    let unscaled5 = header.unscaled_5_prime.first().copied().unwrap_or_default();
    let unscaled3 = header.unscaled_3_prime.first().copied().unwrap_or_default();

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
        let labels = vec!["genes".to_string(), "peaks".to_string()];
        let counts = vec![100, 50];
        let line = build_values_header_line1(&labels, &counts);
        assert_eq!(line, "#genes:100\tpeaks:50\n");
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
