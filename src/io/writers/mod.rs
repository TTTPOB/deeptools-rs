use std::fs::File;
use std::io::{self, BufReader, BufWriter, Seek, SeekFrom, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};
use flate2::write::GzEncoder;
use flate2::{Compression, GzBuilder};
use tempfile::{NamedTempFile, TempPath};

use crate::config::IoOptions;
use crate::pipeline::matrix::{MatrixData, MatrixHeader, MatrixRow};

const STREAMING_CELL_THRESHOLD: usize = 100_000;
const RESERVED_HEADER_COMPRESSED: usize = 8192;
const RESERVED_HEADER_PAYLOAD: usize = RESERVED_HEADER_COMPRESSED - 23;

/// Persist all requested outputs derived from the matrix computation.
pub fn write_outputs(mut matrix: MatrixData, io: &IoOptions) -> Result<()> {
    if should_use_streaming(&matrix, io) {
        write_matrix_gz_streaming(&io.matrix_output, &mut matrix)?;
        return Ok(());
    }

    write_matrix_gz(&io.matrix_output, &matrix)?;

    if let Some(path) = &io.matrix_values_output {
        write_matrix_values(path, &matrix)?;
    }

    if let Some(path) = &io.sorted_regions_output {
        write_sorted_regions(path, &matrix)?;
    }

    Ok(())
}

fn should_use_streaming(matrix: &MatrixData, io: &IoOptions) -> bool {
    should_use_streaming_for_plan(
        matrix.rows.len(),
        matrix.sample_count,
        matrix.bin_count,
        matrix.header.sort_regions == "keep",
        io,
    )
}

pub fn should_use_streaming_for_plan(
    row_count: usize,
    sample_count: usize,
    bin_count: usize,
    sort_is_keep: bool,
    io: &IoOptions,
) -> bool {
    if io.matrix_values_output.is_some() || io.sorted_regions_output.is_some() {
        return false;
    }

    if !sort_is_keep {
        return false;
    }

    if row_count == 0 {
        return false;
    }

    let cell_count = row_count
        .saturating_mul(sample_count)
        .saturating_mul(bin_count);

    cell_count >= STREAMING_CELL_THRESHOLD
}

fn write_matrix_gz(path: &Path, matrix: &MatrixData) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Failed to create matrix file '{}'", path.display()))?;
    let mut encoder = GzEncoder::new(file, Compression::default());

    write_header_line(&mut encoder, &matrix.header)?;

    for row in &matrix.rows {
        write_matrix_row(&mut encoder, row)?;
    }

    encoder.finish()?;
    Ok(())
}

fn write_matrix_gz_streaming(path: &Path, matrix: &mut MatrixData) -> Result<()> {
    let final_header_payload = match build_padded_header_payload(&matrix.header) {
        Ok(payload) => payload,
        Err(_) => {
            return write_matrix_gz(path, matrix);
        }
    };

    if matrix.rows.is_empty() {
        let file = File::create(path)
            .with_context(|| format!("Failed to create matrix file '{}'", path.display()))?;
        write_header_member(file, &final_header_payload)?;
        return Ok(());
    }

    let rows = std::mem::take(&mut matrix.rows);
    let spool_path = spool_rows(rows)?;

    let placeholder_payload =
        placeholder_header_payload().expect("placeholder header payload should fit reserved size");

    let mut file = File::create(path)
        .with_context(|| format!("Failed to create matrix file '{}'", path.display()))?;

    file.seek(SeekFrom::Current(0))
        .context("Streaming output requires a seekable destination")?;

    file = write_header_member(file, &placeholder_payload)?;

    {
        let spool_file = File::open(spool_path.as_ref() as &Path)
            .context("Failed to reopen temporary matrix stream")?;
        let mut reader = BufReader::new(spool_file);
        let builder = GzBuilder::new().mtime(0);
        let mut encoder = builder.write(file, Compression::default());
        io::copy(&mut reader, &mut encoder)
            .context("Failed to stream matrix rows into gzip writer")?;
        file = encoder
            .finish()
            .context("Failed to finalise streamed matrix member")?;
    }

    let _ = rewrite_header_member(file, &final_header_payload)?;

    spool_path.close().ok();

    Ok(())
}

fn build_header_line_from_header(header: &MatrixHeader) -> Result<Vec<u8>> {
    let header = serde_json::to_string(header)?;
    let mut line = Vec::with_capacity(header.len() + 2);
    line.push(b'@');
    line.extend_from_slice(header.as_bytes());
    line.push(b'\n');
    Ok(line)
}

fn pad_header_payload(mut data: Vec<u8>) -> Result<Vec<u8>> {
    if RESERVED_HEADER_PAYLOAD == 0 {
        bail!("reserved header payload must be greater than zero");
    }

    let Some(last) = data.pop() else {
        bail!("header payload must include a trailing newline");
    };

    if last != b'\n' {
        bail!("header payload must end with a newline");
    }

    if data.len() + 1 > RESERVED_HEADER_PAYLOAD {
        bail!(
            "header payload of {} bytes exceeds reserved capacity of {} bytes",
            data.len() + 1,
            RESERVED_HEADER_PAYLOAD
        );
    }

    data.resize(RESERVED_HEADER_PAYLOAD - 1, b' ');
    data.push(b'\n');
    Ok(data)
}

fn placeholder_header_payload() -> Result<Vec<u8>> {
    pad_header_payload(b"@{}\n".to_vec())
}

pub fn build_padded_header_payload(header: &MatrixHeader) -> Result<Vec<u8>> {
    let line = build_header_line_from_header(header)?;
    pad_header_payload(line)
}

pub fn ensure_streaming_header_capacity(header: &MatrixHeader) -> Result<()> {
    let _ = build_padded_header_payload(header)?;
    Ok(())
}

fn write_header_member(file: File, payload: &[u8]) -> Result<File> {
    let builder = GzBuilder::new().mtime(0);
    let mut encoder = builder.write(file, Compression::none());
    encoder
        .write_all(payload)
        .context("Failed to write header member payload")?;
    encoder
        .finish()
        .context("Failed to finalise header member stream")
}

fn rewrite_header_member(mut file: File, payload: &[u8]) -> Result<File> {
    file.seek(SeekFrom::Start(0))
        .context("Failed to seek to start of matrix file for header rewrite")?;

    let builder = GzBuilder::new().mtime(0);
    let mut encoder = builder.write(file, Compression::none());
    encoder
        .write_all(payload)
        .context("Failed to rewrite header member payload")?;
    let mut file = encoder
        .finish()
        .context("Failed to finalise rewritten header member")?;

    file.seek(SeekFrom::End(0))
        .context("Failed to restore file cursor after header rewrite")?;

    Ok(file)
}

fn write_header_line<W: Write>(writer: &mut W, header: &MatrixHeader) -> Result<()> {
    let line = build_header_line_from_header(header)?;
    writer
        .write_all(&line)
        .context("Failed to write matrix header line")?;
    Ok(())
}

pub struct StreamingMatrixWriter {
    encoder: GzEncoder<File>,
}

impl StreamingMatrixWriter {
    pub fn start(path: &Path) -> Result<Self> {
        let placeholder_payload = placeholder_header_payload()
            .expect("placeholder header payload should fit reserved size");

        let mut file = File::create(path)
            .with_context(|| format!("Failed to create matrix file '{}'", path.display()))?;

        file.seek(SeekFrom::Current(0))
            .context("Streaming output requires a seekable destination")?;

        file = write_header_member(file, &placeholder_payload)?;

        let builder = GzBuilder::new().mtime(0);
        let encoder = builder.write(file, Compression::default());

        Ok(Self { encoder })
    }

    pub fn write_row(&mut self, row: &MatrixRow) -> Result<()> {
        write_matrix_row(&mut self.encoder, row)
    }

    pub fn finish(self, header: &MatrixHeader) -> Result<()> {
        let final_payload = build_padded_header_payload(header)?;
        let file = self
            .encoder
            .finish()
            .context("Failed to finalise streamed matrix member")?;
        let _ = rewrite_header_member(file, &final_payload)?;
        Ok(())
    }
}

fn spool_rows(rows: Vec<MatrixRow>) -> Result<TempPath> {
    let mut temp = NamedTempFile::new().context("Failed to allocate temporary matrix stream")?;
    {
        let mut writer = BufWriter::new(temp.as_file_mut());
        for row in &rows {
            write_matrix_row(&mut writer, row)?;
        }
        writer
            .flush()
            .context("Failed to flush temporary matrix stream")?;
    }

    drop(rows);

    Ok(temp.into_temp_path())
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
