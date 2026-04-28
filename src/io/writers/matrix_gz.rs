use std::fs::File;
use std::io::{self, BufReader, BufWriter, Seek, SeekFrom, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};
use flate2::write::GzEncoder;
use flate2::{Compression, GzBuilder};
use tempfile::{NamedTempFile, TempPath};

use crate::pipeline::matrix::{MatrixData, MatrixHeader, MatrixRow};
use super::formatting::write_matrix_row;

pub const STREAMING_CELL_THRESHOLD: usize = 100_000;
const RESERVED_HEADER_COMPRESSED: usize = 8192;
const RESERVED_HEADER_PAYLOAD: usize = RESERVED_HEADER_COMPRESSED - 23;

pub fn write_matrix_gz(path: &Path, matrix: &MatrixData) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Failed to create matrix file '{}'", path.display()))?;
    let mut encoder = GzEncoder::new(
        BufWriter::with_capacity(131_072, file),
        Compression::fast(),
    );

    write_header_line(&mut encoder, &matrix.header)?;

    for row in &matrix.rows {
        write_matrix_row(&mut encoder, row)?;
    }

    let writer = encoder
        .finish()
        .context("Failed to finalise matrix gzip stream")?;
    writer
        .into_inner()
        .map_err(|err| err.into_error())
        .context("Failed to flush buffered matrix gzip stream")?;
    Ok(())
}

pub fn write_matrix_gz_streaming(path: &Path, matrix: &mut MatrixData) -> Result<()> {
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
        let mut encoder = builder.write(
            BufWriter::with_capacity(131_072, file),
            Compression::fast(),
        );
        io::copy(&mut reader, &mut encoder)
            .context("Failed to stream matrix rows into gzip writer")?;
        let writer = encoder
            .finish()
            .context("Failed to finalise streamed matrix member")?;
        file = writer
            .into_inner()
            .map_err(|err| err.into_error())
            .context("Failed to flush buffered streamed matrix member")?;
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
    encoder: GzEncoder<BufWriter<File>>,
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
        let encoder = builder.write(
            BufWriter::with_capacity(131_072, file),
            Compression::fast(),
        );

        Ok(Self { encoder })
    }

    pub fn write_row(&mut self, row: &MatrixRow) -> Result<()> {
        write_matrix_row(&mut self.encoder, row)
    }

    pub fn finish(self, header: &MatrixHeader) -> Result<()> {
        let final_payload = build_padded_header_payload(header)?;
        let writer = self
            .encoder
            .finish()
            .context("Failed to finalise streamed matrix member")?;
        let file = writer
            .into_inner()
            .map_err(|err| err.into_error())
            .context("Failed to flush buffered streamed matrix member")?;
        let _ = rewrite_header_member(file, &final_payload)?;
        Ok(())
    }

    /// Discard the in-progress writer without finalising.  The output file is
    /// left in a corrupt/partial state; callers are responsible for removing it
    /// if desired.  Dropping the encoder closes the underlying file handle.
    pub fn abort(self) {
        drop(self.encoder);
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
