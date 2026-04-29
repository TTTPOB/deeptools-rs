use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};
use flate2::write::GzEncoder;
use flate2::{Compression, GzBuilder};

use crate::pipeline::matrix::{MatrixHeader, MatrixRow};
use super::formatting::write_matrix_row;

const RESERVED_HEADER_COMPRESSED: usize = 8192;
const RESERVED_HEADER_PAYLOAD: usize = RESERVED_HEADER_COMPRESSED - 23;

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
