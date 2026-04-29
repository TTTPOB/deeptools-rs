use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use flate2::write::GzEncoder;
use flate2::{Compression, GzBuilder};

use super::formatting::write_matrix_row;
use crate::pipeline::matrix::{MatrixHeader, MatrixRow};

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
    /// Wrapped in `Option` so `finish()` can take ownership of the encoder
    /// without violating the `Drop` implementation's partial-move restriction.
    encoder: Option<GzEncoder<BufWriter<File>>>,
    /// Path to the output file.  `None` after a successful `finish()`, so
    /// `Drop` only removes the file when the writer is abandoned.
    path: Option<PathBuf>,
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
        let encoder = builder.write(BufWriter::with_capacity(131_072, file), Compression::fast());

        Ok(Self {
            encoder: Some(encoder),
            path: Some(path.to_path_buf()),
        })
    }

    pub fn write_row(&mut self, row: &MatrixRow) -> Result<()> {
        write_matrix_row(self.encoder.as_mut().expect("encoder present"), row)
    }

    pub fn finish(mut self, header: &MatrixHeader) -> Result<()> {
        let final_payload = build_padded_header_payload(header)?;

        let encoder = self.encoder.take().expect("encoder present");
        let writer = encoder
            .finish()
            .context("Failed to finalise streamed matrix member")?;
        let file = writer
            .into_inner()
            .map_err(|err| err.into_error())
            .context("Failed to flush buffered streamed matrix member")?;
        let _ = rewrite_header_member(file, &final_payload)?;

        // Disarm the Drop guard only after all fallible steps succeed,
        // so that a partial matrix file is still cleaned up on failure.
        let _path = self.path.take();
        Ok(())
    }

    /// Discard the in-progress writer without finalising.  The output file is
    /// removed from disk.
    pub fn abort(self) {
        // Cleanup happens in Drop impl.
        drop(self);
    }
}

impl Drop for StreamingMatrixWriter {
    fn drop(&mut self) {
        if let Some(ref path) = self.path {
            let _ = std::fs::remove_file(path);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pipeline::matrix::MatrixHeader;

    // --- pad_header_payload ---

    #[test]
    fn pad_header_payload_normal_case() {
        // Input ends with \n; result must be exactly RESERVED_HEADER_PAYLOAD bytes.
        let input = b"@{}\n".to_vec();
        let result = pad_header_payload(input).expect("should succeed");
        assert_eq!(result.len(), RESERVED_HEADER_PAYLOAD);
    }

    #[test]
    fn pad_header_payload_result_ends_with_newline() {
        let input = b"@{}\n".to_vec();
        let result = pad_header_payload(input).expect("should succeed");
        assert_eq!(*result.last().unwrap(), b'\n');
    }

    #[test]
    fn pad_header_payload_padding_uses_spaces() {
        let input = b"@{}\n".to_vec();
        let result = pad_header_payload(input).expect("should succeed");
        // All bytes between the content and the trailing newline must be spaces.
        let content_end = 3; // b"@{}" is 3 bytes
        let padding = &result[content_end..result.len() - 1];
        assert!(padding.iter().all(|&b| b == b' '));
    }

    #[test]
    fn pad_header_payload_without_trailing_newline_is_error() {
        let input = b"@{}".to_vec(); // no trailing newline
        let err = pad_header_payload(input).unwrap_err();
        assert!(err.to_string().contains("newline"));
    }

    #[test]
    fn pad_header_payload_empty_input_is_error() {
        let input = vec![];
        let err = pad_header_payload(input).unwrap_err();
        assert!(err.to_string().contains("newline"));
    }

    #[test]
    fn pad_header_payload_oversized_input_is_error() {
        // Build a payload that is one byte larger than the reserved capacity.
        let oversized: Vec<u8> = std::iter::repeat(b'x')
            .take(RESERVED_HEADER_PAYLOAD) // RESERVED_HEADER_PAYLOAD - 1 'x' + '\n' = too big
            .chain(std::iter::once(b'\n'))
            .collect();
        let err = pad_header_payload(oversized).unwrap_err();
        assert!(err.to_string().contains("exceeds reserved capacity"));
    }

    // --- placeholder_header_payload ---

    #[test]
    fn placeholder_header_payload_correct_length() {
        let result = placeholder_header_payload().expect("should succeed");
        assert_eq!(result.len(), RESERVED_HEADER_PAYLOAD);
    }

    #[test]
    fn placeholder_header_payload_starts_with_at_brace() {
        let result = placeholder_header_payload().expect("should succeed");
        assert!(result.starts_with(b"@{}"));
    }

    // --- build_padded_header_payload ---

    #[test]
    fn build_padded_header_payload_correct_length() {
        let header = MatrixHeader::default_for_test(vec![3]);
        let result = build_padded_header_payload(&header).expect("should succeed");
        assert_eq!(result.len(), RESERVED_HEADER_PAYLOAD);
    }

    #[test]
    fn build_padded_header_payload_ends_with_newline() {
        let header = MatrixHeader::default_for_test(vec![3]);
        let result = build_padded_header_payload(&header).expect("should succeed");
        assert_eq!(*result.last().unwrap(), b'\n');
    }

    // --- ensure_streaming_header_capacity ---

    #[test]
    fn ensure_streaming_header_capacity_normal_header_ok() {
        let header = MatrixHeader::default_for_test(vec![3]);
        ensure_streaming_header_capacity(&header).expect("normal header should fit");
    }
}
