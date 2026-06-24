use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use serde_json::Value;

/// A single row parsed from the matrix file (BED fields + numeric values).
#[derive(Debug, Clone)]
pub struct MatrixFileRow {
    pub chrom: String,
    pub start: String,
    pub end: String,
    pub name: String,
    pub score: String,
    pub strand: String,
    pub values: Vec<f64>,
}

impl MatrixFileRow {
    /// Returns a tuple key suitable for row identity checks.
    pub fn key(&self) -> String {
        format!(
            "{}:{}:{}:{}:{}:{}",
            self.chrom, self.start, self.end, self.name, self.score, self.strand
        )
    }
}

/// In-memory representation of a parsed matrix file.
#[derive(Debug)]
pub struct Matrix {
    pub header_json: Value,
    pub rows: Vec<MatrixFileRow>,
}

/// Load a matrix file from `path`. Detects plain text vs gzip (including
/// multi-member gzip produced by the streaming writer).
pub fn load_matrix(path: &Path) -> Result<Matrix> {
    let mut file =
        File::open(path).with_context(|| format!("Failed to open '{}'", path.display()))?;

    let mut buf = [0u8; 2];
    let n = file.read(&mut buf).context("Failed to read magic bytes")?;
    file.seek(SeekFrom::Start(0))
        .context("Failed to seek back to start")?;

    let reader: Box<dyn BufRead> = if n == 2 && buf == [0x1f, 0x8b] {
        // gzip (possibly multi-member — MultiGzDecoder handles both)
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        // plain text
        Box::new(BufReader::new(file))
    };

    parse_matrix_from_reader(reader)
}

fn parse_matrix_from_reader(mut reader: Box<dyn BufRead>) -> Result<Matrix> {
    // First line must be the header: @{JSON...}
    let mut header_line = String::new();
    reader
        .read_line(&mut header_line)
        .context("Failed to read header line")?;

    let header_line = header_line.trim_end_matches('\n').trim_end_matches('\r');
    if !header_line.starts_with('@') {
        bail!("Matrix file does not start with '@' header line");
    }

    let json_str = &header_line[1..];
    // Handle padded header: strip trailing spaces before the newline was removed
    let json_str = json_str.trim_end();
    let header_json: Value =
        serde_json::from_str(json_str).context("Failed to parse header JSON")?;

    let mut rows = Vec::new();
    let mut line_buf = String::new();
    let mut line_number = 1usize;

    loop {
        line_buf.clear();
        let n = reader
            .read_line(&mut line_buf)
            .with_context(|| format!("Failed to read line {}", line_number + 1))?;
        if n == 0 {
            break; // EOF
        }
        line_number += 1;

        let line = line_buf.trim_end_matches('\n').trim_end_matches('\r');
        if line.is_empty() {
            continue;
        }

        let row = parse_row(line, line_number)?;
        rows.push(row);
    }

    Ok(Matrix { header_json, rows })
}

fn parse_row(line: &str, line_number: usize) -> Result<MatrixFileRow> {
    let mut fields = line.splitn(7, '\t');

    macro_rules! next_field {
        ($name:expr) => {
            fields
                .next()
                .with_context(|| format!("Missing field '{}' on line {}", $name, line_number))?
                .to_owned()
        };
    }

    let chrom = next_field!("chrom");
    let start = next_field!("start");
    let end = next_field!("end");
    let name = next_field!("name");
    let score = next_field!("score");
    let strand = next_field!("strand");
    // Seventh capture is everything remaining (all value columns)
    let values_str = fields.next().unwrap_or("").trim_end();

    let values: Vec<f64> = if values_str.is_empty() {
        Vec::new()
    } else {
        values_str
            .split('\t')
            .enumerate()
            .map(|(col_idx, s)| {
                parse_value(s).with_context(|| {
                    format!(
                        "Failed to parse value column {} ('{}') on line {}",
                        col_idx + 1,
                        s,
                        line_number
                    )
                })
            })
            .collect::<Result<Vec<_>>>()?
    };

    Ok(MatrixFileRow {
        chrom,
        start,
        end,
        name,
        score,
        strand,
        values,
    })
}

/// Parse a single value token: numeric, "nan", "inf", "-inf".
fn parse_value(s: &str) -> Result<f64> {
    let s = s.trim();
    match s {
        "nan" | "NaN" | "NAN" => Ok(f64::NAN),
        "inf" | "Inf" | "INF" => Ok(f64::INFINITY),
        "-inf" | "-Inf" | "-INF" => Ok(f64::NEG_INFINITY),
        other => other
            .parse::<f64>()
            .with_context(|| format!("Cannot parse '{}' as a floating-point number", other)),
    }
}
