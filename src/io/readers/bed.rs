use std::cell::RefCell;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::sync::Arc;

use thiserror::Error;

thread_local! {
    static CHROM_INTERNER: RefCell<HashMap<String, Arc<str>>> = RefCell::new(HashMap::new());
}

pub(crate) fn intern_chrom(s: String) -> Arc<str> {
    CHROM_INTERNER.with(|map| {
        let mut map = map.borrow_mut();
        if let Some(existing) = map.get(&s) {
            Arc::clone(existing)
        } else {
            let arc: Arc<str> = Arc::from(s.as_str());
            map.insert(s, arc.clone());
            arc
        }
    })
}

fn parse_comma_separated_u32(value: &str, expected: usize) -> Option<Vec<u32>> {
    let mut numbers = Vec::with_capacity(expected);
    for part in value.split(',') {
        if part.is_empty() {
            continue;
        }
        let parsed = part.trim().parse::<u32>().ok()?;
        numbers.push(parsed);
    }
    if numbers.len() != expected {
        return None;
    }
    Some(numbers)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Positive,
    Negative,
    Unstranded,
}

impl Strand {
    pub fn from_symbol(symbol: &str) -> Option<Self> {
        match symbol {
            "+" => Some(Strand::Positive),
            "-" => Some(Strand::Negative),
            "." | "?" | "*" => Some(Strand::Unstranded),
            _ => None,
        }
    }

    pub fn as_char(self) -> char {
        match self {
            Strand::Positive => '+',
            Strand::Negative => '-',
            Strand::Unstranded => '.',
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct BedRecord {
    pub chrom: Arc<str>,
    pub start: u32,
    pub end: u32,
    pub name: Option<String>,
    pub score: Option<f32>,
    pub score_raw: Option<String>,
    pub strand: Strand,
    pub strand_raw: Option<String>,
    pub extra_fields: Vec<String>,
}

impl BedRecord {
    pub fn parse(line: &str) -> Result<Self, String> {
        let mut fields: Vec<&str> = if line.contains('\t') {
            line.split('\t').collect()
        } else {
            line.split_whitespace().collect()
        };

        // remove potential empty trailing fields from trailing tabs
        if let Some(last) = fields.last() {
            if last.is_empty() {
                fields.pop();
            }
        }

        if fields.len() < 3 {
            return Err("BED line requires at least 3 columns".to_string());
        }

        let chrom = intern_chrom(fields[0].to_string());
        let start = fields[1]
            .parse::<u32>()
            .map_err(|_| "BED start column must be an unsigned integer".to_string())?;
        let end = fields[2]
            .parse::<u32>()
            .map_err(|_| "BED end column must be an unsigned integer".to_string())?;

        if end < start {
            return Err("BED end column must be greater than or equal to start".to_string());
        }

        let name = fields
            .get(3)
            .map(|value| value.to_string())
            .filter(|v| !v.is_empty());
        let (score, score_raw) = if let Some(raw) = fields.get(4) {
            if raw.is_empty() || *raw == "." {
                (None, None)
            } else if let Ok(parsed) = raw.parse::<f32>() {
                (Some(parsed), None)
            } else {
                (None, Some(raw.to_string()))
            }
        } else {
            (None, None)
        };

        let (strand, strand_raw) = if let Some(raw) = fields.get(5) {
            if let Some(parsed) = Strand::from_symbol(raw) {
                (parsed, None)
            } else {
                (Strand::Unstranded, Some(raw.to_string()))
            }
        } else {
            (Strand::Unstranded, None)
        };

        let extra_fields = if fields.len() > 6 {
            fields[6..].iter().map(|value| value.to_string()).collect()
        } else {
            Vec::new()
        };

        Ok(Self {
            chrom,
            start,
            end,
            name,
            score,
            score_raw,
            strand,
            strand_raw,
            extra_fields,
        })
    }

    pub fn length(&self) -> u32 {
        self.end - self.start
    }

    pub fn exons(&self) -> Option<Vec<(u32, u32)>> {
        if self.extra_fields.len() < 6 {
            return None;
        }

        let block_count = self.extra_fields[3].parse::<usize>().ok()?;
        if block_count == 0 {
            return Some(Vec::new());
        }

        let block_sizes = parse_comma_separated_u32(&self.extra_fields[4], block_count)?;
        let block_starts = parse_comma_separated_u32(&self.extra_fields[5], block_count)?;

        let mut exons = Vec::with_capacity(block_count);
        for (size, start_offset) in block_sizes.into_iter().zip(block_starts.into_iter()) {
            let exon_start = self.start.checked_add(start_offset)?;
            let exon_end = exon_start.checked_add(size)?;
            if exon_end <= exon_start {
                return None;
            }
            exons.push((exon_start, exon_end));
        }

        Some(exons)
    }
}

#[derive(Debug, Error)]
pub enum BedReadError {
    #[error("I/O error while reading BED file: {0}")]
    Io(#[from] io::Error),
    #[error("line {line_number}: {message}\n  \u{2514} raw: {line}")]
    Parse {
        line_number: usize,
        message: String,
        line: String,
    },
}

pub struct BedReader<R: BufRead> {
    lines: std::iter::Enumerate<io::Lines<R>>,
}

impl BedReader<BufReader<File>> {
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self, BedReadError> {
        let file = File::open(path)?;
        Ok(Self::new(BufReader::new(file)))
    }
}

impl<R: BufRead> BedReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            lines: reader.lines().enumerate(),
        }
    }

    pub fn read_all(mut self) -> Result<Vec<BedRecord>, BedReadError> {
        let mut records = Vec::new();
        while let Some(record) = self.next() {
            records.push(record?);
        }
        Ok(records)
    }
}

impl<R: BufRead> Iterator for BedReader<R> {
    type Item = Result<BedRecord, BedReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some((index, line)) = self.lines.next() {
            let line_number = index + 1;
            let line = match line {
                Ok(value) => value,
                Err(err) => return Some(Err(BedReadError::Io(err))),
            };

            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            match BedRecord::parse(trimmed) {
                Ok(record) => return Some(Ok(record)),
                Err(message) => {
                    return Some(Err(BedReadError::Parse {
                        line_number,
                        message,
                        line,
                    }));
                }
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufReader, Cursor};

    #[test]
    fn extracts_bed12_exons() {
        let record =
            BedRecord::parse("chr1\t100\t500\tname\t5.0\t+\t100\t500\t0\t2\t150,200,\t0,200,")
                .expect("should parse");
        let exons = record.exons().expect("exons");
        assert_eq!(exons.len(), 2);
        assert_eq!(exons[0], (100, 250));
        assert_eq!(exons[1], (300, 500));
    }

    #[test]
    fn parses_basic_bed_line() {
        let record = BedRecord::parse("chr1\t10\t20\tname\t5.0\t+").expect("should parse");
        assert_eq!(&*record.chrom, "chr1");
        assert_eq!(record.start, 10);
        assert_eq!(record.end, 20);
        assert_eq!(record.name.as_deref(), Some("name"));
        assert_eq!(record.score, Some(5.0));
        assert!(record.score_raw.is_none());
        assert_eq!(record.strand, Strand::Positive);
        assert!(record.strand_raw.is_none());
        assert!(record.extra_fields.is_empty());
    }

    #[test]
    fn reader_skips_comments_and_blank_lines() {
        let data = b"# comment\n\nchr2\t0\t50\n";
        let reader = BedReader::new(BufReader::new(Cursor::new(&data[..])));
        let records: Vec<_> = reader.collect::<Result<_, _>>().expect("valid records");
        assert_eq!(records.len(), 1);
        assert_eq!(&*records[0].chrom, "chr2");
    }

    #[test]
    fn keeps_invalid_score_and_strand_as_raw_strings() {
        let record =
            BedRecord::parse("chr1\t10\t20\tname\tabc\tstrandx").expect("should parse");
        assert!(record.score.is_none());
        assert_eq!(record.score_raw.as_deref(), Some("abc"));
        assert_eq!(record.strand, Strand::Unstranded);
        assert_eq!(record.strand_raw.as_deref(), Some("strandx"));
    }
}
