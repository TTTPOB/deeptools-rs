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
    pub bed_field_count: Option<usize>,
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
        let start_raw = fields[1]
            .parse::<i64>()
            .map_err(|_| "BED start column must be an integer".to_string())?;
        let start = if start_raw < 0 {
            0
        } else {
            u32::try_from(start_raw).map_err(|_| "BED start column overflowed u32".to_string())?
        };
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
            bed_field_count: Some(fields.len()),
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
    #[error("no data records found in BED file")]
    EmptyFile,
}

/// A named group of BED records, typically delimited by a `#Label` line.
#[derive(Debug, Clone)]
pub struct Group {
    pub label: String,
    pub records: Vec<BedRecord>,
}

/// Reads a BED file and yields [`Group`]s delimited by `#` comment lines.
///
/// A `#Label` line terminates and names the **previous** accumulated group
/// (trailing-delimiter style).  Groups that reach EOF without a preceding `#`
/// line receive `default_label`.  When a `#` line provides an empty label the
/// `default_label` is also used as a fallback.
///
/// Labels are emitted **raw** — no cross-file deduplication is performed.
/// Callers should apply deduplication when merging groups from multiple files.
pub struct GroupedBedReader<R: BufRead> {
    lines: std::iter::Enumerate<io::Lines<R>>,
    default_label: String,
    /// Records accumulated since the last group boundary.
    current_records: Vec<BedRecord>,
    /// Whether at least one group has been yielded (for `EmptyFile` detection).
    yielded_any: bool,
    /// Whether the underlying line iterator has been exhausted.
    done: bool,
}

impl GroupedBedReader<BufReader<File>> {
    /// Open a BED file at `path` and return a reader that yields groups.
    pub fn open(path: impl AsRef<Path>, default_label: String) -> Result<Self, BedReadError> {
        let file = File::open(path)?;
        Ok(Self::new(BufReader::new(file), default_label))
    }
}

impl<R: BufRead> GroupedBedReader<R> {
    /// Create a new reader from any buffered reader.
    ///
    /// Useful for unit-testing grouping behavior without filesystem setup.
    pub fn new(reader: R, default_label: String) -> Self {
        Self {
            lines: reader.lines().enumerate(),
            default_label,
            current_records: Vec::new(),
            yielded_any: false,
            done: false,
        }
    }
}

impl<R: BufRead> Iterator for GroupedBedReader<R> {
    type Item = Result<Group, BedReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.done {
                return None;
            }

            match self.lines.next() {
                Some((idx, Ok(line))) => {
                    let line_number = idx + 1;
                    let trimmed = line.trim();
                    if trimmed.is_empty() {
                        continue;
                    }
                    if trimmed.starts_with('#') {
                        let raw_label = trimmed.strip_prefix('#').unwrap_or("").trim();
                        if !self.current_records.is_empty() {
                            let records = std::mem::take(&mut self.current_records);
                            let label = if raw_label.is_empty() {
                                self.default_label.clone()
                            } else {
                                raw_label.to_string()
                            };
                            self.yielded_any = true;
                            return Some(Ok(Group { label, records }));
                        }
                        // No preceding records: discard label, keep
                        // accumulating.
                        continue;
                    }
                    // Parse as a BED data line.
                    match BedRecord::parse(trimmed) {
                        Ok(record) => self.current_records.push(record),
                        Err(message) => {
                            self.done = true;
                            return Some(Err(BedReadError::Parse {
                                line_number,
                                message,
                                line,
                            }));
                        }
                    }
                }
                Some((_idx, Err(err))) => {
                    self.done = true;
                    return Some(Err(BedReadError::Io(err)));
                }
                None => {
                    self.done = true;
                    // Emit the final group if records remain.
                    if !self.current_records.is_empty() {
                        let records = std::mem::take(&mut self.current_records);
                        self.yielded_any = true;
                        return Some(Ok(Group {
                            label: self.default_label.clone(),
                            records,
                        }));
                    }
                    if !self.yielded_any {
                        return Some(Err(BedReadError::EmptyFile));
                    }
                    return None;
                }
            }
        }
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
    fn clamps_negative_bed_start_to_zero() {
        let record = BedRecord::parse("chr1\t-5\t20\tname\t5.0\t+").expect("should parse");
        assert_eq!(record.start, 0);
        assert_eq!(record.end, 20);
    }

    // --- GroupedBedReader tests ---

    fn make_reader(data: &[u8]) -> GroupedBedReader<BufReader<Cursor<&[u8]>>> {
        GroupedBedReader::new(BufReader::new(Cursor::new(data)), "default".to_string())
    }

    #[test]
    fn groups_with_trailing_delimiter() {
        // `#Group 1` terminates the first three records; `#Group 2`
        // terminates the next three.  This matches the test2.bed fixture.
        let data = b"\
ch1\t100\t150\tCG11023\t0\t+\n\
ch2\t150\t175\tcda5\t0\t-\n\
ch3\t100\t125\tcda8\t0\t+\n\
#Group 1\n\
ch1\t75\t125\tC11023\t0\t+\n\
ch2\t125\t150\tca5\t0\t-\n\
ch3\t75\t100\tca8\t0\t+\n\
#Group 2\n\
";
        let groups: Vec<_> = make_reader(data)
            .collect::<Result<_, _>>()
            .expect("valid groups");
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].label, "Group 1");
        assert_eq!(groups[0].records.len(), 3);
        assert_eq!(groups[1].label, "Group 2");
        assert_eq!(groups[1].records.len(), 3);
    }

    #[test]
    fn single_group_no_delimiters_gets_default_label() {
        let data = b"chr1\t100\t200\nchr1\t300\t400\n";
        let groups: Vec<_> = make_reader(data)
            .collect::<Result<_, _>>()
            .expect("valid groups");
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].label, "default");
        assert_eq!(groups[0].records.len(), 2);
    }

    #[test]
    fn eof_group_gets_default_label() {
        let data = b"chr1\t100\t200\n#Label1\nchr2\t300\t400\n";
        let groups: Vec<_> = make_reader(data)
            .collect::<Result<_, _>>()
            .expect("valid groups");
        assert_eq!(groups.len(), 2);
        // First record(s) terminated by #Label1.
        assert_eq!(groups[0].label, "Label1");
        // Remaining record(s) at EOF get default_label.
        assert_eq!(groups[1].label, "default");
    }

    #[test]
    fn hash_with_no_preceding_records_discards_label() {
        // File starts with `#` — no records precede it, so label is discarded.
        let data = b"#Orphan\nchr1\t100\t200\n#Real\nchr2\t300\t400\n";
        let groups: Vec<_> = make_reader(data)
            .collect::<Result<_, _>>()
            .expect("valid groups");
        assert_eq!(groups.len(), 2);
        // First record batch is terminated by #Real (NOT #Orphan).
        assert_eq!(groups[0].label, "Real");
        // Final batch at EOF.
        assert_eq!(groups[1].label, "default");
    }

    #[test]
    fn consecutive_hash_lines_no_intervening_records() {
        let data = b"chr1\t100\t200\n#A\n#B\nchr2\t300\t400\n";
        let groups: Vec<_> = make_reader(data)
            .collect::<Result<_, _>>()
            .expect("valid groups");
        assert_eq!(groups.len(), 2);
        // chr1 is terminated by #A. #B has no preceding records so is
        // discarded.
        assert_eq!(groups[0].label, "A");
        // chr2 is the final group.
        assert_eq!(groups[1].label, "default");
    }

    #[test]
    fn empty_hash_line_uses_default_label() {
        let data = b"chr1\t100\t200\n#\n";
        let groups: Vec<_> = make_reader(data)
            .collect::<Result<_, _>>()
            .expect("valid groups");
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].label, "default");
    }

    #[test]
    fn empty_file_yields_empty_file_error() {
        let data = b"";
        let result: Result<Vec<_>, _> = make_reader(data).collect();
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(matches!(err, BedReadError::EmptyFile));
    }

    #[test]
    fn comment_only_file_yields_empty_file_error() {
        let data = b"# just a comment\n\n# another comment\n";
        let result: Result<Vec<_>, _> = make_reader(data).collect();
        assert!(matches!(result.unwrap_err(), BedReadError::EmptyFile));
    }

    #[test]
    fn parse_error_includes_line_number_and_content() {
        let data = b"chr1\n"; // only 1 column – will fail BedRecord::parse
        let result: Result<Vec<_>, _> = make_reader(data).collect();
        let err = result.unwrap_err();
        assert!(matches!(err, BedReadError::Parse { .. }));
        if let BedReadError::Parse {
            line_number,
            message: _,
            line,
        } = &err
        {
            assert_eq!(*line_number, 1);
            assert_eq!(line, "chr1");
        } else {
            unreachable!();
        }
    }

    #[test]
    fn whitespace_only_lines_are_skipped() {
        let data = b"  \t  \nchr1\t100\t200\n";
        let groups: Vec<_> = make_reader(data)
            .collect::<Result<_, _>>()
            .expect("valid groups");
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].records.len(), 1);
    }

    #[test]
    fn raw_labels_not_deduplicated_within_file() {
        // Two groups with the same `#` label — both keep the raw label.
        let data = b"chr1\t100\t200\n#dup\nchr2\t300\t400\n#dup\n";
        let groups: Vec<_> = make_reader(data)
            .collect::<Result<_, _>>()
            .expect("valid groups");
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].label, "dup");
        assert_eq!(groups[1].label, "dup");
    }

    #[test]
    fn hash_line_with_whitespace_label() {
        // `#  ` → raw_label is empty after trim, so fall back to default.
        let data = b"chr1\t100\t200\n#  \n";
        let groups: Vec<_> = make_reader(data)
            .collect::<Result<_, _>>()
            .expect("valid groups");
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].label, "default");
    }

    #[test]
    fn keeps_invalid_score_and_strand_as_raw_strings() {
        let record = BedRecord::parse("chr1\t10\t20\tname\tabc\tstrandx").expect("should parse");
        assert!(record.score.is_none());
        assert_eq!(record.score_raw.as_deref(), Some("abc"));
        assert_eq!(record.strand, Strand::Unstranded);
        assert_eq!(record.strand_raw.as_deref(), Some("strandx"));
    }
}
