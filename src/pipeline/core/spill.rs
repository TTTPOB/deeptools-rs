use std::collections::HashMap;
use std::io::Write;
use std::path::PathBuf;
use std::sync::Arc;
use std::thread::JoinHandle;

use anyhow::{Context, Result};

use crate::io::readers::bed::{BedRecord, Strand};
use crate::pipeline::matrix::MatrixRow;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Default capacity for the reusable serialization buffer.
pub(crate) const SPILL_BUF_CAPACITY: usize = 1_048_576; // 1 MB

/// Memory threshold (bytes) above which a bucket should spill rows to disk.
pub(crate) const DEFAULT_MEMORY_SPILL_THRESHOLD: usize = 1024 * 1024 * 1024; // 1 GB

// ---------------------------------------------------------------------------
// SpillIndex
// ---------------------------------------------------------------------------

/// Index entry kept in memory for each row that has been spilled to disk.
/// Contains everything the merge/sort step needs without reading the row data.
#[derive(Debug, Clone)]
pub(crate) struct SpillIndex {
    /// Original input order index (for keep-order emit).
    pub(crate) orig_idx: usize,
    /// Which group bucket this row belongs to.
    pub(crate) group_index: usize,
    /// Pre-computed sort key (e.g. mean of values for sort-by-mean).
    pub(crate) sort_key: f64,
    /// Per-bucket push order, used as a stable tie-break when sort keys are equal.
    pub(crate) insertion_seq: u32,
    /// Byte offset within the spill file where this row's data begins.
    pub(crate) file_offset: u64,
    /// Length of the serialized row in bytes.
    pub(crate) row_byte_len: u32,
}

// ---------------------------------------------------------------------------
// ChromTable — bidirectional chrom <-> u16 mapping
// ---------------------------------------------------------------------------

/// Compact bidirectional mapping between chromosome names (`Arc<str>`) and
/// `u16` identifiers used in the on-disk spill format.
#[derive(Debug, Clone, Default)]
pub(crate) struct ChromTable {
    to_id: HashMap<Arc<str>, u16>,
    to_name: Vec<Arc<str>>,
}

impl ChromTable {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// Intern a chromosome name, returning a compact `u16` id. Re-uses an
    /// existing id when the same name is interned again.
    pub(crate) fn intern(&mut self, chrom: &Arc<str>) -> u16 {
        if let Some(&id) = self.to_id.get(chrom) {
            return id;
        }
        let id = self.to_name.len() as u16;
        self.to_name.push(Arc::clone(chrom));
        self.to_id.insert(Arc::clone(chrom), id);
        id
    }

    /// Resolve a u16 id back to the interned chromosome name.
    pub(crate) fn resolve(&self, id: u16) -> &Arc<str> {
        &self.to_name[id as usize]
    }
}

// ---------------------------------------------------------------------------
// Flag bits for per-row serialization
// ---------------------------------------------------------------------------

const FLAG_HAS_NAME: u8 = 0b0000_0001;
const FLAG_HAS_SCORE_RAW: u8 = 0b0000_0010;
const FLAG_HAS_STRAND_RAW: u8 = 0b0000_0100;
const FLAG_HAS_EXON_COORDS: u8 = 0b0000_1000;

// ---------------------------------------------------------------------------
// Serialization helpers
// ---------------------------------------------------------------------------

/// Write a length-prefixed UTF-8 string (u16 LE length + bytes).
fn write_len_prefixed_str(buf: &mut Vec<u8>, s: &str) -> Result<()> {
    let len = u16::try_from(s.len()).context("field too large for spill format")?;
    buf.extend_from_slice(&len.to_le_bytes());
    buf.extend_from_slice(s.as_bytes());
    Ok(())
}

/// Read a length-prefixed UTF-8 string, advancing `pos`.
fn read_len_prefixed_str(data: &[u8], pos: &mut usize) -> Result<String> {
    let len = read_u16(data, pos)? as usize;
    if *pos + len > data.len() {
        anyhow::bail!("unexpected end of data reading string of length {len}");
    }
    let s = std::str::from_utf8(&data[*pos..*pos + len])
        .context("invalid UTF-8 in spill data")?
        .to_owned();
    *pos += len;
    Ok(s)
}

fn read_u16(data: &[u8], pos: &mut usize) -> Result<u16> {
    if *pos + 2 > data.len() {
        anyhow::bail!("unexpected end of data reading u16 at offset {}", *pos);
    }
    let val = u16::from_le_bytes([data[*pos], data[*pos + 1]]);
    *pos += 2;
    Ok(val)
}

fn read_u32(data: &[u8], pos: &mut usize) -> Result<u32> {
    if *pos + 4 > data.len() {
        anyhow::bail!("unexpected end of data reading u32 at offset {}", *pos);
    }
    let val = u32::from_le_bytes([
        data[*pos],
        data[*pos + 1],
        data[*pos + 2],
        data[*pos + 3],
    ]);
    *pos += 4;
    Ok(val)
}

fn read_f64(data: &[u8], pos: &mut usize) -> Result<f64> {
    if *pos + 8 > data.len() {
        anyhow::bail!("unexpected end of data reading f64 at offset {}", *pos);
    }
    let val = f64::from_le_bytes([
        data[*pos],
        data[*pos + 1],
        data[*pos + 2],
        data[*pos + 3],
        data[*pos + 4],
        data[*pos + 5],
        data[*pos + 6],
        data[*pos + 7],
    ]);
    *pos += 8;
    Ok(val)
}

// Sentinel value indicating score is None (when score_raw is also None).
const SCORE_NONE_SENTINEL: u32 = 0xFFFF_FFFF;

// ---------------------------------------------------------------------------
// serialize_row
// ---------------------------------------------------------------------------

/// Serialize a `MatrixRow` into `buf`, returning the byte length written.
///
/// The caller is responsible for clearing `buf` before calling this function.
/// The format written is described in the module-level documentation.
pub(crate) fn serialize_row(
    buf: &mut Vec<u8>,
    row: &MatrixRow,
    chrom_table: &mut ChromTable,
) -> Result<u32> {
    buf.clear();

    let rec = &row.record;

    // chrom_id (2 bytes)
    let chrom_id = chrom_table.intern(&rec.chrom);
    buf.extend_from_slice(&chrom_id.to_le_bytes());

    // start (4 bytes)
    buf.extend_from_slice(&rec.start.to_le_bytes());

    // end (4 bytes)
    buf.extend_from_slice(&rec.end.to_le_bytes());

    // strand (1 byte)
    let strand_byte: u8 = match rec.strand {
        Strand::Positive => 0,
        Strand::Negative => 1,
        Strand::Unstranded => 2,
    };
    buf.push(strand_byte);

    // flags (1 byte) — reserve position, fill in after
    let flags_pos = buf.len();
    buf.push(0u8); // placeholder

    let mut flags: u8 = 0;

    // name (optional)
    if let Some(ref name) = rec.name {
        flags |= FLAG_HAS_NAME;
        write_len_prefixed_str(buf, name)?;
    }

    // score — three cases:
    //   1) score_raw is Some => flag bit 1 set, write len-prefixed raw string
    //   2) score_raw is None, score is Some => write f32 LE
    //   3) both None => write sentinel 0xFFFFFFFF
    if let Some(ref raw) = rec.score_raw {
        flags |= FLAG_HAS_SCORE_RAW;
        write_len_prefixed_str(buf, raw)?;
    } else if let Some(score) = rec.score {
        buf.extend_from_slice(&score.to_le_bytes());
    } else {
        buf.extend_from_slice(&SCORE_NONE_SENTINEL.to_le_bytes());
    }

    // strand_raw (optional)
    if let Some(ref strand_raw) = rec.strand_raw {
        flags |= FLAG_HAS_STRAND_RAW;
        write_len_prefixed_str(buf, strand_raw)?;
    }

    // extra_fields
    let extra_count =
        u16::try_from(rec.extra_fields.len()).context("too many extra fields for spill format")?;
    buf.extend_from_slice(&extra_count.to_le_bytes());
    for field in &rec.extra_fields {
        write_len_prefixed_str(buf, field)?;
    }

    // exon_coords (optional)
    if let Some(ref coords) = row.exon_coords {
        flags |= FLAG_HAS_EXON_COORDS;
        let count =
            u16::try_from(coords.len()).context("too many exon coords for spill format")?;
        buf.extend_from_slice(&count.to_le_bytes());
        for &(start, end) in coords {
            buf.extend_from_slice(&start.to_le_bytes());
            buf.extend_from_slice(&end.to_le_bytes());
        }
    }

    // sample_count, bin_count (4+4 bytes)
    buf.extend_from_slice(&(row.sample_count as u32).to_le_bytes());
    buf.extend_from_slice(&(row.bin_count as u32).to_le_bytes());

    // values (N × 8 bytes)
    for &val in &row.values {
        buf.extend_from_slice(&val.to_le_bytes());
    }

    // patch flags byte
    buf[flags_pos] = flags;

    let row_byte_len =
        u32::try_from(buf.len()).context("serialized row exceeds u32::MAX bytes")?;
    Ok(row_byte_len)
}

// ---------------------------------------------------------------------------
// deserialize_row
// ---------------------------------------------------------------------------

/// Deserialize a `MatrixRow` from a byte slice produced by `serialize_row`.
pub(crate) fn deserialize_row(data: &[u8], chrom_table: &ChromTable) -> Result<MatrixRow> {
    let mut pos: usize = 0;

    // chrom_id
    let chrom_id = read_u16(data, &mut pos)?;
    let chrom = Arc::clone(chrom_table.resolve(chrom_id));

    // start, end
    let start = read_u32(data, &mut pos)?;
    let end = read_u32(data, &mut pos)?;

    // strand
    if pos >= data.len() {
        anyhow::bail!("unexpected end of data reading strand byte");
    }
    let strand = match data[pos] {
        0 => Strand::Positive,
        1 => Strand::Negative,
        2 => Strand::Unstranded,
        other => anyhow::bail!("invalid strand byte: {other}"),
    };
    pos += 1;

    // flags
    if pos >= data.len() {
        anyhow::bail!("unexpected end of data reading flags byte");
    }
    let flags = data[pos];
    pos += 1;

    // name
    let name = if flags & FLAG_HAS_NAME != 0 {
        Some(read_len_prefixed_str(data, &mut pos)?)
    } else {
        None
    };

    // score
    let (score, score_raw) = if flags & FLAG_HAS_SCORE_RAW != 0 {
        let raw = read_len_prefixed_str(data, &mut pos)?;
        (None, Some(raw))
    } else {
        let raw_bits = read_u32(data, &mut pos)?;
        if raw_bits == SCORE_NONE_SENTINEL {
            (None, None)
        } else {
            (Some(f32::from_le_bytes(raw_bits.to_le_bytes())), None)
        }
    };

    // strand_raw
    let strand_raw = if flags & FLAG_HAS_STRAND_RAW != 0 {
        Some(read_len_prefixed_str(data, &mut pos)?)
    } else {
        None
    };

    // extra_fields
    let extra_count = read_u16(data, &mut pos)? as usize;
    let mut extra_fields = Vec::with_capacity(extra_count);
    for _ in 0..extra_count {
        extra_fields.push(read_len_prefixed_str(data, &mut pos)?);
    }

    // exon_coords
    let exon_coords = if flags & FLAG_HAS_EXON_COORDS != 0 {
        let count = read_u16(data, &mut pos)? as usize;
        let mut coords = Vec::with_capacity(count);
        for _ in 0..count {
            let s = read_u32(data, &mut pos)?;
            let e = read_u32(data, &mut pos)?;
            coords.push((s, e));
        }
        Some(coords)
    } else {
        None
    };

    // sample_count, bin_count
    let sample_count = read_u32(data, &mut pos)? as usize;
    let bin_count = read_u32(data, &mut pos)? as usize;

    // values
    let value_count = sample_count * bin_count;
    let mut values = Vec::with_capacity(value_count);
    for _ in 0..value_count {
        values.push(read_f64(data, &mut pos)?);
    }

    let record = BedRecord {
        chrom,
        start,
        end,
        name,
        score,
        score_raw,
        strand,
        strand_raw,
        extra_fields,
    };

    Ok(MatrixRow {
        record,
        values,
        sample_count,
        bin_count,
        exon_coords,
    })
}

// ---------------------------------------------------------------------------
// FlushResult
// ---------------------------------------------------------------------------

/// Result returned by a background flush thread.
struct FlushResult {
    /// Spill indices for the rows written.
    indices: Vec<SpillIndex>,
    /// The original Vec with data cleared but capacity retained, for reuse.
    returned_buf: Vec<(usize, f64, u32, MatrixRow)>,
    /// Chromosome name <-> id mapping built during serialization.
    chrom_table: ChromTable,
    /// Path to the temporary file holding the serialized rows.
    temp_path: PathBuf,
}

// ---------------------------------------------------------------------------
// flush_chunk — runs on a spawned thread
// ---------------------------------------------------------------------------

/// Serialize all rows to a new temporary file.
///
/// Each row's `(orig_idx, sort_key, insertion_seq)` is recorded in the
/// returned `SpillIndex` vector. The input Vec is cleared but its capacity
/// is retained so the caller can reuse it as a spare buffer.
fn flush_chunk(
    mut rows: Vec<(usize, f64, u32, MatrixRow)>,
    group_index: usize,
) -> Result<FlushResult> {
    let mut chrom_table = ChromTable::new();
    let mut indices = Vec::with_capacity(rows.len());
    let mut ser_buf = Vec::with_capacity(SPILL_BUF_CAPACITY);

    // Create a temporary file that persists (we keep the path for later reads).
    let tmp = tempfile::NamedTempFile::new().context("failed to create spill temp file")?;
    let temp_path = tmp.path().to_path_buf();
    let mut writer = std::io::BufWriter::new(tmp);

    let mut file_offset: u64 = 0;

    for &(orig_idx, sort_key, insertion_seq, ref row) in &rows {
        let row_byte_len = serialize_row(&mut ser_buf, row, &mut chrom_table)?;

        writer
            .write_all(&ser_buf)
            .context("failed to write spill data")?;

        indices.push(SpillIndex {
            orig_idx,
            group_index,
            sort_key,
            insertion_seq,
            file_offset,
            row_byte_len,
        });

        file_offset += row_byte_len as u64;
    }

    writer.flush().context("failed to flush spill writer")?;
    // Keep the underlying NamedTempFile alive by persisting it —
    // into_temp_path() would delete on drop, but persist() keeps the file.
    // Actually, NamedTempFile deletes on drop, so we need to persist it.
    let inner = writer.into_inner().context("failed to unwrap BufWriter")?;
    // persist without a target — keeps the file at the original temp path.
    inner
        .persist(&temp_path)
        .context("failed to persist temp file")?;

    // Clear the Vec but keep its capacity for reuse.
    rows.clear();

    Ok(FlushResult {
        indices,
        returned_buf: rows,
        chrom_table,
        temp_path,
    })
}

// ---------------------------------------------------------------------------
// CollectorBucket — per-group bucket with double-buffer spilling
// ---------------------------------------------------------------------------

/// A per-group bucket that accumulates rows in memory and flushes to disk
/// when the estimated byte usage exceeds the threshold.
struct CollectorBucket {
    /// Rows currently being accumulated.
    active: Vec<(usize, f64, u32, MatrixRow)>,
    /// Spare Vec with pre-allocated capacity, returned from a completed flush.
    spare: Option<Vec<(usize, f64, u32, MatrixRow)>>,
    /// Approximate memory usage of rows in `active`.
    estimated_bytes: usize,
    /// The group index this bucket belongs to.
    group_index: usize,
    /// In-flight flush thread handles.
    flush_handles: Vec<JoinHandle<Result<FlushResult>>>,
    /// Completed flush results (with returned_buf already drained).
    completed_flushes: Vec<FlushResult>,
    /// Per-bucket monotonically incrementing insertion counter.
    next_insertion_seq: u32,
}

impl CollectorBucket {
    fn new(group_index: usize) -> Self {
        Self {
            active: Vec::new(),
            spare: None,
            estimated_bytes: 0,
            group_index,
            flush_handles: Vec::new(),
            completed_flushes: Vec::new(),
            next_insertion_seq: 0,
        }
    }

    /// Push a row into the active buffer and trigger a flush if the estimated
    /// byte usage exceeds the threshold.
    fn push(
        &mut self,
        row: MatrixRow,
        orig_idx: usize,
        sort_key: f64,
        row_estimated_bytes: usize,
        threshold: usize,
    ) -> Result<()> {
        let seq = self.next_insertion_seq;
        self.next_insertion_seq = seq.wrapping_add(1);

        self.active.push((orig_idx, sort_key, seq, row));
        self.estimated_bytes += row_estimated_bytes;

        if self.estimated_bytes > threshold {
            self.trigger_flush()?;
        }

        Ok(())
    }

    /// Flush the active buffer to disk on a spawned thread.
    fn trigger_flush(&mut self) -> Result<()> {
        // Determine the replacement buffer for `active`.
        let replacement = if let Some(spare) = self.spare.take() {
            // Best case: reuse spare with pre-allocated capacity.
            spare
        } else if !self.flush_handles.is_empty() {
            // Back-pressure: join the oldest handle to get its returned_buf.
            let oldest = self.flush_handles.remove(0);
            let mut result = oldest
                .join()
                .map_err(|_| anyhow::anyhow!("flush thread panicked"))?
                .context("flush thread returned error")?;
            let spare = std::mem::take(&mut result.returned_buf);
            self.completed_flushes.push(result);
            spare
        } else {
            // First flush ever — allocate fresh.
            Vec::new()
        };

        let full_buf = std::mem::replace(&mut self.active, replacement);
        self.estimated_bytes = 0;

        let group_index = self.group_index;
        let handle = std::thread::spawn(move || flush_chunk(full_buf, group_index));
        self.flush_handles.push(handle);

        Ok(())
    }

    /// Join all remaining in-flight flush handles and collect their results.
    fn join_all(&mut self) -> Result<()> {
        for handle in self.flush_handles.drain(..) {
            let mut result = handle
                .join()
                .map_err(|_| anyhow::anyhow!("flush thread panicked"))?
                .context("flush thread returned error")?;
            // Release capacity — data is on disk, we don't need the returned_buf.
            result.returned_buf = Vec::new();
            self.completed_flushes.push(result);
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// HybridBucketCollector
// ---------------------------------------------------------------------------

/// Collector that keeps rows in per-group buckets and spills to disk when
/// memory usage exceeds a configurable threshold. Uses a double-buffer scheme
/// to reuse Vec capacity across flushes, with back-pressure when the spare
/// buffer is unavailable.
pub(crate) struct HybridBucketCollector {
    buckets: Vec<CollectorBucket>,
    sample_count: usize,
    bin_count: usize,
    threshold: usize,
    row_estimated_bytes: usize,
}

impl HybridBucketCollector {
    /// Create a new collector with the default memory threshold.
    pub(crate) fn new(group_count: usize, sample_count: usize, bin_count: usize) -> Self {
        Self::with_threshold(
            group_count,
            sample_count,
            bin_count,
            DEFAULT_MEMORY_SPILL_THRESHOLD,
        )
    }

    /// Create a new collector with a custom memory threshold per bucket.
    pub(crate) fn with_threshold(
        group_count: usize,
        sample_count: usize,
        bin_count: usize,
        threshold: usize,
    ) -> Self {
        let buckets = (0..group_count)
            .map(|i| CollectorBucket::new(i))
            .collect();
        // Approximate bytes per row: f64 values + overhead for BedRecord, etc.
        let row_estimated_bytes = sample_count * bin_count * 8 + 100;
        Self {
            buckets,
            sample_count,
            bin_count,
            threshold,
            row_estimated_bytes,
        }
    }

    /// Push a row into the appropriate group bucket.
    pub(crate) fn push(
        &mut self,
        row: MatrixRow,
        orig_idx: usize,
        group_index: usize,
        sort_key: f64,
    ) -> Result<()> {
        let bucket = self
            .buckets
            .get_mut(group_index)
            .context("group_index out of range")?;
        bucket.push(row, orig_idx, sort_key, self.row_estimated_bytes, self.threshold)
    }

    /// Join all in-flight flushes across all buckets.
    pub(crate) fn join_all(&mut self) -> Result<()> {
        for bucket in &mut self.buckets {
            bucket.join_all()?;
        }
        Ok(())
    }

    /// Access the buckets (for reading results after join_all).
    #[cfg(test)]
    fn buckets(&self) -> &[CollectorBucket] {
        &self.buckets
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a fully-populated MatrixRow for round-trip testing.
    fn make_full_row() -> MatrixRow {
        let record = BedRecord {
            chrom: Arc::from("chr1"),
            start: 1000,
            end: 2000,
            name: Some("gene_A".to_string()),
            score: Some(42.5),
            score_raw: None,
            strand: Strand::Negative,
            strand_raw: None,
            extra_fields: vec!["field1".to_string(), "field2".to_string()],
        };
        MatrixRow {
            record,
            values: vec![1.0, f64::NAN, 3.5, -0.0, f64::INFINITY],
            sample_count: 1,
            bin_count: 5,
            exon_coords: Some(vec![(1000, 1200), (1500, 2000)]),
        }
    }

    /// Build a minimal MatrixRow (no optional fields).
    fn make_minimal_row() -> MatrixRow {
        let record = BedRecord {
            chrom: Arc::from("chrX"),
            start: 0,
            end: 100,
            name: None,
            score: None,
            score_raw: None,
            strand: Strand::Unstranded,
            strand_raw: None,
            extra_fields: vec![],
        };
        MatrixRow {
            record,
            values: vec![],
            sample_count: 0,
            bin_count: 0,
            exon_coords: None,
        }
    }

    #[test]
    fn round_trip_full_row() {
        let row = make_full_row();
        let mut ct = ChromTable::new();
        let mut buf = Vec::with_capacity(SPILL_BUF_CAPACITY);

        let byte_len = serialize_row(&mut buf, &row, &mut ct).unwrap();
        assert_eq!(byte_len as usize, buf.len());

        let restored = deserialize_row(&buf, &ct).unwrap();

        assert_eq!(&*restored.record.chrom, "chr1");
        assert_eq!(restored.record.start, 1000);
        assert_eq!(restored.record.end, 2000);
        assert_eq!(restored.record.name.as_deref(), Some("gene_A"));
        assert_eq!(restored.record.score, Some(42.5));
        assert!(restored.record.score_raw.is_none());
        assert_eq!(restored.record.strand, Strand::Negative);
        assert!(restored.record.strand_raw.is_none());
        assert_eq!(restored.record.extra_fields, vec!["field1", "field2"]);
        assert_eq!(restored.sample_count, 1);
        assert_eq!(restored.bin_count, 5);

        // Check values — NaN needs special handling
        assert_eq!(restored.values.len(), 5);
        assert_eq!(restored.values[0], 1.0);
        assert!(restored.values[1].is_nan());
        assert_eq!(restored.values[2], 3.5);
        assert_eq!(restored.values[3], 0.0); // -0.0 bit pattern preserved via f64 LE
        assert_eq!(restored.values[4], f64::INFINITY);

        // exon_coords
        let coords = restored.exon_coords.unwrap();
        assert_eq!(coords, vec![(1000, 1200), (1500, 2000)]);
    }

    #[test]
    fn round_trip_minimal_row() {
        let row = make_minimal_row();
        let mut ct = ChromTable::new();
        let mut buf = Vec::with_capacity(SPILL_BUF_CAPACITY);

        let byte_len = serialize_row(&mut buf, &row, &mut ct).unwrap();
        assert_eq!(byte_len as usize, buf.len());

        let restored = deserialize_row(&buf, &ct).unwrap();

        assert_eq!(&*restored.record.chrom, "chrX");
        assert_eq!(restored.record.start, 0);
        assert_eq!(restored.record.end, 100);
        assert!(restored.record.name.is_none());
        assert!(restored.record.score.is_none());
        assert!(restored.record.score_raw.is_none());
        assert_eq!(restored.record.strand, Strand::Unstranded);
        assert!(restored.record.strand_raw.is_none());
        assert!(restored.record.extra_fields.is_empty());
        assert_eq!(restored.sample_count, 0);
        assert_eq!(restored.bin_count, 0);
        assert!(restored.values.is_empty());
        assert!(restored.exon_coords.is_none());
    }

    #[test]
    fn chrom_table_interning() {
        let mut ct = ChromTable::new();

        let chr1: Arc<str> = Arc::from("chr1");
        let chr2: Arc<str> = Arc::from("chr2");

        let id1 = ct.intern(&chr1);
        let id2 = ct.intern(&chr2);
        let id1_again = ct.intern(&chr1);

        // Same name returns same id
        assert_eq!(id1, id1_again);
        // Different names get different ids
        assert_ne!(id1, id2);

        // Resolve round-trips
        assert_eq!(&**ct.resolve(id1), "chr1");
        assert_eq!(&**ct.resolve(id2), "chr2");
    }

    #[test]
    fn round_trip_score_raw_and_strand_raw() {
        let record = BedRecord {
            chrom: Arc::from("chr3"),
            start: 500,
            end: 600,
            name: Some("item".to_string()),
            score: None,
            score_raw: Some("abc".to_string()),
            strand: Strand::Unstranded,
            strand_raw: Some("strandx".to_string()),
            extra_fields: vec![],
        };
        let row = MatrixRow {
            record,
            values: vec![9.9],
            sample_count: 1,
            bin_count: 1,
            exon_coords: None,
        };

        let mut ct = ChromTable::new();
        let mut buf = Vec::new();
        let _len = serialize_row(&mut buf, &row, &mut ct).unwrap();

        let restored = deserialize_row(&buf, &ct).unwrap();

        assert!(restored.record.score.is_none());
        assert_eq!(restored.record.score_raw.as_deref(), Some("abc"));
        assert_eq!(restored.record.strand_raw.as_deref(), Some("strandx"));
        assert_eq!(restored.values, vec![9.9]);
    }

    // -----------------------------------------------------------------------
    // HybridBucketCollector tests
    // -----------------------------------------------------------------------

    /// Helper: create a simple MatrixRow with given chrom/start and N values.
    fn make_test_row(chrom: &str, start: u32, n_values: usize) -> MatrixRow {
        let record = BedRecord {
            chrom: Arc::from(chrom),
            start,
            end: start + 100,
            name: None,
            score: None,
            score_raw: None,
            strand: Strand::Unstranded,
            strand_raw: None,
            extra_fields: vec![],
        };
        MatrixRow {
            record,
            values: vec![1.0; n_values],
            sample_count: 1,
            bin_count: n_values,
            exon_coords: None,
        }
    }

    #[test]
    fn hybrid_collector_small_data_stays_in_memory() {
        // With a high threshold, small data should never spill.
        let mut collector = HybridBucketCollector::with_threshold(2, 1, 3, 1_000_000);

        // Push a few rows into group 0 and group 1.
        for i in 0..5 {
            let row = make_test_row("chr1", i * 100, 3);
            collector.push(row, i as usize, 0, i as f64).unwrap();
        }
        for i in 0..3 {
            let row = make_test_row("chr2", i * 200, 3);
            collector.push(row, (5 + i) as usize, 1, i as f64).unwrap();
        }

        collector.join_all().unwrap();

        // No flushes should have occurred.
        assert!(collector.buckets()[0].completed_flushes.is_empty());
        assert!(collector.buckets()[1].completed_flushes.is_empty());

        // All rows should still be in the active buffer.
        assert_eq!(collector.buckets()[0].active.len(), 5);
        assert_eq!(collector.buckets()[1].active.len(), 3);

        // Verify insertion_seq is monotonically incrementing per bucket.
        for (i, entry) in collector.buckets()[0].active.iter().enumerate() {
            assert_eq!(entry.2, i as u32); // insertion_seq
        }
        for (i, entry) in collector.buckets()[1].active.iter().enumerate() {
            assert_eq!(entry.2, i as u32);
        }
    }

    #[test]
    fn hybrid_collector_low_threshold_triggers_spill() {
        // With threshold of 100 bytes and rows that are ~124 bytes each
        // (1 * 3 * 8 + 100 = 124), even a single row should trigger a flush.
        // Actually, the first push makes estimated_bytes=124 > 100, so flush
        // is triggered after the first push.
        let mut collector = HybridBucketCollector::with_threshold(1, 1, 3, 100);

        let n_rows = 10;
        for i in 0..n_rows {
            let row = make_test_row("chr1", i * 100, 3);
            collector.push(row, i as usize, 0, i as f64).unwrap();
        }

        collector.join_all().unwrap();

        let bucket = &collector.buckets()[0];
        // Multiple flushes should have occurred.
        assert!(
            !bucket.completed_flushes.is_empty(),
            "expected at least one flush with low threshold"
        );

        // Verify that all spill indices have the correct group_index.
        for flush in &bucket.completed_flushes {
            for idx in &flush.indices {
                assert_eq!(idx.group_index, 0);
            }
        }

        // Count total rows: flushed + still in active buffer.
        let flushed_count: usize = bucket
            .completed_flushes
            .iter()
            .map(|f| f.indices.len())
            .sum();
        let total = flushed_count + bucket.active.len();
        assert_eq!(total, n_rows as usize, "all rows must be accounted for");

        // Verify temp files exist on disk.
        for flush in &bucket.completed_flushes {
            assert!(
                flush.temp_path.exists(),
                "temp file should exist: {:?}",
                flush.temp_path
            );
        }

        // Verify that returned_buf has been cleared (capacity released) in
        // completed flushes.
        for flush in &bucket.completed_flushes {
            assert!(
                flush.returned_buf.is_empty(),
                "returned_buf should be empty after join_all"
            );
            assert_eq!(
                flush.returned_buf.capacity(),
                0,
                "returned_buf capacity should be 0 after join_all"
            );
        }

        // Verify insertion_seq values are unique and correct across all rows.
        let mut all_seqs: Vec<u32> = bucket
            .completed_flushes
            .iter()
            .flat_map(|f| f.indices.iter().map(|idx| idx.insertion_seq))
            .chain(bucket.active.iter().map(|entry| entry.2))
            .collect();
        all_seqs.sort();
        let expected: Vec<u32> = (0..n_rows as u32).collect();
        assert_eq!(all_seqs, expected, "insertion_seq must be 0..N contiguous");
    }

    #[test]
    fn hybrid_collector_spill_data_is_deserializable() {
        // Verify that data written to temp files can actually be read back.
        let mut collector = HybridBucketCollector::with_threshold(1, 1, 2, 100);

        // Push enough rows to trigger at least one flush.
        for i in 0..5 {
            let row = make_test_row("chr1", i * 100, 2);
            collector.push(row, i as usize, 0, (i as f64) * 1.5).unwrap();
        }

        collector.join_all().unwrap();

        let bucket = &collector.buckets()[0];
        assert!(!bucket.completed_flushes.is_empty());

        // Read back from the first flush and verify the data.
        let flush = &bucket.completed_flushes[0];
        let file_data = std::fs::read(&flush.temp_path).unwrap();

        for idx in &flush.indices {
            let offset = idx.file_offset as usize;
            let len = idx.row_byte_len as usize;
            let slice = &file_data[offset..offset + len];
            let restored = deserialize_row(slice, &flush.chrom_table).unwrap();

            assert_eq!(&*restored.record.chrom, "chr1");
            assert_eq!(restored.sample_count, 1);
            assert_eq!(restored.bin_count, 2);
            assert_eq!(restored.values.len(), 2);
        }
    }

    #[test]
    fn hybrid_collector_multiple_groups_independent() {
        // Each group has its own insertion_seq counter and flush state.
        let mut collector = HybridBucketCollector::with_threshold(3, 1, 2, 100);

        // Push rows to different groups.
        for i in 0..4 {
            let row = make_test_row("chr1", i * 100, 2);
            collector.push(row, i as usize, 0, 0.0).unwrap();
        }
        for i in 0..4 {
            let row = make_test_row("chr2", i * 100, 2);
            collector.push(row, (10 + i) as usize, 1, 0.0).unwrap();
        }
        // Group 2 gets no rows — should be fine.

        collector.join_all().unwrap();

        // Group 2 should have no flushes and empty active.
        assert!(collector.buckets()[2].completed_flushes.is_empty());
        assert!(collector.buckets()[2].active.is_empty());

        // Each of groups 0 and 1 should have some rows (flushed + active).
        for g in 0..2 {
            let bucket = &collector.buckets()[g];
            let flushed: usize = bucket
                .completed_flushes
                .iter()
                .map(|f| f.indices.len())
                .sum();
            let total = flushed + bucket.active.len();
            assert_eq!(total, 4, "group {g} should have 4 rows total");
        }
    }
}
