use std::collections::HashMap;
use std::io::Write;
use std::path::PathBuf;
use std::sync::Arc;
use std::thread::JoinHandle;

use anyhow::{Context, Result};
use memmap2::Mmap;

use crate::io::readers::bed::{BedRecord, Strand};
use crate::pipeline::matrix::{MatrixHeader, MatrixRow, compare_ascending};

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
    #[allow(dead_code)]
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
    #[allow(dead_code)]
    sample_count: usize,
    #[allow(dead_code)]
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

    /// Finalize the collector in sorted mode.
    ///
    /// Joins all in-flight flushes, then for each bucket (group) emits rows in
    /// sorted order by `(sort_key, insertion_seq)`. When `sort_ascending` is
    /// false the final order is reversed (matching the `sort_by(ascending) +
    /// reverse()` behaviour of the previous in-memory sort).
    ///
    /// The `header_builder` closure receives the final per-group row counts and
    /// returns a `MatrixHeader`. The `emit` closure receives `(group_index, row)`.
    pub(crate) fn finalize_sorted<F>(
        mut self,
        sort_ascending: bool,
        header_builder: impl FnOnce(Vec<usize>) -> Result<MatrixHeader>,
        mut emit: F,
    ) -> Result<MatrixHeader>
    where
        F: FnMut(usize, MatrixRow) -> Result<()>,
    {
        // 1. Join all in-flight flushes.
        self.join_all()?;

        // 2. Compute group counts for the header.
        let group_counts: Vec<usize> = self
            .buckets
            .iter()
            .map(|b| {
                let flushed: usize = b.completed_flushes.iter().map(|f| f.indices.len()).sum();
                flushed + b.active.len()
            })
            .collect();
        let header = header_builder(group_counts)?;

        // 3. For each bucket, merge spilled + in-memory rows and emit in sorted order.
        for bucket in self.buckets.drain(..) {
            let group_index = bucket.group_index;

            // mmap all spill files for this bucket.
            let mut mmaps: Vec<Mmap> = Vec::with_capacity(bucket.completed_flushes.len());
            let mut chrom_tables: Vec<&ChromTable> =
                Vec::with_capacity(bucket.completed_flushes.len());
            let mut temp_paths: Vec<PathBuf> = Vec::new();

            // We need to hold references to completed_flushes, so collect everything first.
            // Build a unified list of emit entries.

            // An entry that can be either a spill reference or an owned in-memory row.
            enum EmitEntry {
                Spill {
                    sort_key: f64,
                    insertion_seq: u32,
                    mmap_idx: usize,
                    file_offset: u64,
                    row_byte_len: u32,
                },
                InMemory {
                    sort_key: f64,
                    insertion_seq: u32,
                    row: MatrixRow,
                },
            }

            impl EmitEntry {
                fn sort_key(&self) -> f64 {
                    match self {
                        EmitEntry::Spill { sort_key, .. } => *sort_key,
                        EmitEntry::InMemory { sort_key, .. } => *sort_key,
                    }
                }
                fn insertion_seq(&self) -> u32 {
                    match self {
                        EmitEntry::Spill { insertion_seq, .. } => *insertion_seq,
                        EmitEntry::InMemory { insertion_seq, .. } => *insertion_seq,
                    }
                }
            }

            // Open mmaps for all spill files.
            for flush in &bucket.completed_flushes {
                let file = std::fs::File::open(&flush.temp_path)
                    .with_context(|| format!("failed to open spill file {:?}", flush.temp_path))?;
                // SAFETY: the file is complete and no longer being written to.
                let mmap = unsafe { Mmap::map(&file) }
                    .with_context(|| format!("failed to mmap spill file {:?}", flush.temp_path))?;
                mmaps.push(mmap);
                chrom_tables.push(&flush.chrom_table);
                temp_paths.push(flush.temp_path.clone());
            }

            let mut all_entries: Vec<EmitEntry> = Vec::new();

            // Collect spill references.
            for (flush_idx, flush) in bucket.completed_flushes.iter().enumerate() {
                for idx in &flush.indices {
                    all_entries.push(EmitEntry::Spill {
                        sort_key: idx.sort_key,
                        insertion_seq: idx.insertion_seq,
                        mmap_idx: flush_idx,
                        file_offset: idx.file_offset,
                        row_byte_len: idx.row_byte_len,
                    });
                }
            }

            // Collect in-memory rows.
            for (_orig_idx, sort_key, insertion_seq, row) in bucket.active {
                all_entries.push(EmitEntry::InMemory {
                    sort_key,
                    insertion_seq,
                    row,
                });
            }

            // Sort ascending by (sort_key, insertion_seq).
            all_entries.sort_by(|a, b| {
                compare_ascending(a.sort_key(), b.sort_key())
                    .then(a.insertion_seq().cmp(&b.insertion_seq()))
            });
            if !sort_ascending {
                all_entries.reverse();
            }

            // Emit rows.
            for entry in all_entries {
                let row = match entry {
                    EmitEntry::Spill {
                        mmap_idx,
                        file_offset,
                        row_byte_len,
                        ..
                    } => {
                        let offset = file_offset as usize;
                        let len = row_byte_len as usize;
                        let data = &mmaps[mmap_idx][offset..offset + len];
                        deserialize_row(data, chrom_tables[mmap_idx])?
                    }
                    EmitEntry::InMemory { row, .. } => row,
                };
                emit(group_index, row)?;
            }

            // Cleanup: drop mmaps then remove temp files.
            drop(mmaps);
            for path in &temp_paths {
                let _ = std::fs::remove_file(path);
            }
        }

        Ok(header)
    }

    /// Finalize the collector in keep-order mode.
    ///
    /// Joins all in-flight flushes, then places every row into a flat placement
    /// array indexed by `orig_idx`. A linear scan emits rows in their original
    /// input order. Slots at filtered-row indices remain `None` and are skipped.
    ///
    /// `task_count` is the total number of input tasks (regions) — it determines
    /// the size of the placement array.
    pub(crate) fn finalize_keep_order<F>(
        mut self,
        task_count: usize,
        header_builder: impl FnOnce(Vec<usize>) -> Result<MatrixHeader>,
        mut emit: F,
    ) -> Result<MatrixHeader>
    where
        F: FnMut(usize, MatrixRow) -> Result<()>,
    {
        // 1. Join all in-flight flushes.
        self.join_all()?;

        // Slot for placement array: group_index + either a spill reference or owned row.
        enum RowSlot {
            Spill {
                group_index: usize,
                mmap_idx: usize,
                file_offset: u64,
                row_byte_len: u32,
            },
            InMemory {
                group_index: usize,
                row: MatrixRow,
            },
        }

        impl RowSlot {
            fn group_index(&self) -> usize {
                match self {
                    RowSlot::Spill { group_index, .. } => *group_index,
                    RowSlot::InMemory { group_index, .. } => *group_index,
                }
            }
        }

        let mut slots: Vec<Option<RowSlot>> = (0..task_count).map(|_| None).collect();

        // Global mmap / chrom_table storage keyed by a global flush index.
        let mut mmaps: Vec<Mmap> = Vec::new();
        let mut chrom_tables: Vec<ChromTable> = Vec::new();
        let mut temp_paths: Vec<PathBuf> = Vec::new();

        // 2. Populate slots from all buckets.
        let group_count = self.buckets.len();
        for bucket in self.buckets.drain(..) {
            let group_index = bucket.group_index;

            // Process spilled flushes.
            for flush in bucket.completed_flushes {
                let mmap_idx = mmaps.len();
                let file = std::fs::File::open(&flush.temp_path).with_context(|| {
                    format!("failed to open spill file {:?}", flush.temp_path)
                })?;
                let mmap = unsafe { Mmap::map(&file) }.with_context(|| {
                    format!("failed to mmap spill file {:?}", flush.temp_path)
                })?;
                mmaps.push(mmap);
                chrom_tables.push(flush.chrom_table);
                temp_paths.push(flush.temp_path);

                for idx in &flush.indices {
                    if idx.orig_idx < task_count {
                        slots[idx.orig_idx] = Some(RowSlot::Spill {
                            group_index,
                            mmap_idx,
                            file_offset: idx.file_offset,
                            row_byte_len: idx.row_byte_len,
                        });
                    }
                }
            }

            // Process in-memory rows.
            for (orig_idx, _sort_key, _insertion_seq, row) in bucket.active {
                if orig_idx < task_count {
                    slots[orig_idx] = Some(RowSlot::InMemory {
                        group_index,
                        row,
                    });
                }
            }
        }

        // 3. Count rows per group.
        let mut group_counts = vec![0usize; group_count];
        for slot in &slots {
            if let Some(s) = slot {
                let gi = s.group_index();
                if gi < group_counts.len() {
                    group_counts[gi] += 1;
                }
            }
        }

        let header = header_builder(group_counts)?;

        // 4. Linear scan — emit in original order, skip None slots (filtered rows).
        for slot in slots {
            if let Some(s) = slot {
                let gi = s.group_index();
                let row = match s {
                    RowSlot::Spill {
                        mmap_idx,
                        file_offset,
                        row_byte_len,
                        ..
                    } => {
                        let offset = file_offset as usize;
                        let len = row_byte_len as usize;
                        let data = &mmaps[mmap_idx][offset..offset + len];
                        deserialize_row(data, &chrom_tables[mmap_idx])?
                    }
                    RowSlot::InMemory { row, .. } => row,
                };
                emit(gi, row)?;
            }
        }

        // 5. Cleanup: drop mmaps then remove temp files.
        drop(mmaps);
        for path in &temp_paths {
            let _ = std::fs::remove_file(path);
        }

        Ok(header)
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

    // -----------------------------------------------------------------------
    // finalize_sorted tests
    // -----------------------------------------------------------------------

    /// Helper: build a trivial MatrixHeader from group counts (test only).
    fn test_header_builder(group_counts: Vec<usize>) -> Result<MatrixHeader> {
        Ok(MatrixHeader::default_for_test(group_counts))
    }

    /// Helper: create a row whose first value encodes a tag so we can track it.
    fn make_tagged_row(chrom: &str, start: u32, tag: f64) -> MatrixRow {
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
            values: vec![tag],
            sample_count: 1,
            bin_count: 1,
            exon_coords: None,
        }
    }

    #[test]
    fn finalize_sorted_ascending_in_memory_only() {
        // All rows stay in memory (high threshold). Verify ascending sort order.
        let mut collector = HybridBucketCollector::with_threshold(1, 1, 1, 1_000_000);

        // Push rows with sort keys: 3.0, 1.0, 2.0
        collector
            .push(make_tagged_row("chr1", 300, 30.0), 0, 0, 3.0)
            .unwrap();
        collector
            .push(make_tagged_row("chr1", 100, 10.0), 1, 0, 1.0)
            .unwrap();
        collector
            .push(make_tagged_row("chr1", 200, 20.0), 2, 0, 2.0)
            .unwrap();

        let mut emitted: Vec<(usize, f64)> = Vec::new();
        let header = collector
            .finalize_sorted(true, test_header_builder, |gi, row| {
                emitted.push((gi, row.values[0]));
                Ok(())
            })
            .unwrap();

        // Ascending: sort_key 1.0, 2.0, 3.0 → tags 10.0, 20.0, 30.0
        assert_eq!(emitted.len(), 3);
        assert_eq!(emitted[0], (0, 10.0));
        assert_eq!(emitted[1], (0, 20.0));
        assert_eq!(emitted[2], (0, 30.0));

        // Header should reflect 3 rows in group 0.
        assert_eq!(header.group_boundaries, vec![0, 3]);
    }

    #[test]
    fn finalize_sorted_descending_in_memory_only() {
        let mut collector = HybridBucketCollector::with_threshold(1, 1, 1, 1_000_000);

        collector
            .push(make_tagged_row("chr1", 300, 30.0), 0, 0, 3.0)
            .unwrap();
        collector
            .push(make_tagged_row("chr1", 100, 10.0), 1, 0, 1.0)
            .unwrap();
        collector
            .push(make_tagged_row("chr1", 200, 20.0), 2, 0, 2.0)
            .unwrap();

        let mut emitted: Vec<f64> = Vec::new();
        let _header = collector
            .finalize_sorted(false, test_header_builder, |_gi, row| {
                emitted.push(row.values[0]);
                Ok(())
            })
            .unwrap();

        // Descending: sort_key 3.0, 2.0, 1.0 → tags 30.0, 20.0, 10.0
        assert_eq!(emitted, vec![30.0, 20.0, 10.0]);
    }

    #[test]
    fn finalize_sorted_with_spill_ascending() {
        // Use a very low threshold so every row triggers a spill.
        let mut collector = HybridBucketCollector::with_threshold(1, 1, 1, 1);

        // Push rows with sort keys: 5.0, 1.0, 3.0, 2.0, 4.0
        let keys = [5.0, 1.0, 3.0, 2.0, 4.0];
        for (i, &key) in keys.iter().enumerate() {
            let tag = key * 10.0;
            collector
                .push(make_tagged_row("chr1", i as u32 * 100, tag), i, 0, key)
                .unwrap();
        }

        let mut emitted: Vec<f64> = Vec::new();
        let header = collector
            .finalize_sorted(true, test_header_builder, |_gi, row| {
                emitted.push(row.values[0]);
                Ok(())
            })
            .unwrap();

        // Ascending: sort_key 1.0, 2.0, 3.0, 4.0, 5.0 → tags 10.0, 20.0, 30.0, 40.0, 50.0
        assert_eq!(emitted, vec![10.0, 20.0, 30.0, 40.0, 50.0]);
        assert_eq!(header.group_boundaries, vec![0, 5]);
    }

    #[test]
    fn finalize_sorted_with_spill_descending() {
        let mut collector = HybridBucketCollector::with_threshold(1, 1, 1, 1);

        let keys = [5.0, 1.0, 3.0, 2.0, 4.0];
        for (i, &key) in keys.iter().enumerate() {
            let tag = key * 10.0;
            collector
                .push(make_tagged_row("chr1", i as u32 * 100, tag), i, 0, key)
                .unwrap();
        }

        let mut emitted: Vec<f64> = Vec::new();
        let _header = collector
            .finalize_sorted(false, test_header_builder, |_gi, row| {
                emitted.push(row.values[0]);
                Ok(())
            })
            .unwrap();

        // Descending: 50.0, 40.0, 30.0, 20.0, 10.0
        assert_eq!(emitted, vec![50.0, 40.0, 30.0, 20.0, 10.0]);
    }

    #[test]
    fn finalize_sorted_stable_tiebreak() {
        // When sort keys are equal, insertion_seq determines order.
        // Ascending: equal keys keep insertion order.
        // Descending: equal keys appear in reverse insertion order.
        let mut collector = HybridBucketCollector::with_threshold(1, 1, 1, 1);

        // All sort keys are 1.0. Tags encode insertion order: 10, 20, 30.
        collector
            .push(make_tagged_row("chr1", 0, 10.0), 0, 0, 1.0)
            .unwrap();
        collector
            .push(make_tagged_row("chr1", 100, 20.0), 1, 0, 1.0)
            .unwrap();
        collector
            .push(make_tagged_row("chr1", 200, 30.0), 2, 0, 1.0)
            .unwrap();

        // Ascending: insertion order preserved.
        let mut emitted_asc: Vec<f64> = Vec::new();
        let collector_asc = HybridBucketCollector::with_threshold(1, 1, 1, 1);
        // Need a fresh collector for ascending test.
        let mut c = HybridBucketCollector::with_threshold(1, 1, 1, 1);
        c.push(make_tagged_row("chr1", 0, 10.0), 0, 0, 1.0)
            .unwrap();
        c.push(make_tagged_row("chr1", 100, 20.0), 1, 0, 1.0)
            .unwrap();
        c.push(make_tagged_row("chr1", 200, 30.0), 2, 0, 1.0)
            .unwrap();
        let _ = c
            .finalize_sorted(true, test_header_builder, |_gi, row| {
                emitted_asc.push(row.values[0]);
                Ok(())
            })
            .unwrap();
        assert_eq!(emitted_asc, vec![10.0, 20.0, 30.0]);

        // Descending: reversed insertion order (sort ascending then .reverse()).
        let mut emitted_desc: Vec<f64> = Vec::new();
        let _header = collector
            .finalize_sorted(false, test_header_builder, |_gi, row| {
                emitted_desc.push(row.values[0]);
                Ok(())
            })
            .unwrap();
        assert_eq!(emitted_desc, vec![30.0, 20.0, 10.0]);

        drop(collector_asc); // suppress unused warning
    }

    #[test]
    fn finalize_sorted_multiple_groups() {
        // Two groups, each with their own sort order.
        let mut collector = HybridBucketCollector::with_threshold(2, 1, 1, 1);

        // Group 0: keys 3.0, 1.0
        collector
            .push(make_tagged_row("chr1", 0, 30.0), 0, 0, 3.0)
            .unwrap();
        collector
            .push(make_tagged_row("chr1", 100, 10.0), 1, 0, 1.0)
            .unwrap();

        // Group 1: keys 2.0, 4.0
        collector
            .push(make_tagged_row("chr2", 0, 20.0), 2, 1, 2.0)
            .unwrap();
        collector
            .push(make_tagged_row("chr2", 100, 40.0), 3, 1, 4.0)
            .unwrap();

        let mut emitted: Vec<(usize, f64)> = Vec::new();
        let header = collector
            .finalize_sorted(true, test_header_builder, |gi, row| {
                emitted.push((gi, row.values[0]));
                Ok(())
            })
            .unwrap();

        // Group 0 ascending: 10.0, 30.0
        // Group 1 ascending: 20.0, 40.0
        assert_eq!(emitted.len(), 4);
        assert_eq!(emitted[0], (0, 10.0));
        assert_eq!(emitted[1], (0, 30.0));
        assert_eq!(emitted[2], (1, 20.0));
        assert_eq!(emitted[3], (1, 40.0));

        assert_eq!(header.group_boundaries, vec![0, 2, 4]);
    }

    #[test]
    fn finalize_sorted_temp_files_cleaned_up() {
        let mut collector = HybridBucketCollector::with_threshold(1, 1, 1, 1);

        for i in 0..5 {
            collector
                .push(make_tagged_row("chr1", i * 100, i as f64), i as usize, 0, i as f64)
                .unwrap();
        }

        // Join to get temp paths before finalize consumes the collector.
        collector.join_all().unwrap();
        let temp_paths: Vec<PathBuf> = collector
            .buckets()
            .iter()
            .flat_map(|b| b.completed_flushes.iter().map(|f| f.temp_path.clone()))
            .collect();
        assert!(!temp_paths.is_empty(), "should have at least one spill file");

        // Re-create collector to run finalize (since join_all already consumed handles).
        // Actually, we already joined — finalize_sorted calls join_all again which is a no-op.
        let _header = collector
            .finalize_sorted(true, test_header_builder, |_gi, _row| Ok(()))
            .unwrap();

        // All temp files should have been removed.
        for path in &temp_paths {
            assert!(
                !path.exists(),
                "temp file should have been removed: {:?}",
                path
            );
        }
    }

    // -----------------------------------------------------------------------
    // finalize_keep_order tests
    // -----------------------------------------------------------------------

    #[test]
    fn finalize_keep_order_preserves_original_indices() {
        let mut collector = HybridBucketCollector::with_threshold(1, 1, 1, 1_000_000);

        // Push rows with orig_idx: 2, 0, 4 (out of task_count=5).
        // Indices 1 and 3 are "filtered" (absent).
        collector
            .push(make_tagged_row("chr1", 200, 20.0), 2, 0, 0.0)
            .unwrap();
        collector
            .push(make_tagged_row("chr1", 0, 0.0), 0, 0, 0.0)
            .unwrap();
        collector
            .push(make_tagged_row("chr1", 400, 40.0), 4, 0, 0.0)
            .unwrap();

        let mut emitted: Vec<(usize, f64)> = Vec::new();
        let header = collector
            .finalize_keep_order(5, test_header_builder, |gi, row| {
                emitted.push((gi, row.values[0]));
                Ok(())
            })
            .unwrap();

        // Should emit in orig_idx order: 0, 2, 4 (skipping 1 and 3).
        assert_eq!(emitted.len(), 3);
        assert_eq!(emitted[0], (0, 0.0));   // orig_idx=0
        assert_eq!(emitted[1], (0, 20.0));  // orig_idx=2
        assert_eq!(emitted[2], (0, 40.0));  // orig_idx=4

        assert_eq!(header.group_boundaries, vec![0, 3]);
    }

    #[test]
    fn finalize_keep_order_with_spill() {
        // Very low threshold forces spilling.
        let mut collector = HybridBucketCollector::with_threshold(1, 1, 1, 1);

        // task_count=6, push at indices 5, 3, 1 (gaps at 0, 2, 4).
        collector
            .push(make_tagged_row("chr1", 500, 50.0), 5, 0, 0.0)
            .unwrap();
        collector
            .push(make_tagged_row("chr1", 300, 30.0), 3, 0, 0.0)
            .unwrap();
        collector
            .push(make_tagged_row("chr1", 100, 10.0), 1, 0, 0.0)
            .unwrap();

        let mut emitted: Vec<f64> = Vec::new();
        let _header = collector
            .finalize_keep_order(6, test_header_builder, |_gi, row| {
                emitted.push(row.values[0]);
                Ok(())
            })
            .unwrap();

        // Emitted in original index order: idx 1, 3, 5 → tags 10.0, 30.0, 50.0
        assert_eq!(emitted, vec![10.0, 30.0, 50.0]);
    }

    #[test]
    fn finalize_keep_order_multiple_groups() {
        let mut collector = HybridBucketCollector::with_threshold(2, 1, 1, 1_000_000);

        // task_count=4. Group 0 at indices 0, 2. Group 1 at index 3. Index 1 is filtered.
        collector
            .push(make_tagged_row("chr1", 0, 0.0), 0, 0, 0.0)
            .unwrap();
        collector
            .push(make_tagged_row("chr1", 200, 20.0), 2, 0, 0.0)
            .unwrap();
        collector
            .push(make_tagged_row("chr2", 300, 30.0), 3, 1, 0.0)
            .unwrap();

        let mut emitted: Vec<(usize, f64)> = Vec::new();
        let header = collector
            .finalize_keep_order(4, test_header_builder, |gi, row| {
                emitted.push((gi, row.values[0]));
                Ok(())
            })
            .unwrap();

        // Original order: idx 0 (g0), skip 1, idx 2 (g0), idx 3 (g1)
        assert_eq!(emitted.len(), 3);
        assert_eq!(emitted[0], (0, 0.0));   // orig_idx=0
        assert_eq!(emitted[1], (0, 20.0));  // orig_idx=2
        assert_eq!(emitted[2], (1, 30.0));  // orig_idx=3

        // Group 0 has 2 rows, group 1 has 1 row.
        assert_eq!(header.group_boundaries, vec![0, 2, 3]);
    }

    #[test]
    fn finalize_keep_order_temp_files_cleaned_up() {
        let mut collector = HybridBucketCollector::with_threshold(1, 1, 1, 1);

        for i in 0..5u32 {
            collector
                .push(make_tagged_row("chr1", i * 100, i as f64), i as usize, 0, 0.0)
                .unwrap();
        }

        // Join to capture temp paths.
        collector.join_all().unwrap();
        let temp_paths: Vec<PathBuf> = collector
            .buckets()
            .iter()
            .flat_map(|b| b.completed_flushes.iter().map(|f| f.temp_path.clone()))
            .collect();
        assert!(!temp_paths.is_empty());

        let _header = collector
            .finalize_keep_order(5, test_header_builder, |_gi, _row| Ok(()))
            .unwrap();

        for path in &temp_paths {
            assert!(!path.exists(), "temp file should be removed: {:?}", path);
        }
    }

    // -----------------------------------------------------------------------
    // Integration tests — full pipeline with spilling
    // -----------------------------------------------------------------------

    /// Build a tagged row with a specific value at index 0 (for tracking identity
    /// through the round-trip) and additional values to bulk up the row size.
    fn make_large_tagged_row(chrom: &str, start: u32, tag: f64, n_extra: usize) -> MatrixRow {
        let mut values = vec![tag];
        values.extend(std::iter::repeat(tag * 0.1).take(n_extra));
        let bin_count = 1 + n_extra;
        let record = BedRecord {
            chrom: Arc::from(chrom),
            start,
            end: start + 200,
            name: Some(format!("region_{start}")),
            score: Some(1.0),
            score_raw: None,
            strand: Strand::Unstranded,
            strand_raw: None,
            extra_fields: vec!["extra1".to_string()],
        };
        MatrixRow {
            record,
            values,
            sample_count: 1,
            bin_count,
            exon_coords: None,
        }
    }

    /// Build a tagged row that contains NaN and Inf values to test special-value
    /// preservation through serialize → spill file → mmap → deserialize.
    fn make_special_values_row(chrom: &str, start: u32, tag: f64) -> MatrixRow {
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
            values: vec![tag, f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.0],
            sample_count: 1,
            bin_count: 5,
            exon_coords: None,
        }
    }

    /// Integration test: 24 rows across 3 groups with threshold=100 bytes.
    /// Verifies:
    ///   - Multiple spills are triggered.
    ///   - finalize_sorted (ascending) emits each group sorted by sort_key.
    ///   - finalize_sorted (descending) emits each group in reverse sorted order.
    ///   - Tie-break: rows with equal sort_key appear in insertion_seq order for
    ///     ascending and reverse insertion_seq order for descending.
    #[test]
    fn integration_sorted_ascending_and_descending_with_spill() {
        // Parameters: 3 groups, 8 rows each, threshold small enough to force spills.
        // row_estimated_bytes = 1 * 4 * 8 + 100 = 132, threshold=100 → every push spills.
        let n_groups = 3usize;
        let n_rows_per_group = 8usize;

        // Pre-build expected data: sort_key and tag for each (group, row_within_group).
        // Groups 0..2, sort keys cycle through [5, 3, 7, 1, 9, 2, 8, 4] per group.
        let sort_keys = [5.0f64, 3.0, 7.0, 1.0, 9.0, 2.0, 8.0, 4.0];

        // Build ascending — push rows interleaved across groups.
        let mut collector_asc =
            HybridBucketCollector::with_threshold(n_groups, 1, 4, 100);
        let mut orig_idx = 0usize;
        for row_i in 0..n_rows_per_group {
            for group_i in 0..n_groups {
                let key = sort_keys[row_i];
                // tag encodes (group_i * 100 + row_i) so we can verify identity.
                let tag = (group_i * 100 + row_i) as f64;
                let row = make_large_tagged_row("chr1", orig_idx as u32 * 200, tag, 3);
                collector_asc.push(row, orig_idx, group_i, key).unwrap();
                orig_idx += 1;
            }
        }

        let mut emitted_asc: Vec<(usize, f64, f64)> = Vec::new(); // (group_index, sort_key, tag)
        // We need the sort key in emitted — encode as tag = group*100+row_i, and retrieve
        // sort_key from values[0] being the tag, but we need another field. Instead, we
        // re-derive sort_key from expected order after emission. We'll verify by checking
        // that within each group the emitted tags are in expected ascending sort_key order.
        let _header = collector_asc
            .finalize_sorted(true, test_header_builder, |gi, row| {
                emitted_asc.push((gi, row.values[0], 0.0)); // tag is values[0]
                Ok(())
            })
            .unwrap();

        // Verify total count.
        assert_eq!(
            emitted_asc.len(),
            n_groups * n_rows_per_group,
            "ascending: total row count mismatch"
        );

        // Verify group-contiguity and per-group ascending order.
        // Expected ascending tag order for each group: row_i sorted by sort_keys.
        // sort_keys[0..8] = [5,3,7,1,9,2,8,4], sorted indices = [3,5,1,7,0,2,6,4]
        // → row_i in ascending sort_key order: [1,5,3,7,0,2,8,4] (by values 1,2,3,4,5,7,8,9)
        let sorted_row_indices: Vec<usize> = {
            let mut pairs: Vec<(f64, usize)> = sort_keys
                .iter()
                .copied()
                .enumerate()
                .map(|(i, k)| (k, i))
                .collect();
            pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            pairs.iter().map(|&(_, i)| i).collect()
        };

        // Verify group 0 comes before group 1, group 1 before group 2.
        let group_ranges: Vec<std::ops::Range<usize>> = {
            let mut ranges = Vec::new();
            let mut start = 0;
            for g in 0..n_groups {
                let count = emitted_asc.iter().filter(|&&(gi, _, _)| gi == g).count();
                ranges.push(start..start + count);
                start += count;
            }
            ranges
        };
        for g in 0..n_groups {
            assert_eq!(
                group_ranges[g].len(),
                n_rows_per_group,
                "group {g} row count mismatch in ascending"
            );
        }

        // Check ascending sort order within each group.
        for g in 0..n_groups {
            let group_tags: Vec<f64> = emitted_asc[group_ranges[g].clone()]
                .iter()
                .map(|&(_, tag, _)| tag)
                .collect();
            let expected_tags: Vec<f64> = sorted_row_indices
                .iter()
                .map(|&row_i| (g * 100 + row_i) as f64)
                .collect();
            assert_eq!(
                group_tags, expected_tags,
                "group {g} ascending order mismatch"
            );
        }

        // Descending test — same data but sort_ascending=false.
        let mut collector_desc =
            HybridBucketCollector::with_threshold(n_groups, 1, 4, 100);
        orig_idx = 0;
        for row_i in 0..n_rows_per_group {
            for group_i in 0..n_groups {
                let key = sort_keys[row_i];
                let tag = (group_i * 100 + row_i) as f64;
                let row = make_large_tagged_row("chr1", orig_idx as u32 * 200, tag, 3);
                collector_desc.push(row, orig_idx, group_i, key).unwrap();
                orig_idx += 1;
            }
        }

        let mut emitted_desc: Vec<(usize, f64)> = Vec::new();
        let _header = collector_desc
            .finalize_sorted(false, test_header_builder, |gi, row| {
                emitted_desc.push((gi, row.values[0]));
                Ok(())
            })
            .unwrap();

        assert_eq!(emitted_desc.len(), n_groups * n_rows_per_group);

        // Descending expected: reverse of ascending.
        let desc_row_indices: Vec<usize> = sorted_row_indices.iter().rev().copied().collect();
        let group_ranges_desc: Vec<std::ops::Range<usize>> = {
            let mut ranges = Vec::new();
            let mut start = 0;
            for g in 0..n_groups {
                let count = emitted_desc.iter().filter(|&&(gi, _)| gi == g).count();
                ranges.push(start..start + count);
                start += count;
            }
            ranges
        };
        for g in 0..n_groups {
            let group_tags: Vec<f64> = emitted_desc[group_ranges_desc[g].clone()]
                .iter()
                .map(|&(_, tag)| tag)
                .collect();
            let expected_tags: Vec<f64> = desc_row_indices
                .iter()
                .map(|&row_i| (g * 100 + row_i) as f64)
                .collect();
            assert_eq!(
                group_tags, expected_tags,
                "group {g} descending order mismatch"
            );
        }
    }

    /// Integration test: keep-order mode with gaps (filtered rows) and spilling.
    /// 20 tasks with alternating group assignment; every other task is "filtered"
    /// (not pushed). Verifies emitted rows appear in original input order and
    /// the group label matches the expected assignment for each row.
    #[test]
    fn integration_keep_order_with_gaps_and_spill() {
        // threshold=100, row_estimated_bytes = 1*3*8+100 = 124 > 100 → every push spills.
        let task_count = 20usize;
        let n_groups = 2usize;
        let mut collector = HybridBucketCollector::with_threshold(n_groups, 1, 3, 100);

        // Push only even-indexed tasks (odd indices are "filtered" gaps).
        // Even tasks alternate groups: task 0 → group 0, task 2 → group 1, task 4 → group 0, …
        for task_i in (0..task_count).step_by(2) {
            let group_i = (task_i / 2) % n_groups;
            let tag = task_i as f64;
            let row = make_large_tagged_row("chr1", task_i as u32 * 100, tag, 2);
            collector.push(row, task_i, group_i, 0.0).unwrap();
        }

        let mut emitted: Vec<(usize, f64)> = Vec::new(); // (group_index, tag)
        let _header = collector
            .finalize_keep_order(task_count, test_header_builder, |gi, row| {
                emitted.push((gi, row.values[0]));
                Ok(())
            })
            .unwrap();

        // Only even-indexed tasks were pushed: 0, 2, 4, … 18 → 10 rows.
        assert_eq!(emitted.len(), 10, "expected 10 non-filtered rows");

        // Verify original input order and correct group assignment.
        for (emit_i, &(gi, tag)) in emitted.iter().enumerate() {
            let expected_task_i = emit_i * 2; // 0, 2, 4, …
            let expected_group = (emit_i) % n_groups;
            assert_eq!(
                tag, expected_task_i as f64,
                "row {emit_i}: tag mismatch (expected task {expected_task_i})"
            );
            assert_eq!(
                gi, expected_group,
                "row {emit_i}: group mismatch for task {expected_task_i}"
            );
        }
    }

    /// Integration test: special float values (NaN, Inf, -Inf, -0.0) survive the
    /// full round-trip through serialize → temp file → mmap → deserialize → emit.
    #[test]
    fn integration_special_values_survive_spill_round_trip() {
        // threshold=1 forces every push to spill immediately.
        let mut collector = HybridBucketCollector::with_threshold(1, 1, 5, 1);

        // Push 6 rows; each contains the special values pattern.
        for i in 0..6u32 {
            let row = make_special_values_row("chrM", i * 50, i as f64 * 10.0);
            collector.push(row, i as usize, 0, i as f64).unwrap();
        }

        let mut emitted_rows: Vec<MatrixRow> = Vec::new();
        let _header = collector
            .finalize_sorted(true, test_header_builder, |_gi, row| {
                emitted_rows.push(row);
                Ok(())
            })
            .unwrap();

        assert_eq!(emitted_rows.len(), 6);

        // Verify ascending sort order and special-value preservation.
        for (i, row) in emitted_rows.iter().enumerate() {
            let expected_tag = i as f64 * 10.0;
            assert_eq!(
                row.values[0], expected_tag,
                "row {i}: tag value mismatch"
            );
            assert!(row.values[1].is_nan(), "row {i}: NaN not preserved");
            assert_eq!(
                row.values[2],
                f64::INFINITY,
                "row {i}: +Inf not preserved"
            );
            assert_eq!(
                row.values[3],
                f64::NEG_INFINITY,
                "row {i}: -Inf not preserved"
            );
            // -0.0 == 0.0 in IEEE 754, so check the bit pattern.
            assert_eq!(
                row.values[4].to_bits(),
                (-0.0f64).to_bits(),
                "row {i}: -0.0 bit pattern not preserved"
            );
        }
    }

    /// Integration test: tie-break behaviour with spilling.
    /// All rows have the same sort_key=1.0; insertion_seq should determine order.
    /// Ascending: 0,1,2,3,4,… Descending: …4,3,2,1,0.
    #[test]
    fn integration_tiebreak_insertion_seq_with_spill() {
        let n_rows = 10usize;

        // Ascending tie-break.
        let mut c_asc = HybridBucketCollector::with_threshold(1, 1, 1, 100);
        for i in 0..n_rows {
            // tag encodes insertion order.
            let row = make_tagged_row("chr1", i as u32 * 100, i as f64);
            c_asc.push(row, i, 0, 1.0).unwrap(); // all same sort_key
        }
        let mut tags_asc: Vec<f64> = Vec::new();
        c_asc
            .finalize_sorted(true, test_header_builder, |_gi, row| {
                tags_asc.push(row.values[0]);
                Ok(())
            })
            .unwrap();
        let expected_asc: Vec<f64> = (0..n_rows).map(|i| i as f64).collect();
        assert_eq!(
            tags_asc, expected_asc,
            "ascending tie-break: expected insertion order"
        );

        // Descending tie-break.
        let mut c_desc = HybridBucketCollector::with_threshold(1, 1, 1, 100);
        for i in 0..n_rows {
            let row = make_tagged_row("chr1", i as u32 * 100, i as f64);
            c_desc.push(row, i, 0, 1.0).unwrap();
        }
        let mut tags_desc: Vec<f64> = Vec::new();
        c_desc
            .finalize_sorted(false, test_header_builder, |_gi, row| {
                tags_desc.push(row.values[0]);
                Ok(())
            })
            .unwrap();
        let expected_desc: Vec<f64> = (0..n_rows).rev().map(|i| i as f64).collect();
        assert_eq!(
            tags_desc, expected_desc,
            "descending tie-break: expected reverse insertion order"
        );
    }

    // -----------------------------------------------------------------------
    // sort=No group-contiguous regression test
    // -----------------------------------------------------------------------

    /// Regression test: when sort=No, the executor places rows by group (all
    /// group-0 rows before group-1 rows etc.) using finalize_keep_order.
    /// This test simulates multi-group inputs on shared chromosomes where the
    /// I/O arrival order could be interleaved between groups, and verifies that
    /// the emitted output is group-contiguous.
    ///
    /// Concretely: task indices are not grouped — row for group 1 can arrive
    /// before a later row for group 0. finalize_keep_order must still emit
    /// all group-0 rows before group-1 rows because the placement array is
    /// indexed by orig_idx and the group assignment comes from the slot, not
    /// from the emit order.
    ///
    /// NOTE: finalize_keep_order emits in orig_idx order, so true "group
    /// contiguity" requires that the caller assigns orig_idx values that are
    /// already block-separated by group (as the executor does for sort=No).
    /// This test verifies that the collector faithfully preserves whatever
    /// assignment the caller made — it does not re-sort by group itself.
    #[test]
    fn sort_no_group_contiguous_regression() {
        // Simulate the executor layout for sort=No:
        //   tasks for group 0 get orig_idx 0..n_per_group
        //   tasks for group 1 get orig_idx n_per_group..2*n_per_group
        // Rows are pushed in interleaved I/O order (group-1 row can arrive
        // before group-0 row), but finalize_keep_order must restore the
        // block layout.
        let n_per_group = 6usize;
        let n_groups = 3usize;
        let task_count = n_per_group * n_groups;

        // threshold=100, row_estimated_bytes=124 → every push spills.
        let mut collector = HybridBucketCollector::with_threshold(n_groups, 1, 3, 100);

        // Push in deliberately scrambled order: cycle through groups within
        // each "position" to simulate async I/O reordering.
        for pos in 0..n_per_group {
            for gi in (0..n_groups).rev() {
                // orig_idx is block-separated: gi * n_per_group + pos
                let orig_idx = gi * n_per_group + pos;
                let tag = orig_idx as f64;
                let row = make_large_tagged_row(
                    &format!("chr{}", gi + 1),
                    pos as u32 * 100,
                    tag,
                    2,
                );
                collector.push(row, orig_idx, gi, 0.0).unwrap();
            }
        }

        let mut emitted: Vec<(usize, f64)> = Vec::new(); // (group_index, tag)
        let header = collector
            .finalize_keep_order(task_count, test_header_builder, |gi, row| {
                emitted.push((gi, row.values[0]));
                Ok(())
            })
            .unwrap();

        assert_eq!(
            emitted.len(),
            task_count,
            "all rows should be emitted"
        );

        // Verify group-contiguity: all group-0 rows first, then group-1, then group-2.
        let mut current_group = 0usize;
        for (emit_i, &(gi, _tag)) in emitted.iter().enumerate() {
            if gi > current_group {
                // Verify we've seen exactly n_per_group rows for the previous group.
                let prev_group_count = emitted[..emit_i]
                    .iter()
                    .filter(|&&(g, _)| g == current_group)
                    .count();
                assert_eq!(
                    prev_group_count, n_per_group,
                    "group {current_group} should have exactly {n_per_group} rows before group {gi} starts"
                );
                current_group = gi;
            }
            assert!(
                gi >= current_group,
                "group index went backwards at emit position {emit_i}: got {gi}, expected >= {current_group}"
            );
        }

        // Verify the last group's count.
        let last_group_count = emitted
            .iter()
            .filter(|&&(g, _)| g == n_groups - 1)
            .count();
        assert_eq!(
            last_group_count, n_per_group,
            "last group should have exactly {n_per_group} rows"
        );

        // Verify each group's rows have the correct tags (orig_idx = gi*n_per_group + pos).
        for gi in 0..n_groups {
            let group_tags: Vec<f64> = emitted
                .iter()
                .filter(|&&(g, _)| g == gi)
                .map(|&(_, tag)| tag)
                .collect();
            let expected_tags: Vec<f64> = (0..n_per_group)
                .map(|pos| (gi * n_per_group + pos) as f64)
                .collect();
            assert_eq!(
                group_tags, expected_tags,
                "group {gi}: tags do not match expected orig_idx-derived values"
            );
        }

        // Verify the header group boundaries reflect n_per_group rows per group.
        let expected_boundaries: Vec<usize> = (0..=n_groups).map(|g| g * n_per_group).collect();
        assert_eq!(
            header.group_boundaries, expected_boundaries,
            "header group_boundaries mismatch"
        );
    }
}
