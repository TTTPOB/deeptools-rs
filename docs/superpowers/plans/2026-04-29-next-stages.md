# Next Stages Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Five improvements: remove `num_cpus` dep, add file-spilling for large sorted matrices, build a matrix comparison dev binary, migrate Python regression tests to Rust integration tests, and improve the profiling bench script.

**Architecture:** File spilling uses a hybrid per-group bucket design (in-memory below 4 GB, chunked-flush to temp file above) to eliminate OOM on large sorted matrices. Each bucket accumulates rows in memory up to ~4 GB, then flushes the entire chunk to a temp file via a dedicated writer thread (`std::thread::spawn` + move). A double-buffer scheme reuses the flushed Vec's capacity; back-pressure (join oldest flush handle) prevents unbounded memory growth on slow I/O (HDD). At finalize, the temp file is mmap'd (`memmap2`) for zero-copy sorted readback. A new dev-only `compare_matrix` binary provides subcommand-based matrix comparison (header, values, full diff) and serves as test infrastructure. The `StreamOrdered` path is extended to support auxiliary outputs (`--outFileNameMatrix`, `--outFileSortedRegions`), eliminating the forced in-memory fallback.

**Tech Stack:** Rust (edition 2024), `memmap2` for mmap, `clap` for the comparison binary CLI, existing `flate2`/`serde_json` for matrix I/O.

---

## File Map

### Task 1 (num_cpus removal)
- Modify: `src/config.rs` (replace `num_cpus::get()` with `std::thread::available_parallelism()`)
- Modify: `Cargo.toml` (remove `num_cpus` dependency)

### Task 2 (file spilling)
- Create: `src/pipeline/core/spill.rs` (HybridBucket, SpillIndex, serialization/deserialization, mmap readback)
- Modify: `src/pipeline/core/mod.rs` (re-export spill module)
- Modify: `src/pipeline/core/collector.rs` (replace `GroupBucketCollector` with `HybridBucketCollector`)
- Modify: `src/pipeline/core/executor.rs` (new `Spilling` output strategy, integrate HybridBucketCollector)
- Modify: `src/pipeline/matrix.rs` (`compute_sort_metric` made `pub(crate)`)
- Modify: `src/pipeline/run.rs` (remove forced in-memory fallback for auxiliary outputs)
- Modify: `src/io/writers/mod.rs` (remove `matrix_values_output`/`sorted_regions_output` streaming guard, add streaming auxiliary output support)
- Modify: `src/io/writers/auxiliary.rs` (add per-row streaming write functions)
- Modify: `Cargo.toml` (add `memmap2` dependency)

### Task 3 (matrix comparison binary)
- Create: `src/bin/compare_matrix.rs` (CLI entry point with subcommands)
- Create: `src/bin/compare_matrix/` directory with modules:
  - `src/bin/compare_matrix/parse.rs` (matrix loading: gzip decompress + JSON header + tab-separated rows)
  - `src/bin/compare_matrix/header.rs` (header comparison logic)
  - `src/bin/compare_matrix/values.rs` (numerical comparison with tolerance)
  - `src/bin/compare_matrix/diff.rs` (detailed per-row diff output)
- Modify: `Cargo.toml` (add `[[bin]]` target for `compare_matrix`, mark dev-only)

### Task 4 (Rust integration tests)
- Create: `tests/python_compatibility.rs` (integration tests using `compare_matrix` binary against vendored reference matrices)

### Task 5 (profile_bench.sh improvement)
- Modify: `scripts/profile_bench.sh` (add warm-cache run before profiling runs)

---

### Task 1: Remove `num_cpus` Dependency

**Files:**
- Modify: `src/config.rs:237-241`
- Modify: `Cargo.toml:10`

- [ ] **Step 1: Write a test for processor resolution**

```rust
// src/config.rs — add inside existing #[cfg(test)] or at bottom
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resolve_max_returns_at_least_one() {
        let count = ProcessorRequest::Max.resolve();
        assert!(count >= 1);
    }

    #[test]
    fn resolve_max_half_returns_at_least_one() {
        let count = ProcessorRequest::MaxHalf.resolve();
        assert!(count >= 1);
    }

    #[test]
    fn resolve_fixed_zero_clamps_to_one() {
        let count = ProcessorRequest::Fixed(0).resolve();
        assert_eq!(count, 1);
    }

    #[test]
    fn resolve_fixed_value() {
        let count = ProcessorRequest::Fixed(4).resolve();
        assert_eq!(count, 4);
    }
}
```

- [ ] **Step 2: Run test to verify it passes with current implementation**

Run: `cargo test --lib -- config::tests`
Expected: PASS (4 tests)

- [ ] **Step 3: Replace `num_cpus` with `std::thread::available_parallelism`**

Replace `available_cpus()` in `src/config.rs:237-241`:

```rust
fn available_cpus() -> u32 {
    let count = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    u32::try_from(count).unwrap_or(u32::MAX)
}
```

Remove from `Cargo.toml`:

```toml
# Remove this line:
num_cpus = "1.16"
```

- [ ] **Step 4: Run tests and build**

Run: `cargo test --lib -- config::tests && cargo build --release`
Expected: PASS, no warnings about unused `num_cpus`

- [ ] **Step 5: Commit**

```bash
git add src/config.rs Cargo.toml Cargo.lock
git commit -m "refactor: remove num_cpus dep, use std::thread::available_parallelism"
```

---

### Task 2: File Spilling for Large Sorted Matrices

This is the largest task. It is split into sub-tasks 2a–2f for manageable commits.

#### Task 2a: Add `memmap2` Dependency and Spill Module Skeleton

**Files:**
- Modify: `Cargo.toml`
- Create: `src/pipeline/core/spill.rs`
- Modify: `src/pipeline/core/mod.rs`

- [ ] **Step 1: Add `memmap2` to Cargo.toml**

```toml
memmap2 = "0.9"
```

- [ ] **Step 2: Create `src/pipeline/core/spill.rs` with types and serialization**

```rust
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::Arc;

use anyhow::{Context, Result};
use memmap2::Mmap;
use tempfile::NamedTempFile;

use crate::io::readers::bed::{BedRecord, Strand, intern_chrom};
use crate::pipeline::matrix::MatrixRow;

const SPILL_BUF_CAPACITY: usize = 1_048_576; // 1 MB — large sequential writes

/// Threshold in bytes per bucket before spilling to disk (~4 GB).
const MEMORY_SPILL_THRESHOLD: usize = 4 * 1024 * 1024 * 1024;

/// Lightweight index entry kept in memory for each spilled row.
#[derive(Clone)]
pub(crate) struct SpillIndex {
    pub(crate) orig_idx: usize,
    pub(crate) group_index: usize,
    pub(crate) sort_key: f64,
    pub(crate) file_offset: u64,
    pub(crate) row_byte_len: u32,
}

/// Intern table mapping chrom strings to compact IDs for the spill file format.
struct ChromTable {
    to_id: HashMap<Arc<str>, u16>,
    to_str: Vec<Arc<str>>,
}

impl ChromTable {
    fn new() -> Self {
        Self {
            to_id: HashMap::new(),
            to_str: Vec::new(),
        }
    }

    fn get_or_insert(&mut self, chrom: &Arc<str>) -> u16 {
        if let Some(&id) = self.to_id.get(chrom) {
            return id;
        }
        let id = self.to_str.len() as u16;
        self.to_str.push(Arc::clone(chrom));
        self.to_id.insert(Arc::clone(chrom), id);
        id
    }

    fn resolve(&self, id: u16) -> Arc<str> {
        Arc::clone(&self.to_str[id as usize])
    }
}

/// Per-row binary format written to the spill file.
/// No outer length prefix — row byte length is tracked in `SpillIndex` in memory.
///
/// ```text
/// [2 bytes: chrom_id (u16 LE)]
/// [4 bytes: start (u32 LE)]
/// [4 bytes: end (u32 LE)]
/// [1 byte:  strand (0=Positive, 1=Negative, 2=Unstranded)]
/// [1 byte:  flags (bit 0: has_name, bit 1: has_score_raw, bit 2: has_strand_raw,
///                   bit 3: has_exon_coords)]
/// [variable: name (2-byte len prefix + UTF-8, if has_name)]
/// [variable: score — if has_score_raw: 2-byte len + UTF-8, else 4 bytes f32 or 0xFF×4 for None]
/// [variable: strand_raw (2-byte len + UTF-8, if has_strand_raw)]
/// [variable: extra_fields (2-byte count, then each: 2-byte len + UTF-8)]
/// [variable: exon_coords (2-byte count, then each: 4-byte start + 4-byte end, if has_exon_coords)]
/// [2 bytes:  sample_count (u16 LE)]
/// [2 bytes:  bin_count (u16 LE)]
/// [N × 8 bytes: values (f64 LE, N = sample_count × bin_count)]
/// ```
pub(crate) fn serialize_row(buf: &mut Vec<u8>, row: &MatrixRow, chrom_table: &mut ChromTable) {
    buf.clear();

    let chrom_id = chrom_table.get_or_insert(&row.record.chrom);
    buf.extend_from_slice(&chrom_id.to_le_bytes());
    buf.extend_from_slice(&row.record.start.to_le_bytes());
    buf.extend_from_slice(&row.record.end.to_le_bytes());

    let strand_byte: u8 = match row.record.strand {
        Strand::Positive => 0,
        Strand::Negative => 1,
        Strand::Unstranded => 2,
    };
    buf.push(strand_byte);

    let has_name = row.record.name.is_some();
    let has_score_raw = row.record.score_raw.is_some();
    let has_strand_raw = row.record.strand_raw.is_some();
    let has_exon_coords = row.exon_coords.is_some();
    let flags: u8 = (has_name as u8)
        | ((has_score_raw as u8) << 1)
        | ((has_strand_raw as u8) << 2)
        | ((has_exon_coords as u8) << 3);
    buf.push(flags);

    if let Some(ref name) = row.record.name {
        let bytes = name.as_bytes();
        buf.extend_from_slice(&(bytes.len() as u16).to_le_bytes());
        buf.extend_from_slice(bytes);
    }

    if let Some(ref raw) = row.record.score_raw {
        let bytes = raw.as_bytes();
        buf.extend_from_slice(&(bytes.len() as u16).to_le_bytes());
        buf.extend_from_slice(bytes);
    } else if let Some(score) = row.record.score {
        buf.extend_from_slice(&score.to_le_bytes());
    } else {
        buf.extend_from_slice(&[0xFF, 0xFF, 0xFF, 0xFF]);
    }

    if let Some(ref raw) = row.record.strand_raw {
        let bytes = raw.as_bytes();
        buf.extend_from_slice(&(bytes.len() as u16).to_le_bytes());
        buf.extend_from_slice(bytes);
    }

    let extra = &row.record.extra_fields;
    buf.extend_from_slice(&(extra.len() as u16).to_le_bytes());
    for field in extra {
        let bytes = field.as_bytes();
        buf.extend_from_slice(&(bytes.len() as u16).to_le_bytes());
        buf.extend_from_slice(bytes);
    }

    if let Some(ref coords) = row.exon_coords {
        buf.extend_from_slice(&(coords.len() as u16).to_le_bytes());
        for &(start, end) in coords {
            buf.extend_from_slice(&start.to_le_bytes());
            buf.extend_from_slice(&end.to_le_bytes());
        }
    }

    buf.extend_from_slice(&(row.sample_count as u16).to_le_bytes());
    buf.extend_from_slice(&(row.bin_count as u16).to_le_bytes());
    for &value in &row.values {
        buf.extend_from_slice(&value.to_le_bytes());
    }
}

/// Deserialize a row from a byte slice (typically an mmap region).
pub(crate) fn deserialize_row(data: &[u8], chrom_table: &ChromTable) -> MatrixRow {
    let mut cursor = 0usize;

    let chrom_id = u16::from_le_bytes([data[cursor], data[cursor + 1]]);
    cursor += 2;
    let chrom = chrom_table.resolve(chrom_id);

    let start = u32::from_le_bytes(data[cursor..cursor + 4].try_into().unwrap());
    cursor += 4;
    let end = u32::from_le_bytes(data[cursor..cursor + 4].try_into().unwrap());
    cursor += 4;

    let strand = match data[cursor] {
        0 => Strand::Positive,
        1 => Strand::Negative,
        _ => Strand::Unstranded,
    };
    cursor += 1;

    let flags = data[cursor];
    cursor += 1;
    let has_name = flags & 1 != 0;
    let has_score_raw = flags & 2 != 0;
    let has_strand_raw = flags & 4 != 0;
    let has_exon_coords = flags & 8 != 0;

    let name = if has_name {
        let len = u16::from_le_bytes([data[cursor], data[cursor + 1]]) as usize;
        cursor += 2;
        let s = std::str::from_utf8(&data[cursor..cursor + len]).unwrap().to_string();
        cursor += len;
        Some(s)
    } else {
        None
    };

    let (score, score_raw) = if has_score_raw {
        let len = u16::from_le_bytes([data[cursor], data[cursor + 1]]) as usize;
        cursor += 2;
        let s = std::str::from_utf8(&data[cursor..cursor + len]).unwrap().to_string();
        cursor += len;
        (None, Some(s))
    } else {
        let bytes: [u8; 4] = data[cursor..cursor + 4].try_into().unwrap();
        cursor += 4;
        if bytes == [0xFF, 0xFF, 0xFF, 0xFF] {
            (None, None)
        } else {
            (Some(f32::from_le_bytes(bytes)), None)
        }
    };

    let strand_raw = if has_strand_raw {
        let len = u16::from_le_bytes([data[cursor], data[cursor + 1]]) as usize;
        cursor += 2;
        let s = std::str::from_utf8(&data[cursor..cursor + len]).unwrap().to_string();
        cursor += len;
        Some(s)
    } else {
        None
    };

    let extra_count = u16::from_le_bytes([data[cursor], data[cursor + 1]]) as usize;
    cursor += 2;
    let mut extra_fields = Vec::with_capacity(extra_count);
    for _ in 0..extra_count {
        let len = u16::from_le_bytes([data[cursor], data[cursor + 1]]) as usize;
        cursor += 2;
        extra_fields.push(std::str::from_utf8(&data[cursor..cursor + len]).unwrap().to_string());
        cursor += len;
    }

    let exon_coords = if has_exon_coords {
        let count = u16::from_le_bytes([data[cursor], data[cursor + 1]]) as usize;
        cursor += 2;
        let mut coords = Vec::with_capacity(count);
        for _ in 0..count {
            let s = u32::from_le_bytes(data[cursor..cursor + 4].try_into().unwrap());
            cursor += 4;
            let e = u32::from_le_bytes(data[cursor..cursor + 4].try_into().unwrap());
            cursor += 4;
            coords.push((s, e));
        }
        Some(coords)
    } else {
        None
    };

    let sample_count = u16::from_le_bytes([data[cursor], data[cursor + 1]]) as usize;
    cursor += 2;
    let bin_count = u16::from_le_bytes([data[cursor], data[cursor + 1]]) as usize;
    cursor += 2;

    let value_count = sample_count * bin_count;
    let mut values = Vec::with_capacity(value_count);
    for _ in 0..value_count {
        let v = f64::from_le_bytes(data[cursor..cursor + 8].try_into().unwrap());
        cursor += 8;
        values.push(v);
    }

    MatrixRow {
        record: BedRecord {
            chrom,
            start,
            end,
            name,
            score,
            score_raw,
            strand,
            strand_raw,
            extra_fields,
        },
        values,
        sample_count,
        bin_count,
        exon_coords,
    }
}

/// Result returned by a writer thread after flushing a chunk to disk.
struct FlushResult {
    indices: Vec<SpillIndex>,
    returned_buf: Vec<(usize, f64, MatrixRow)>,
    chrom_table: ChromTable,
    temp_path: std::path::PathBuf,
}

/// Flush a chunk of rows to a new temp file. Runs on a dedicated thread.
fn flush_chunk(
    rows: Vec<(usize, f64, MatrixRow)>,
    group_index: usize,
) -> Result<FlushResult> {
    let temp = NamedTempFile::new().context("Failed to create spill temp file")?;
    let temp_path = temp.path().to_path_buf();
    let file = temp.into_file();
    let mut writer = BufWriter::with_capacity(SPILL_BUF_CAPACITY, file);
    let mut chrom_table = ChromTable::new();
    let mut index = Vec::with_capacity(rows.len());
    let mut serialize_buf = Vec::with_capacity(32768);
    let mut offset: u64 = 0;

    for &(orig_idx, sort_key, ref row) in &rows {
        serialize_row(&mut serialize_buf, row, &mut chrom_table);
        writer
            .write_all(&serialize_buf)
            .context("Failed to write row to spill file")?;
        index.push(SpillIndex {
            orig_idx,
            group_index,
            sort_key,
            file_offset: offset,
            row_byte_len: serialize_buf.len() as u32,
        });
        offset += serialize_buf.len() as u64;
    }

    writer.flush().context("Failed to flush spill file")?;

    // Clear the Vec but keep its capacity for reuse
    let mut returned_buf = rows;
    returned_buf.clear();

    Ok(FlushResult {
        indices: index,
        returned_buf,
        chrom_table,
        temp_path,
    })
}

/// A single per-group bucket using chunked flush with double-buffering.
///
/// Accumulates rows in memory up to ~4 GB, then flushes the entire chunk
/// to a temp file on a dedicated writer thread. A spare buffer is reused
/// to avoid re-allocation. Back-pressure is applied when the spare buffer
/// is not yet returned (join oldest flush handle).
struct CollectorBucket {
    active: Vec<(usize, f64, MatrixRow)>,
    spare: Option<Vec<(usize, f64, MatrixRow)>>,
    estimated_bytes: usize,
    group_index: usize,
    flush_handles: Vec<std::thread::JoinHandle<Result<FlushResult>>>,
    completed_flushes: Vec<FlushResult>,
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
        }
    }

    fn push(
        &mut self,
        row: MatrixRow,
        orig_idx: usize,
        sort_key: f64,
    ) -> Result<()> {
        let row_bytes = row.values.len() * 8 + 200;
        self.estimated_bytes += row_bytes;

        if self.estimated_bytes > MEMORY_SPILL_THRESHOLD && !self.active.is_empty() {
            self.flush_active()?;
            self.estimated_bytes = row_bytes;
        }

        self.active.push((orig_idx, sort_key, row));
        Ok(())
    }

    fn flush_active(&mut self) -> Result<()> {
        // Get spare buffer (reuse capacity) or apply back-pressure
        let spare = match self.spare.take() {
            Some(buf) => buf,
            None if !self.flush_handles.is_empty() => {
                // Back-pressure: wait for oldest flush to complete
                let handle = self.flush_handles.remove(0);
                let result = handle.join().unwrap()?;
                let buf = result.returned_buf;
                self.completed_flushes.push(FlushResult {
                    returned_buf: Vec::new(), // capacity already taken
                    ..result
                });
                buf
            }
            None => Vec::new(), // first flush, no spare yet
        };

        let full = std::mem::replace(&mut self.active, spare);
        let group_index = self.group_index;
        let handle = std::thread::spawn(move || flush_chunk(full, group_index));
        self.flush_handles.push(handle);
        Ok(())
    }

    /// Join all pending flush handles and return total row count
    /// (in-memory + all flushed chunks).
    fn join_all(&mut self) -> Result<usize> {
        for handle in self.flush_handles.drain(..) {
            let result = handle.join().unwrap()?;
            self.completed_flushes.push(result);
        }
        let spilled: usize = self.completed_flushes.iter().map(|f| f.indices.len()).sum();
        Ok(self.active.len() + spilled)
    }
}

/// Collector that manages multiple hybrid buckets — one per group for
/// ascend/descend sorting, or a single bucket for sort=keep reordering.
///
/// Each bucket accumulates rows in memory up to ~4 GB, then flushes the
/// chunk to a temp file on a writer thread. Double-buffering reuses Vec
/// capacity; back-pressure prevents unbounded memory growth on slow I/O.
pub(crate) struct HybridBucketCollector {
    buckets: Vec<CollectorBucket>,
    sample_count: usize,
    bin_count: usize,
}

impl HybridBucketCollector {
    pub(crate) fn new(group_count: usize, sample_count: usize, bin_count: usize) -> Self {
        Self {
            buckets: (0..group_count)
                .map(|i| CollectorBucket::new(i))
                .collect(),
            sample_count,
            bin_count,
        }
    }

    pub(crate) fn push(
        &mut self,
        row: MatrixRow,
        orig_idx: usize,
        group_index: usize,
        sort_key: f64,
    ) -> Result<()> {
        self.buckets[group_index].push(row, orig_idx, sort_key)
    }
}
```

This file will be extended in subsequent sub-tasks with the finalize/readback logic.

- [ ] **Step 3: Add module to `src/pipeline/core/mod.rs`**

Add after `mod coalesce;`:

```rust
pub(crate) mod spill;
```

And add to the pub use block:

```rust
pub use spill::HybridBucketCollector;
```

- [ ] **Step 4: Verify it compiles**

Run: `cargo check`
Expected: PASS (warnings about dead code are acceptable at this stage)

- [ ] **Step 5: Write round-trip serialization test**

Add to the bottom of `src/pipeline/core/spill.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::readers::bed::{BedRecord, Strand};
    use std::sync::Arc;

    fn sample_row() -> MatrixRow {
        MatrixRow {
            record: BedRecord {
                chrom: Arc::from("chr1"),
                start: 1000,
                end: 2000,
                name: Some("gene_A".to_string()),
                score: Some(3.14),
                score_raw: None,
                strand: Strand::Positive,
                strand_raw: None,
                extra_fields: vec!["extra1".to_string(), "extra2".to_string()],
            },
            values: vec![1.0, 2.5, f64::NAN, 0.0, -1.0, 99.999999],
            sample_count: 2,
            bin_count: 3,
            exon_coords: Some(vec![(100, 200), (300, 400)]),
        }
    }

    #[test]
    fn serialize_deserialize_round_trip() {
        let row = sample_row();
        let mut chrom_table = ChromTable::new();
        let mut buf = Vec::new();
        serialize_row(&mut buf, &row, &mut chrom_table);

        let restored = deserialize_row(&buf, &chrom_table);

        assert_eq!(restored.record.chrom.as_ref(), "chr1");
        assert_eq!(restored.record.start, 1000);
        assert_eq!(restored.record.end, 2000);
        assert_eq!(restored.record.name.as_deref(), Some("gene_A"));
        assert_eq!(restored.record.score, Some(3.14));
        assert_eq!(restored.record.score_raw, None);
        assert!(matches!(restored.record.strand, Strand::Positive));
        assert_eq!(restored.record.extra_fields, vec!["extra1", "extra2"]);
        assert_eq!(restored.sample_count, 2);
        assert_eq!(restored.bin_count, 3);
        assert_eq!(restored.exon_coords, Some(vec![(100, 200), (300, 400)]));

        // Check values (NaN needs special comparison)
        assert_eq!(restored.values.len(), 6);
        assert_eq!(restored.values[0], 1.0);
        assert_eq!(restored.values[1], 2.5);
        assert!(restored.values[2].is_nan());
        assert_eq!(restored.values[3], 0.0);
        assert_eq!(restored.values[4], -1.0);
        assert_eq!(restored.values[5], 99.999999);
    }

    #[test]
    fn serialize_minimal_row() {
        let row = MatrixRow {
            record: BedRecord {
                chrom: Arc::from("chrX"),
                start: 0,
                end: 100,
                name: None,
                score: None,
                score_raw: None,
                strand: Strand::Unstranded,
                strand_raw: None,
                extra_fields: vec![],
            },
            values: vec![42.0],
            sample_count: 1,
            bin_count: 1,
            exon_coords: None,
        };

        let mut chrom_table = ChromTable::new();
        let mut buf = Vec::new();
        serialize_row(&mut buf, &row, &mut chrom_table);

        let restored = deserialize_row(&buf, &chrom_table);

        assert_eq!(restored.record.chrom.as_ref(), "chrX");
        assert_eq!(restored.record.name, None);
        assert_eq!(restored.record.score, None);
        assert_eq!(restored.exon_coords, None);
        assert_eq!(restored.values, vec![42.0]);
    }

    #[test]
    fn chrom_table_interning() {
        let mut table = ChromTable::new();
        let chr1: Arc<str> = Arc::from("chr1");
        let chr2: Arc<str> = Arc::from("chr2");

        assert_eq!(table.get_or_insert(&chr1), 0);
        assert_eq!(table.get_or_insert(&chr2), 1);
        assert_eq!(table.get_or_insert(&chr1), 0); // same ID

        assert_eq!(table.resolve(0).as_ref(), "chr1");
        assert_eq!(table.resolve(1).as_ref(), "chr2");
    }
}
```

- [ ] **Step 6: Run tests**

Run: `cargo test --lib -- spill::tests`
Expected: PASS (3 tests)

- [ ] **Step 7: Commit**

```bash
git add Cargo.toml Cargo.lock src/pipeline/core/spill.rs src/pipeline/core/mod.rs
git commit -m "feat: add spill module with hybrid bucket serialization and chrom interning"
```

#### Task 2b: Finalize and Readback Logic for HybridBucketCollector

**Files:**
- Modify: `src/pipeline/core/spill.rs` (add finalize methods with mmap readback)

- [ ] **Step 1: Add finalize methods to `HybridBucketCollector`**

Append to the `impl HybridBucketCollector` block in `spill.rs`:

```rust
    /// Finalize for ascend/descend: join all flush handles, sort each
    /// bucket by sort_key, emit rows via the provided callback in group order.
    pub(crate) fn finalize_sorted<F>(
        mut self,
        sort_ascending: bool,
        header_builder: impl FnOnce(Vec<usize>) -> Result<MatrixHeader>,
        mut emit: F,
    ) -> Result<MatrixHeader>
    where
        F: FnMut(MatrixRow) -> Result<()>,
    {
        // Join all pending writer threads
        let mut group_counts = Vec::with_capacity(self.buckets.len());
        for bucket in &mut self.buckets {
            group_counts.push(bucket.join_all()?);
        }
        let header = header_builder(group_counts)?;

        for bucket in self.buckets {
            Self::emit_bucket_sorted(bucket, sort_ascending, &mut emit)?;
        }
        Ok(header)
    }

    /// Finalize for sort=keep: join all flush handles, emit all rows
    /// across all buckets in orig_idx order.
    pub(crate) fn finalize_keep_order<F>(
        mut self,
        header_builder: impl FnOnce(Vec<usize>) -> Result<MatrixHeader>,
        mut emit: F,
    ) -> Result<MatrixHeader>
    where
        F: FnMut(MatrixRow) -> Result<()>,
    {
        let mut group_counts = Vec::with_capacity(self.buckets.len());
        for bucket in &mut self.buckets {
            group_counts.push(bucket.join_all()?);
        }
        let header = header_builder(group_counts)?;
        let total: usize = group_counts.iter().sum();

        // Full implementation in Task 2c — uses pre-allocated placement array
        // indexed by orig_idx for O(n) emit with zero sorting overhead.
        todo!("Complete keep-order finalize — see Task 2c")
    }

    fn emit_bucket_sorted<F>(
        bucket: CollectorBucket,
        sort_ascending: bool,
        emit: &mut F,
    ) -> Result<()>
    where
        F: FnMut(MatrixRow) -> Result<()>,
    {
        // Wrapper pairing a SpillIndex with its mmap index
        struct SpillRef {
            sort_key: f64,
            mmap_idx: usize,
            file_offset: u64,
            row_byte_len: u32,
        }

        // Collect all spill refs + mmap all chunk files
        let mut all_refs: Vec<SpillRef> = Vec::new();
        let mut mmaps: Vec<(Mmap, ChromTable)> = Vec::new();
        let mut temp_paths: Vec<std::path::PathBuf> = Vec::new();

        for flush in bucket.completed_flushes {
            let file = File::open(&flush.temp_path)
                .context("Failed to reopen spill file for mmap")?;
            let mmap = unsafe { Mmap::map(&file) }
                .context("Failed to mmap spill file")?;
            let mmap_idx = mmaps.len();

            for entry in flush.indices {
                all_refs.push(SpillRef {
                    sort_key: entry.sort_key,
                    mmap_idx,
                    file_offset: entry.file_offset,
                    row_byte_len: entry.row_byte_len,
                });
            }
            mmaps.push((mmap, flush.chrom_table));
            temp_paths.push(flush.temp_path);
        }

        // Merge in-memory rows: assign sort_key, add to a combined list
        // In-memory rows don't need mmap — emit directly after sorting
        let mut in_memory_rows = bucket.active;

        // Sort everything by sort_key
        all_refs.sort_by(|a, b| {
            let cmp = crate::pipeline::matrix::compare_ascending(a.sort_key, b.sort_key);
            if sort_ascending { cmp } else { cmp.reverse() }
        });
        in_memory_rows.sort_by(|a, b| {
            let cmp = crate::pipeline::matrix::compare_ascending(a.1, b.1);
            if sort_ascending { cmp } else { cmp.reverse() }
        });

        // Merge-emit spilled and in-memory rows in sorted order
        let mut spill_iter = all_refs.iter().peekable();
        let mut mem_iter = in_memory_rows.into_iter().peekable();

        loop {
            let pick_spill = match (spill_iter.peek(), mem_iter.peek()) {
                (Some(s), Some(m)) => {
                    let cmp = crate::pipeline::matrix::compare_ascending(s.sort_key, m.1);
                    if sort_ascending { cmp.is_le() } else { cmp.is_ge() }
                }
                (Some(_), None) => true,
                (None, Some(_)) => false,
                (None, None) => break,
            };
            if pick_spill {
                let entry = spill_iter.next().unwrap();
                let (ref mmap, ref ct) = mmaps[entry.mmap_idx];
                let start = entry.file_offset as usize;
                let end = start + entry.row_byte_len as usize;
                emit(deserialize_row(&mmap[start..end], ct))?;
            } else {
                let (_, _, row) = mem_iter.next().unwrap();
                emit(row)?;
            }
        }

        // Cleanup temp files
        drop(mmaps);
        for path in temp_paths {
            let _ = std::fs::remove_file(path);
        }

        Ok(())
    }
```

Note: `compare_ascending` in `src/pipeline/matrix.rs` needs to be made `pub(crate)`:

```rust
// Change from:
fn compare_ascending(left: f64, right: f64) -> Ordering {
// To:
pub(crate) fn compare_ascending(left: f64, right: f64) -> Ordering {
```

- [ ] **Step 2: Run tests**

Run: `cargo check`
Expected: PASS (the `todo!()` in `finalize_keep_order` is acceptable — it's completed in Task 2c)

- [ ] **Step 3: Commit**

```bash
git add src/pipeline/core/spill.rs src/pipeline/matrix.rs
git commit -m "feat: add finalize_sorted with mmap readback for hybrid buckets"
```

#### Task 2c: Complete Keep-Order Finalize and Integration with Executor

**Files:**
- Modify: `src/pipeline/core/spill.rs` (complete `finalize_keep_order`)
- Modify: `src/pipeline/core/executor.rs` (replace `InMemoryKeep` and `InMemoryGroupBucket` with spilling strategy)
- Modify: `src/pipeline/matrix.rs` (make `compute_sort_metric` pub(crate))
- Modify: `src/pipeline/core/collector.rs` (remove `GroupBucketCollector`)

- [ ] **Step 1: Complete `finalize_keep_order` in `spill.rs`**

Replace the `todo!()` in `finalize_keep_order` with the full implementation using the pre-allocated placement approach. This uses the new chunked-flush design — `bucket.active` for in-memory rows, `bucket.completed_flushes` for on-disk chunks (already joined by `join_all()` in the earlier block):

```rust
    pub(crate) fn finalize_keep_order<F>(
        mut self,
        header_builder: impl FnOnce(Vec<usize>) -> Result<MatrixHeader>,
        mut emit: F,
    ) -> Result<MatrixHeader>
    where
        F: FnMut(MatrixRow) -> Result<()>,
    {
        // join_all() already called above — group_counts and header built
        let mut group_counts = Vec::with_capacity(self.buckets.len());
        for bucket in &mut self.buckets {
            group_counts.push(bucket.join_all()?);
        }
        let header = header_builder(group_counts)?;
        let total: usize = group_counts.iter().sum();

        // Pre-allocate placement array indexed by orig_idx
        let mut slots: Vec<Option<RowSlot>> = (0..total).map(|_| None).collect();

        // Mmap all chunk files, record slot references
        let mut mmaps: Vec<(Mmap, ChromTable)> = Vec::new();
        let mut temp_paths: Vec<std::path::PathBuf> = Vec::new();
        let mut in_memory_rows: Vec<Option<MatrixRow>> = Vec::with_capacity(total);
        in_memory_rows.resize_with(total, || None);

        for bucket in self.buckets {
            // Place spilled rows from completed flush chunks
            for flush in bucket.completed_flushes {
                let file = File::open(&flush.temp_path)
                    .context("Failed to reopen spill file for mmap")?;
                let mmap = unsafe { Mmap::map(&file) }
                    .context("Failed to mmap spill file")?;
                let mmap_idx = mmaps.len();

                for entry in &flush.indices {
                    slots[entry.orig_idx] = Some(RowSlot::Spilled {
                        mmap_idx,
                        offset: entry.file_offset,
                        len: entry.row_byte_len,
                    });
                }
                mmaps.push((mmap, flush.chrom_table));
                temp_paths.push(flush.temp_path);
            }

            // Place in-memory rows (from the last unflushed active buffer)
            for (orig_idx, _sort_key, row) in bucket.active {
                in_memory_rows[orig_idx] = Some(row);
            }
        }

        // Emit in orig_idx order — O(n) scan, no sort needed
        for idx in 0..total {
            match &slots[idx] {
                Some(RowSlot::Spilled { mmap_idx, offset, len }) => {
                    let (ref mmap, ref chrom_table) = mmaps[*mmap_idx];
                    let start = *offset as usize;
                    let end = start + *len as usize;
                    let row = deserialize_row(&mmap[start..end], chrom_table);
                    emit(row)?;
                }
                None => {
                    if let Some(row) = in_memory_rows[idx].take() {
                        emit(row)?;
                    }
                }
            }
        }

        // Cleanup: drop mmaps first, then delete temp files
        drop(mmaps);
        for path in temp_paths {
            let _ = std::fs::remove_file(path);
        }

        Ok(header)
    }
```

Add the `RowSlot` enum:

```rust
enum RowSlot {
    Spilled {
        mmap_idx: usize,
        offset: u64,
        len: u32,
    },
}
```

- [ ] **Step 2: Make `compute_sort_metric` pub(crate) in `src/pipeline/matrix.rs`**

```rust
// Change from:
fn compute_sort_metric(
// To:
pub(crate) fn compute_sort_metric(
```

Also make `collect_values` pub(crate):

```rust
pub(crate) fn collect_values(row: &MatrixRow, sample_list: Option<&[usize]>) -> Vec<f64> {
```

- [ ] **Step 3: Integrate into `executor.rs`**

Replace `OutputStrategy::InMemoryKeep` and `OutputStrategy::InMemoryGroupBucket` branches in `execute_mode`. The new `Spilling` variant uses `HybridBucketCollector`:

In `OutputStrategy` enum, replace:

```rust
enum OutputStrategy {
    StreamOrdered,
    Spilling,    // replaces both InMemoryKeep and InMemoryGroupBucket
}
```

Update the strategy selection (line 102-107):

```rust
    let output_strategy = match general.sort_regions {
        SortRegions::Keep if already_sorted => OutputStrategy::StreamOrdered,
        SortRegions::No => OutputStrategy::StreamOrdered,
        _ => OutputStrategy::Spilling,
    };
```

Replace the `InMemoryKeep` and `InMemoryGroupBucket` match arms with a single `Spilling` arm that:
1. Creates `HybridBucketCollector::new(group_count, sample_count, total_bins)`
2. In the compute loop, for each result `(orig_idx, group_index, row)`:
   - Computes `sort_key` via `compute_sort_metric(&row, sort_using, sort_sample_list)` (for ascend/descend) or uses `0.0` placeholder (for keep)
   - Calls `hybrid_collector.push(row, orig_idx, group_index, sort_key)`
3. After the compute loop, calls either `finalize_sorted` or `finalize_keep_order` depending on `sort_regions`
4. The emit callback writes to the passed-in `collector`

- [ ] **Step 4: Remove `GroupBucketCollector` from `collector.rs`**

Remove the `GroupBucketCollector` struct and its `impl` block. Update `mod.rs` to remove the re-export.

- [ ] **Step 5: Run all tests**

Run: `cargo test`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add src/pipeline/core/spill.rs src/pipeline/core/executor.rs \
       src/pipeline/core/collector.rs src/pipeline/core/mod.rs \
       src/pipeline/matrix.rs
git commit -m "feat: integrate hybrid bucket collector into executor, replace InMemoryKeep and GroupBucket"
```

#### Task 2d: Extend StreamOrdered to Support Auxiliary Outputs

**Files:**
- Modify: `src/io/writers/mod.rs` (remove streaming guard, add auxiliary streaming support)
- Modify: `src/io/writers/auxiliary.rs` (add per-row streaming functions)
- Modify: `src/pipeline/core/collector.rs` (extend `FileCollector` with auxiliary writers)
- Modify: `src/pipeline/run.rs` (remove forced in-memory fallback)

- [ ] **Step 1: Add per-row streaming functions to `auxiliary.rs`**

```rust
use crate::pipeline::matrix::MatrixRow;

/// Write a single sorted-region BED line for streaming output.
pub fn write_sorted_region_row<W: Write>(
    writer: &mut W,
    row: &MatrixRow,
    group_label: &str,
) -> Result<()> {
    let name = row.record.name.as_deref().unwrap_or(".");
    let score = row
        .record
        .score_raw
        .as_deref()
        .map(|s| s.to_string())
        .or_else(|| row.record.score.map(|s| format!("{s:.6}")))
        .unwrap_or_else(|| ".".to_string());
    let strand = row
        .record
        .strand_raw
        .as_deref()
        .map(|s| s.to_string())
        .unwrap_or_else(|| row.record.strand.as_char().to_string());
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
        row.record.chrom, row.record.start, row.record.end, name, score, strand, group_label
    )?;
    Ok(())
}
```

- [ ] **Step 2: Extend `FileCollector` to optionally write auxiliary outputs**

In `collector.rs`, modify `FileCollector` to hold optional auxiliary writers:

```rust
pub struct FileCollector {
    writer: StreamingMatrixWriter,
    values_writer: Option<BufWriter<File>>,
    regions_writer: Option<BufWriter<File>>,
    group_labels: Vec<String>,
    group_boundaries_builder: Vec<usize>, // running group counts
    current_row_index: usize,
}
```

The `on_row` method writes to all active writers simultaneously.

- [ ] **Step 3: Remove the streaming guard in `src/io/writers/mod.rs`**

Remove lines 52-54 from `should_use_streaming_for_plan`:

```rust
// Remove:
if io.matrix_values_output.is_some() || io.sorted_regions_output.is_some() {
    return false;
}
```

- [ ] **Step 4: Update `run.rs` to pass auxiliary paths to FileCollector**

When creating `FileCollector`, also pass `io.matrix_values_output` and `io.sorted_regions_output` so it can open those writers.

- [ ] **Step 5: Run tests**

Run: `cargo test`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add src/io/writers/mod.rs src/io/writers/auxiliary.rs \
       src/pipeline/core/collector.rs src/pipeline/run.rs
git commit -m "feat: extend streaming path to support auxiliary outputs (outFileNameMatrix, outFileSortedRegions)"
```

#### Task 2e: Wire Spilling Path Auxiliary Outputs

**Files:**
- Modify: `src/pipeline/core/spill.rs` (emit callback for auxiliary outputs)
- Modify: `src/pipeline/core/executor.rs` (pass auxiliary writers through the spilling finalize)

- [ ] **Step 1: Update executor's Spilling branch**

In the `Spilling` match arm, after calling `finalize_sorted` or `finalize_keep_order`, the emit callback should write to all outputs (main gzip + auxiliary). This is already handled by emitting into the `collector` which now has auxiliary writers from Task 2d.

Verify the emit callback in the Spilling arm calls `collector.on_row(row)` which writes to all outputs.

- [ ] **Step 2: Run full test suite**

Run: `cargo test`
Expected: PASS

- [ ] **Step 3: Commit**

```bash
git add src/pipeline/core/executor.rs src/pipeline/core/spill.rs
git commit -m "feat: wire spilling path through auxiliary-aware collector"
```

#### Task 2f: Integration Test for File Spilling

**Files:**
- Modify: `src/pipeline/core/spill.rs` (add integration-level tests)

- [ ] **Step 1: Add test for spilling threshold transition**

```rust
#[cfg(test)]
mod integration_tests {
    use super::*;
    use crate::io::readers::bed::{BedRecord, Strand};
    use std::sync::Arc;

    fn make_row(chrom: &str, start: u32, values: Vec<f64>) -> MatrixRow {
        MatrixRow {
            record: BedRecord {
                chrom: Arc::from(chrom),
                start,
                end: start + 100,
                name: Some(format!("region_{start}")),
                score: None,
                score_raw: None,
                strand: Strand::Unstranded,
                strand_raw: None,
                extra_fields: vec![],
            },
            values,
            sample_count: 1,
            bin_count: 3,
            exon_coords: None,
        }
    }

    #[test]
    fn hybrid_bucket_stays_in_memory_for_small_data() {
        let mut collector = HybridBucketCollector::new(2, 1, 3);
        collector
            .push(make_row("chr1", 100, vec![1.0, 2.0, 3.0]), 0, 0, 2.0)
            .unwrap();
        collector
            .push(make_row("chr1", 200, vec![4.0, 5.0, 6.0]), 1, 1, 5.0)
            .unwrap();

        let mut emitted = Vec::new();
        let header = collector
            .finalize_sorted(true, |counts| {
                Ok(crate::pipeline::matrix::MatrixHeader::default_for_test(counts))
            }, |row| {
                emitted.push(row);
                Ok(())
            })
            .unwrap();

        assert_eq!(emitted.len(), 2);
        // Group 0 first (1 row), then group 1 (1 row)
        assert_eq!(emitted[0].record.start, 100);
        assert_eq!(emitted[1].record.start, 200);
    }
}
```

- [ ] **Step 2: Run tests**

Run: `cargo test --lib -- spill`
Expected: PASS

- [ ] **Step 3: Commit**

```bash
git add src/pipeline/core/spill.rs
git commit -m "test: add integration tests for hybrid bucket collector"
```

---

### Task 3: Matrix Comparison Dev Binary

**Files:**
- Create: `src/bin/compare_matrix.rs` (main entry point)
- Create: `src/bin/compare_matrix/parse.rs` (matrix loading)
- Create: `src/bin/compare_matrix/header.rs` (header comparison)
- Create: `src/bin/compare_matrix/values.rs` (numerical comparison)
- Create: `src/bin/compare_matrix/diff.rs` (detailed diff)
- Modify: `Cargo.toml` (add binary target)

- [ ] **Step 1: Add binary target to Cargo.toml**

```toml
[[bin]]
name = "compare_matrix"
path = "src/bin/compare_matrix.rs"
```

- [ ] **Step 2: Create `src/bin/compare_matrix/parse.rs`**

Matrix loading module — reads `.mat.gz` files:

```rust
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result, bail};
use flate2::read::GzDecoder;

pub struct Matrix {
    pub header: serde_json::Value,
    pub rows: Vec<MatrixFileRow>,
}

pub struct MatrixFileRow {
    pub bed_fields: Vec<String>, // chrom, start, end, name, score, strand
    pub values: Vec<f64>,
}

pub fn load_matrix(path: &Path) -> Result<Matrix> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open matrix file '{}'", path.display()))?;
    let decoder = GzDecoder::new(BufReader::new(file));
    let reader = BufReader::new(decoder);

    let mut lines = reader.lines();

    // First line: @{JSON header}
    let header_line = lines
        .next()
        .context("Matrix file is empty")?
        .context("Failed to read header line")?;
    let header_json = header_line
        .strip_prefix('@')
        .context("Header line must start with '@'")?
        .trim();
    let header: serde_json::Value =
        serde_json::from_str(header_json).context("Failed to parse matrix header JSON")?;

    let mut rows = Vec::new();
    for line in lines {
        let line = line.context("Failed to read matrix row")?;
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 7 {
            bail!("Matrix row has fewer than 7 fields: {}", line);
        }
        let bed_fields: Vec<String> = fields[..6].iter().map(|s| s.to_string()).collect();
        let values: Vec<f64> = fields[6..]
            .iter()
            .map(|s| {
                if *s == "nan" {
                    f64::NAN
                } else if *s == "inf" {
                    f64::INFINITY
                } else if *s == "-inf" {
                    f64::NEG_INFINITY
                } else {
                    s.parse::<f64>().unwrap_or(f64::NAN)
                }
            })
            .collect();
        rows.push(MatrixFileRow { bed_fields, values });
    }

    Ok(Matrix { header, rows })
}
```

- [ ] **Step 3: Create `src/bin/compare_matrix/header.rs`**

```rust
use anyhow::Result;

pub struct HeaderDiff {
    pub key: String,
    pub left: String,
    pub right: String,
}

pub fn compare_headers(
    left: &serde_json::Value,
    right: &serde_json::Value,
    ignore_keys: &[&str],
) -> Vec<HeaderDiff> {
    let mut diffs = Vec::new();
    let left_obj = left.as_object();
    let right_obj = right.as_object();

    let (Some(left_obj), Some(right_obj)) = (left_obj, right_obj) else {
        if left != right {
            diffs.push(HeaderDiff {
                key: "(root)".to_string(),
                left: left.to_string(),
                right: right.to_string(),
            });
        }
        return diffs;
    };

    let mut all_keys: Vec<&String> = left_obj.keys().chain(right_obj.keys()).collect();
    all_keys.sort();
    all_keys.dedup();

    for key in all_keys {
        if ignore_keys.contains(&key.as_str()) {
            continue;
        }
        let lv = left_obj.get(key);
        let rv = right_obj.get(key);
        match (lv, rv) {
            (Some(l), Some(r)) if l != r => {
                diffs.push(HeaderDiff {
                    key: key.clone(),
                    left: l.to_string(),
                    right: r.to_string(),
                });
            }
            (None, Some(r)) => {
                diffs.push(HeaderDiff {
                    key: key.clone(),
                    left: "(missing)".to_string(),
                    right: r.to_string(),
                });
            }
            (Some(l), None) => {
                diffs.push(HeaderDiff {
                    key: key.clone(),
                    left: l.to_string(),
                    right: "(missing)".to_string(),
                });
            }
            _ => {}
        }
    }

    diffs
}
```

- [ ] **Step 4: Create `src/bin/compare_matrix/values.rs`**

```rust
pub struct ValueDiff {
    pub row: usize,
    pub col: usize,
    pub left: f64,
    pub right: f64,
    pub abs_diff: f64,
}

pub struct ComparisonResult {
    pub total_cells: usize,
    pub differing_cells: usize,
    pub max_abs_diff: f64,
    pub worst_diffs: Vec<ValueDiff>, // top N worst
    pub row_count_match: bool,
    pub col_count_match: bool,
    pub left_rows: usize,
    pub right_rows: usize,
}

pub fn compare_values(
    left: &[Vec<f64>],
    right: &[Vec<f64>],
    tolerance: f64,
    max_diffs_reported: usize,
) -> ComparisonResult {
    let row_count_match = left.len() == right.len();
    let col_count_match = left
        .first()
        .zip(right.first())
        .map_or(true, |(l, r)| l.len() == r.len());

    let compare_rows = left.len().min(right.len());
    let mut total_cells = 0usize;
    let mut differing_cells = 0usize;
    let mut max_abs_diff = 0.0f64;
    let mut worst_diffs: Vec<ValueDiff> = Vec::new();

    for row_idx in 0..compare_rows {
        let lrow = &left[row_idx];
        let rrow = &right[row_idx];
        let cols = lrow.len().min(rrow.len());
        for col_idx in 0..cols {
            total_cells += 1;
            let lv = lrow[col_idx];
            let rv = rrow[col_idx];

            if lv.is_nan() && rv.is_nan() {
                continue;
            }
            if lv.is_nan() != rv.is_nan() {
                differing_cells += 1;
                let diff = f64::INFINITY;
                if worst_diffs.len() < max_diffs_reported {
                    worst_diffs.push(ValueDiff {
                        row: row_idx,
                        col: col_idx,
                        left: lv,
                        right: rv,
                        abs_diff: diff,
                    });
                }
                max_abs_diff = f64::INFINITY;
                continue;
            }

            let abs_diff = (lv - rv).abs();
            if abs_diff > tolerance {
                differing_cells += 1;
                if abs_diff > max_abs_diff {
                    max_abs_diff = abs_diff;
                }
                if worst_diffs.len() < max_diffs_reported {
                    worst_diffs.push(ValueDiff {
                        row: row_idx,
                        col: col_idx,
                        left: lv,
                        right: rv,
                        abs_diff,
                    });
                }
            }
        }
    }

    ComparisonResult {
        total_cells,
        differing_cells,
        max_abs_diff,
        worst_diffs,
        row_count_match,
        col_count_match,
        left_rows: left.len(),
        right_rows: right.len(),
    }
}
```

- [ ] **Step 5: Create `src/bin/compare_matrix/diff.rs`**

```rust
use super::parse::Matrix;
use super::header::compare_headers;
use super::values::compare_values;

use anyhow::Result;

pub fn full_diff(
    left: &Matrix,
    right: &Matrix,
    tolerance: f64,
    ignore_header_keys: &[&str],
) -> Result<String> {
    let mut output = String::new();

    // Header diff
    let header_diffs = compare_headers(&left.header, &right.header, ignore_header_keys);
    if header_diffs.is_empty() {
        output.push_str("Header: MATCH\n");
    } else {
        output.push_str(&format!("Header: {} difference(s)\n", header_diffs.len()));
        for d in &header_diffs {
            output.push_str(&format!("  {}: {} vs {}\n", d.key, d.left, d.right));
        }
    }

    // Value diff
    let left_values: Vec<Vec<f64>> = left.rows.iter().map(|r| r.values.clone()).collect();
    let right_values: Vec<Vec<f64>> = right.rows.iter().map(|r| r.values.clone()).collect();
    let result = compare_values(&left_values, &right_values, tolerance, 20);

    if !result.row_count_match {
        output.push_str(&format!(
            "Row count: MISMATCH (left={}, right={})\n",
            result.left_rows, result.right_rows
        ));
    }
    if !result.col_count_match {
        output.push_str("Column count: MISMATCH\n");
    }

    output.push_str(&format!(
        "Values: {}/{} cells differ (tolerance={:.0e}), max abs diff={:.6e}\n",
        result.differing_cells, result.total_cells, tolerance, result.max_abs_diff
    ));

    if !result.worst_diffs.is_empty() {
        output.push_str("Worst differences:\n");
        for d in &result.worst_diffs {
            output.push_str(&format!(
                "  row={} col={}: {:.6} vs {:.6} (diff={:.6e})\n",
                d.row, d.col, d.left, d.right, d.abs_diff
            ));
        }
    }

    Ok(output)
}
```

- [ ] **Step 6: Create `src/bin/compare_matrix.rs`**

```rust
mod compare_matrix {
    pub mod parse;
    pub mod header;
    pub mod values;
    pub mod diff;
}

use std::path::PathBuf;
use std::process::ExitCode;

use clap::{Parser, Subcommand};

use compare_matrix::{parse, header, values, diff};

#[derive(Parser)]
#[command(name = "compare_matrix", about = "Compare two computeMatrix .mat.gz files")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Compare only the JSON headers
    Header {
        left: PathBuf,
        right: PathBuf,
        /// Header keys to ignore (e.g. "proc number")
        #[arg(long, value_delimiter = ',')]
        ignore: Vec<String>,
    },
    /// Compare only numerical values
    Values {
        left: PathBuf,
        right: PathBuf,
        /// Tolerance for absolute difference
        #[arg(long, default_value = "5e-6")]
        tolerance: f64,
    },
    /// Full comparison (header + values + detailed diff)
    Diff {
        left: PathBuf,
        right: PathBuf,
        /// Tolerance for absolute difference
        #[arg(long, default_value = "5e-6")]
        tolerance: f64,
        /// Header keys to ignore
        #[arg(long, value_delimiter = ',')]
        ignore: Vec<String>,
    },
}

fn main() -> ExitCode {
    let cli = Cli::parse();

    let result = match cli.command {
        Commands::Header { left, right, ignore } => run_header(&left, &right, &ignore),
        Commands::Values { left, right, tolerance } => run_values(&left, &right, tolerance),
        Commands::Diff { left, right, tolerance, ignore } => {
            run_diff(&left, &right, tolerance, &ignore)
        }
    };

    match result {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => ExitCode::from(1),
        Err(e) => {
            eprintln!("Error: {e:#}");
            ExitCode::from(2)
        }
    }
}

fn run_header(left: &PathBuf, right: &PathBuf, ignore: &[String]) -> anyhow::Result<bool> {
    let lm = parse::load_matrix(left)?;
    let rm = parse::load_matrix(right)?;
    let ignore_refs: Vec<&str> = ignore.iter().map(|s| s.as_str()).collect();
    let diffs = header::compare_headers(&lm.header, &rm.header, &ignore_refs);
    if diffs.is_empty() {
        println!("Headers match.");
        Ok(true)
    } else {
        println!("{} header difference(s):", diffs.len());
        for d in &diffs {
            println!("  {}: {} vs {}", d.key, d.left, d.right);
        }
        Ok(false)
    }
}

fn run_values(left: &PathBuf, right: &PathBuf, tolerance: f64) -> anyhow::Result<bool> {
    let lm = parse::load_matrix(left)?;
    let rm = parse::load_matrix(right)?;
    let lv: Vec<Vec<f64>> = lm.rows.iter().map(|r| r.values.clone()).collect();
    let rv: Vec<Vec<f64>> = rm.rows.iter().map(|r| r.values.clone()).collect();
    let result = values::compare_values(&lv, &rv, tolerance, 20);
    println!(
        "{}/{} cells differ (tolerance={:.0e}), max abs diff={:.6e}",
        result.differing_cells, result.total_cells, tolerance, result.max_abs_diff
    );
    if !result.worst_diffs.is_empty() {
        for d in &result.worst_diffs {
            println!(
                "  row={} col={}: {:.6} vs {:.6} (diff={:.6e})",
                d.row, d.col, d.left, d.right, d.abs_diff
            );
        }
    }
    Ok(result.differing_cells == 0 && result.row_count_match)
}

fn run_diff(
    left: &PathBuf,
    right: &PathBuf,
    tolerance: f64,
    ignore: &[String],
) -> anyhow::Result<bool> {
    let lm = parse::load_matrix(left)?;
    let rm = parse::load_matrix(right)?;
    let ignore_refs: Vec<&str> = ignore.iter().map(|s| s.as_str()).collect();
    let report = diff::full_diff(&lm, &rm, tolerance, &ignore_refs)?;
    print!("{report}");

    let header_ok = header::compare_headers(&lm.header, &rm.header, &ignore_refs).is_empty();
    let lv: Vec<Vec<f64>> = lm.rows.iter().map(|r| r.values.clone()).collect();
    let rv: Vec<Vec<f64>> = rm.rows.iter().map(|r| r.values.clone()).collect();
    let values_ok = values::compare_values(&lv, &rv, tolerance, 0).differing_cells == 0;

    Ok(header_ok && values_ok)
}
```

Note: The module structure uses `src/bin/compare_matrix.rs` as the entry point and inline module declarations. Alternatively, create `src/bin/compare_matrix/main.rs` with sibling modules — use whichever Cargo supports for the binary target. The `[[bin]]` path must point to the file containing `fn main()`.

- [ ] **Step 7: Verify it compiles and runs**

Run: `cargo build --release --bin compare_matrix && target/release/compare_matrix --help`
Expected: Shows subcommands `header`, `values`, `diff`

- [ ] **Step 8: Commit**

```bash
git add Cargo.toml src/bin/
git commit -m "feat: add compare_matrix dev binary with header/values/diff subcommands"
```

---

### Task 4: Rust Integration Tests

**Files:**
- Create: `tests/python_compatibility.rs`

- [ ] **Step 1: Create the integration test file**

```rust
use std::path::{Path, PathBuf};
use std::process::Command;

fn project_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn data_root() -> PathBuf {
    project_root().join("deeptools/deeptools/test/test_heatmapper")
}

fn test_data_root() -> PathBuf {
    project_root().join("deeptools/deeptools/test/test_data")
}

fn compute_matrix_bin() -> PathBuf {
    project_root().join("target/release/compute_matrix_rs")
}

fn compare_matrix_bin() -> PathBuf {
    project_root().join("target/release/compare_matrix")
}

fn run_compute_and_compare(
    args: &[&str],
    reference_matrix: &str,
    tolerance: f64,
) {
    let bin = compute_matrix_bin();
    assert!(
        bin.exists(),
        "Binary not found at {:?}. Run `cargo build --release` first.",
        bin
    );

    let temp_dir = tempfile::tempdir().unwrap();
    let output_path = temp_dir.path().join("output.mat.gz");

    let mut cmd = Command::new(&bin);
    for arg in args {
        let expanded = arg
            .replace("{data_root}", data_root().to_str().unwrap())
            .replace("{test_data_root}", test_data_root().to_str().unwrap());
        cmd.arg(expanded);
    }
    cmd.arg("-o").arg(&output_path);

    let status = cmd.status().expect("Failed to execute compute_matrix_rs");
    assert!(status.success(), "compute_matrix_rs exited with {status}");
    assert!(output_path.exists(), "Output matrix not created");

    let ref_path = data_root().join(reference_matrix);
    assert!(
        ref_path.exists(),
        "Reference matrix not found at {:?}",
        ref_path
    );

    let cmp_bin = compare_matrix_bin();
    let cmp_status = Command::new(&cmp_bin)
        .arg("values")
        .arg(&ref_path)
        .arg(&output_path)
        .arg("--tolerance")
        .arg(format!("{tolerance:.0e}"))
        .status()
        .expect("Failed to execute compare_matrix");

    assert!(
        cmp_status.success(),
        "Matrix comparison failed for {reference_matrix} (tolerance={tolerance:.0e})"
    );
}

#[test]
fn reference_point_basic() {
    run_compute_and_compare(
        &[
            "reference-point",
            "-R", "{data_root}/test2.bed",
            "-S", "{data_root}/test.bw",
            "-b", "100", "-a", "100",
            "--bs", "1", "-p", "1",
        ],
        "master.mat",
        5e-6,
    );
}

#[test]
fn reference_point_center() {
    run_compute_and_compare(
        &[
            "reference-point",
            "-R", "{data_root}/test2.bed",
            "-S", "{data_root}/test.bw",
            "-b", "100", "-a", "100",
            "--referencePoint", "center",
            "--bs", "1", "-p", "1",
        ],
        "master_center.mat",
        5e-6,
    );
}

#[test]
fn reference_point_tes() {
    run_compute_and_compare(
        &[
            "reference-point",
            "-R", "{data_root}/test2.bed",
            "-S", "{data_root}/test.bw",
            "-b", "100", "-a", "100",
            "--referencePoint", "TES",
            "--bs", "1", "-p", "1",
        ],
        "master_TES.mat",
        5e-6,
    );
}

#[test]
fn reference_point_missing_data_as_zero() {
    run_compute_and_compare(
        &[
            "reference-point",
            "-R", "{data_root}/test2.bed",
            "-S", "{data_root}/test.bw",
            "-b", "100", "-a", "100",
            "--bs", "1", "-p", "1",
            "--missingDataAsZero",
        ],
        "master_nan_to_zero.mat",
        5e-6,
    );
}

#[test]
fn scale_regions_basic() {
    run_compute_and_compare(
        &[
            "scale-regions",
            "-R", "{data_root}/test2.bed",
            "-S", "{data_root}/test.bw",
            "-b", "100", "-a", "100",
            "-m", "100",
            "--bs", "1", "-p", "1",
        ],
        "master_scale_reg.mat",
        5e-6,
    );
}

#[test]
fn multiple_bed() {
    run_compute_and_compare(
        &[
            "reference-point",
            "-R", "{data_root}/test2.bed", "{data_root}/test2.bed",
            "-S", "{data_root}/test.bw",
            "-b", "100", "-a", "100",
            "--bs", "1", "-p", "1",
        ],
        "master_multi.mat",
        5e-6,
    );
}

#[test]
fn region_extend_beyond_chr() {
    run_compute_and_compare(
        &[
            "reference-point",
            "-R", "{data_root}/test2.bed",
            "-S", "{data_root}/test.bw",
            "-b", "100", "-a", "500",
            "--bs", "1", "-p", "1",
        ],
        "master_extend_beyond.mat",
        5e-6,
    );
}

#[test]
fn scale_regions_unscaled() {
    run_compute_and_compare(
        &[
            "scale-regions",
            "-R", "{data_root}/test2.bed",
            "-S", "{data_root}/test.bw",
            "-b", "100", "-a", "100",
            "-m", "100",
            "--unscaled5prime", "20",
            "--unscaled3prime", "20",
            "--bs", "1", "-p", "1",
        ],
        "master_scale_reg_unscaled.mat",
        5e-6,
    );
}

#[test]
fn gtf_input() {
    run_compute_and_compare(
        &[
            "scale-regions",
            "-R", "{test_data_root}/test.gtf",
            "-S", "{data_root}/test.bw",
            "-b", "100", "-a", "100",
            "-m", "100",
            "--bs", "10", "-p", "1",
        ],
        "master_gtf.mat",
        5e-6,
    );
}

#[test]
fn metagene() {
    run_compute_and_compare(
        &[
            "scale-regions",
            "-R", "{test_data_root}/test.gtf",
            "-S", "{data_root}/test.bw",
            "-b", "100", "-a", "100",
            "-m", "100",
            "--bs", "10", "-p", "1",
            "--metagene",
        ],
        "master_metagene.mat",
        5e-6,
    );
}
```

- [ ] **Step 2: Build release binaries and run integration tests**

Run: `cargo build --release && cargo test --test python_compatibility -- --test-threads=1`
Expected: 10/10 PASS

- [ ] **Step 3: Commit**

```bash
git add tests/python_compatibility.rs
git commit -m "test: add Rust integration tests for Python compatibility (10 scenarios)"
```

---

### Task 5: Improve `profile_bench.sh` with Warm Cache

**Files:**
- Modify: `scripts/profile_bench.sh`

- [ ] **Step 1: Add warm cache run before profiling**

Insert after the `mkdir -p bench_reports` line (line 28) and before "Run 1/4":

```bash
echo "=== Warm-up run (populating page cache) ===" >&2
"$@" > /dev/null 2>&1 || true
echo "  warm-up complete" >&2
echo "" >&2
```

This runs the command once silently, discarding output. The `|| true` ensures the script continues even if the warm-up run fails (e.g., output file already exists).

- [ ] **Step 2: Verify the script still works**

Run: `bash scripts/profile_bench.sh test-warmup "test" "test" -- /usr/bin/echo hello`
Expected: Report generated in `bench_reports/`, warm-up line visible in stderr

- [ ] **Step 3: Commit**

```bash
git add scripts/profile_bench.sh
git commit -m "feat: add warm-cache run to profile_bench.sh before profiling"
```

---

## Execution Order

Tasks 1 and 5 are independent and trivial. Tasks 2a–2f must be sequential. Task 3 is independent of Task 2. Task 4 depends on Task 3.

Recommended order: **1 → 5 → 3 → 2a → 2b → 2c → 2d → 2e → 2f → 4**

Tasks 1, 5, and 3 can be parallelized if using subagent-driven development.
