# Next Stages Implementation Plan (v2)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Five improvements: remove `num_cpus` dep, add file-spilling for large sorted matrices, build a matrix comparison dev binary, migrate Python regression tests to Rust integration tests, and improve the profiling bench script.

**Architecture:** `run.rs` always creates a `FileCollector` (with optional auxiliary writers) and passes it to `executor`. The executor decides the output strategy internally:
1. **StreamOrdered**: sort=no (with per-group I/O sort), or sort=keep + already compute-sorted → rows piped directly to `FileCollector`.
2. **HybridBucket**: any reordering needed (keep+coalesced / ascend / descend), any matrix size → `HybridBucketCollector` accumulates rows with chunked ~1 GB flush to temp files via writer threads. Small matrices stay entirely in memory (never trigger spill). At finalize, sorted rows are emitted through the same `FileCollector`.

`InMemoryCollector`, `GroupBucketCollector`, `MatrixData` (the struct bundling header + rows), `write_outputs()`, `RunOutcome::Matrix`, and `spawn_writer_thread` become dead code and are removed. `MatrixHeader` and `MatrixRow` are retained (used by all paths). All output is written directly through `FileCollector`.

Each HybridBucket accumulates rows in memory up to ~1 GB (default, injectable for testing), then flushes the chunk to a temp file via `std::thread::spawn` + move. A double-buffer scheme reuses the flushed Vec's capacity; back-pressure (join oldest flush handle) prevents unbounded memory growth on slow I/O (HDD). At finalize, temp files are mmap'd (`memmap2`) for zero-copy sorted readback.

Auxiliary outputs (`--outFileNameMatrix`, `--outFileSortedRegions`) are handled by `FileCollector`:
- **`outFileSortedRegions`**: Static BED header written first, then row-by-row. Group label is **passed explicitly** alongside each row — the emit interface is `on_row(group_index, row)`. Although all paths now produce group-contiguous output (sort=No uses per-group I/O sort, keep/ascend/descend restore group order), explicit group_index is more robust than deriving from boundaries and simplifies the code.
- **`outFileNameMatrix`**: Three-line header format. Only line 1 (`#Group1:N\tGroup2:M`) depends on final group counts. Lines 2-3 are parameter-derived and fixed. Placeholder trick:
  1. At creation: compute line 1 using `group_capacity` counts, record its byte length as `reserved_line1_len`. Write line 1 + lines 2-3. **Flush the BufWriter** before writing any data rows.
  2. Data rows are appended after the header.
  3. At `finish()`: **flush the BufWriter first**. Compute final line 1 with actual group counts. Assert `final_len <= reserved_line1_len` (initial counts ≥ final counts since filtering only removes rows). Pad final line 1 with trailing spaces before `\n` to exactly `reserved_line1_len` bytes. **Seek to offset 0**, write the padded line 1. If `final_len > reserved_line1_len`, bail with error (should never happen, but safety check prevents overwriting line 2).
  4. Lines 2-3 are never rewritten.
Both StreamOrdered and HybridBucket paths emit through the same `FileCollector`. `on_row(group_index, row)` writes to all active writers (main gzip + optional auxiliary).

**Tech Stack:** Rust (edition 2024), `memmap2` for mmap, `clap` for the comparison binary CLI, existing `flate2`/`serde_json` for matrix I/O.

---

## Key Design Decisions & Constraints

### No MatrixData — rows flow through FileCollector
Both paths emit `MatrixRow` through `FileCollector`. The `MatrixData` struct (which bundled header + all rows in memory) is deleted. `MatrixHeader` is built at finalize time from group counts; `MatrixRow` remains the per-row type.

### orig_idx is task index, not row index
`orig_idx` is assigned as `task.index` in `run.rs` (0..task_count). Filtered rows (skipZeros, threshold) produce `None` from the worker — no MatrixRow is emitted for them. For `finalize_keep_order`, the placement array must be sized by `task_count` (not `group_counts.sum()`), with empty slots for filtered rows.

### Sort tie-break uses insertion sequence, not orig_idx
Current `sort_groups()` does stable `sort_by(ascending)` on rows in their GroupBucket push order (which is compute order — I/O locality sorted, NOT orig_idx order), then `indices.reverse()` for descend. Equal-key rows preserve compute order in ascending, and reverse compute order in descending.

The HybridBucketCollector must replicate this exactly. Each row gets an `insertion_seq: u32` (monotonically increasing counter per bucket) at push time. Sort comparator: `(sort_key ascending, insertion_seq ascending)`, then `.reverse()` the entire sequence for descend. `orig_idx` is only used for `sort=keep` (placement array), not for sort tie-breaking.

### sortUsingSamples normalization
The HybridBucket path computes sort_key via `compute_sort_metric`. The 1-based `--sortUsingSamples` indices must be normalized to 0-based via `normalize_sort_sample_indices()` in the executor before computing sort keys.

### group_index travels with each row
The emit interface is `on_row(group_index, row)` — not derived from boundaries. All paths now produce group-contiguous output (sort=No uses per-group I/O sort; keep/ascend/descend restore group order via HybridBucket). Explicit group_index is still used because it's more robust and simpler than deriving from boundaries. The group_index is available from the executor's batch result tuple `(orig_idx, group_index, Option<MatrixRow>)` throughout the pipeline.

### sort=No must preserve group-contiguous output
The current executor globally sorts all work_items by chrom/query for I/O locality, even for sort=No. This interleaves groups and makes `group_boundaries` inaccurate — breaking plotHeatmap/plotProfile.

**Fix**: For sort=No, sort work_items **within each group** by I/O locality, but keep groups in their original order. This guarantees group-contiguous output while preserving intra-group I/O locality.

```rust
// In executor, Phase 2:
if matches!(general.sort_regions, SortRegions::No) && !already_sorted {
    // Sort within each group, keeping groups contiguous
    let mut start = 0;
    for &count in &group_item_counts {
        let end = start + count;
        work_items[start..end].sort_by(|a, b| {
            a.record.chrom.cmp(&b.record.chrom)
                .then(a.query_start.cmp(&b.query_start))
                .then(a.query_end.cmp(&b.query_end))
        });
        start = end;
    }
} else if !already_sorted {
    // Global sort for keep/ascend/descend (cross-group I/O locality OK)
    work_items.sort_by(|a, b| { ... });
}
```

Cross-group coalescing is lost, but the performance impact is minimal (groups typically have different regions). Header correctness is non-negotiable.

### ChromTable per chunk file
Each flush produces its own `ChromTable`. Deserialization uses the matching chunk's table. Tables are tiny (~25 entries) and stay in memory alongside the mmap.

### Overflow protection
`serialize_row` uses `u16` for string lengths and field counts (name, extra_fields, exon_coords), `u32` for `sample_count`, `bin_count`, and `row_byte_len`. `bin_count` as u16 would reject legitimate large matrices (e.g., long window + `--bs 1` easily exceeds 65535 bins), so `sample_count` and `bin_count` use `u32`. String lengths and collection counts stay `u16` (65535 is a reasonable limit for a gene name or extra fields). Add `u16::try_from().context()` / `u32::try_from().context()` at serialization time.

### Spilling threshold is injectable for testing
`DEFAULT_MEMORY_SPILL_THRESHOLD` is a `const` (1 GB). `HybridBucketCollector::new()` uses the default. `HybridBucketCollector::with_threshold(threshold)` accepts an override for testing (e.g., 100 bytes to trigger spilling in unit tests). The threshold is stored as an instance field, not a global.

---

## File Map

### Task 1 (num_cpus removal)
- Modify: `src/config.rs`
- Modify: `Cargo.toml`

### Task 2 (file spilling) — sub-tasks 2a through 2e
- Create: `src/pipeline/core/spill.rs` (SpillIndex, serialize/deserialize, ChromTable, HybridBucketCollector, FlushResult)
- Modify: `src/pipeline/core/mod.rs` (re-export)
- Modify: `src/pipeline/core/executor.rs` (OutputStrategy: StreamOrdered + HybridBucket; executor decides internally)
- Modify: `src/pipeline/core/collector.rs` (remove `InMemoryCollector` + `GroupBucketCollector`; extend `FileCollector` with AuxValuesWriter + regions writer; change `RowCollector::on_row` to include `group_index`)
- Modify: `src/pipeline/matrix.rs` (`compute_sort_metric`, `compare_ascending` → `pub(crate)`; remove `MatrixData` struct, `sort_groups()`, `prune_zero_rows()`, `GroupStats`; keep `MatrixHeader`, `MatrixRow`, `MatrixHeaderBuilder`, `LayoutVectors`)
- Modify: `src/pipeline/run.rs` (always create FileCollector, always call execute_mode, remove InMemory path)
- Modify: `src/pipeline/mod.rs` (remove `RunOutcome::Matrix`, `spawn_writer_thread`)
- Modify: `src/io/writers/mod.rs` (remove `write_outputs()`, `should_use_streaming()`, auxiliary guard; re-export `STREAMING_CELL_THRESHOLD`)
- Modify: `src/io/writers/auxiliary.rs` (add per-row streaming functions)
- Modify: `Cargo.toml` (add `memmap2`)

### Task 3 (matrix comparison binary)
- Create: `src/bin/compare_matrix/main.rs` (CLI entry with subcommands)
- Create: `src/bin/compare_matrix/parse.rs` (supports plain text + gzip + multi-member gzip via `MultiGzDecoder` with fallback)
- Create: `src/bin/compare_matrix/header.rs`
- Create: `src/bin/compare_matrix/values.rs`
- Create: `src/bin/compare_matrix/diff.rs`
- Modify: `Cargo.toml` (`[[bin]]` target)

### Task 4 (Rust integration tests)
- Create: `tests/python_compatibility.rs` (uses `env!("CARGO_BIN_EXE_compute_matrix_rs")` and `env!("CARGO_BIN_EXE_compare_matrix")`)

### Task 5 (profile_bench.sh improvement)
- Modify: `scripts/profile_bench.sh`

---

### Task 1: Remove `num_cpus` Dependency

**Files:**
- Modify: `src/config.rs:237-241`
- Modify: `Cargo.toml:10`

- [ ] **Step 1: Write a test for processor resolution**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resolve_max_returns_at_least_one() {
        assert!(ProcessorRequest::Max.resolve() >= 1);
    }

    #[test]
    fn resolve_max_half_returns_at_least_one() {
        assert!(ProcessorRequest::MaxHalf.resolve() >= 1);
    }

    #[test]
    fn resolve_fixed_zero_clamps_to_one() {
        assert_eq!(ProcessorRequest::Fixed(0).resolve(), 1);
    }

    #[test]
    fn resolve_fixed_value() {
        assert_eq!(ProcessorRequest::Fixed(4).resolve(), 4);
    }
}
```

- [ ] **Step 2: Run test to verify it passes**

Run: `cargo test --lib -- config::tests`
Expected: PASS (4 tests)

- [ ] **Step 3: Replace `num_cpus` with `std::thread::available_parallelism`**

In `src/config.rs`, replace `available_cpus()`:

```rust
fn available_cpus() -> u32 {
    let count = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    u32::try_from(count).unwrap_or(u32::MAX)
}
```

Remove `num_cpus = "1.16"` from `Cargo.toml`.

- [ ] **Step 4: Run tests and build**

Run: `cargo test --lib -- config::tests && cargo build --release`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/config.rs Cargo.toml && git add Cargo.lock 2>/dev/null; true
git commit -m "refactor: remove num_cpus dep, use std::thread::available_parallelism"
```

---

### Task 2: File Spilling for Large Sorted Matrices

Split into sub-tasks 2a–2e.

#### Task 2a: Spill Module — Types, Serialization, ChromTable

**Files:**
- Modify: `Cargo.toml` (add `memmap2 = "0.9"`)
- Create: `src/pipeline/core/spill.rs`
- Modify: `src/pipeline/core/mod.rs`

- [ ] **Step 1: Create `spill.rs` with core types and serialization**

Key types:

```rust
const SPILL_BUF_CAPACITY: usize = 1_048_576; // 1 MB
const DEFAULT_MEMORY_SPILL_THRESHOLD: usize = 1024 * 1024 * 1024; // 1 GB per bucket

pub(crate) struct SpillIndex {
    pub(crate) orig_idx: usize,
    pub(crate) group_index: usize,
    pub(crate) sort_key: f64,
    pub(crate) insertion_seq: u32,  // per-bucket push order, for stable sort tie-break
    pub(crate) file_offset: u64,
    pub(crate) row_byte_len: u32,
}
```

On-disk format per row (no outer length prefix — `row_byte_len` is in SpillIndex):

```text
[2 bytes: chrom_id (u16 LE)]
[4 bytes: start (u32 LE)]
[4 bytes: end (u32 LE)]
[1 byte:  strand (0=Positive, 1=Negative, 2=Unstranded)]
[1 byte:  flags (bit 0: has_name, bit 1: has_score_raw, bit 2: has_strand_raw, bit 3: has_exon_coords)]
[variable: name (u16 len + UTF-8, if bit 0)]
[variable: score (if bit 1: u16 len + UTF-8 raw; else: 4 bytes f32 LE, or 0xFFFFFFFF if None)]
[variable: strand_raw (u16 len + UTF-8, if bit 2)]
[variable: extra_fields (u16 count, then each: u16 len + UTF-8)]
[variable: exon_coords (u16 count, then each: u32 start + u32 end, if bit 3)]
[4+4 bytes: sample_count, bin_count (u32 LE)]
[N × 8 bytes: values (f64 LE)]
```

`serialize_row`: clears buf, writes fields, returns buf.len() as row_byte_len. Uses `u16::try_from(len).context("field too large for spill format")?` for overflow protection (return `Result` instead of silent truncation).

`deserialize_row`: reads from `&[u8]` slice + `ChromTable`.

`ChromTable`: `HashMap<Arc<str>, u16>` + `Vec<Arc<str>>` for bidirectional lookup.

- [ ] **Step 2: Write round-trip serialization tests**

Three tests: full round-trip, minimal row, ChromTable interning. Deserialize directly from `&buf` (no LP to skip).

- [ ] **Step 3: Add module to `mod.rs`, verify compilation**

Run: `cargo check`

- [ ] **Step 4: Commit**

```bash
git commit -m "feat: add spill module with row serialization and chrom interning"
```

#### Task 2b: HybridBucketCollector — Chunked Flush with Double Buffer

**Files:**
- Modify: `src/pipeline/core/spill.rs`

- [ ] **Step 1: Add FlushResult, flush_chunk, CollectorBucket, HybridBucketCollector**

```rust
struct FlushResult {
    indices: Vec<SpillIndex>,
    returned_buf: Vec<(usize, f64, u32, MatrixRow)>,
    chrom_table: ChromTable,
    temp_path: PathBuf,
}
```

`flush_chunk(rows: Vec<(usize, f64, u32, MatrixRow)>, group_index) -> Result<FlushResult>`: runs on spawned thread, serializes all rows to a new temp file. Each row's `(orig_idx, sort_key, insertion_seq)` is recorded in the SpillIndex. Returns indices + cleared-but-capacity-retained Vec.

`CollectorBucket`:
```rust
struct CollectorBucket {
    active: Vec<(usize, f64, u32, MatrixRow)>,  // (orig_idx, sort_key, insertion_seq, row)
    spare: Option<Vec<(usize, f64, u32, MatrixRow)>>,
    estimated_bytes: usize,
    group_index: usize,
    flush_handles: Vec<JoinHandle<Result<FlushResult>>>,
    completed_flushes: Vec<FlushResult>,
}
```

`push()`: accumulates in `active`. When `estimated_bytes > threshold`:
1. Try `spare.take()` for the replacement buffer.
2. If spare unavailable + flush handles exist: **back-pressure** — join oldest handle, take its `returned_buf`, store FlushResult **without** `returned_buf` (set to empty Vec to release memory).
3. If no spare and no handles: `Vec::new()` (first flush).
4. `std::mem::replace(&mut self.active, spare)` → spawn writer thread with the full Vec.

`join_all()`: joins remaining handles. **Important**: after extracting `returned_buf` for spare reuse, clear the `returned_buf` field in the stored `FlushResult` to avoid holding dead capacity:
```rust
for handle in self.flush_handles.drain(..) {
    let mut result = handle.join().unwrap()?;
    result.returned_buf = Vec::new(); // release capacity — data is on disk
    self.completed_flushes.push(result);
}
```

`HybridBucketCollector`:
```rust
pub(crate) struct HybridBucketCollector { ... }

impl HybridBucketCollector {
    pub(crate) fn new(group_count, sample_count, bin_count) -> Self;
    pub(crate) fn with_threshold(group_count, sample_count, bin_count, threshold: usize) -> Self;
    pub(crate) fn push(row, orig_idx, group_index, sort_key) -> Result<()>;
}
```

- [ ] **Step 2: Add small-data in-memory test + low-threshold spilling test**

Use `with_threshold(100)` so even a few rows trigger spilling.

- [ ] **Step 3: Commit**

```bash
git commit -m "feat: add HybridBucketCollector with chunked flush and back-pressure"
```

#### Task 2c: Finalize — Sorted Emit with Merge + Keep-Order Emit

**Files:**
- Modify: `src/pipeline/core/spill.rs`
- Modify: `src/pipeline/matrix.rs` (`compare_ascending`, `compute_sort_metric` → `pub(crate)`)

- [ ] **Step 1: Add `finalize_sorted`**

Joins all handles. For each bucket: mmaps all chunk files, collects `SpillRef` from completed flushes, sorts in-memory `active` rows. Merge-emits spilled + in-memory using peekable iterators.

```rust
struct SpillRef {
    sort_key: f64,
    insertion_seq: u32,  // for stable sort tie-break (NOT orig_idx)
    mmap_idx: usize,
    file_offset: u64,
    row_byte_len: u32,
}
```

**Sort replicates current `sort_groups()` exactly**: always sort ascending by `(sort_key, insertion_seq)`, collect into a Vec, then `.reverse()` for descend:

```rust
// Merge all spill_refs + in-memory into one Vec<EmitEntry>
// (EmitEntry: SpillRef or owned MatrixRow, each carrying sort_key + insertion_seq)
// In-memory rows carry their insertion_seq from the active buffer tuple.
all_entries.sort_by(|a, b| {
    compare_ascending(a.sort_key, b.sort_key)
        .then(a.insertion_seq.cmp(&b.insertion_seq))
});
if !sort_ascending {
    all_entries.reverse();
}
for entry in all_entries { emit(entry)?; }
```

This means for descend, equal-key rows appear in reverse insertion order — matching the existing `stable sort_by(ascending)` + `reverse()` behavior exactly.

Note: **`normalize_sort_sample_indices`** must be called before computing sort keys. The executor receives `general.sort_using_samples` (1-based from CLI), normalizes to 0-based, and passes to `compute_sort_metric`.

Cleanup: drop mmaps, remove temp files.

- [ ] **Step 2: Add `finalize_keep_order`**

Joins all handles. Uses pre-allocated placement array sized by `task_count` (passed through from `run.rs`, NOT `group_counts.sum()`):

```rust
pub(crate) fn finalize_keep_order<F>(
    mut self,
    task_count: usize,
    header_builder: impl FnOnce(Vec<usize>) -> Result<MatrixHeader>,
    mut emit: F,
) -> Result<MatrixHeader>
```

`slots: Vec<Option<RowSlot>>` sized `task_count`. Slots at filtered-row indices remain `None` and are silently skipped during the O(n) scan. No sort needed.

- [ ] **Step 3: Tests**

Test finalize_sorted with low threshold, verify sort order. Test finalize_keep_order with gaps in orig_idx (simulating filtered rows).

- [ ] **Step 4: Commit**

```bash
git commit -m "feat: add finalize_sorted/keep_order with mmap readback and stable merge"
```

#### Task 2d: Integrate into executor.rs + run.rs (two paths only)

**Files:**
- Modify: `src/pipeline/core/executor.rs`
- Modify: `src/pipeline/run.rs`
- Modify: `src/pipeline/mod.rs` (remove `RunOutcome::Matrix`, `spawn_writer_thread`)
- Modify: `src/io/writers/mod.rs` (remove auxiliary guard, remove `write_outputs()`, re-export `STREAMING_CELL_THRESHOLD`)
- Modify: `src/io/writers/auxiliary.rs` (add per-row streaming functions)
- Modify: `src/pipeline/core/collector.rs` (remove `InMemoryCollector` + `GroupBucketCollector`; extend `FileCollector`; change `RowCollector::on_row` signature)
- Modify: `src/pipeline/matrix.rs` (remove `MatrixData`, `GroupStats`, `sort_groups()`, `prune_zero_rows()`)

- [ ] **Step 1: Simplify executor to two strategies**

```rust
enum OutputStrategy {
    StreamOrdered,   // sort=no (with per-group I/O sort), or sort=keep+already_sorted
    HybridBucket,    // everything else (keep+coalesced, ascend, descend) — any matrix size
}
```

For sort=No, **change I/O sort to per-group** (not global):
```rust
let output_strategy = match general.sort_regions {
    SortRegions::Keep if already_sorted => OutputStrategy::StreamOrdered,
    SortRegions::No => OutputStrategy::StreamOrdered,
    _ => OutputStrategy::HybridBucket,
};

// I/O locality sort
match general.sort_regions {
    SortRegions::No if !already_sorted => {
        // Per-group sort: keeps group-contiguous order, header stays correct
        let mut start = 0;
        for &count in &group_item_counts {
            let end = start + count;
            work_items[start..end].sort_by(|a, b| {
                a.record.chrom.cmp(&b.record.chrom)
                    .then(a.query_start.cmp(&b.query_start))
                    .then(a.query_end.cmp(&b.query_end))
            });
            start = end;
        }
    }
    _ if !already_sorted => {
        // Global sort (OK for keep/ascend/descend — HybridBucket restores order)
        work_items.sort_by(|a, b| { /* chrom, query_start, query_end */ });
    }
    _ => {}
}
```

HybridBucket branch:
1. Calls `normalize_sort_sample_indices(general.sort_using_samples, sample_count)` for sort key computation.
2. Creates `HybridBucketCollector::new(group_count, sample_count, total_bins)`.
3. Compute loop: for each `(orig_idx, group_index, Some(row))`, compute `sort_key` and `push(row, orig_idx, group_index, sort_key)`.
4. After compute: call `finalize_sorted` or `finalize_keep_order`, emitting through the `FileCollector` passed from `run.rs`.

StreamOrdered branch: same as current, but passes `group_index` to `collector.on_row(group_index, row)`.

- [ ] **Step 2: Change `RowCollector::on_row` to include `group_index`**

```rust
pub trait RowCollector: Send {
    type Output: Send;
    fn on_row(&mut self, group_index: usize, row: MatrixRow) -> Result<()>;
    fn finalize(self, header: MatrixHeader) -> Result<Self::Output>;
    fn abort(self) where Self: Sized {}
}
```

- [ ] **Step 3: Extend `FileCollector` with auxiliary writers**

```rust
pub struct FileCollector {
    writer: StreamingMatrixWriter,
    values_writer: Option<AuxValuesWriter>,
    regions_writer: Option<BufWriter<File>>,
    group_labels: Vec<String>,
}
```

`FileCollector::new(writer, group_labels, group_capacity, header_estimate, io)`:
- For outFileSortedRegions: opens file, writes static BED header.
- For outFileNameMatrix: creates `AuxValuesWriter` — writes placeholder line 1 (using group_capacity), then lines 2-3 (from header_estimate). **Flushes before returning** so rows start at known offset.

`on_row(group_index, row)`: writes main gzip row + optional plain values row + optional sorted region BED line (using `group_labels[group_index]`).

`finalize(header)`: for AuxValuesWriter, **flush**, seek to 0, overwrite line 1 with actual group counts (space-padded to match initial byte length), assert no overflow. Finalize main gzip.

- [ ] **Step 4: Simplify `run.rs` — always create FileCollector, pass to executor**

```rust
let header_estimate = mode.build_header(&general, metadata.as_ref(), &sample_labels,
    &group_labels, &group_capacity, thread_count, sample_count);
writers::ensure_streaming_header_capacity(&header_estimate)?;

let writer = writers::StreamingMatrixWriter::start(&io.matrix_output)?;
let collector = FileCollector::new(
    writer, &group_labels, &group_capacity, &header_estimate, io,
)?;

core::execute_mode(
    tasks, general, sample_paths, collector, thread_count,
    &mode, metadata, header_builder, group_labels.len(), task_count,
)?;

Ok(RunOutcome::Streamed)
```

`run.rs` no longer decides strategy — that's the executor's job. `RunOutcome::Matrix` is removed. `write_outputs()` is removed. `spawn_writer_thread` is removed. All output is written through `FileCollector`.

- [ ] **Step 5: Remove dead code**

- `InMemoryCollector` from `collector.rs`
- `GroupBucketCollector` from `collector.rs`
- `MatrixData`, `GroupStats`, `sort_groups()`, `prune_zero_rows()` from `matrix.rs`
- `write_outputs()`, `should_use_streaming()` from `io/writers/mod.rs`
- `RunOutcome::Matrix`, `spawn_writer_thread` from `pipeline/mod.rs`
- Auxiliary streaming guard from `should_use_streaming_for_plan`

- [ ] **Step 6: Add per-row streaming functions to `auxiliary.rs`**

```rust
pub fn write_sorted_region_row<W: Write>(writer: &mut W, row: &MatrixRow, group_label: &str) -> Result<()>
pub fn write_plain_values_row<W: Write>(writer: &mut W, row: &MatrixRow) -> Result<()>
```

- [ ] **Step 7: Add AuxValuesWriter unit tests**

Test cases:
- Group counts decrease after filtering → line 1 rewrite shorter, space-padded correctly
- 0 rows in a group → count is 0, line 1 still valid
- Long group labels → verify no overflow into line 2
- Multiple groups → tab-separated label:count format preserved
- **Re-parse rewritten header**: after finish(), read file back, parse line 1, verify group names and counts are correct (trailing spaces don't corrupt parsing)
- Header byte length exactly equals reserved → no padding needed, still works

- [ ] **Step 8: Run all tests**

Run: `cargo test`

- [ ] **Step 9: Commit**

```bash
git commit -m "feat: unify to two execution paths (StreamOrdered + HybridBucket), remove InMemory path"
```

#### Task 2e: Integration Tests

- [ ] **Step 1: Add spilling integration test with injectable threshold**

Use `with_threshold(100)` to trigger spilling on small test data. Verify:
- Sorted output matches expected order (ascend + descend).
- Keep-order output with gaps (filtered rows) matches expected order.
- Round-trip through serialize → temp file → mmap → deserialize → emit produces correct values.
- Tie-break: rows with equal sort key appear in insertion_seq order (ascending) or reverse (descending).

- [ ] **Step 2: Add sort=No group-contiguous regression test**

Multi-group test with regions on shared chromosomes (forcing I/O reordering). Verify:
- Output `.mat.gz` header `group_boundaries` correctly describes contiguous group segments.
- `outFileSortedRegions` group labels match actual group assignment for each row.

- [ ] **Step 3: Commit**

```bash
git commit -m "test: add integration tests for spilling with injectable threshold"
```

---

### Task 3: Matrix Comparison Dev Binary

**Files:**
- Create: `src/bin/compare_matrix/main.rs`
- Create: `src/bin/compare_matrix/parse.rs`
- Create: `src/bin/compare_matrix/header.rs`
- Create: `src/bin/compare_matrix/values.rs`
- Create: `src/bin/compare_matrix/diff.rs`
- Modify: `Cargo.toml`

- [ ] **Step 1: Add `[[bin]]` target**

```toml
[[bin]]
name = "compare_matrix"
path = "src/bin/compare_matrix/main.rs"
```

- [ ] **Step 2: Create `parse.rs` — supports plain text + gzip + multi-member gzip**

```rust
pub fn load_matrix(path: &Path) -> Result<Matrix> {
    let file = File::open(path)?;
    let mut buf = [0u8; 2];
    file.read_exact(&mut buf)?;
    file.seek(SeekFrom::Start(0))?;

    let reader: Box<dyn BufRead> = if buf == [0x1f, 0x8b] {
        // gzip (possibly multi-member)
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        // plain text
        Box::new(BufReader::new(file))
    };
    parse_matrix_from_reader(reader)
}
```

Uses `flate2::read::MultiGzDecoder` for multi-member gzip (produced by streaming writer's separate header + body gzip members).

- [ ] **Step 3: Create `header.rs`**

Same as before. Compare JSON objects key-by-key, skip `ignore_keys`.

- [ ] **Step 4: Create `values.rs`**

`compare_values(left, right, tolerance, max_diffs_reported) -> ComparisonResult`

Fix: check ALL rows for column count mismatch, not just first row. Include `col_count_match` in the success condition.

- [ ] **Step 5: Create `diff.rs`**

Full diff output: header diffs + value diffs + BED field comparison (detect row reordering by comparing chrom:start:end:name tuples).

- [ ] **Step 6: Create `main.rs`**

Subcommands: `header`, `values`, `diff`. Exit code 0 = match, 1 = mismatch, 2 = error.

- [ ] **Step 7: Build and test**

Run: `cargo build --release --bin compare_matrix && target/release/compare_matrix --help`

- [ ] **Step 8: Commit**

```bash
git commit -m "feat: add compare_matrix binary with plain/gzip/multi-member support"
```

---

### Task 4: Rust Integration Tests

**Files:**
- Create: `tests/python_compatibility.rs`

- [ ] **Step 1: Create integration test file**

Uses `env!("CARGO_BIN_EXE_compute_matrix_rs")` and `env!("CARGO_BIN_EXE_compare_matrix")` — no manual `cargo build --release` needed. Cargo builds the binaries automatically for integration tests.

Test parameters match `scripts/config/python_compatibility.yaml` exactly:

| Test | Reference | Key differences from v1 plan |
|---|---|---|
| reference_point_basic | master.mat | unchanged |
| reference_point_center | master_center.mat | unchanged |
| reference_point_tes | master_TES.mat | unchanged |
| reference_point_missing_data_as_zero | master_nan_to_zero.mat | unchanged |
| scale_regions_basic | master_scale_reg.mat | unchanged |
| multiple_bed | **master_multibed.mat** | uses **group1.bed + group2.bed** |
| region_extend_beyond_chr | **master_extend_beyond_chr_size.mat** | uses **group1.bed + group2.bed**, -a 500 |
| scale_regions_unscaled | **master_unscaled.mat** | uses **unscaled.bed + unscaled.bigWig**, -a 300 -b 500, --unscaled5prime 100 --unscaled3prime 50, --bs 10 |
| gtf_input | master_gtf.mat | uses **{test_data_root}/test1.bw.bw**, -a 300 -b 500, --unscaled5prime 20 --unscaled3prime 50, --bs 10 |
| metagene | master_metagene.mat | same as gtf_input + --metagene |

```rust
fn compute_matrix_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_compute_matrix_rs"))
}

fn compare_matrix_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_compare_matrix"))
}
```

`run_compute_and_compare` uses `compare_matrix diff` (not just `values`) to also catch header and BED field mismatches. The `--ignore "proc number"` flag is passed to ignore thread count differences.

- [ ] **Step 2: Run integration tests**

Run: `cargo test --test python_compatibility -- --test-threads=1`
Expected: 10/10 PASS

- [ ] **Step 3: Commit**

```bash
git commit -m "test: add Rust integration tests for Python compatibility (10 scenarios)"
```

---

### Task 5: Improve `profile_bench.sh` with Warm Cache

**Files:**
- Modify: `scripts/profile_bench.sh`

- [ ] **Step 1: Add warm cache run with isolated output**

Insert after `mkdir -p bench_reports` and before "Run 1/4":

```bash
echo "=== Warm-up run (populating page cache) ===" >&2
WARMUP_DIR=$(mktemp -d)
# Redirect all output-producing flags to temp dir
WARMUP_MODIFIED=()
SKIP_NEXT=false
for arg in "$@"; do
    if $SKIP_NEXT; then
        WARMUP_MODIFIED+=("$WARMUP_DIR/warmup_output")
        SKIP_NEXT=false
    elif [ "$arg" = "-o" ] || [ "$arg" = "--outFileName" ] \
      || [ "$arg" = "--outFileNameMatrix" ] \
      || [ "$arg" = "--outFileSortedRegions" ]; then
        WARMUP_MODIFIED+=("$arg")
        SKIP_NEXT=true
    else
        WARMUP_MODIFIED+=("$arg")
    fi
done
"${WARMUP_MODIFIED[@]}" > /dev/null 2>&1 || true
rm -rf "$WARMUP_DIR"
echo "  warm-up complete" >&2
echo "" >&2
```

Redirects `-o`, `--outFileName`, `--outFileNameMatrix`, and `--outFileSortedRegions` to temp dir. Handles both `-o path` and `--outFileName path` (space-separated) forms. Does not handle `--outFileName=path` (equals form) — acceptable since our CLI uses space-separated args.

- [ ] **Step 2: Test**

Run: `bash scripts/profile_bench.sh test-warmup "test" "test" -- /usr/bin/echo hello`

- [ ] **Step 3: Commit**

```bash
git commit -m "feat: add warm-cache run to profile_bench.sh with isolated output"
```

---

## Execution Order

```
Task 1 (trivial, independent)
Task 5 (trivial, independent)
Task 3 (independent of Task 2)
Task 2a → 2b → 2c → 2d → 2e (sequential)
Task 4 (depends on Task 3 + Task 2d for full coverage)
```

Tasks 1, 5, 3 can be parallelized. Task 2 subtasks are sequential. Task 4 should run last.

---

## Known Risks & Mitigations

| Risk | Mitigation |
|---|---|
| Many groups × 1 GB per bucket → OOM | Default 1 GB per bucket. Future: global memory budget. Current: acceptable for typical use (1-10 groups, 1-10 GB total before spilling). |
| Spill temp files leak on abort/panic | `HybridBucketCollector` implements `Drop` to join handles and remove temp files. `CollectorBucket::drop` removes `completed_flushes[*].temp_path`. |
| `compare_matrix` built with normal `[[bin]]` even in production | Don't feature-gate. Binary is tiny (shares all deps with main binary), adds no overhead. `cargo install` users can `--bin compute_matrix_rs` to skip it. Simpler than feature gates that complicate `cargo test`. |

---

## Appendix: Issues Addressed in v2

| # | Issue | Fix |
|---|---|---|
| 1 | Spilling rows re-collected into InMemoryCollector → OOM not solved | HybridBucket emits through `FileCollector`, bypasses `MatrixData` entirely |
| 2 | `finalize_keep_order` placement array sized by group_counts.sum(), filtered rows cause panic | Array sized by `task_count` (passed from run.rs), filtered slots are `None` |
| 3 | `join_all` retains FlushResult.returned_buf capacity | Clear `returned_buf` to `Vec::new()` after extracting for spare |
| 4 | Auxiliary output header needs group counts not yet known | `outFileSortedRegions`: static header + row-by-row with explicit group_index. `outFileNameMatrix`: placeholder line 1 (group_capacity counts) + seek-overwrite at `finish()` |
| 5 | `load_matrix` only reads gzip, references are plain text, streaming produces multi-member gzip | Magic-byte detection: `[0x1f, 0x8b]` → `MultiGzDecoder`, else plain `BufReader` |
| 6 | Task 4 test params/filenames wrong vs YAML | All tests now match `python_compatibility.yaml` exactly |
| 7 | No stable tie-break in spilling sort | `SpillRef` carries `insertion_seq`, sort comparator uses `.then(insertion_seq.cmp)` |
| 8 | u16/u32 overflow in serialization | `u16::try_from().context()` in `serialize_row` |
| 9 | Spilling test doesn't trigger spilling | `with_threshold(threshold)` constructor for injectable threshold |
| 10 | `compare_values` only checks first row columns, `col_count_match` not in success | Check all rows, include in success condition |
| 11 | Task 4 requires manual `cargo build --release` | `env!("CARGO_BIN_EXE_...")` auto-builds |
| 12 | `profile_bench.sh` warm-up may clobber output | Warm-up uses temp dir for output |
| 13 | StreamOrdered path missing FileCollector auxiliary writer impl | FileCollector extended with AuxValuesWriter + BufWriter for regions |
| 14 | outFileSortedRegions group label wrong if rows not group-contiguous | emit interface passes `group_index` explicitly — no boundary derivation needed. |
| 15 | sortUsingSamples not normalized in HybridBucket path | Executor calls `normalize_sort_sample_indices` before computing sort keys |
| 16 | Descend tie-break differs between old and new paths | HybridBucket replicates exact behavior: ascending sort, then `.reverse()` for descend |
| 17 | u16 bin_count rejects legitimate large matrices | bin_count/sample_count use u32; string lengths/counts stay u16 |
| 18 | outFileNameMatrix "seek prepend" impossible on plain file | Placeholder header + seek-overwrite (same as gzip path, simpler for plain text) |
| 19 | Per-bucket threshold may OOM with many groups | Default lowered to 1 GB; future: global memory budget |
| 20 | Spill temp file leak on abort | HybridBucketCollector Drop impl cleans up |
| 21 | compare_matrix built in production | Not gated — shares deps, no overhead, simplifies `cargo test` |
| 22 | group label derived from boundaries + row_index is wrong | emit interface passes `group_index` explicitly: `on_row(group_index, row)` |
| 23 | "emit order always group-contiguous" is false for sort=No | Confirmed: I/O coalescing interleaves groups. Group label must come from explicit group_index, not boundaries. |
| 24 | outFileNameMatrix placeholder header format unclear | Only line 1 has variable counts. Initial counts (group_capacity) ≥ final counts → space-pad final to match. Lines 2-3 are fixed. |
| 25 | write_outputs() early return skips auxiliary | `write_outputs()` removed entirely; all output through FileCollector |
| 26 | Feature gate + cargo test contradiction | Dropped feature gate; binary shares all deps, no overhead. |
| 27 | run.rs can't know if keep+already_sorted | Executor decides strategy internally; run.rs always passes FileCollector. StreamOrdered preserved for keep+already_sorted. |
| 28 | Small matrix InMemory path conflicts with GroupBucketCollector deletion | HybridBucketCollector serves all sizes; InMemoryCollector deleted. Small matrices stay in-memory (never trigger spill). |
| 29 | Spilling sort tie-break uses orig_idx, doesn't match compute order | Changed to `insertion_seq` (per-bucket push counter) — matches current stable sort behavior exactly. |
| 30 | AuxValuesWriter can't write lines 2-3 without header info | Receives `header_estimate: &MatrixHeader` (already computed in run.rs) for fixed header fields. |
| 31 | STREAMING_CELL_THRESHOLD not visible to run.rs | Re-exported from `io::writers::mod.rs`. |
| 32 | outFileNameMatrix seek/flush sequence unclear | Explicit: flush before rows, flush before seek, assert final_len <= reserved, pad to exact byte length. |
| 33 | sort=No group interleaving in .mat.gz | Fixed: sort=No uses per-group I/O sort (not global), preserving group-contiguous output. |
| 34 | Execution paths self-contradictory (two vs three) | run.rs always creates FileCollector; executor decides StreamOrdered vs HybridBucket internally. |
| 35 | Appendix has stale orig_idx / group-contiguous claims | Fixed: insertion_seq for tie-break, group_index passthrough for labels. |
| 36 | Default 4GB per bucket too high | Lowered to 1 GB. |
| 37 | Task 1 assumes Cargo.lock always changes | Changed to "if changed" |
| 38 | Task 2d step numbering duplicate | Fixed |
| 39 | AuxValuesWriter needs unit tests | Added: re-parse test, exact-length test |
| 40 | keep+already_sorted vs run.rs HybridBucket contradiction | Executor decides; run.rs is strategy-agnostic |
| 41 | Remaining InMemory references in plan text | All removed |
| 42 | outFileNameMatrix padding may corrupt header parsing | Added re-parse unit test to verify trailing spaces don't break consumers |
| 43 | sort=No group-contiguous needs regression test | Added to Task 2e: multi-group + shared chromosomes test |
| 44 | Warm-up doesn't handle --outFileNameMatrix/--outFileSortedRegions | All output flags redirected to temp dir |
| 45 | "Spilling" vs "HybridBucket" naming inconsistent | Unified to "HybridBucket" throughout |
| 46 | MEMORY_SPILL_THRESHOLD naming unclear | Renamed to DEFAULT_MEMORY_SPILL_THRESHOLD; instance field for override |
| 47 | Flags byte layout undefined | Explicit bit definitions in on-disk format doc |
| 48 | Task 2 says "2a-2f" but only 2a-2e exist | Fixed to "2a-2e" |
