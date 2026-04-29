# Next Stages Implementation Plan (v2)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Five improvements: remove `num_cpus` dep, add file-spilling for large sorted matrices, build a matrix comparison dev binary, migrate Python regression tests to Rust integration tests, and improve the profiling bench script.

**Architecture:** File spilling introduces a third execution path in `run.rs`:
1. **StreamOrdered** (existing): small/large matrix + no reordering → `FileCollector` streams directly.
2. **InMemory** (existing): small matrix + any sort → `InMemoryCollector` → `MatrixData` → `sort_groups()` → `write_outputs()`.
3. **Spilling** (new): large matrix + reordering needed (keep+coalesced / ascend / descend) → `HybridBucketCollector` accumulates rows with chunked 4 GB flush to temp files via writer threads. At finalize, sorted rows are emitted **directly to output writers** (main gzip + optional auxiliary), bypassing `MatrixData` entirely.

Each spilling bucket accumulates rows in memory up to ~4 GB, then flushes the chunk to a temp file via `std::thread::spawn` + move. A double-buffer scheme reuses the flushed Vec's capacity; back-pressure (join oldest flush handle) prevents unbounded memory growth on slow I/O (HDD). At finalize, temp files are mmap'd (`memmap2`) for zero-copy sorted readback.

Auxiliary outputs (`--outFileNameMatrix`, `--outFileSortedRegions`) are supported in both StreamOrdered and Spilling paths. `outFileSortedRegions` writes a static header then row-by-row (group label derived from group boundaries). `outFileNameMatrix` writes rows first, then seeks back to prepend the header (or uses a temp-file-then-copy approach, since it's uncompressed and small overhead).

**Tech Stack:** Rust (edition 2024), `memmap2` for mmap, `clap` for the comparison binary CLI, existing `flate2`/`serde_json` for matrix I/O.

---

## Key Design Decisions & Constraints

### Spilling path writes output directly — no MatrixData intermediate
The spilling path's `finalize_sorted`/`finalize_keep_order` emits rows directly to a `SpillingOutputWriter` that wraps the main gzip writer + optional auxiliary writers. `run.rs` returns `RunOutcome::Streamed` for this path. `MatrixData` is never constructed for large sorted matrices.

### orig_idx is task index, not row index
`orig_idx` is assigned as `task.index` in `run.rs` (0..task_count). Filtered rows (skipZeros, threshold) produce `None` from the worker — no MatrixRow is emitted for them. For `finalize_keep_order`, the placement array must be sized by `task_count` (not `group_counts.sum()`), with empty slots for filtered rows.

### Spilling sort must be stable
`MatrixData::sort_groups()` uses `sort_by` (stable). The spilling path's merge-emit of spilled chunks + in-memory rows must also be stable — equal sort keys preserve insertion order. `SpillRef` carries `orig_idx` for tie-breaking.

### ChromTable per chunk file
Each flush produces its own `ChromTable`. Deserialization uses the matching chunk's table. Tables are tiny (~25 entries) and stay in memory alongside the mmap.

### Overflow protection
`serialize_row` uses `u16` for string lengths and counts, `u32` for row_byte_len. Add `checked` conversions or assertions at serialization time; bail with a clear error if a field exceeds the limit.

### Spilling threshold is injectable for testing
`MEMORY_SPILL_THRESHOLD` is a `const` but tests can use a separate `HybridBucketCollector::with_threshold(threshold)` constructor to set a tiny threshold (e.g., 100 bytes) and actually trigger spilling in unit tests.

---

## File Map

### Task 1 (num_cpus removal)
- Modify: `src/config.rs`
- Modify: `Cargo.toml`

### Task 2 (file spilling)
- Create: `src/pipeline/core/spill.rs` (SpillIndex, serialize/deserialize, ChromTable, HybridBucketCollector, FlushResult, SpillingOutputWriter)
- Modify: `src/pipeline/core/mod.rs` (re-export)
- Modify: `src/pipeline/core/executor.rs` (new `Spilling` output strategy that writes directly to SpillingOutputWriter)
- Modify: `src/pipeline/matrix.rs` (`compute_sort_metric`, `compare_ascending` → `pub(crate)`)
- Modify: `src/pipeline/run.rs` (add third Spilling code path, pass `task_count` to executor)
- Modify: `src/io/writers/mod.rs` (remove auxiliary streaming guard)
- Modify: `src/io/writers/auxiliary.rs` (add per-row streaming functions with group_label param)
- Modify: `src/pipeline/core/collector.rs` (remove `GroupBucketCollector`; optionally extend `FileCollector` for auxiliary outputs in StreamOrdered path)
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
git add src/config.rs Cargo.toml Cargo.lock
git commit -m "refactor: remove num_cpus dep, use std::thread::available_parallelism"
```

---

### Task 2: File Spilling for Large Sorted Matrices

Split into sub-tasks 2a–2f.

#### Task 2a: Spill Module — Types, Serialization, ChromTable

**Files:**
- Modify: `Cargo.toml` (add `memmap2 = "0.9"`)
- Create: `src/pipeline/core/spill.rs`
- Modify: `src/pipeline/core/mod.rs`

- [ ] **Step 1: Create `spill.rs` with core types and serialization**

Key types:

```rust
const SPILL_BUF_CAPACITY: usize = 1_048_576; // 1 MB
const DEFAULT_MEMORY_SPILL_THRESHOLD: usize = 4 * 1024 * 1024 * 1024;

pub(crate) struct SpillIndex {
    pub(crate) orig_idx: usize,
    pub(crate) group_index: usize,
    pub(crate) sort_key: f64,
    pub(crate) file_offset: u64,
    pub(crate) row_byte_len: u32,
}
```

On-disk format per row (no outer length prefix — `row_byte_len` is in SpillIndex):

```text
[2 bytes: chrom_id (u16 LE)]
[4 bytes: start (u32 LE)]
[4 bytes: end (u32 LE)]
[1 byte:  strand]
[1 byte:  flags]
[variable: name, score, strand_raw, extra_fields, exon_coords — each with own length prefix]
[2+2 bytes: sample_count, bin_count (u16 LE)]
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
    returned_buf: Vec<(usize, f64, MatrixRow)>,
    chrom_table: ChromTable,
    temp_path: PathBuf,
}
```

`flush_chunk(rows, group_index) -> Result<FlushResult>`: runs on spawned thread, serializes all rows to a new temp file, returns indices + cleared-but-capacity-retained Vec.

`CollectorBucket`:
```rust
struct CollectorBucket {
    active: Vec<(usize, f64, MatrixRow)>,  // (orig_idx, sort_key, row)
    spare: Option<Vec<(usize, f64, MatrixRow)>>,
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

- [ ] **Step 1: Add `SpillingOutputWriter`**

A trait or struct that wraps final output destinations:

```rust
pub(crate) struct SpillingOutputWriter {
    matrix_writer: StreamingMatrixWriter,
    values_writer: Option<BufWriter<File>>,
    regions_writer: Option<BufWriter<File>>,
    group_labels: Vec<String>,
    group_boundaries: Vec<usize>,
}
```

`emit_row(&mut self, row: &MatrixRow)`: writes to main gzip + optional auxiliary.
- `outFileSortedRegions`: derives group label from `group_boundaries` + running row count.
- `outFileNameMatrix`: writes plain values row (header written at `finish()`).

`finish(self, header: &MatrixHeader)`: finalizes main gzip (rewrite header), writes `outFileNameMatrix` header (seek to start or use temp-file-then-copy approach), flushes all.

- [ ] **Step 2: Add `finalize_sorted`**

Joins all handles. For each bucket: mmaps all chunk files, collects `SpillRef { sort_key, orig_idx, mmap_idx, file_offset, row_byte_len }` from completed flushes, sorts in-memory `active` rows. Merge-emits spilled + in-memory using peekable iterators. **Stable tie-break**: when sort keys are equal, compare by `orig_idx` ascending.

```rust
struct SpillRef {
    sort_key: f64,
    orig_idx: usize,  // for stable tie-break
    mmap_idx: usize,
    file_offset: u64,
    row_byte_len: u32,
}
```

Sort comparator:
```rust
all_refs.sort_by(|a, b| {
    let cmp = compare_ascending(a.sort_key, b.sort_key);
    let cmp = if sort_ascending { cmp } else { cmp.reverse() };
    cmp.then(a.orig_idx.cmp(&b.orig_idx))
});
```

Cleanup: drop mmaps, remove temp files.

- [ ] **Step 3: Add `finalize_keep_order`**

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

- [ ] **Step 4: Tests**

Test finalize_sorted with low threshold, verify sort order. Test finalize_keep_order with gaps in orig_idx (simulating filtered rows).

- [ ] **Step 5: Commit**

```bash
git commit -m "feat: add finalize_sorted/keep_order with mmap readback and stable merge"
```

#### Task 2d: Integrate into executor.rs + run.rs

**Files:**
- Modify: `src/pipeline/core/executor.rs`
- Modify: `src/pipeline/run.rs`
- Modify: `src/io/writers/mod.rs`
- Modify: `src/io/writers/auxiliary.rs`
- Modify: `src/pipeline/core/collector.rs` (remove `GroupBucketCollector`)

- [ ] **Step 1: Add `OutputStrategy::Spilling` to executor**

```rust
enum OutputStrategy {
    StreamOrdered,
    Spilling,  // replaces InMemoryKeep + InMemoryGroupBucket for large matrices
}
```

The Spilling branch:
1. Creates `HybridBucketCollector`.
2. Compute loop: for each result `(orig_idx, group_index, Option<MatrixRow>)`, if `Some(row)`, compute sort_key and push.
3. After compute: calls `finalize_sorted` or `finalize_keep_order` depending on `sort_regions`.
4. The emit callback writes to a passed-in `SpillingOutputWriter`.

New executor function signature for spilling path:
```rust
pub fn execute_mode_spilling<M>(
    tasks: Vec<RegionTask>,
    general: &GeneralOptions,
    // ... same params ...
    output_writer: SpillingOutputWriter,
    task_count: usize,
) -> Result<()>
```

- [ ] **Step 2: Add third code path in `run.rs`**

```rust
fn needs_spilling(row_count: usize, sample_count: usize, total_bins: usize, sort_regions: SortRegions) -> bool {
    if matches!(sort_regions, SortRegions::No) {
        return false;
    }
    if matches!(sort_regions, SortRegions::Keep) {
        // Keep might go StreamOrdered if already sorted — executor decides.
        // But for large matrices that WILL be coalesced, estimate spilling.
        // Conservative: if cell_count is large AND sort != No, flag for spilling.
    }
    let cell_count = row_count.saturating_mul(sample_count).saturating_mul(total_bins);
    cell_count >= STREAMING_CELL_THRESHOLD
}
```

When `needs_spilling` is true:
```rust
let writers = SpillingOutputWriter::new(io, &group_labels, &group_capacity)?;
core::execute_mode_spilling(tasks, general, ..., writers, task_count)?;
return Ok(RunOutcome::Streamed);
```

Keep the existing InMemory path for small matrices (below threshold), and StreamOrdered for sort=no / sort=keep+already_sorted.

- [ ] **Step 3: Remove `should_use_streaming_for_plan` auxiliary guard**

Remove from `src/io/writers/mod.rs`:
```rust
if io.matrix_values_output.is_some() || io.sorted_regions_output.is_some() {
    return false;
}
```

- [ ] **Step 4: Add per-row streaming to `auxiliary.rs`**

```rust
pub fn write_sorted_region_row<W: Write>(writer: &mut W, row: &MatrixRow, group_label: &str) -> Result<()>
pub fn write_plain_values_row<W: Write>(writer: &mut W, row: &MatrixRow) -> Result<()>
```

- [ ] **Step 5: Remove `GroupBucketCollector` from `collector.rs`**

- [ ] **Step 6: Run all tests**

Run: `cargo test`

- [ ] **Step 7: Commit**

```bash
git commit -m "feat: integrate spilling path into executor and run.rs, remove GroupBucketCollector"
```

#### Task 2e: Integration Test for File Spilling

- [ ] **Step 1: Add spilling integration test with injectable threshold**

Use `with_threshold(100)` to trigger spilling on small test data. Verify:
- Sorted output matches expected order.
- Keep-order output with gaps (filtered rows) matches expected order.
- Round-trip through serialize → temp file → mmap → deserialize → emit produces correct values.

- [ ] **Step 2: Commit**

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
WARMUP_ARGS=()
for arg in "$@"; do
    WARMUP_ARGS+=("$arg")
done
# Replace -o argument with temp path to avoid clobbering
WARMUP_MODIFIED=()
SKIP_NEXT=false
for arg in "${WARMUP_ARGS[@]}"; do
    if $SKIP_NEXT; then
        WARMUP_MODIFIED+=("$WARMUP_DIR/warmup_output.mat.gz")
        SKIP_NEXT=false
    elif [ "$arg" = "-o" ] || [ "$arg" = "--outFileName" ]; then
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

This redirects the warm-up run's output to a temp directory, avoiding clobbering the actual profiling output.

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

## Appendix: Issues Addressed in v2

| # | Issue | Fix |
|---|---|---|
| 1 | Spilling rows re-collected into InMemoryCollector → OOM not solved | New `execute_mode_spilling` writes directly to `SpillingOutputWriter`, bypasses `MatrixData` |
| 2 | `finalize_keep_order` placement array sized by group_counts.sum(), filtered rows cause panic | Array sized by `task_count` (passed from run.rs), filtered slots are `None` |
| 3 | `join_all` retains FlushResult.returned_buf capacity | Clear `returned_buf` to `Vec::new()` after extracting for spare |
| 4 | Auxiliary output header needs group counts not yet known | `outFileSortedRegions`: static header + row-by-row with group label. `outFileNameMatrix`: rows first, header prepended at `finish()` |
| 5 | `load_matrix` only reads gzip, references are plain text, streaming produces multi-member gzip | Magic-byte detection: `[0x1f, 0x8b]` → `MultiGzDecoder`, else plain `BufReader` |
| 6 | Task 4 test params/filenames wrong vs YAML | All tests now match `python_compatibility.yaml` exactly |
| 7 | No stable tie-break in spilling sort | `SpillRef` carries `orig_idx`, sort comparator uses `.then(orig_idx.cmp)` |
| 8 | u16/u32 overflow in serialization | `u16::try_from().context()` in `serialize_row` |
| 9 | Spilling test doesn't trigger spilling | `with_threshold(threshold)` constructor for injectable threshold |
| 10 | `compare_values` only checks first row columns, `col_count_match` not in success | Check all rows, include in success condition |
| 11 | Task 4 requires manual `cargo build --release` | `env!("CARGO_BIN_EXE_...")` auto-builds |
| 12 | `profile_bench.sh` warm-up may clobber output | Warm-up uses temp dir for output |
