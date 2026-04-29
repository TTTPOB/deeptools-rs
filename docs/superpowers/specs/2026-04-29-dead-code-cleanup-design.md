# Dead Code Cleanup & Missing Feature Spec

Date: 2026-04-29
Branch: `feat/performance-architecture-improvements`

## Overview

Remove dead code, eliminate stale abstractions, rename ambiguous types, and implement the missing `--blackListFileName` feature. The codebase was audited for dead code and over-abstraction; this spec captures all agreed-upon changes.

## Change 1: Replace BedReader with GroupedBedReader

**Problem:** `BedReader` (iterator over individual `BedRecord`s) is never used in production. The actual BED parsing in `regions.rs::parse_grouped_bed()` is inlined because it needs `#`-delimited group boundary handling, which `BedReader` cannot express.

**Design:**

- Delete `BedReader` from `src/io/readers/bed.rs` (struct, `from_path`, `new`, `read_all`, `Iterator` impl).
- Remove `BedReader` from `pub use` re-exports in `src/io/readers/mod.rs` and `src/io/mod.rs`.
- Create `GroupedBedReader` in `src/io/readers/bed.rs`:
  - Constructor: `GroupedBedReader::open(path, default_label)`.
  - Implements `Iterator<Item = Result<Group, BedReadError>>` where `Group { label: String, records: Vec<BedRecord> }`.
  - Handles `#` lines as group boundaries (current `finalize_group` logic).
  - Group labels emitted are **raw** (not deduplicated). Cross-file label deduplication remains in `load_groups()` via the existing shared `seen_labels: HashSet<String>`. This preserves the invariant that labels are globally unique across all input BED/GTF files.
  - No blacklist filtering — blacklist is applied at the signal level (see Change 2).
- Refactor `regions.rs::parse_grouped_bed()` to use `GroupedBedReader`.
- Move `Group` struct from `regions.rs` to `bed.rs` (it is the iterator's output type).

## Change 2: Implement --blackListFileName

**Problem:** The CLI accepts `--blackListFileName` and stores it in `GeneralOptions.blacklist`, but no pipeline code reads it. Users get silently incorrect results when using this flag.

**Python deeptools behavior:** Python's `mapReduce.blSubtract()` subtracts blacklist intervals from genome chunks before dispatching work. The exact effect on `computeMatrix` output (whether it is region-level filtering, signal-level masking, or a hybrid due to chunk boundaries) depends on chunk size and region/blacklist geometry. Rather than assume the semantics, we generate reference output first and implement to match.

**Implementation approach: reference-driven**

1. **Generate reference output first** (before writing implementation code):
   - Use pixi + deeptools 3.5.6 to run `computeMatrix` with `--blackListFileName` on test data.
   - Use existing test BED/bigWig files from `deeptools/` if they have sufficient blacklist overlap; otherwise create synthetic test data with known overlap patterns.
   - Scenarios to cover:
     - **reference-point, BED partial overlap:** Region partially overlaps a blacklist interval.
     - **reference-point, BED full overlap:** Region is entirely inside a blacklist interval.
     - **scale-regions, BED12/GTF metagene:** Blacklist falls in an intron of a multi-exon region.
     - **--missingDataAsZero interaction:** Same scenarios with `--missingDataAsZero`.
   - Store reference output as test fixtures alongside existing `tests/` data.

2. **Analyze reference output** to determine exact Python behavior:
   - Are partially-overlapping regions present in output? If yes, what values do blacklisted bins have?
   - Are fully-overlapping regions present or dropped?
   - How does `--missingDataAsZero` interact with blacklisted positions?

3. **Implement to match observed behavior.** The most likely implementation:

   - **Loading blacklist:** In `pipeline/mod.rs::execute()` (or `run_pipeline`), if `general.blacklist` is `Some`:
     - Read the blacklist BED file using `GroupedBedReader::open(path, "blacklist")`.
     - Flatten all groups, discard labels, collect into `Vec<(Arc<str>, u32, u32)>` of (chrom, start, end).
     - Sort by (chrom, start) and merge overlapping intervals.
     - Store as `Arc<Vec<...>>` and pass through to the worker layer.

   - **Applying blacklist:** In `worker.rs`, mask blacklisted positions in the signal. The exact mechanism (coverage buffer zeroing, query range splitting, or region-level filtering) is determined by the reference analysis in step 2. All three worker paths must be covered:
     - `compute_sample_bins` (per-item fallback path)
     - `process_batch` coverage buffer path (Median/Std/Min/Max)
     - `process_batch` direct aggregation path (Mean/Sum)

   - **Threading blacklist to workers:** Add blacklist as a parameter to:
     - `execute_mode()` — receives the loaded blacklist.
     - `process_batch()` and `compute_row()` / `compute_sample_bins()` — receive a reference to the blacklist slice.
     - The blacklist is shared via `Arc`, no per-worker copies.

   - **Overlap helper:** Add `fn blacklist_intervals_for_chrom(blacklist, chrom) -> &[(u32, u32)]` that binary-searches the sorted blacklist for the relevant chrom's intervals.

4. **Validate** with parity tests (see Testing Strategy).

## Change 3: Rename BigWigReader / SharedBigWigReader

**Problem:** `SharedBigWigReader` sounds like "a shared reader" but it is actually the opened file with metadata. `BigWigReader` sounds like the primary reader but is actually a per-worker cursor.

**Design:**

- Rename `SharedBigWigReader` → `BigWigFile` everywhere.
- Keep `BigWigReader` name (it is the per-worker reader, and this name is now clearer paired with `BigWigFile`).
- Update all re-exports in `src/io/mod.rs` and `src/io/readers/mod.rs`.
- Update type alias `SharedBigWigReader` in `src/io/mod.rs` re-export (remove the old name).

## Change 4: Fold RowCollector trait into FileCollector

**Problem:** `RowCollector` trait has only one implementor (`FileCollector`). The in-memory collector was removed. The trait adds indirection with no polymorphism benefit.

**Design:**

- Delete the `RowCollector` trait from `src/pipeline/core/collector.rs`.
- Convert `on_row`, `finalize`, and `abort` to inherent methods on `FileCollector`.
- Update call sites in `executor.rs` to call methods directly on `FileCollector`.

## Change 5: Delete SpillIndex::group_index

**Problem:** `SpillIndex::group_index` is written during serialization but never read. Each `CollectorBucket` already tracks its own `group_index`, making this field redundant.

**Design:**

- Remove `group_index` field from `SpillIndex`.
- Remove `group_index` parameter from `flush_chunk()`.
- Remove `#[allow(dead_code)]` annotation.

## Change 6: Delete remaining dead code

The following items are deleted with no replacement:

| Item | Location | Reason |
|------|----------|--------|
| `MatrixRow::flattened_values()` | `matrix.rs:221-223` | Never called; trivial `.values.clone()` wrapper |
| `sample_boundaries_from_counts()` | `matrix.rs:248-257` | Never called; `sample_boundaries_uniform()` used instead |
| `WorkerSamples::new()` | `samples.rs:72-88` | Old path; production uses `from_shared()` |
| `WorkerSamples::samples()` | `samples.rs:114-119` | Never called; `samples_and_bufs()` used instead |
| `open_samples()` | `samples.rs:132-138` | Only called by dead `WorkerSamples::new()` |
| `Sample::open()` | `samples.rs:14-22` | Only called by dead `open_samples()` |
| `Sample::reader()` (immutable) | `samples.rs:28-30` | Only called by dead `WorkerSamples::new()` for buf size |
| `BigWigReader::chroms()` | `bwig.rs:394-396` | Never called in production |
| `SharedBlockCache::new()` | `block_cache.rs:45-47` | Never called; `with_capacity()` used instead |
| `SharedBigWigReader::block_cache()` | `bwig.rs:191-193` | Never called (accessor for the block cache field) |
| `SharedBigWigReader::open()` (no cache param) | `bwig.rs:131-135` | Only used by dead `BigWigReader::open()` |
| `BigWigReader::open()` | `bwig.rs:376-380` | Only used by dead `Sample::open()` |

Note: the table above uses pre-rename names. After Change 3 renames `SharedBigWigReader` → `BigWigFile`, the deletions apply to the renamed type.

Also remove:
- `#[allow(unused_imports)]` on `use std::sync::Arc` in `zones/mod.rs` if Arc is no longer needed after cleanup.

## Items explicitly kept

- `BigWigReaderStats` — will be used in the future for diagnostics.
- `HybridBucketCollector` sample_count/bin_count stored fields with `#[allow(dead_code)]` — acceptable as-is (user confirmed).

## Testing strategy

- `cargo check` — zero warnings (no dead_code, no unused imports).
- `cargo test` — all existing tests pass. Tests that reference deleted APIs are updated or removed.
- Unit tests:
  - `GroupedBedReader`: grouping, `#` boundary handling, raw label emission.
  - Cross-file label deduplication in `load_groups()` (two files with same `#` label → `_1` suffix).
  - Blacklist overlap helper functions.
- **Blacklist parity tests (mandatory, not deferred):**
  - Reuse existing heatmapper test data (`test.bw` with ch1/ch2/ch3 @ 400bp, `test.bed`/`test2.bed`).
  - Create a synthetic blacklist BED file `test_blacklist.bed` with two intervals:
    - `ch1 110 130` — partial overlap with the ch1 region (100-150), masks part of the signal at 100-125.
    - `ch2 140 180` — full overlap with the ch2 region (150-175), masks all signal for that region.
  - Generate reference output using pixi + deeptools 3.5.6:
    1. `computeMatrix reference-point -S test.bw -R test2.bed -b 100 -a 100 -bs 1 --blackListFileName test_blacklist.bed` → `master_blacklist.mat`
    2. `computeMatrix scale-regions -S test.bw -R test2.bed -b 100 -a 100 -m 100 -bs 1 --blackListFileName test_blacklist.bed` → `master_scale_reg_blacklist.mat`
  - Compare Rust output vs deeptools reference using `compare_matrix diff`.
  - These tests are added to `tests/python_compatibility.rs`.
  - The reference output also serves as ground truth for resolving any ambiguity about Python's exact blacklist behavior (see Change 2, step 2).
- Existing integration tests in `tests/python_compatibility.rs` continue to pass.
