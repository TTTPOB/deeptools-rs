# Streaming `matrix.gz` Investigation

## Current Behaviour Snapshot
- `reference_point::run` returns a `MatrixData` that owns `header`, full `Vec<MatrixRow>`, and counters (`src/pipeline/reference_point.rs:160`).
- `write_outputs` serialises the matrix synchronously on the caller thread; the gz encoder is fed by iterating the full `rows` vector (`src/io/writers/mod.rs:13`–`55`).
- Sorting (`MatrixData::sort_groups`) and zero-row pruning mutate the in-memory `Vec<MatrixRow>` before any IO happens (`src/pipeline/matrix.rs:114` onward).
- Future clustering hooks (`MatrixData::group_stats`, `hmcluster_placeholder`) assume random read access to all rows.

Peak memory therefore scales with `rows.len() * samples * bins`, matching Python’s behaviour but preventing us from handling million-region workloads without resorting to swap.

## Constraints That Block Naïve Streaming
- **Header upfront**: the deepTools header embeds `group_boundaries` and `sample_boundaries`; these depend on the final row counts after skip-zero/threshold pruning. We cannot write the header (and thus start the gzip stream) until every region has been evaluated.
- **Ordering guarantees**: `sortRegions ascend|descend` reorders rows within each group. Streaming rows as soon as their worker finishes would break this invariant unless we buffer an entire group until its sort metric is computed.
- **Shared consumers**: optional outputs (`matrixValues`, `sortedRegions`) and upcoming clustering features need the same row data in the final sorted order. Any streaming solution must either broadcast rows to multiple sinks or provide a backing store that can be re-read.

## Proposed Architecture

### 1. Split matrix assembly from output sinks
Introduce a `MatrixStream` abstraction that owns:
- the immutable `MatrixHeader`,
- a `RowSource` enum with two modes:
  1. `InMemory(Vec<MatrixRow>)` (status quo, useful for tests),
  2. `Spool(SpoolIndex)` pointing at a temporary on-disk buffer.

`MatrixData::into_stream()` would convert the current struct into this streaming form by moving the rows and (optionally) spooling them to disk. This isolates streaming logic from computation and keeps clustering entry points talking to a single type.

### 2. Dedicated writer thread via row sinks
Define a `RowSink` trait with `fn write_rows(&mut self, header: &MatrixHeader, source: RowSource) -> Result<()>` plus a helper to move rows through bounded channels. Implementations:
- `GzipSink` — spawns a writer thread that owns `GzEncoder<File>` and drains a `crossbeam_channel` receiver. Recommended settings:
  - Bounded channel (~4–8 groups) to cap memory if we stay in-memory.
  - `RowCommand::Header(MatrixHeader)` followed by `RowCommand::Row(MatrixRow)` messages, ending with `RowCommand::Finish`.
- `PlainSink` and `SortedRegionsSink` reuse the same plumbing for other outputs, keeping serialization off the compute thread.

`write_outputs` becomes a coordinator that:
1. Chooses a `RowSource` (in-memory vs spool) based on CLI flag (`--stream-output`?) or heuristic.
2. Spins up sinks as needed, cloning the source for each consumer (spool handles multi-reader; in-memory copies the vector).
3. Waits for all writer threads to finish.

### 3. Getting rid of the in-memory vector
Two viable strategies:

**Option A – Spool-once**
- While workers are running we accumulate rows in a `tempfile` (e.g. bincode or Arrow). Each entry stores: group index, sort metric, and the raw values.
- After computation completes we know `group_counts`; we sort per group by reading from the spool into chunked vectors (one group at a time), then immediately feed them to the sinks.
- Memory bound drops to “largest group × sample bins” instead of entire matrix. Disk footprint equals final output but uncompressed.

**Option B – Dual-pass header discovery**
- First pass computes skip-zero decisions + sort metrics but discards per-bin values, only tallying counts. (We can reuse `compute_sample_bins` but drop the vector once thresholds evaluated.)
- Second pass recomputes the bins and streams directly to the writer thread, now that the header is known. This doubles BigWig IO but keeps implementation simple and allows true streaming to the final gzip file.

Both options keep the gzip writer fully streaming (constant RAM) once the header is available. Option A favours single-pass IO with extra disk; Option B favours zero extra disk at the cost of CPU/IO.

### 4. Clustering extensibility
- The `RowSource::Spool` variant can expose an iterator that yields `MatrixRow` chunks; clustering algorithms (k-means, hclust) can reuse this iterator without forcing everything back into RAM.
- Provide a trait `RowView` that abstracts over `Iterator<Item = MatrixRow>` vs `ChunkedReader`. Clustering routines can accept any `RowView`, allowing them to operate on streaming chunks or in-memory vectors interchangeably.
- For algorithms requiring random access (e.g. repeated centroid updates), back the spool with `memmap2` so chunks can be mmapped into process space without copying.

## Implementation Outline
1. **Refactor `MatrixData`**
   - Add `MatrixData::into_spool(spooler: &mut SpoolWriter) -> Result<MatrixStream>` that writes rows into a temp file and returns indexing metadata.
   - Add `MatrixStream::for_each_row<F>(&self, f: F)` to abstract iteration source.
2. **Introduce concurrency plumbing**
   - Add `RowCommand` enum and `RowSink` trait under `io::writers`.
   - Implement `GzipSink`, `PlainSink`, `SortedRegionsSink` using Rayon’s thread pool or dedicated `std::thread::spawn`.
3. **Update CLI / config**
   - Add optional flag `--stream-output` (default `auto`) mapped to an enum `RowStorageMode` controlling whether we stay in-memory, spool, or dual-pass.
4. **Refine reference-point pipeline**
   - Replace the `Vec<MatrixRow>` return type with `MatrixStream`.
   - During computation, capture per-group counts and sort metrics alongside spool offsets so that sorting can be performed group-by-group without rehydrating the entire matrix.
5. **Adapt downstream consumers**
   - Rewrite `write_outputs` to accept `MatrixStream` and dispatch to sinks.
   - Ensure regression harness still compares byte output by consuming the stream.
6. **Testing hooks**
   - Add property tests that stream through a tiny bounded channel to ensure backpressure works and all rows are written.
   - Extend the pixi regression harness to run in both `in-memory` and `streaming` modes and diff outputs against deepTools.

## Risks & Open Questions
- Tempfile format choice (custom bincode vs Arrow/Parquet) affects clustering reuse; we need to balance implementation effort with tooling.
- Dual-pass mode doubles BigWig reads; need benchmarking to confirm throughput remains acceptable.
- Crossbeam dependency adds ~70 KB but simplifies bounded channels; verify policy allows it.
- Mmap-based clustering will require careful lifetime management to avoid holding stale pointers after the spool is deleted.

## Recommendation
Start with Option A (spool-once) because it preserves single-pass BigWig reads and gives us a durable intermediate representation that clustering can reuse. Hide the complexity behind `RowSource` so we can fall back to the current in-memory path when desired. Once the streaming sinks are stabilised, revisit Option B if profiles show spare IO budget.
