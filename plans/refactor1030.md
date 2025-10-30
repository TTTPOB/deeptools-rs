# Refactor Roadmap — 2025-10-30

## Code Path Inventory
- `src/main.rs` → `src/cli.rs` builds a `Config` and hands control to `pipeline::execute` in `src/pipeline/mod.rs`.
- `pipeline::execute` fans out into two mode drivers:
  - Reference-point: `src/pipeline/reference_point.rs::run`.
  - Scale-regions: `src/pipeline/scale_regions.rs::run`.
- Each driver currently supports two execution surfaces:
  1. **In-memory matrix** – build `RegionTask` list, process with Rayon, collect `MatrixRow`s into `MatrixData`, then spawn a writer thread via `io::writers::write_outputs`.
  2. **Streaming matrix** – same task scheduling, but rows are sent over an `mpsc` channel to a `StreamingMatrixWriter` (see `stream_*_rows` in each mode) that buffers, counts groups, and rewrites the header footer in-place.
- Shared binning happens inside `pipeline::core::compute_row`, which consumes a mode-specific `RegionPlan` (`zones::ReferencePointPlan` or `zones::ScaleRegionsPlan`) before post-processing (`nan_after_end` masking or strand reversal) in the mode driver.
- Header construction and final sorting/pruning live in `pipeline::matrix`, invoked once per driver after computation.

## Redundant / Fragile Path Segments
- **Task scaffolding duplication** – Both drivers define identical `RegionTask`/`RegionResult` structs, group-capacity bookkeeping, Rayon pool creation, and `map_init` closures (cf. `reference_point.rs` lines ~24–180 and `scale_regions.rs` lines ~24–190).
- **Streaming aggregation clone** – `run_streaming_*` and `stream_*_rows` functions are near copies differing only in header builder and row post-processing (`nan_after_end` vs strand reverse). The `BTreeMap` buffering, channel wiring, and completion checks are duplicated.
- **Validation helpers repeated** – `ensure_positive` and `ensure_multiple` exist in both mode files with identical bodies.
- **Header assembly forks** – `build_reference_point_header` and `build_scale_regions_header` share the same skeleton with only field payload differences; future options will exacerbate drift if maintained separately.
- **Row post-processing asymmetry** – Reference-point injects `nan_after_end`, scale-regions reverses bins on negative strand. Both paths currently bake this logic inline rather than plugging into a shared post-step hook, making additional behaviours (e.g. upcoming GTF exon expansion) harder to compose.

### Streaming vs Buffered Behaviour Check
- Verified that zero-pruning and threshold checks already occur before rows reach the writer: `core::compute_row` (via `should_skip_row`) returns `None` for all-zero/threshold-violating rows when `general.skip_zeros` or min/max thresholds are set. Both streaming and buffered runners drop `None` rows, so the streaming path does not emit zero-only records today.
- The buffered path performs a second guard via `MatrixData::prune_zero_rows`, which re-checks after sorting; the streaming variant bypasses that pass. Once the collector abstraction lands, ensure the shared row pipeline keeps the pre-write pruning hook so both sinks remain consistent even if sort-dependent filtering grows more complex.

## Abstraction Direction (Rust Style)
- Introduce a `PipelineMode` trait (or generic struct) encapsulating the mode-specific pieces: argument validation, plan construction per record, per-row post-processing, header assembly, and any static metadata (e.g. `nan_after_end` default). Each mode implements the trait, enabling a single generic driver to orchestrate work.
- Build a reusable `TaskRunner` in `pipeline::core` that:
  1. Constructs tasks from grouped BED records.
  2. Drives Rayon computation with shared `WorkerSamples`.
  3. Chooses between in-memory vs streaming via a pluggable strategy, supplying callbacks for “row ready” and “header finalize”.
  This aligns with idiomatic Rust patterns—generic functions over traits plus small structs carrying lifetimes—while keeping compile-time dispatch.
- Replace ad-hoc header builders with a `MatrixHeaderBuilder` struct that receives per-mode invariants (upstream/downstream lengths, unscaled zones, ref-point labels) and shared knobs (sort metadata, thread count). This reduces duplicate JSON assembly and enforces consistent defaults.
- Centralize validation helpers into `pipeline::core::validation` (or similar) to avoid drift and enable richer error contexts with mode tags.
- For strand-sensitive transforms, define a lightweight `RowPostProcessor` enum (e.g. `None`, `ReverseBins`, `MaskBeyondRegion`) so the generic runner can apply them uniformly.

## Stepwise Refactor Plan
1. **Extract shared task driver** – ✅ 2025-10-30: introduced a mode-agnostic row aggregator (`spawn_row_aggregator`) so both pipelines share the same Rayon loop while swapping only the row collection strategy.
2. **Introduce `PipelineMode` trait** – ✅ 2025-10-30: added a shared trait (`pipeline::core::PipelineMode`) plus generic `execute_mode` runner; both reference-point and scale-regions now implement the trait and delegate scheduling/header wiring through it.
3. **Unify row aggregation targets** – ✅ 2025-10-30: replaced the stop-gap row sink with a dedicated `RowCollector` trait (implemented by `FileCollector` and `InMemoryCollector`) so `execute_mode` now streams rows directly to gzip or yields a ready-to-sort `MatrixData` without extra conversion.
4. **Consolidate header construction** – ✅ 2025-10-30: introduced `MatrixHeaderBuilder` + `LayoutVectors`, rewired both modes to share it, and removed the bespoke header helpers now that the builder emits identical metadata for streaming and buffered runs.
5. **Generalize validation utilities** – ✅ 2025-10-30: lifted the `ensure_*` helpers into `pipeline::core`, added mode-tagged error messages, and pointed both drivers at the shared checks for consistent CLI diagnostics.
6. **Regression + performance guard rails** –✅ confirmed, correct. Re-run pixi-backed comparisons for both modes, and measure streaming/non-streaming throughput to ensure refactor does not regress hotspots documented in `plans/write_performance.md`.
