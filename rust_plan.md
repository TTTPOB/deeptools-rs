# Rust `computeMatrix` Reimplementation Plan

## Objectives
- Keep the Rust CLI interface 1:1 compatible with the Python `computeMatrix` command (flag names, defaults, mutually-exclusive groups, help text semantics). If a perfect match is infeasible, document the gap and choose the closest behavior.
- Match the Python pipeline in `deeptools/deeptools/computeMatrix.py` and `deeptools/deeptools/heatmapper.py` feature-for-feature for both `reference-point` and `scale-regions` modes, prioritizing identical output over new optimizations.
- Preserve on-disk formats produced by `heatmapper.heatmapper.save_matrix()` and `save_matrix_values()` so downstream DeepTools consumers remain compatible.
- Leverage the `bigtools` crate for bigWig IO and `polars` for tabular/columnar manipulation where aggregate statistics are needed.
- Use existing fixtures under `deeptools/deeptools/test/test_data/` (e.g. `computeMatrixOperations.bed`, `testA.bw`, etc.) as regression inputs to confirm byte-identical gzipped matrix output compared to the Python CLI.

## Python Implementation Notes (Ground Truth)
- CLI orchestration lives in `computeMatrix.py` (argument parsing, mode selection, sorting, file emission). Heavy lifting delegates to `heatmapper.heatmapper`, which:
  - Loads bigWig metadata via `getScorePerBigWigBin.getChromSizes()` and dispatches work through `mapReduce.mapReduce` (chunked by chromosome windows).
  - Uses `compute_sub_matrix_worker()` to construct each region’s binned values, handling strand-aware upstream/downstream windows, optional unscaled flanks, NaN padding, and reference-point logic (`TSS`, `TES`, `center`).
  - Wraps the final numpy array in a `_matrix` helper that tracks group/sample boundaries, labels, and optional sorting/skipping rules before calling `save_matrix`.
- Output file header is JSON followed by BED-like lines: `chrom\tcomma_joined_starts\tcomma_joined_ends\tname\tscore\tstrand\t<tab-separated values>`.
- Tests for expected matrix content reside in `deeptools/deeptools/test/test_computeMatrixOperations.py` and related fixtures in `test_data/`.

## Planned Work Breakdown
1. **Spec Extraction & Validation (Python)**
   - Catalogue required CLI options by reviewing `computeMatrixRequiredArgs`, `computeMatrixOptArgs`, and dependent helpers (gtf options, filtering, sorting flags). Record expected defaults, aliases, and help text to mirror in Clap.
   - Trace parameter normalization inside `heatmapper.computeMatrix()` to understand derived fields (bin counts, group boundaries, labels).
   - Document edge cases from `compute_sub_matrix_worker()` (short regions, NaN handling, strand flips, `nan_after_end`, unscaled 5'/3' sections, padding behavior).
   - Identify minimum viable subset of options needed to reproduce regression fixtures; mark advanced options (clustering, silhouette) for later phases.

2. **Rust Project Architecture**
   - Introduce workspace structure with modules: `cli` (clap-based argument parsing mirroring Python flags/aliases exactly), `input` (BED/GTF parsing, grouping, blacklist support), `bigwig` (wrappers around `bigtools::bigwig` readers, caching chrom sizes), `matrix` (data model akin to `_matrix` with group/sample boundaries), and `pipeline` (mode-specific binning engines).
   - Establish shared configuration struct analogous to Python `parameters` dict, typed and validated on construction.
   - Plan concurrency model (likely Rayon) to parallelize region chunks similar to `mapReduce.mapReduce`.

3. **Core Pipeline Implementation**
   - BED/GTF ingestion: parse regions into exon lists, keep group segmentation indicated by `#`, support strand and score metadata; cross-check with Python expectations from `parserCommon`.
   - Chromosome metadata: use `bigtools` to read chromosome sizes once and validate requested windows.
   - Region binning: implement Rust equivalent of `compute_sub_matrix_worker`, covering:
     - Computation of upstream/downstream/unscaled/body zones with strand-aware bin counts.
     - Handling of `reference-point` modes (`TSS`, `TES`, `center`) and `scale-regions`.
     - NaN vs zero padding decisions toggled by flags like `missingDataAsZero`, `nanAfterEnd`.
     - Aggregation of bigWig signal across potentially gapped exons (respecting masking/padding).
   - Matrix assembly: stitch per-region vectors into contiguous matrices, maintain group/sample boundaries, optional skip-zero filtering, and sorting hooks.

4. **Output & Interop**
   - Serialize header JSON matching Python key casing and list semantics (convert scalar params to per-sample lists where required).
   - Emit gzipped matrix and optional plain tabular output (`save_matrix_values` behavior) using `flate2`/`xz2` as needed.
   - Ensure command-line UX mirrors DeepTools (help text, version, quiet/verbose flags).

5. **Validation Strategy**
   - Build thin Python harness to generate expected outputs via upstream script for fixtures in `deeptools/deeptools/test/test_data/`.
   - Create Rust integration tests comparing gzip payloads line-by-line (after normalizing floating-point formatting to `%f` with 6 decimals, like Python).
   - Add sanity tests for edge cases: reverse-strand regions, short bodies, `scale-regions` with unscaled flanks, `nanAfterEnd`, blacklist filtering.
   - Set up CI workflow (cargo test) to guard regressions once implementation lands.

6. **Performance Follow-up**
   - Defer profiling until functional parity is confirmed; placeholder tasks include benchmarking against Python on large synthetic inputs and investigating memory layout improvements via Polars (e.g., columnar storage for batched operations) or Arrow buffers.

## Existing Test Assets to Reuse
- `deeptools/deeptools/test/test_data/computeMatrixOperations.bed` plus the bundled bigWig files (`testA.bw`, `testA_offset*.bw`, `testB.bw`, etc.) cover multi-sample, multi-group scenarios.
- `deeptools/deeptools/test/test_computeMatrixOperations.py` outlines expected group/sample arrangements and provides reference command lines to replicate.
- Additional edge-case inputs (e.g., `test.gtf`, `othergenes.txt.gz`) support GTF parsing and complex exon layouts.

## Immediate Next Steps
- Flesh out CLI spec doc & config schema based on Python arg groups, ensuring flag names, aliases, defaults, and mutual exclusivity mirror the Python CLI.
- Prototype bigWig reader module calling `bigtools` to load a narrow window and compare values against `pyBigWig` on fixture data.
- Draft data structures for region representation and matrix buffering, ensuring we can round-trip through `save_matrix` serialization tests.
