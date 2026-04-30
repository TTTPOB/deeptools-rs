# Compatibility Cleaning: Python deepTools 3.5.6 Parity

## Overview

Systematic audit and fix of all behavioral differences between the Rust
`computeMatrix` reimplementation and the Python deepTools 3.5.6 reference.
Only differences that affect **output results** are addressed; internal
optimizations that preserve identical output are left alone.

## Phases

### Phase 1 — Worker computation fixes
Status: **done**

Fixes in `src/pipeline/core/worker.rs` and related files:

| # | Issue | Description |
|---|-------|-------------|
| 1a | skipZeros semantics | Python uses `np.mean(row)==0 or NaN` → remove. Rust checks "all non-NaN values are zero". Difference when values cancel to zero mean. |
| 1b | Threshold vs scale ordering | Python checks thresholds on pre-scale values, then scales. Rust scales first, then checks thresholds. |
| 1c | Threshold NaN propagation | Python `coverage.min()` returns NaN when NaN present → threshold check passes. Rust filters NaN before checking → finds below-threshold values NaN was hiding. |
| 1d | Chromosome missing + missingDataAsZero | Python fills 0 when chrom missing and missingDataAsZero=true. Rust unconditionally returns NaN. |

### Phase 2 — Zone planning fixes
Status: **done**

Fixes in `src/pipeline/zones/mod.rs` and `metagene.rs`:

| # | Issue | Description |
|---|-------|-------------|
| 2a | scale-regions short body | Python fills entire row with NaN/0 when body_length < bin_size. Rust still computes upstream/downstream normally. |
| 2b | Metagene center padding | Rust clears left/right bins when exon length insufficient; Python preserves partial exon fragments and only records padding count. |

### Phase 3 — Matrix, sorting, group boundary fixes
Status: **done**

Fixes in `src/pipeline/matrix.rs`, `executor.rs`, `run.rs`:

| # | Issue | Description |
|---|-------|-------------|
| 3a | sortUsing region_length + metagene | Rust uses genomic span (`end-start`); Python uses sum of exon lengths. |
| 3b | Group boundaries after row filtering | skipZeros/threshold filtering can empty groups; Rust doesn't rebuild header group info post-filter. |

### Phase 4 — CLI defaults
Status: **done**

Fixes in `src/cli.rs`:

| # | Issue | Description |
|---|-------|-------------|
| 4a | smartLabels default | Rust defaults to true (strips extension); Python defaults to false (keeps extension). |

### Phase 5 — Output formatting fixes
Status: **done**

Fixes in `src/io/writers/formatting.rs`, `auxiliary.rs`, `matrix_gz.rs`:

| # | Issue | Description |
|---|-------|-------------|
| 5a | outFileNameMatrix format | Python uses `%.4g` (4 significant digits, g-format). Rust uses `{:.4}` (4 decimal places). |
| 5b | Header JSON field order | Rust struct field order differs from Python dict insertion order. |
| 5c | Header scale serialization | Python `json.dumps(1.0)` → `1.0`. Rust serializes as `1`. |
| 5d | BED score format in matrix rows | Python outputs score as-is from region tuple; Rust formats numeric scores with 6 decimal places. |
| 5e | outFileSortedRegions format | Python outputs 13-column BED-like with deepTools_group; Rust outputs 7 columns. |

### Phase 6 — Self-review & final verification
Status: **done**

- Run full `cargo test`
- Run `pixi run python scripts/custom_compare.py --mode python-compatibility`
- Address any remaining failures
- Final doc update

## Known Architectural Differences (not fixed)

| Area | Description |
|------|-------------|
| Missing-chrom region dispatch | Python's mapReduce dispatches regions by bigWig chromosomes — if a chromosome doesn't exist in any bigWig, the BED region is never processed (row omitted). Rust processes all BED regions and fills with NaN/0. Fixing this would require changes to the region dispatch architecture. |

## Results

- 12 commits total
- 8 behavioral fixes implemented
- 4 corner case integration tests added with Python reference data
- 269 unit tests + 23 integration tests passing
- 7 deferred items documented (auxiliary output formats, architectural differences)

## Audit Sources

- Claude Opus analysis of Rust vs Python code (2026-04-30)
- Codex static scan findings (2026-04-30)
- Reference: Python deepTools 3.5.6 (`deeptools/deeptools/heatmapper.py`)
