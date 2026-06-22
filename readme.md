compute-matrix-rs
=================

Rust reimplementation of deepTools `computeMatrix` (reference-point and scale-regions) targeting output compatibility with deepTools 3.5.6 at ≤5e-6 absolute tolerance.

Status
- Production-ready: CLI parity (including metagene), streaming output, and regression harness checked in.
- Manifest-driven compatibility suite passes against vendored deepTools fixtures.
- Test datasets are bundled under `deeptools/deeptools/test/`; larger ENCODE fixtures download into `target/compute-matrix-datasets/` on demand.
- **CLI caveat:** deepTools multi-letter short flags `-bs` (bin size) and `-bl` (blacklist) are not available; use `--bs/--binSize` and `--bl/--blackListFileName` instead because rust cli parsing library `clap` only supports single-letter short flags.

Quick Start
- Build/inspect the CLI: `cargo run -- --help`
- Run committed compatibility cases:
  `cargo test --test python_compatibility -- --test-threads=1`
- Run the unified harness:
  `pixi run compat`
- Regenerate and verify Python reference artifacts:
  `pixi run regen-artifacts` and `pixi run verify-artifacts`
- Prepare ENCODE benchmark data:
  `pixi run prepare-data encode_k562_atac`
- Local performance smoke run:
  `pixi run bench-smoke`

Latest Compatibility & Performance
- Command: `cargo test --test python_compatibility -- --test-threads=1`
- Result: all manifest `compat` cases pass with tolerance ≤5e-6.
- ENCODE benchmark cases live in `scripts/config/compute_matrix_cases.json` and run through `pixi run encode`.

Known Behavior Differences from Python deepTools
- **Header `scale` field type:** Python serializes the `scale` value as int when it is not explicitly passed (the default `1` is a Python int) and as float when explicitly passed (`1.0`, `2.0`, …). Rust always emits a JSON float. This is harmless — no downstream consumer (`plotHeatmap`, `plotProfile`, `computeMatrixOperations`) reads the `scale` header field back after matrix generation; it is only used during `computeMatrix` execution itself. Integration tests skip this field.
- **Blacklist + empty region group:** When `--blackListFileName` removes all regions in a group, Python deepTools only errors under `--sortRegions keep` (via `computeMatrixOperations.sortMatrix()`) and silently drops the empty group from the header for `no`/`ascend`/`descend`. We match this behavior for output parity, but note that silently dropping a group is arguably a design gap in deepTools — an empty group usually signals a blacklist/BED mismatch rather than intended behavior.

ENCODE K562 ATAC Benchmark (4 cores)
- Dataset files are prepared by `pixi run prepare-data encode_k562_atac`.
- Reference-point (`center`, ±100 bp, bin 10): Python 171.35s vs Rust 17.90s → **9.57× faster**.
- Scale-regions (body 200, ±100 bp, unscaled 50/50, bin 10): Python 346.56s vs Rust 18.64s → **18.59× faster**.
