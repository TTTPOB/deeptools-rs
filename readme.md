compute-matrix-rs
=================

Rust reimplementation of deepTools `computeMatrix` (reference-point and scale-regions) targeting output compatibility with deepTools 3.5.6 at ≤5e-6 absolute tolerance.

Status
- Production-ready: CLI parity (including metagene), streaming output, and regression harness checked in.
- Python compatibility suite (10 scenarios from `test_heatmapper.py`) passes against vendored fixtures.
- Test datasets are bundled under `deeptools/deeptools/test/`; larger ENCODE fixtures download into `target/compute-matrix-datasets/` on demand.
- **CLI caveat:** deepTools multi-letter short flags `-bs` (bin size) and `-bl` (blacklist) are not available; use `--bs/--binSize` and `--bl/--blackListFileName` instead because rust cli parsing library `clap` only supports single-letter short flags.

Quick Start
- Build/inspect the CLI: `cargo run -- --help`
- Run the compatibility suite (uses pixi with deepTools 3.5.6):
  `pixi run python scripts/custom_compare.py --mode python-compatibility`
- Ad-hoc regression/perf run (reference-point example):
  `pixi run python scripts/custom_compare.py --mode reference-point --reference-point center --upstream 100 --downstream 100 --bin-size 10`
- Local performance smoke run:
  `scripts/perf_smoke.sh`

Latest Compatibility & Performance (2025-12-01)
- Command: `pixi run python scripts/custom_compare.py --mode python-compatibility`
- Result: 10/10 tests passed (tolerance ≤5e-6)
- Performance: total Rust time 11.48s (0 cached, 10 fresh); slowest case `reference_point_basic` at 9.09s on the deepTools test corpus.

ENCODE K562 ATAC Benchmark (4 cores)
- Dataset/command hash cached under `target/*-regression/.cache`.
- Reference-point (`center`, ±100 bp, bin 10): Python 171.35s vs Rust 17.90s → **9.57× faster**.
- Scale-regions (body 200, ±100 bp, unscaled 50/50, bin 10): Python 346.56s vs Rust 18.64s → **18.59× faster**.
