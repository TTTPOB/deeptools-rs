#!/usr/bin/env bash
set -euo pipefail

cargo build --release

mkdir -p target/perf-smoke

run_case() {
  local name="$1"
  shift
  local log="target/perf-smoke/${name}.log"
  local out="target/perf-smoke/${name}.mat.gz"

  rm -f "$out" "$log"
  echo "[perf-smoke] running ${name}"
  /usr/bin/time -f "elapsed=%e user=%U sys=%S max_rss_kb=%M" \
    -o "$log" \
    target/release/compute_matrix_rs "$@" --outFileName "$out"
  cat "$log"
}

run_case reference_point_basic \
  reference-point \
  --referencePoint center \
  -R deeptools/deeptools/test/test_data/genes.bed \
  -S deeptools/deeptools/test/test_data/test.bw \
  --beforeRegionStartLength 100 \
  --afterRegionStartLength 100 \
  --binSize 10 \
  --numberOfProcessors 4

run_case scale_regions_basic \
  scale-regions \
  -R deeptools/deeptools/test/test_data/genes.bed \
  -S deeptools/deeptools/test/test_data/test.bw \
  --regionBodyLength 200 \
  --beforeRegionStartLength 100 \
  --afterRegionStartLength 100 \
  --binSize 10 \
  --numberOfProcessors 4
