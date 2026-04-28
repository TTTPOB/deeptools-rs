# Performance Optimization Handoff

## Baseline (commit 5acc6ab)

| Metric | Reference-point | Scale-regions |
|--------|----------------|---------------|
| Wall time | 8.59s | 13.70s |
| User CPU | 24.00s | 39.99s |
| Max RSS | 416 MB | 415 MB |
| Parity | 10/10 PASS | 10/10 PASS |

ENCODE: 269,800 regions, 3 bigWig samples. Reference-point: 20 bins/sample. Scale-regions: 50 bins/sample.

## Architecture Overview

Hot path is `compute_sample_bins()` in `src/pipeline/core/mod.rs` — called for every (region, sample) pair. For scale-regions encode: 809,400 calls.

### Call chain:
1. `execute_mode()` — rayon parallel over RegionTasks
2. `compute_row()` — per-region, iterates samples
3. `compute_sample_bins()` — per-sample: allocates coverage Vec, fetches bigWig intervals, fills coverage, aggregates per-bin
4. `aggregate_slice()` — per-bin: allocates Vec<f64>, collects valid values, computes statistic
5. `write_matrix_row()` — per-row: thread-local buffer, formats values

## Optimization Candidates (priority order)

### 1. aggregate_slice — single-pass computation (HIGH impact, LOW risk)
**File**: `src/pipeline/core/mod.rs:408-456`
**Problem**: Every bin aggregates by allocating `Vec<f64>`, collecting + filtering, then computing. For scale-regions encode: 269,800×3×50 = 40M allocations.
**Fix**: Compute Mean/Sum/Min/Max in a single pass over the slice without collecting. Only collect for Median/Std (rare modes).
**Memory saving**: Eliminates ~40M temporary Vec<f64> per encode run.
**Speed**: Avoids allocation + collect overhead on the hottest inner loop.

### 2. Coverage buffer reuse (HIGH impact, LOW risk)
**File**: `src/pipeline/core/mod.rs:309`
**Problem**: Each `compute_sample_bins` call allocates `vec![default_fill; window_len]`. For encode: ~800K allocations.
**Fix**: Add a thread-local `Vec<f32>` buffer that resizes to max window_len and gets reused (cleared/reset each call).
**Memory saving**: ~800K fewer allocations, peak memory unchanged.
**Speed**: Reduced allocator pressure in hot loop.

### 3. MatrixRow value flattening (MEDIUM impact, MEDIUM risk)
**File**: `src/pipeline/matrix.rs:199-216`
**Problem**: Values stored as `Vec<Vec<f32>>` — N+1 allocations per row, poor cache locality.
**Fix**: Store as single `Vec<f32>` with `sample_boundaries` for indexing. For encode: 269,800 rows each with 3 sample vecs = ~1M fewer allocations.
**Risk**: Touches all row value readers (writers, sort, prune). Must preserve sample-major order.

### 4. sort_groups clone elimination (MEDIUM impact, LOW risk)
**File**: `src/pipeline/matrix.rs:309`
**Problem**: `let mut reordered = self.rows.clone()` clones all 269,800 rows (deep clone of all value vecs).
**Fix**: Build the reordered Vec by draining/moving from original, or use indices-based reorder.
**Memory saving**: Avoids duplicating entire matrix during sort (~100MB for encode).

### 5. write_matrix_value micro-optimizations (LOW impact, LOW risk)
**File**: `src/io/writers/mod.rs:489-538`
**Problem**: `round_ties_even` for 6-decimal formatting — Python uses standard rounding.
**Fix**: Use `libm::round` or standard rounding. Pre-size ROW_BUFFER more aggressively.
**Speed**: Small constant factor improvement in output path.

## Safety Rules

1. **Parity tests must pass after every change** — run `pixi run python scripts/custom_compare.py --mode python-compatibility`
2. **Encode benchmark values must match ≤ 5e-6** — run encode tests with `--no-cache` flag
3. **No overfitting** — never optimize for encode dataset specifically. All changes must be general.
4. **Memory constraint** — total RAM available ~8GB. Keep max RSS under 2GB.
5. **Release builds** — always `cargo build --release` before benchmark.
6. **Atomic commits** — one commit per optimization with before/after metrics.

## Benchmark Commands

```bash
# Parity tests (fast, run after every change)
pixi run python scripts/custom_compare.py --mode python-compatibility

# Encode reference-point (skip reference: --keep-ref)
pixi run python scripts/custom_compare.py --mode reference-point \
  --reference-point center --upstream 100 --downstream 100 --bin-size 10 \
  --keep-ref --no-cache

# Encode scale-regions
pixi run python scripts/custom_compare.py --mode scale-regions \
  --region-body-length 200 --upstream 100 --downstream 100 \
  --unscaled-5-prime 50 --unscaled-3-prime 50 --bin-size 10 \
  --keep-ref --no-cache

# Memory measurement
/usr/bin/time -v ./target/release/compute_matrix_rs <mode> [args] -o /tmp/test.mat.gz
```

## Iterator Order

For maximum impact with minimum risk:
1. aggregate_slice single-pass
2. Coverage buffer reuse
3. sort_groups clone elimination
4. MatrixRow value flattening
5. write_matrix_value micro-optimizations
