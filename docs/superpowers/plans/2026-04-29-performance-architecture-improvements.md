# Performance Architecture Improvements Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Improve bigWig cache correctness/fairness, add a repeatable performance smoke gate, and prepare focused maintainability/performance work without changing user-visible matrix semantics.

**Architecture:** Keep the current two-path executor architecture (`StreamOrdered` and `HybridBucket`) intact. Replace the single global bigWig block cache with per-bigWig caches whose capacities are statically split from one total entry cap, then add a smoke benchmark, extract aggregation helpers, prototype guarded direct mean/sum aggregation, and split spill serialization from collector logic.

**Tech Stack:** Rust 2024, existing `quick_cache`, `flate2`/`zlib-rs`, `rayon`, `memmap2`, existing Pixi-based Python/deepTools regression harness, existing Rust unit/integration tests.

---

## Scope

This plan intentionally does **not** include the original repeated `plan_for` optimization, row value buffer reuse, or runtime tuning centralization tasks. Those were removed after discussion.

The remaining work is split into independent tasks:

1. Per-bigWig block caches with a strict total entry cap.
2. Local performance smoke script.
3. Aggregation helper extraction from `worker.rs`.
4. Guarded direct interval-to-bin aggregation for `Mean` and `Sum`. (**Depends on Task 3**: `aggregation.rs` must exist first.)
5. Spill serialization format extraction from `spill.rs`.

Each task should be implemented, tested, and committed separately.

---

## File Map

### Task 1: Per-bigWig block cache with total cap

- Modify: `src/io/readers/block_cache.rs`
  - Replace the single fixed-capacity shared cache API with a capacity-aware per-file block cache.
  - Add deterministic capacity splitting helpers.
  - Support capacity `0` as no-cache mode.
- Modify: `src/io/readers/bwig.rs`
  - Give each `SharedBigWigReader` its own block cache.
  - Keep `SharedBigWigReader::open()` working with the default total capacity for single-file use.
  - Add/open a capacity-aware constructor used by the executor.
- Modify: `src/pipeline/core/executor.rs`
  - Compute each sample's cache capacity from `sample_paths.len()`.
  - Open each bigWig with its assigned per-file cache capacity.
- Modify: `plans/overall_status.md`
  - Record the cache isolation/fairness change after implementation.

### Task 2: Performance smoke script

- Create: `scripts/perf_smoke.sh`
  - Run fixed local reference-point and scale-regions examples.
  - Record elapsed/user/sys/max RSS.
- Modify: `readme.md`
  - Document the smoke command.
- Modify: `plans/overall_status.md`
  - Record the new smoke benchmark entry.

### Task 3: Aggregation helper extraction

- Create: `src/pipeline/core/aggregation.rs`
  - Move `aggregate_slice`, `index_from_coordinate`, and their tests out of `worker.rs`.
- Modify: `src/pipeline/core/mod.rs`
  - Register the new module.
- Modify: `src/pipeline/core/worker.rs`
  - Import aggregation helpers from the new module.

### Task 4: Direct mean/sum interval-to-bin aggregation

- Modify: `src/pipeline/core/aggregation.rs`
  - Add direct aggregation helpers and tests.
- Modify: `src/pipeline/core/worker.rs`
  - Use the direct path only under explicit semantic guards.

### Task 5: Spill serialization split

- Create: `src/pipeline/core/spill_format.rs`
  - Move spill row serialization/deserialization and chromosome table code.
- Modify: `src/pipeline/core/mod.rs`
  - Register the new module.
- Modify: `src/pipeline/core/spill.rs`
  - Import moved format helpers and keep collector/finalization logic.

---

### Task 1: Per-bigWig Block Cache With Total Cap

**Files:**

- Modify: `src/io/readers/block_cache.rs`
- Modify: `src/io/readers/bwig.rs`
- Modify: `src/pipeline/core/executor.rs`
- Modify: `plans/overall_status.md`

**Design:**

Use one independent block cache per bigWig file. Split the existing total cap of `500` entries across all sample files:

```text
total_cap = 500
sample_count = N
base = total_cap / N
remainder = total_cap % N

sample_index < remainder -> base + 1 entries
sample_index >= remainder -> base entries
```

Examples:

```text
1 file  -> [500]
2 files -> [250, 250]
3 files -> [167, 167, 166]
8 files -> [63, 63, 63, 63, 62, 62, 62, 62]
800 files with total_cap=500 -> first 500 files get 1 entry, remaining 300 files get 0 entries
```

Capacity `0` means the bigWig reader does not cache blocks. `get()` always returns `None`; `insert()` is a no-op.

- [ ] **Step 1: Add failing capacity split and zero-cap tests**

Add these tests to `src/io/readers/block_cache.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    #[test]
    fn split_cache_capacity_divides_total_with_remainder() {
        let caps: Vec<_> = (0..3)
            .map(|sample_index| split_block_cache_capacity(500, 3, sample_index))
            .collect();
        assert_eq!(caps, vec![167, 167, 166]);
        assert_eq!(caps.iter().sum::<usize>(), 500);
    }

    #[test]
    fn split_cache_capacity_allows_zero_capacity_when_samples_exceed_total() {
        assert_eq!(split_block_cache_capacity(3, 5, 0), 1);
        assert_eq!(split_block_cache_capacity(3, 5, 1), 1);
        assert_eq!(split_block_cache_capacity(3, 5, 2), 1);
        assert_eq!(split_block_cache_capacity(3, 5, 3), 0);
        assert_eq!(split_block_cache_capacity(3, 5, 4), 0);
    }

    #[test]
    fn zero_capacity_cache_never_stores_entries() {
        let cache = SharedBlockCache::with_capacity(0);
        let key = (128, 64);
        cache.insert(key, Arc::from(&b"payload"[..]));
        assert!(cache.get(&key).is_none());
    }

    #[test]
    fn independent_per_file_caches_do_not_share_offsets() {
        let cache_a = SharedBlockCache::with_capacity(2);
        let cache_b = SharedBlockCache::with_capacity(2);
        let key = (128, 64);

        cache_a.insert(key, Arc::from(&b"file-a"[..]));
        cache_b.insert(key, Arc::from(&b"file-b"[..]));

        assert_eq!(cache_a.get(&key).unwrap().as_ref(), b"file-a");
        assert_eq!(cache_b.get(&key).unwrap().as_ref(), b"file-b");
    }
}
```

- [ ] **Step 2: Run focused tests and verify they fail**

Run:

```bash
cargo test io::readers::block_cache
```

Expected: FAIL because `SharedBlockCache::with_capacity` and `split_block_cache_capacity` do not exist yet.

- [ ] **Step 3: Implement capacity-aware per-file cache**

Rename the existing constant `MAX_SHARED_BLOCK_CACHE_ENTRIES` → `DEFAULT_TOTAL_BLOCK_CACHE_ENTRIES` to reflect its new role as a configurable total budget rather than a fixed per-cache maximum. Replace `src/io/readers/block_cache.rs` with this shape:

```rust
use std::sync::Arc;

use quick_cache::sync::Cache;

pub(crate) const DEFAULT_TOTAL_BLOCK_CACHE_ENTRIES: usize = 500;

type BlockCacheKey = (u64, u64);

pub fn split_block_cache_capacity(
    total_capacity: usize,
    sample_count: usize,
    sample_index: usize,
) -> usize {
    if sample_count == 0 || sample_index >= sample_count {
        return 0;
    }

    let base = total_capacity / sample_count;
    let remainder = total_capacity % sample_count;
    if sample_index < remainder {
        base + 1
    } else {
        base
    }
}

pub struct SharedBlockCache {
    cache: Option<Cache<BlockCacheKey, Arc<[u8]>>>,
}

impl SharedBlockCache {
    pub fn new() -> Self {
        Self::with_capacity(DEFAULT_TOTAL_BLOCK_CACHE_ENTRIES)
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            cache: (capacity > 0).then(|| Cache::new(capacity)),
        }
    }

    pub fn get(&self, key: &BlockCacheKey) -> Option<Arc<[u8]>> {
        self.cache.as_ref().and_then(|cache| cache.get(key))
    }

    pub fn insert(&self, key: BlockCacheKey, value: Arc<[u8]>) {
        if let Some(cache) = &self.cache {
            cache.insert(key, value);
        }
    }
}
```

- [ ] **Step 4: Run focused cache tests**

Run:

```bash
cargo test io::readers::block_cache
```

Expected: PASS.

- [ ] **Step 5: Add capacity-aware bigWig open constructor**

Modify `src/io/readers/bwig.rs` so `SharedBigWigReader` still stores `Arc<SharedBlockCache>`, but each reader gets its own cache instance.

Keep `open()` as single-file default:

```rust
pub fn open(path: impl AsRef<Path>) -> Result<Self, BigWigReadError> {
    Self::open_with_block_cache_capacity(
        path,
        super::block_cache::DEFAULT_TOTAL_BLOCK_CACHE_ENTRIES,
    )
}
```

Add a capacity-aware constructor:

```rust
pub fn open_with_block_cache_capacity(
    path: impl AsRef<Path>,
    block_cache_capacity: usize,
) -> Result<Self, BigWigReadError> {
    Self::open_with_cache(
        path,
        Arc::new(SharedBlockCache::with_capacity(block_cache_capacity)),
    )
}
```

Keep `open_with_cache(path, block_cache)` for tests and internal callers that need explicit cache injection. Do not add file IDs to the key because caches are now per-file.

- [ ] **Step 6: Assign per-sample capacities in executor**

Modify `src/pipeline/core/executor.rs` imports:

```rust
use crate::io::readers::block_cache::{
    DEFAULT_TOTAL_BLOCK_CACHE_ENTRIES,
    split_block_cache_capacity,
};
```

Then replace the shared cache construction/opening block with per-file cache capacity assignment:

```rust
let sample_count_for_cache = sample_paths.len();
let shared_readers = Arc::new(
    sample_paths
        .iter()
        .enumerate()
        .map(|(sample_index, path)| {
            let cache_capacity = split_block_cache_capacity(
                DEFAULT_TOTAL_BLOCK_CACHE_ENTRIES,
                sample_count_for_cache,
                sample_index,
            );
            SharedBigWigReader::open_with_block_cache_capacity(path, cache_capacity)
                .map(Arc::new)
                .with_context(|| {
                    format!("Failed to open bigWig file '{}'", path.display())
                })
        })
        .collect::<Result<Vec<_>>>()?,
);
```

Remove the old executor-local `Arc::new(SharedBlockCache::new())` shared across all sample files.

- [ ] **Step 7: Run cache and full tests**

Run:

```bash
cargo test io::readers::block_cache
cargo test
pixi run python scripts/custom_compare.py --mode python-compatibility
```

Expected:

- `cargo test io::readers::block_cache`: PASS
- `cargo test`: PASS
- Python compatibility: 10/10 scenarios pass within ≤5e-6 tolerance

- [ ] **Step 8: Update project status**

Add this entry to `plans/overall_status.md` under recent fixes or v0.3.0 changes:

```markdown
- ✅ **2026-04-29**: Changed bigWig block caching to per-file caches with a strict total entry cap, preventing cross-sample block collisions while keeping multi-sample cache usage bounded.
```

- [ ] **Step 9: Commit**

```bash
git add src/io/readers/block_cache.rs src/io/readers/bwig.rs src/pipeline/core/executor.rs plans/overall_status.md
git commit -m "perf: isolate bigwig block caches per file with total cap"
```

---

### Task 2: Add a Repeatable Performance Smoke Script

> **Note:** The existing `scripts/profile_bench.sh` is a full profiling harness (requires `perf`, `heaptrack`, writes detailed reports). This new script is a lightweight smoke test — no profiling tools required, just `/usr/bin/time` — intended for quick regression checks during development.

**Files:**

- Create: `scripts/perf_smoke.sh`
- Modify: `readme.md`
- Modify: `plans/overall_status.md`

- [ ] **Step 1: Write the smoke script**

Create `scripts/perf_smoke.sh`:

```bash
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
```

- [ ] **Step 2: Make the script executable**

Run:

```bash
chmod +x scripts/perf_smoke.sh
```

- [ ] **Step 3: Run the script**

Run:

```bash
scripts/perf_smoke.sh
```

Expected: both cases complete and write timing logs under `target/perf-smoke/`.

- [ ] **Step 4: Document the command**

Add this to `readme.md` under Quick Start:

```markdown
- Local performance smoke run:
  `scripts/perf_smoke.sh`
```

- [ ] **Step 5: Update project status**

Add this to `plans/overall_status.md`:

```markdown
- ✅ **2026-04-29**: Added `scripts/perf_smoke.sh` for repeatable local performance smoke timing of reference-point and scale-regions cases.
```

- [ ] **Step 6: Commit**

```bash
git add scripts/perf_smoke.sh readme.md plans/overall_status.md
git commit -m "test: add performance smoke script"
```

---

### Task 3: Split Aggregation Helpers Out of `worker.rs`

**Files:**

- Create: `src/pipeline/core/aggregation.rs`
- Modify: `src/pipeline/core/mod.rs`
- Modify: `src/pipeline/core/worker.rs`

- [ ] **Step 1: Create aggregation module with moved helper code**

Create `src/pipeline/core/aggregation.rs` with the current `index_from_coordinate` and `aggregate_slice` implementations from `src/pipeline/core/worker.rs`:

```rust
use crate::config::AverageTypeBins;

pub(super) fn index_from_coordinate(value: i64, base: i64, window_len: usize) -> usize {
    if value <= base {
        return 0;
    }
    let diff = value - base;
    let idx = usize::try_from(diff).unwrap_or(window_len);
    idx.min(window_len)
}

pub(super) fn aggregate_slice(slice: &[f64], average_type: AverageTypeBins) -> Option<f64> {
    let len = slice.len();
    if len == 0 {
        return None;
    }

    match average_type {
        AverageTypeBins::Mean => {
            let mut sum = 0.0f64;
            let mut count = 0u32;
            for &value in slice {
                if !value.is_nan() {
                    sum += value;
                    count += 1;
                }
            }
            if count == 0 { None } else { Some(sum / count as f64) }
        }
        AverageTypeBins::Sum => {
            let mut sum = 0.0f64;
            let mut found = false;
            for &value in slice {
                if !value.is_nan() {
                    sum += value;
                    found = true;
                }
            }
            if found { Some(sum) } else { None }
        }
        AverageTypeBins::Min => {
            let mut min = f64::INFINITY;
            let mut found = false;
            for &value in slice {
                if !value.is_nan() {
                    min = min.min(value);
                    found = true;
                }
            }
            if found { Some(min) } else { None }
        }
        AverageTypeBins::Max => {
            let mut max = f64::NEG_INFINITY;
            let mut found = false;
            for &value in slice {
                if !value.is_nan() {
                    max = max.max(value);
                    found = true;
                }
            }
            if found { Some(max) } else { None }
        }
        AverageTypeBins::Std => {
            let mut sum = 0.0f64;
            let mut count = 0u32;
            for &value in slice {
                if !value.is_nan() {
                    sum += value;
                    count += 1;
                }
            }
            if count == 0 {
                return None;
            }
            let mean = sum / count as f64;
            let mut variance_sum = 0.0f64;
            for &value in slice {
                if !value.is_nan() {
                    let delta = value - mean;
                    variance_sum += delta * delta;
                }
            }
            Some((variance_sum / count as f64).sqrt())
        }
        AverageTypeBins::Median => {
            let mut values: Vec<f64> = slice
                .iter()
                .copied()
                .filter(|v| !v.is_nan())
                .collect();
            if values.is_empty() {
                return None;
            }
            values.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let mid = values.len() / 2;
            if values.len() % 2 == 0 {
                Some((values[mid - 1] + values[mid]) / 2.0)
            } else {
                Some(values[mid])
            }
        }
    }
}
```

Move the existing aggregation and index unit tests from `src/pipeline/core/worker.rs` into this new module without changing expected values.

- [ ] **Step 2: Register the module**

Modify `src/pipeline/core/mod.rs`:

```rust
mod aggregation;
```

- [ ] **Step 3: Update worker imports and remove moved functions**

Modify `src/pipeline/core/worker.rs`:

```rust
use super::aggregation::{aggregate_slice, index_from_coordinate};
```

Remove the original `aggregate_slice` and `index_from_coordinate` definitions from `worker.rs` after the import compiles.

- [ ] **Step 4: Run focused tests**

Run:

```bash
cargo test pipeline::core::aggregation
cargo test pipeline::core::worker
```

Expected: PASS.

- [ ] **Step 5: Run full tests**

Run:

```bash
cargo test
pixi run python scripts/custom_compare.py --mode python-compatibility
```

Expected: all tests pass and compatibility remains within ≤5e-6.

- [ ] **Step 6: Commit**

```bash
git add src/pipeline/core/aggregation.rs src/pipeline/core/mod.rs src/pipeline/core/worker.rs
git commit -m "refactor: extract aggregation helpers from worker"
```

---

### Task 4: Prototype Direct Interval-to-Bin Aggregation for Mean and Sum

**Files:**

- Modify: `src/pipeline/core/aggregation.rs`
- Modify: `src/pipeline/core/worker.rs`

**Design:**

The direct path must preserve current coverage-buffer semantics. In particular, when `missing_data_as_zero` is enabled, uncovered bases count as valid zeroes for `Mean`, because the existing buffer path fills uncovered positions with `0.0` before calling `aggregate_slice`.

Do not enable the direct path for `Median`, `Std`, `Min`, `Max`, or metagene/intron-masked rows in this task.

- [ ] **Step 1: Add direct aggregation tests**

Add these tests to `src/pipeline/core/aggregation.rs`:

```rust
#[test]
fn direct_mean_weights_interval_overlap_by_base_count() {
    let bins = vec![(0, 10), (10, 20)];
    let intervals = vec![(0, 5, 2.0), (5, 20, 4.0)];
    let values = direct_mean_bins(&bins, &intervals, false);
    assert_eq!(values, vec![Some(3.0), Some(4.0)]);
}

#[test]
fn direct_mean_counts_uncovered_bases_as_zero_when_missing_data_as_zero() {
    let bins = vec![(0, 10)];
    let intervals = vec![(0, 5, 2.0)];
    let values = direct_mean_bins(&bins, &intervals, true);
    assert_eq!(values, vec![Some(1.0)]);
}

#[test]
fn direct_mean_ignores_uncovered_bases_when_missing_data_stays_nan() {
    let bins = vec![(0, 10)];
    let intervals = vec![(0, 5, 2.0)];
    let values = direct_mean_bins(&bins, &intervals, false);
    assert_eq!(values, vec![Some(2.0)]);
}

#[test]
fn direct_mean_returns_none_for_fully_uncovered_nan_bin() {
    let bins = vec![(0, 10)];
    let intervals = Vec::new();
    let values = direct_mean_bins(&bins, &intervals, false);
    assert_eq!(values, vec![None]);
}

#[test]
fn direct_sum_weights_interval_overlap_by_base_count() {
    let bins = vec![(0, 10), (10, 20)];
    let intervals = vec![(0, 5, 2.0), (5, 20, 4.0)];
    let values = direct_sum_bins(&bins, &intervals, false);
    assert_eq!(values, vec![Some(30.0), Some(40.0)]);
}

#[test]
fn direct_sum_returns_zero_for_uncovered_bin_when_missing_data_as_zero() {
    let bins = vec![(0, 10)];
    let intervals = Vec::new();
    let values = direct_sum_bins(&bins, &intervals, true);
    assert_eq!(values, vec![Some(0.0)]);
}
```

- [ ] **Step 2: Run tests and verify they fail**

Run:

```bash
cargo test pipeline::core::aggregation::direct_
```

Expected: FAIL because `direct_mean_bins` and `direct_sum_bins` do not exist yet.

- [ ] **Step 3: Implement direct mean/sum helpers**

Add helpers to `src/pipeline/core/aggregation.rs`:

```rust
pub(super) fn direct_mean_bins(
    bins: &[(i64, i64)],
    intervals: &[(i64, i64, f64)],
    missing_data_as_zero: bool,
) -> Vec<Option<f64>> {
    bins.iter()
        .map(|&(bin_start, bin_end)| {
            if bin_start >= bin_end {
                return if missing_data_as_zero { Some(0.0) } else { None };
            }

            let mut weighted_sum = 0.0;
            let mut covered = 0_i64;
            for &(start, end, value) in intervals {
                if value.is_nan() {
                    continue;
                }
                let overlap_start = start.max(bin_start);
                let overlap_end = end.min(bin_end);
                if overlap_start < overlap_end {
                    let width = overlap_end - overlap_start;
                    weighted_sum += value * width as f64;
                    covered += width;
                }
            }

            if missing_data_as_zero {
                Some(weighted_sum / (bin_end - bin_start) as f64)
            } else if covered == 0 {
                None
            } else {
                Some(weighted_sum / covered as f64)
            }
        })
        .collect()
}

pub(super) fn direct_sum_bins(
    bins: &[(i64, i64)],
    intervals: &[(i64, i64, f64)],
    missing_data_as_zero: bool,
) -> Vec<Option<f64>> {
    bins.iter()
        .map(|&(bin_start, bin_end)| {
            let mut sum = 0.0;
            let mut found = false;
            for &(start, end, value) in intervals {
                if value.is_nan() {
                    continue;
                }
                let overlap_start = start.max(bin_start);
                let overlap_end = end.min(bin_end);
                if overlap_start < overlap_end {
                    sum += value * (overlap_end - overlap_start) as f64;
                    found = true;
                }
            }

            if found || missing_data_as_zero {
                Some(sum)
            } else {
                None
            }
        })
        .collect()
}
```

- [ ] **Step 4: Run direct aggregation tests**

Run:

```bash
cargo test pipeline::core::aggregation::direct_
```

Expected: PASS.

- [ ] **Step 5: Integrate only behind exact semantic guard**

In `src/pipeline/core/worker.rs`, use direct aggregation only when all of these are true:

- `general.average_type_bins` is `AverageTypeBins::Mean` or `AverageTypeBins::Sum`.
- `plan.included_intervals().is_none()`.
- The batch path can reuse intervals already fetched for the merged query window without an extra bigWig read.
- Unit tests and Python compatibility prove parity.

**Interval data source:** The `process_batch` function (worker.rs:355) already fetches `BigWigValue` intervals per sample via `sample.reader_mut().values(chrom, fetch_start, fetch_end)` (worker.rs:414-424). Each `BigWigValue` has `start: u32`, `end: u32`, `value: f32` (bwig.rs:35-39). Currently these are expanded into a per-base coverage buffer (`cov[rs..re].fill(value)`) then sliced per bin via `aggregate_slice`. The direct path should instead convert the fetched `&[BigWigValue]` to `Vec<(i64, i64, f64)>` and pass them to `direct_mean_bins`/`direct_sum_bins`, skipping the coverage buffer entirely. This avoids the O(window_span) fill and gives O(intervals × bins) complexity instead. The conversion is:
```rust
let intervals: Vec<(i64, i64, f64)> = raw_intervals
    .iter()
    .map(|v| (i64::from(v.start), i64::from(v.end), f64::from(v.value)))
    .collect();
```
The direct path applies only in `process_batch`'s per-bin loop (worker.rs:460-494), **not** in the single-item `compute_sample_bins` path. The single-item path already has per-base semantics that are harder to shortcut safely.

For every direct helper result, apply the same post-processing currently applied after `aggregate_slice`. Note: the `missing_data_as_zero` → `Some(0.0)` fallback is already handled inside `direct_mean_bins`/`direct_sum_bins` (they take it as a parameter), so do **not** duplicate it in the post-processing:

```rust
// missing_data_as_zero is already handled inside direct_*_bins — do NOT add it here
let mut value = value_option.unwrap_or(f64::NAN);

if nan_after_end && bin.beyond_region() {
    value = f64::NAN;
}

if value.is_finite() {
    value *= general.scale_factor;
}
```

- [ ] **Step 6: Run full compatibility and perf smoke**

Run:

```bash
cargo test
pixi run python scripts/custom_compare.py --mode python-compatibility
scripts/perf_smoke.sh
```

Expected:

- Rust tests pass.
- Python compatibility remains 10/10 within ≤5e-6.
- Perf smoke shows no regression.

- [ ] **Step 7: Commit**

```bash
git add src/pipeline/core/aggregation.rs src/pipeline/core/worker.rs
git commit -m "perf: add direct mean and sum bin aggregation"
```

---

### Task 5: Split Spill Serialization from Collector Logic

**Files:**

- Create: `src/pipeline/core/spill_format.rs`
- Modify: `src/pipeline/core/mod.rs`
- Modify: `src/pipeline/core/spill.rs`

- [ ] **Step 1: Move format-only code**

Move these items from `src/pipeline/core/spill.rs` into `src/pipeline/core/spill_format.rs`:

- `SpillIndex` struct (spill.rs:30-44)
- `ChromTable` struct and its `impl` block (spill.rs:53-78)
- Flag constants: `FLAG_HAS_NAME`, `FLAG_HAS_SCORE_RAW`, `FLAG_HAS_STRAND_RAW`, `FLAG_HAS_EXON_COORDS` (spill.rs:84-87)
- `SCORE_NONE_SENTINEL` constant (spill.rs:156)
- `write_len_prefixed_str` (spill.rs:94-99)
- `read_len_prefixed_str` (spill.rs:102-112)
- `read_u16` (spill.rs:114-121)
- `read_u32` (spill.rs:123-135)
- `read_f64` (spill.rs:137-153)
- `serialize_row` (spill.rs:166-~)
- `deserialize_row` (spill.rs:~-374)

Keep collector, bucket, flush lifecycle, mmap finalization, sorting, and keep-order placement in `spill.rs`.

- [ ] **Step 2: Register the module**

Modify `src/pipeline/core/mod.rs`:

```rust
pub(crate) mod spill_format;
```

- [ ] **Step 3: Update imports in `spill.rs`**

Import moved items from `spill_format`. The exact import list should match the moved item names; the target shape is:

```rust
use super::spill_format::{
    deserialize_row, serialize_row, ChromTable, SpillIndex,
    FLAG_HAS_NAME, FLAG_HAS_SCORE_RAW, FLAG_HAS_STRAND_RAW, FLAG_HAS_EXON_COORDS,
};
```

Note: the `FLAG_*` constants and helper functions (`write_len_prefixed_str`, etc.) only need to be re-imported if `spill.rs` still references them directly. If they are exclusively used inside `serialize_row`/`deserialize_row`, they can remain private to `spill_format.rs` and no re-import is needed. Check after the move and trim imports accordingly.

- [ ] **Step 4: Keep serialization tests attached to the format module**

Move serialization/deserialization-focused tests with the format code. Keep integration tests that exercise `HybridBucketCollector` in `spill.rs`.

- [ ] **Step 5: Run spill tests**

Run:

```bash
cargo test pipeline::core::spill
cargo test pipeline::core::spill_format
```

Expected: PASS with no behavior changes.

- [ ] **Step 6: Run full tests**

Run:

```bash
cargo test
pixi run python scripts/custom_compare.py --mode python-compatibility
```

Expected: all tests pass and compatibility remains within ≤5e-6.

- [ ] **Step 7: Commit**

```bash
git add src/pipeline/core/spill_format.rs src/pipeline/core/mod.rs src/pipeline/core/spill.rs
git commit -m "refactor: extract spill serialization format"
```

---

## Verification Checklist

Before claiming the full plan complete, run:

```bash
cargo fmt --check
cargo test
pixi run python scripts/custom_compare.py --mode python-compatibility
scripts/perf_smoke.sh
```

Expected final state:

- Formatting passes.
- Rust test suite passes.
- Python compatibility suite reports 10/10 passing within ≤5e-6 tolerance.
- Performance smoke script completes both reference-point and scale-regions cases.
- BigWig block caches are isolated per file and the total configured entry count remains bounded.
- No user-visible CLI or matrix output semantics change.

---

## Self-Review

Spec coverage:

- Per-bigWig cache isolation and strict total cap are covered by Task 1.
- Benchmark repeatability is covered by Task 2.
- Original Task 5, aggregation helper extraction, is retained as Task 3.
- Direct mean/sum aggregation is covered by Task 4 and guarded to avoid semantic drift.
- Spill serialization maintainability is covered by Task 5.
- Original repeated planning, row buffer reuse, and runtime tuning tasks were intentionally removed.

Placeholder scan:

- The plan contains no `TBD` or `TODO` markers.
- Deferred work is explicitly out of scope, not left as incomplete implementation text.

Type consistency:

- `SharedBlockCache::with_capacity` and `split_block_cache_capacity` are introduced before executor integration.
- `aggregation.rs` is introduced before direct aggregation helpers are added.
- `direct_mean_bins` and `direct_sum_bins` return `Vec<Option<f64>>`, matching uncovered-bin semantics before conversion to `f64::NAN` in `worker.rs`.
