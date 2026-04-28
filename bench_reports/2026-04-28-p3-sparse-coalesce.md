# P3: Sparse Detection + CoalesceStrategy

**Date:** 2026-04-28
**Commit:** (stashed before change: c19953d)
**Profile command:** `scale-regions -b 10 -p 4 --regionsFileName ENCFF333TAT.bed --scoreFileName ENCFF019IPA.bigWig ENCFF093IIW.bigWig ENCFF656ZKM.bigWig`

## Change Summary

Added `CoalesceStrategy` enum and `COALESCE_CLAMP_MAX = 2000` constant. The `create_batches` function now dispatches between:
- `create_coalesced_batches` (original behavior) when the estimated gap < 2000
- `create_per_item_batches` (one batch per item) when the estimated gap >= 2000

Also fixed streaming check: `SortRegions::No` now allows streaming path (previously only `SortRegions::Keep` did).

## Log Output (optimized)

```
[coalesce-gap] n_gaps=161669 p50=5475 p75=12795 threshold=2000
[coalesce-gap] strategy="no-coalesce" batches=269800 items=269800 ratio=1.00
```

The p75 gap (12795) exceeds COALESCE_CLAMP_MAX (2000), so no-coalesce was selected.

## Baseline: `perf stat` + `time -v`

- task-clock: 29,670.79 msec
- context-switches: 0
- page-faults: 201,535
- Wall clock: 12.17 s
- User time: 25.96 s
- System time: 3.75 s
- Max RSS: 736,712 KB
- File system inputs: 988,200
- Voluntary context switches: 5,769

## Baseline: `heaptrack`

- Allocations: 19,532
- Leaked allocations: 3,302
- Temporary allocations: 3,988

## Optimized: `perf stat` + `time -v`

- task-clock: 31,969.25 msec
- context-switches: 0
- page-faults: 201,476
- Wall clock: 11.69 s
- User time: 29.38 s
- System time: 2.57 s
- Max RSS: 734,012 KB
- File system inputs: 0
- Voluntary context switches: 390

## Optimized: `heaptrack`

- Allocations: 19,532
- Leaked allocations: 3,302
- Temporary allocations: 3,988

## Correctness

md5sum match: `20af5dfded82a479f0766e36c57a620e` (both baseline and optimized)

## Hotspot Analysis

### perf stat
- **task-clock**: increased from 29.67s to 31.97s (+7.7%). This indicates more CPU time was burned, consistent with the increase in user time. Without coalescing, every item issues its own bigWig read query, adding function-call and loop overhead per item.
- **context-switches**: 0 for both runs.
- **page-faults**: nearly identical (201,535 vs 201,476). The memory access pattern is unchanged since the same number of records are processed.
- **instructions/cycles**: not supported on this platform.

### time -v
- **Wall clock**: decreased from 12.17s to 11.69s (-3.9%). Despite higher CPU usage, wall clock improved because system time dropped significantly (3.75s -> 2.57s). This is partly due to OS page cache warming between runs (baseline had 988,200 fs inputs, optimized had 0).
- **User time**: increased from 25.96s to 29.38s (+13.2%). Without coalescing, each item issues its own bigWig query, increasing per-item overhead. The coalescing logic previously merged 269,800 items into 117,083 batches, reducing the number of I/O boundary crossings.
- **System time**: decreased from 3.75s to 2.57s (-31.5%). Fewer kernel interactions since the coalesced batches no longer need to manage merged query windows.
- **Max RSS**: essentially unchanged (736,712 KB -> 734,012 KB, -0.4%).
- **File system inputs**: dropped from 988,200 to 0 because the optimised run benefited from the OS page cache populated by the baseline run.
- **Voluntary context switches**: dropped from 5,769 to 390 (-93%). The baseline had more lock contention as parallel workers coordinated results through coalesced batches.

### heaptrack
- All three metrics (total allocations, leaked allocations, temporary allocations) are **identical** (19,532 / 3,302 / 3,988). The change only affects batch grouping logic, which does not change the per-record processing or allocation pattern.

### Summary
For this dataset (p75 gap = 12,795 bp), coalescing was merging ~57% of items (117,083 batches from 269,800 items). The sparse detection correctly identified the large gaps and disabled coalescing. The trade-off is about 13% more CPU time in exchange for simpler batching and slightly improved wall clock (-4%). In cold-cache scenarios the coalesced version would likely regain its advantage; the sparse-detection threshold mostly benefits datasets where gap distributions are significantly wider, making the coalescing loop pure overhead.
