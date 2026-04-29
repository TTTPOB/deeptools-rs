# BigWig Reader Memory & Cache Optimization Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reduce per-worker memory waste and improve cache efficiency in the bigWig reader layer through four targeted changes: eliminate redundant chrom_lengths HashMap, share CIR node cache across workers, share decompression buffers within each worker, and implement adaptive block cache sizing.

**Architecture:** The bigWig reader has a two-tier design: `SharedBigWigReader` (one per file, shared via `Arc`) holds immutable state, while `BigWigReader` (one per worker per file) holds mutable working state. This plan moves immutable data up to the shared tier, consolidates mutable buffers that are used serially within a worker, and replaces the fixed 500-entry block cache budget with a per-file default of 200 entries capped at a 2000-entry global hard limit.

**Tech Stack:** Rust, `quick_cache::sync::Cache`, `flate2`, `rayon`

---

## File Map

| File | Changes |
|------|---------|
| `src/io/readers/bwig.rs` | Remove `chrom_lengths` duplication by adding `find_chrom_length()` to `SharedBigWigReader`; move `cir_node_cache` from `BigWigReader` to `SharedBigWigReader`; extract `work_buf`/`decode_buf` from `BigWigReader` into external `&mut` params; add `values_with_bufs()` method |
| `src/io/readers/block_cache.rs` | Replace `split_block_cache_capacity()` with `compute_per_file_block_cache_capacity()` implementing the 200-default / 2000-hard-limit policy |
| `src/pipeline/core/samples.rs` | Remove `chrom_lengths: HashMap<String, u32>` field from `Sample`; delegate `chrom_length()` to `SharedBigWigReader::find_chrom_length()` |
| `src/pipeline/core/executor.rs` | Update block cache capacity computation call site |

---

### Task 1: Eliminate chrom_lengths HashMap in Sample

Each `Sample` clones a full `HashMap<String, u32>` from `SharedBigWigReader.chroms()` on creation. This is unnecessary — `SharedBigWigReader` already has `chroms: Vec<ChromInfo>` sorted by name, which supports O(log n) binary search.

**Files:**
- Modify: `src/io/readers/bwig.rs:109-194` (add `find_chrom_length` to `SharedBigWigReader`)
- Modify: `src/pipeline/core/samples.rs:1-63` (remove HashMap, delegate lookup)
- Test: existing tests + new unit test in `bwig.rs`

- [ ] **Step 1: Add `find_chrom_length` method to `SharedBigWigReader`**

In `src/io/readers/bwig.rs`, add after the existing `find_chrom_id` method (line 194):

```rust
pub fn find_chrom_length(&self, name: &str) -> Option<u32> {
    self.chroms
        .binary_search_by(|c| c.name.as_str().cmp(name))
        .ok()
        .map(|idx| self.chroms[idx].length)
}
```

- [ ] **Step 2: Extract binary search logic into a testable free function**

The `find_chrom_length` method delegates to the chroms `Vec` binary search. To make this directly testable without constructing a full `SharedBigWigReader` (which requires a real file), extract the lookup as a free function in `src/io/readers/bwig.rs`:

```rust
fn binary_search_chrom_length(chroms: &[ChromInfo], name: &str) -> Option<u32> {
    chroms
        .binary_search_by(|c| c.name.as_str().cmp(name))
        .ok()
        .map(|idx| chroms[idx].length)
}
```

Then have `find_chrom_length` delegate to it:

```rust
pub fn find_chrom_length(&self, name: &str) -> Option<u32> {
    binary_search_chrom_length(&self.chroms, name)
}
```

- [ ] **Step 3: Write a unit test that calls the actual helper**

In `src/io/readers/bwig.rs` inside `mod tests`, add:

```rust
#[test]
fn binary_search_chrom_length_found() {
    let chroms = vec![
        ChromInfo { name: "chr1".to_string(), length: 1000 },
        ChromInfo { name: "chr2".to_string(), length: 2000 },
        ChromInfo { name: "chrX".to_string(), length: 3000 },
    ];
    assert_eq!(binary_search_chrom_length(&chroms, "chr1"), Some(1000));
    assert_eq!(binary_search_chrom_length(&chroms, "chr2"), Some(2000));
    assert_eq!(binary_search_chrom_length(&chroms, "chrX"), Some(3000));
}

#[test]
fn binary_search_chrom_length_not_found() {
    let chroms = vec![
        ChromInfo { name: "chr1".to_string(), length: 1000 },
        ChromInfo { name: "chr2".to_string(), length: 2000 },
    ];
    assert_eq!(binary_search_chrom_length(&chroms, "chr3"), None);
    assert_eq!(binary_search_chrom_length(&chroms, ""), None);
}

#[test]
fn binary_search_chrom_length_empty_vec() {
    let chroms: Vec<ChromInfo> = vec![];
    assert_eq!(binary_search_chrom_length(&chroms, "chr1"), None);
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cargo test --lib io::readers::bwig::tests::binary_search_chrom_length`
Expected: PASS

- [ ] **Step 5: Expose `SharedBigWigReader` access in `BigWigReader`**

In `src/io/readers/bwig.rs`, add a method to `BigWigReader` (after `chroms()` at line 374):

```rust
pub fn shared(&self) -> &Arc<SharedBigWigReader> {
    &self.shared
}
```

- [ ] **Step 6: Remove `chrom_lengths` from `Sample` and delegate to shared reader**

Replace the entire contents of `src/pipeline/core/samples.rs` with:

```rust
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context, Result, anyhow};

use crate::io::{BigWigReader, SharedBigWigReader};

pub struct Sample {
    path: PathBuf,
    reader: BigWigReader,
}

impl Sample {
    pub fn open(path: &Path) -> Result<Self> {
        let reader = BigWigReader::open(path)
            .with_context(|| format!("Failed to open bigWig file '{}'", path.display()))?;
        Ok(Self {
            path: path.to_path_buf(),
            reader,
        })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn reader(&self) -> &BigWigReader {
        &self.reader
    }

    pub fn reader_mut(&mut self) -> &mut BigWigReader {
        &mut self.reader
    }

    pub fn from_shared(path: PathBuf, shared: Arc<SharedBigWigReader>) -> Self {
        Self {
            path,
            reader: BigWigReader::from_shared(shared),
        }
    }

    pub fn chrom_length(&self, chrom: &str) -> Option<u32> {
        self.reader.shared().find_chrom_length(chrom)
    }
}

pub struct WorkerSamples {
    samples: Result<Vec<Sample>, String>,
}

impl WorkerSamples {
    pub fn new(paths: Arc<Vec<PathBuf>>) -> Self {
        let samples = open_samples(paths.as_ref()).map_err(|err| err.to_string());
        Self { samples }
    }

    pub fn from_shared(
        paths: Arc<Vec<PathBuf>>,
        shared_readers: Arc<Vec<Arc<SharedBigWigReader>>>,
    ) -> Self {
        let samples = paths
            .iter()
            .zip(shared_readers.iter())
            .map(|(path, shared)| Sample::from_shared(path.clone(), Arc::clone(shared)))
            .collect();
        Self { samples: Ok(samples) }
    }

    pub fn samples(&mut self) -> Result<&mut Vec<Sample>> {
        match &mut self.samples {
            Ok(samples) => Ok(samples),
            Err(message) => Err(anyhow!(message.clone())),
        }
    }
}

fn open_samples(paths: &[PathBuf]) -> Result<Vec<Sample>> {
    let mut samples = Vec::with_capacity(paths.len());
    for path in paths {
        samples.push(Sample::open(path)?);
    }
    Ok(samples)
}
```

- [ ] **Step 7: Run full test suite to verify no regressions**

Run: `cargo test`
Expected: All tests PASS. The `chrom_length()` method now delegates to `SharedBigWigReader::find_chrom_length()` which does the same binary search.

- [ ] **Step 8: Commit**

```bash
git add src/io/readers/bwig.rs src/pipeline/core/samples.rs
git commit -m "refactor: eliminate per-sample chrom_lengths HashMap clone

Delegate Sample::chrom_length() to SharedBigWigReader::find_chrom_length()
which binary-searches the already-sorted chroms Vec. Removes one
HashMap<String, u32> clone per Sample instance."
```

---

### Task 2: Move CIR node cache to SharedBigWigReader

The CIR tree node cache stores immutable `Arc<CachedCirNode>` entries keyed by file offset. Currently each `BigWigReader` instance maintains its own `HashMap<u64, Arc<CachedCirNode>>` with a crude "clear all at 1000 entries" eviction. Since nodes are immutable and identical across all readers of the same file, this cache should live in `SharedBigWigReader` using `quick_cache::sync::Cache` for proper LRU eviction and cross-worker sharing.

**Files:**
- Modify: `src/io/readers/bwig.rs:59-60, 109-116, 338-347, 359-371, 410-472`
- Test: existing tests (CIR tree behavior is covered by integration/regression tests; the shared cache is a transparent backend swap)

- [ ] **Step 1: Add shared CIR cache to `SharedBigWigReader`**

In `src/io/readers/bwig.rs`, update the imports at the top of the file. Add:

```rust
use quick_cache::sync::Cache;
```

Change the constant (line 59):

```rust
const MAX_CIR_CACHE_ENTRIES: usize = 1000;
```

to:

```rust
const SHARED_CIR_CACHE_ENTRIES: usize = 1000;
```

Add a `cir_cache` field to `SharedBigWigReader` (line 109-116):

```rust
pub struct SharedBigWigReader {
    file: File,
    uncompress_buf_size: usize,
    chroms: Vec<ChromInfo>,
    chrom_id_by_name: Vec<(String, u32)>,
    cir_tree_root: u64,
    block_cache: Arc<SharedBlockCache>,
    cir_cache: Cache<u64, Arc<CachedCirNode>>,
}
```

Initialize it in `open_with_cache` (after line 169, before `Ok(Self { ... })`):

```rust
Ok(Self {
    file,
    uncompress_buf_size,
    chroms,
    chrom_id_by_name,
    cir_tree_root,
    block_cache,
    cir_cache: Cache::new(SHARED_CIR_CACHE_ENTRIES),
})
```

Add a lookup method to `SharedBigWigReader`:

```rust
fn get_or_read_cir_node(&self, offset: u64) -> io::Result<Arc<CachedCirNode>> {
    if let Some(node) = self.cir_cache.get(&offset) {
        return Ok(node);
    }
    let parsed = Self::read_cir_node_raw(&self.file, offset)?;
    let arc_node = Arc::new(parsed);
    self.cir_cache.insert(offset, Arc::clone(&arc_node));
    Ok(arc_node)
}
```

- [ ] **Step 2: Remove per-reader CIR cache and update `search_cir_tree`**

In `BigWigReader` struct (line 338-347), remove the `cir_node_cache` field:

```rust
pub struct BigWigReader {
    shared: Arc<SharedBigWigReader>,
    work_buf: Vec<u8>,
    decode_buf: Vec<u8>,
    values_buf: Vec<BigWigValue>,
    blocks_buf: Vec<Block>,
    remaining_buf: VecDeque<u64>,
    pub stats: BigWigReaderStats,
}
```

In `from_shared` (line 359-371), remove the `cir_node_cache` initialization:

```rust
pub fn from_shared(shared: Arc<SharedBigWigReader>) -> Self {
    let uncompress_buf_size = shared.uncompress_buf_size;
    Self {
        shared,
        work_buf: Vec::with_capacity(uncompress_buf_size),
        decode_buf: Vec::new(),
        values_buf: Vec::new(),
        blocks_buf: Vec::new(),
        remaining_buf: VecDeque::new(),
        stats: BigWigReaderStats::default(),
    }
}
```

Replace `search_cir_tree` (line 410-472) to use the shared cache:

```rust
fn search_cir_tree(
    &mut self,
    chrom_ix: u32,
    start: u32,
    end: u32,
) -> io::Result<()> {
    let cir_tree_root = self.shared.cir_tree_root;

    self.blocks_buf.clear();
    self.remaining_buf.clear();
    self.remaining_buf.push_front(cir_tree_root);

    while let Some(node_offset) = self.remaining_buf.pop_front() {
        let node = self.shared.get_or_read_cir_node(node_offset)?;

        for item in &node.items {
            if item.end_chrom_id < chrom_ix || item.start_chrom_id > chrom_ix {
                continue;
            }
            if item.start_chrom_id == item.end_chrom_id {
                if item.end_base <= start || item.start_base >= end {
                    if item.start_chrom_id == chrom_ix {
                        continue;
                    }
                }
            }

            if node.is_leaf {
                self.blocks_buf.push(Block {
                    offset: item.data_offset,
                    size: item.data_size,
                });
            } else {
                self.remaining_buf.push_front(item.data_offset);
            }
        }
    }

    Ok(())
}
```

Note: The `cir_cache_hits`, `cir_cache_misses`, and `cir_cache_clears` stats fields become stale. Either remove them from `BigWigReaderStats` or keep them as zeroes. Since `BigWigReaderStats` is not referenced outside `bwig.rs` (verified earlier), remove those three fields:

```rust
#[derive(Debug, Default)]
pub struct BigWigReaderStats {
    pub values_calls: u64,
    pub values_returned: u64,
    pub blocks_per_query_total: u64,
    pub block_cache_hits: u64,
    pub block_cache_misses: u64,
    pub decoded_bytes: u64,
}
```

Also remove the `use std::collections::HashMap;` import if it's no longer used (the `VecDeque` import stays).

- [ ] **Step 3: Run full test suite**

Run: `cargo test`
Expected: All tests PASS. The CIR tree search produces identical results; only the caching backend changed.

- [ ] **Step 4: Commit**

```bash
git add src/io/readers/bwig.rs
git commit -m "refactor: move CIR node cache to SharedBigWigReader

Replace per-reader HashMap<u64, Arc<CachedCirNode>> with a shared
quick_cache::sync::Cache on SharedBigWigReader. All workers for the
same file now share one LRU cache (1000 entries) instead of each
maintaining a separate copy with crude clear-all eviction."
```

---

### Task 3: Share work_buf / decode_buf within each worker

Within `process_batch`, samples are iterated serially — only one `BigWigReader` is active at a time. Currently each reader carries its own `work_buf` (pre-allocated to `uncompress_buf_size`, typically 256 KB) and `decode_buf`. With T threads × S samples, this wastes (S-1) × T sets of buffers. We share one pair per worker by extracting them from `BigWigReader` and passing `&mut` references.

**Critical design note:** `process_batch` is called once per batch, not once per worker — it is invoked from a `map_init` closure in `executor.rs:227-238`. The `map_init` creates a `WorkerSamples` that is reused across batches on the same rayon thread. The decompression buffers **must live in `WorkerSamples`** (not in `process_batch` locals) so that their `Vec` capacity is retained across batches, avoiding per-batch re-allocation.

**Dependency:** This task depends on Task 1 (for the `shared()` accessor) and Task 2 (for the `BigWigReader` struct changes).

**Files:**
- Modify: `src/io/readers/bwig.rs:338-502` (add `values_with_bufs`, remove buf fields, add `uncompress_buf_size()`)
- Modify: `src/pipeline/core/samples.rs` (add `work_buf`/`decode_buf` to `WorkerSamples`, add buf-forwarding methods)
- Modify: `src/pipeline/core/worker.rs:260-487` (accept `&mut` bufs, pass to reader)
- Modify: `src/pipeline/core/executor.rs:227-238` (pass bufs from `WorkerSamples` to `process_batch`)
- Test: existing tests

- [ ] **Step 1: Add `values_with_bufs` method to `BigWigReader`**

In `src/io/readers/bwig.rs`, add a new public method right after `values()` (after line 408):

```rust
pub fn values_with_bufs(
    &mut self,
    chrom: &str,
    start: u32,
    end: u32,
    work_buf: &mut Vec<u8>,
    decode_buf: &mut Vec<u8>,
) -> Result<&[BigWigValue], BigWigReadError> {
    self.stats.values_calls += 1;

    let chrom_id = self
        .shared
        .find_chrom_id(chrom)
        .ok_or_else(|| BigWigReadError::ChromNotFound(chrom.to_string()))?;

    self.values_buf.clear();
    self.search_cir_tree(chrom_id, start, end)?;

    for i in 0..self.blocks_buf.len() {
        let (offset, size) = (self.blocks_buf[i].offset, self.blocks_buf[i].size);
        let data = self.get_or_cache_block_with_bufs(offset, size, work_buf, decode_buf)?;
        if data.is_empty() {
            continue;
        }
        parse_block_values(&data, start, end, &mut self.values_buf);
    }

    self.stats.values_returned += self.values_buf.len() as u64;
    self.stats.blocks_per_query_total += self.blocks_buf.len() as u64;

    Ok(&self.values_buf)
}

fn get_or_cache_block_with_bufs(
    &mut self,
    offset: u64,
    size: u64,
    work_buf: &mut Vec<u8>,
    decode_buf: &mut Vec<u8>,
) -> io::Result<Arc<[u8]>> {
    let key = (offset, size);
    if let Some(data) = self.shared.block_cache.get(&key) {
        self.stats.block_cache_hits += 1;
        return Ok(data);
    }

    self.stats.block_cache_misses += 1;
    let raw = read_and_decompress(
        &self.shared.file,
        offset,
        size,
        work_buf,
        decode_buf,
    )?;
    self.stats.decoded_bytes += raw.len() as u64;
    let data: Arc<[u8]> = Arc::from(raw);

    if !data.is_empty() {
        self.shared.block_cache.insert(key, Arc::clone(&data));
    }

    Ok(data)
}
```

- [ ] **Step 2: Add `uncompress_buf_size()` accessor to `SharedBigWigReader`**

In `src/io/readers/bwig.rs`, add to the `impl SharedBigWigReader` block:

```rust
pub fn uncompress_buf_size(&self) -> usize {
    self.uncompress_buf_size
}
```

- [ ] **Step 3: Remove `work_buf` and `decode_buf` from `BigWigReader`**

Update the struct definition (should now be, after Task 2):

```rust
pub struct BigWigReader {
    shared: Arc<SharedBigWigReader>,
    values_buf: Vec<BigWigValue>,
    blocks_buf: Vec<Block>,
    remaining_buf: VecDeque<u64>,
    pub stats: BigWigReaderStats,
}
```

Update `from_shared`:

```rust
pub fn from_shared(shared: Arc<SharedBigWigReader>) -> Self {
    Self {
        shared,
        values_buf: Vec::new(),
        blocks_buf: Vec::new(),
        remaining_buf: VecDeque::new(),
        stats: BigWigReaderStats::default(),
    }
}
```

Update the original `values()` method to allocate temporary local buffers (fallback for callers that don't pass external buffers):

```rust
pub fn values(
    &mut self,
    chrom: &str,
    start: u32,
    end: u32,
) -> Result<&[BigWigValue], BigWigReadError> {
    let mut work_buf = Vec::with_capacity(self.shared.uncompress_buf_size);
    let mut decode_buf = Vec::new();
    self.values_with_bufs(chrom, start, end, &mut work_buf, &mut decode_buf)
}
```

Also remove the `get_or_cache_block` method (the old one that used `self.work_buf` / `self.decode_buf`) since it's replaced by `get_or_cache_block_with_bufs`.

- [ ] **Step 4: Add decompression buffers and forwarding method to `WorkerSamples`**

In `src/pipeline/core/samples.rs`, add decompression buffer fields to `WorkerSamples` and a method on `Sample`:

Add to `Sample` impl:

```rust
pub fn values_with_bufs(
    &mut self,
    chrom: &str,
    start: u32,
    end: u32,
    work_buf: &mut Vec<u8>,
    decode_buf: &mut Vec<u8>,
) -> Result<&[crate::io::BigWigValue], anyhow::Error> {
    self.reader
        .values_with_bufs(chrom, start, end, work_buf, decode_buf)
        .map_err(anyhow::Error::new)
}
```

Update `WorkerSamples`:

```rust
pub struct WorkerSamples {
    samples: Result<Vec<Sample>, String>,
    work_buf: Vec<u8>,
    decode_buf: Vec<u8>,
}
```

Update `WorkerSamples::new`:

```rust
pub fn new(paths: Arc<Vec<PathBuf>>) -> Self {
    let samples = open_samples(paths.as_ref()).map_err(|err| err.to_string());
    let max_buf_size = match &samples {
        Ok(s) => s.iter()
            .map(|sample| sample.reader().shared().uncompress_buf_size())
            .max()
            .unwrap_or(0),
        Err(_) => 0,
    };
    Self {
        samples,
        work_buf: Vec::with_capacity(max_buf_size),
        decode_buf: Vec::new(),
    }
}
```

Update `WorkerSamples::from_shared`:

```rust
pub fn from_shared(
    paths: Arc<Vec<PathBuf>>,
    shared_readers: Arc<Vec<Arc<SharedBigWigReader>>>,
) -> Self {
    let max_buf_size = shared_readers.iter()
        .map(|s| s.uncompress_buf_size())
        .max()
        .unwrap_or(0);
    let samples = paths
        .iter()
        .zip(shared_readers.iter())
        .map(|(path, shared)| Sample::from_shared(path.clone(), Arc::clone(shared)))
        .collect();
    Self {
        samples: Ok(samples),
        work_buf: Vec::with_capacity(max_buf_size),
        decode_buf: Vec::new(),
    }
}
```

Note: the `decompression_bufs()` accessor is not needed — Step 6 uses `samples_and_bufs()` instead to avoid borrow checker issues.

- [ ] **Step 5: Update `process_batch` signature and callers**

In `src/pipeline/core/worker.rs`, change `process_batch` to accept buffer references:

```rust
pub(super) fn process_batch<M: PipelineMode>(
    samples: &mut [Sample],
    batch: CoalescedBatch,
    mode: &M,
    general: &GeneralOptions,
    metadata: &M::Metadata,
    work_buf: &mut Vec<u8>,
    decode_buf: &mut Vec<u8>,
) -> Result<Vec<(usize, usize, Option<MatrixRow>)>> {
```

Then update the bigWig read loop (around line 339-349). Replace:

```rust
let intervals = sample
    .reader_mut()
    .values(chrom, fetch_start, fetch_end)
    .map_err(anyhow::Error::new)
    .with_context(|| {
        format!(
            "Failed to read bigWig intervals for '{}' in '{}'",
            chrom,
            sample_paths[si].display()
        )
    })?;
```

With:

```rust
let intervals = sample
    .values_with_bufs(chrom, fetch_start, fetch_end, work_buf, decode_buf)
    .with_context(|| {
        format!(
            "Failed to read bigWig intervals for '{}' in '{}'",
            chrom,
            sample_paths[si].display()
        )
    })?;
```

- [ ] **Step 6: Update the `dispatch_chunk` closure in executor.rs**

In `src/pipeline/core/executor.rs`, update the `map_init` closure (lines 227-238) to pass buffers from `WorkerSamples`. Because `samples()` and `decompression_bufs()` both borrow `&mut self`, calling them separately will fail the borrow checker. Use a single `samples_and_bufs()` method that destructures the struct fields in one match arm:

Add to `WorkerSamples` in `src/pipeline/core/samples.rs` (replace the existing `samples()` method or add alongside it):

```rust
pub fn samples_and_bufs(&mut self) -> Result<(&mut Vec<Sample>, &mut Vec<u8>, &mut Vec<u8>)> {
    match &mut self.samples {
        Ok(samples) => Ok((samples, &mut self.work_buf, &mut self.decode_buf)),
        Err(message) => Err(anyhow!(message.clone())),
    }
}
```

Then the closure in `executor.rs` becomes:

```rust
|worker_samples, batch| {
    let (samples, work_buf, decode_buf) = worker_samples.samples_and_bufs()?;
    process_batch(samples.as_mut_slice(), batch, mode, general, metadata_ref, work_buf, decode_buf)
}
```

- [ ] **Step 7: `compute_sample_bins` fallback path**

In `src/pipeline/core/worker.rs`, `compute_sample_bins` (line 91-227) also calls `sample.reader_mut().values()`. This function is called from `compute_row` which is the per-item fallback path (metagene mode). For this path, the fallback `values()` method already allocates temporary buffers internally. No change needed here — metagene is the rare path and correctness is preserved.

- [ ] **Step 8: Run full test suite**

Run: `cargo test`
Expected: All tests PASS.

- [ ] **Step 9: Commit**

```bash
git add src/io/readers/bwig.rs src/pipeline/core/samples.rs src/pipeline/core/worker.rs src/pipeline/core/executor.rs
git commit -m "refactor: share work_buf/decode_buf within each worker

Extract decompression buffers from BigWigReader into WorkerSamples,
created once per rayon worker via map_init and reused across batches.
Within process_batch, samples are iterated serially so one buffer pair
suffices. Reduces memory from T×S×buf_size to T×1×buf_size."
```

---

### Task 4: Adaptive block cache sizing (200 per file, 2000 global hard limit)

Replace the current fixed 500-entry global budget (split evenly across files) with a smarter policy: each file gets up to 200 entries by default, but the total across all files is hard-capped at 2000. When `file_count * 200 > 2000`, the 2000-entry budget is distributed across files using integer division with remainder: the first `2000 % file_count` files each get `2000 / file_count + 1` entries, the rest get `2000 / file_count`. This ensures `sum(per_file) == min(file_count * 200, 2000)` exactly, and every file gets a fair share. Files that receive 0 entries operate without caching — `SharedBlockCache::with_capacity(0)` already handles this correctly (inserts are no-ops, gets always return `None`).

**Files:**
- Modify: `src/io/readers/block_cache.rs:1-28, 39-48` (new capacity function + constants + update `SharedBlockCache::new()`)
- Modify: `src/pipeline/core/executor.rs:9-11, 164-175` (update call site)
- Modify: `src/io/readers/bwig.rs:120-125` (update `SharedBigWigReader::open` default constant)
- Test: `src/io/readers/block_cache.rs` (update/add tests)

- [ ] **Step 1: Write failing tests for the new capacity function**

In `src/io/readers/block_cache.rs`, replace the existing `split_block_cache_capacity` tests and add new tests at the bottom of the `mod tests` block:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    // ── compute_per_file_block_cache_capacity ──────────────────────────

    #[test]
    fn few_files_get_default_per_file() {
        // 3 files: 3*200=600 <= 2000, each gets 200
        for i in 0..3 {
            assert_eq!(compute_per_file_block_cache_capacity(3, i), 200);
        }
    }

    #[test]
    fn exactly_at_hard_limit() {
        // 10 files: 10*200=2000, each gets 200
        for i in 0..10 {
            assert_eq!(compute_per_file_block_cache_capacity(10, i), 200);
        }
    }

    #[test]
    fn over_hard_limit_distributes_with_remainder() {
        // 11 files: budget=2000, 2000/11=181 base, remainder=9
        // files 0..8 get 182, files 9..10 get 181
        assert_eq!(compute_per_file_block_cache_capacity(11, 0), 182);
        assert_eq!(compute_per_file_block_cache_capacity(11, 8), 182);
        assert_eq!(compute_per_file_block_cache_capacity(11, 9), 181);
        assert_eq!(compute_per_file_block_cache_capacity(11, 10), 181);
        // 20 files: 2000/20=100, remainder=0, all get 100
        for i in 0..20 {
            assert_eq!(compute_per_file_block_cache_capacity(20, i), 100);
        }
    }

    #[test]
    fn many_files_first_2000_get_one_rest_zero() {
        // 3000 files: 2000/3000=0 base, remainder=2000
        // first 2000 files get 1, last 1000 get 0
        assert_eq!(compute_per_file_block_cache_capacity(3000, 0), 1);
        assert_eq!(compute_per_file_block_cache_capacity(3000, 1999), 1);
        assert_eq!(compute_per_file_block_cache_capacity(3000, 2000), 0);
        assert_eq!(compute_per_file_block_cache_capacity(3000, 2999), 0);
    }

    #[test]
    fn total_equals_budget_exactly() {
        for file_count in [1, 3, 10, 11, 20, 50, 100, 500, 2000, 3000] {
            let total: usize = (0..file_count)
                .map(|i| compute_per_file_block_cache_capacity(file_count, i))
                .sum();
            let expected = std::cmp::min(
                file_count * DEFAULT_PER_FILE_BLOCK_CACHE_ENTRIES,
                HARD_LIMIT_TOTAL_BLOCK_CACHE_ENTRIES,
            );
            assert_eq!(
                total, expected,
                "file_count={file_count}: total={total} != expected={expected}"
            );
        }
    }

    #[test]
    fn single_file_gets_default() {
        assert_eq!(compute_per_file_block_cache_capacity(1, 0), 200);
    }

    #[test]
    fn zero_files_returns_zero() {
        assert_eq!(compute_per_file_block_cache_capacity(0, 0), 0);
    }

    // ── SharedBlockCache basic behavior (keep existing) ─────────────────

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

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test --lib io::readers::block_cache::tests`
Expected: FAIL — `compute_per_file_block_cache_capacity` not found.

- [ ] **Step 3: Implement the new capacity function**

In `src/io/readers/block_cache.rs`, replace `DEFAULT_TOTAL_BLOCK_CACHE_ENTRIES` and `split_block_cache_capacity` with:

```rust
pub(crate) const DEFAULT_PER_FILE_BLOCK_CACHE_ENTRIES: usize = 200;
pub(crate) const HARD_LIMIT_TOTAL_BLOCK_CACHE_ENTRIES: usize = 2000;

/// Compute the block-cache capacity for each file.
///
/// Compute the block-cache capacity for file `sample_index` out of `file_count`.
///
/// Each file gets up to `DEFAULT_PER_FILE_BLOCK_CACHE_ENTRIES` (200) entries.
/// If `file_count * 200` exceeds the hard limit (2000), the 2000-entry budget
/// is distributed with remainder: the first `remainder` files each get
/// `base + 1`, the rest get `base` (where `base = 2000 / file_count`).
/// This ensures `sum(per_file) == min(file_count * 200, 2000)` exactly.
pub fn compute_per_file_block_cache_capacity(file_count: usize, sample_index: usize) -> usize {
    if file_count == 0 || sample_index >= file_count {
        return 0;
    }
    let desired = file_count.saturating_mul(DEFAULT_PER_FILE_BLOCK_CACHE_ENTRIES);
    if desired <= HARD_LIMIT_TOTAL_BLOCK_CACHE_ENTRIES {
        DEFAULT_PER_FILE_BLOCK_CACHE_ENTRIES
    } else {
        let base = HARD_LIMIT_TOTAL_BLOCK_CACHE_ENTRIES / file_count;
        let remainder = HARD_LIMIT_TOTAL_BLOCK_CACHE_ENTRIES % file_count;
        if sample_index < remainder { base + 1 } else { base }
    }
}
```

Remove the old `split_block_cache_capacity` function and `DEFAULT_TOTAL_BLOCK_CACHE_ENTRIES` constant entirely.

- [ ] **Step 4: Update `SharedBlockCache::new()` to use the new constant**

In `src/io/readers/block_cache.rs`, the `new()` method (line 40-42) currently references `DEFAULT_TOTAL_BLOCK_CACHE_ENTRIES`. Update it:

```rust
pub fn new() -> Self {
    Self::with_capacity(DEFAULT_PER_FILE_BLOCK_CACHE_ENTRIES)
}
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `cargo test --lib io::readers::block_cache::tests`
Expected: All PASS.

- [ ] **Step 6: Update the executor call site**

In `src/pipeline/core/executor.rs`, update the imports (lines 9-11):

```rust
use crate::io::readers::block_cache::compute_per_file_block_cache_capacity;
```

Remove the old imports `DEFAULT_TOTAL_BLOCK_CACHE_ENTRIES` and `split_block_cache_capacity`.

Replace the cache capacity computation (lines 164-175):

```rust
let sample_count_for_cache = sample_paths.len();
let shared_readers = Arc::new(
    sample_paths
        .iter()
        .enumerate()
        .map(|(sample_index, path)| {
            let cache_capacity = compute_per_file_block_cache_capacity(
                sample_count_for_cache,
                sample_index,
            );
            SharedBigWigReader::open_with_block_cache_capacity(path, cache_capacity)
                .map(Arc::new)
                .with_context(|| format!("Failed to open bigWig file '{}'", path.display()))
        })
        .collect::<Result<Vec<_>>>()?,
);
```

- [ ] **Step 7: Update `SharedBigWigReader::open` default**

In `src/io/readers/bwig.rs`, update the `open` method (line 120-125) to use the new constant:

```rust
pub fn open(path: impl AsRef<Path>) -> Result<Self, BigWigReadError> {
    Self::open_with_block_cache_capacity(
        path,
        super::block_cache::DEFAULT_PER_FILE_BLOCK_CACHE_ENTRIES,
    )
}
```

- [ ] **Step 8: Run full test suite**

Run: `cargo test`
Expected: All tests PASS.

- [ ] **Step 9: Commit**

```bash
git add src/io/readers/block_cache.rs src/pipeline/core/executor.rs src/io/readers/bwig.rs
git commit -m "feat: adaptive block cache sizing (200/file, 2000 global hard cap)

Replace the fixed 500-entry global budget split evenly across files
with a per-file default of 200 entries. When file_count * 200 exceeds
the 2000-entry hard limit, the budget is distributed evenly via integer
division. Files receiving 0 entries operate without caching.
Invariant: sum(per_file) <= 2000 always holds."
```

---

## Post-Implementation: Update overall_status.md

After all four tasks are complete, update `plans/overall_status.md` to record these changes under a new section. Add a summary of what changed (chrom_lengths removal, CIR cache sharing, buffer consolidation, adaptive block cache sizing) and commit:

```bash
git add plans/overall_status.md
git commit -m "docs: update overall_status.md with bigwig reader optimizations"
```

## Verification

After all four tasks are complete:

1. **Unit tests**: `cargo test` — all pass
2. **Regression test**: `scripts/custom_compare.py` with tolerance 5e-6 against Python baseline
3. **Performance benchmark**: `scripts/profile_bench.sh` to verify:
   - RSS decreased (fewer HashMap clones, fewer buffer duplicates)
   - No wall-clock regression (CIR shared cache may slightly improve due to cross-worker warming)
   - Note: block cache hit/miss stats are tracked internally in `BigWigReaderStats` but not yet wired to verbose output. Stats instrumentation will be added separately.

## Task Execution Order

```
Task 1 (chrom_lengths) → Task 2 (CIR cache) → Task 3 (work_buf/decode_buf) → Task 4 (block cache)
```

Although Tasks 1, 2, and 4 are logically independent, they all modify `src/io/readers/bwig.rs`. Running them in parallel via separate workers would cause merge conflicts. **Execute sequentially** in the order above:

- **Task 1 → Task 2**: both modify `BigWigReader` struct and `SharedBigWigReader` impl in bwig.rs
- **Task 2 → Task 3**: Task 3 depends on Task 2's struct changes (removing `cir_node_cache`) and Task 1's `shared()` accessor
- **Task 3 → Task 4**: Task 4 modifies `SharedBigWigReader::open` default in bwig.rs; doing it last avoids conflicts with Task 1/2's changes to the same file
