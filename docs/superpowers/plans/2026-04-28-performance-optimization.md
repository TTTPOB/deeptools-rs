# Performance Optimization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement 12 performance optimizations (CPU, memory, I/O) identified via perf + heaptrack profiling on ENCODE K562 ATAC-seq data.

**Architecture:** Tasks are grouped into 4 phases. Phase 0 sets up profiling infrastructure. Phases 1-2 are independent file-local optimizations. Phase 3 builds on prior phases with structural changes (CoalesceStrategy, channel pipeline). Each task includes a profile step using `scripts/profile_bench.sh` to verify improvement.

**Tech Stack:** Rust 2024 edition, rayon, zlib-rs, zune-inflate, flate2, std::sync::mpsc

**Key files:**
- `src/io/readers/bwig.rs` — BigWig reader (caching, deflate, work_buf) — T4, T5, T7 share this file
- `src/io/readers/bed.rs` — BED record parsing — T3
- `src/io/writers/mod.rs` — Matrix value formatting — T2
- `src/pipeline/core/mod.rs` — Core pipeline (aggregation, batching, execution) — T3, T6, T8, T9, T10, T11 share this file
- `src/pipeline/scale_regions.rs` — Scale-regions mode — T8, T9, T10
- `src/pipeline/reference_point.rs` — Reference-point mode — T8, T9, T10
- `src/pipeline/mod.rs` — Pipeline entry + spawn_writer_thread — T10
- `src/config.rs` — SortRegions enum — T9
- `Cargo.toml` — Dead dependency removal — T1

**File conflict note:** Tasks sharing a file MUST be dispatched sequentially (later task waits for earlier to commit). Orchestrator tracks this.

---

### Task 0: Profile Harness

**Files:**
- Create: `scripts/profile_bench.sh`
- Create: `bench_reports/.gitkeep`

- [ ] **Step 1: Create the profile script**

```bash
#!/bin/bash
set -euo pipefail

if [ $# -lt 4 ]; then
    echo "Usage: $0 <name> <target> <hot-path> -- <command...>"
    echo "Example: $0 p1-i64-div \"eliminate i128 div in write_matrix_value\" \"write_matrix_value -> rint + __divti3\" -- cargo run --release -- scale-regions ..."
    exit 1
fi

NAME="$1"
TARGET="$2"
HOT_PATH="$3"
shift 3

TIMESTAMP=$(date +%Y-%m-%d-%H%M%S)
REPORT="bench_reports/${TIMESTAMP}-${NAME}.md"
mkdir -p bench_reports

cat > "$REPORT" << EOF
# Profile: ${NAME}
Time: $(date -Iseconds)
Command: $@
Target: ${TARGET}
Hot path: ${HOT_PATH}

EOF

# /usr/bin/time -v
echo "## /usr/bin/time -v" >> "$REPORT"
echo '```' >> "$REPORT"
/usr/bin/time -v "$@" 2>> "$REPORT"
echo '```' >> "$REPORT"

# perf stat
echo "## perf stat" >> "$REPORT"
echo '```' >> "$REPORT"
perf stat -d "$@" 2>> "$REPORT"
echo '```' >> "$REPORT"

# heaptrack
echo "## heaptrack" >> "$REPORT"
echo '```' >> "$REPORT"
heaptrack "$@" 2>> "$REPORT"
echo '```' >> "$REPORT"

# Cleanup raw profiler files
rm -f perf.data heaptrack.*.gz 2>/dev/null

echo "Report written to $REPORT"
```

- [ ] **Step 2: Make it executable**

Run: `chmod +x scripts/profile_bench.sh`

- [ ] **Step 3: Create bench_reports/.gitkeep**

Run: `touch bench_reports/.gitkeep`

- [ ] **Step 4: Run baseline profile**

Run:
```bash
./scripts/profile_bench.sh baseline \
  "baseline before any optimizations" \
  "full pipeline" \
  -- cargo run --release -- scale-regions -b 10 -p 4 --unscaled5prime 0 --unscaled3prime 0 --region test_data/encode_k562_atac.bed --scoreFileName test_data/sample1.bw test_data/sample2.bw test_data/sample3.bw -o /tmp/bench_baseline.mat.gz
```

Note: Adjust paths to match your test data. The key is using a fixed, reproducible dataset.

- [ ] **Step 5: Verify report exists**

Run: `ls bench_reports/`
Expected: One `.md` file with the baseline report

- [ ] **Step 6: Commit**

```bash
git add scripts/profile_bench.sh bench_reports/.gitkeep bench_reports/*.md
git commit -m "feat: add profile harness with perf + heaptrack + time reporting

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 1: Remove Dead Dependencies and Code (AR1)

**Files:**
- Delete: `src/io/readers/bigwig.rs`
- Modify: `Cargo.toml:6-19`

- [ ] **Step 1: Delete the dead bigwig.rs file**

Run: `rm src/io/readers/bigwig.rs`

- [ ] **Step 2: Remove bigtools from Cargo.toml**

Edit `Cargo.toml`, remove line 8:
```diff
- bigtools = { version = "0.5.6"}
```

- [ ] **Step 3: Remove crossbeam-channel from Cargo.toml**

Edit `Cargo.toml`, remove line 19:
```diff
- crossbeam-channel = "0.5"
```

- [ ] **Step 4: Verify it compiles**

Run: `cargo build --release 2>&1`
Expected: Build succeeds without bigtools or crossbeam-channel

- [ ] **Step 5: Run profile_bench.sh**

Run:
```bash
./scripts/profile_bench.sh ar1-dead-deps \
  "remove unused bigtools + crossbeam-channel deps" \
  "compile time / binary size" \
  -- cargo run --release -- scale-regions -b 10 -p 4 --unscaled5prime 0 --unscaled3prime 0 --region test_data/encode_k562_atac.bed --scoreFileName test_data/sample1.bw test_data/sample2.bw test_data/sample3.bw -o /tmp/bench_ar1.mat.gz
```

- [ ] **Step 6: Verify correct output**

Run: `scripts/custom_compare.py /tmp/bench_baseline.mat.gz /tmp/bench_ar1.mat.gz`
Expected: tolerance 5e-6, no differences

- [ ] **Step 7: Commit**

```bash
git add Cargo.toml src/io/readers/bigwig.rs bench_reports/
git commit -m "chore: remove unused bigtools and crossbeam-channel dependencies

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 2: Fast i64 Division Path in write_matrix_value (P1)

**Files:**
- Modify: `src/io/writers/mod.rs:494-543`

- [ ] **Step 1: Implement fast path for i64 range values**

Replace the `write_matrix_value` function in `src/io/writers/mod.rs:494-543`:

```rust
fn write_matrix_value<W: Write>(writer: &mut W, value: f32) -> io::Result<()> {
    if value.is_nan() {
        return writer.write_all(b"nan");
    }

    if value.is_infinite() {
        return if value.is_sign_negative() {
            writer.write_all(b"-inf")
        } else {
            writer.write_all(b"inf")
        };
    }

    // Fast path: for values whose scaled integer fits in i64 (abs < ~9e12),
    // avoid the expensive 128-bit division. Genomic data values are almost
    // always in this range.
    if value.is_finite() && value > -1e7 && value < 1e7 {
        let scaled = (value as f64 * 1_000_000.0).round_ties_even();
        return write_scaled_i64(writer, scaled as i64);
    }

    // Slow path: fallback for extreme values that need i128 range
    let scaled = (value as f64 * 1_000_000.0).round_ties_even();
    if !scaled.is_finite() || scaled.abs() > i128::MAX as f64 {
        let fallback = format!("{value:.6}");
        return writer.write_all(fallback.as_bytes());
    }

    write_scaled_i128(writer, scaled as i128)
}

#[inline]
fn write_scaled_i64<W: Write>(writer: &mut W, scaled: i64) -> io::Result<()> {
    let mut buffer = itoa::Buffer::new();

    if scaled == 0 {
        return writer.write_all(b"0.000000");
    }

    if scaled < 0 {
        writer.write_all(b"-")?;
    }
    let abs = scaled.unsigned_abs();
    let integer_part = abs / 1_000_000;
    let fractional_part = (abs % 1_000_000) as u32;

    writer.write_all(buffer.format(integer_part).as_bytes())?;
    writer.write_all(b".")?;

    let mut frac_digits = [b'0'; 6];
    let mut remainder = fractional_part;
    for slot in frac_digits.iter_mut().rev() {
        *slot = b'0' + (remainder % 10) as u8;
        remainder /= 10;
    }
    writer.write_all(&frac_digits)
}

#[inline]
fn write_scaled_i128<W: Write>(writer: &mut W, scaled: i128) -> io::Result<()> {
    let mut buffer = itoa::Buffer::new();
    let sign_negative = scaled < 0;

    if scaled == 0 {
        if sign_negative {
            writer.write_all(b"-")?;
        }
        return writer.write_all(b"0.000000");
    }

    let abs = if sign_negative { -scaled } else { scaled };
    let integer_part = (abs / 1_000_000) as u128;
    let fractional_part = (abs % 1_000_000) as u32;

    if sign_negative {
        writer.write_all(b"-")?;
    }
    writer.write_all(buffer.format(integer_part).as_bytes())?;
    writer.write_all(b".")?;

    let mut frac_digits = [b'0'; 6];
    let mut remainder = fractional_part;
    for slot in frac_digits.iter_mut().rev() {
        *slot = b'0' + (remainder % 10) as u8;
        remainder /= 10;
    }
    writer.write_all(&frac_digits)
}
```

- [ ] **Step 2: Build**

Run: `cargo build --release 2>&1`
Expected: Compiles successfully

- [ ] **Step 3: Run profile_bench.sh**

Run:
```bash
./scripts/profile_bench.sh p1-i64-div \
  "eliminate i128 division in write_matrix_value" \
  "write_matrix_value -> rint + __divti3" \
  -- cargo run --release -- scale-regions -b 10 -p 4 --unscaled5prime 0 --unscaled3prime 0 --region test_data/encode_k562_atac.bed --scoreFileName test_data/sample1.bw test_data/sample2.bw test_data/sample3.bw -o /tmp/bench_p1.mat.gz
```

- [ ] **Step 4: Verify correct output**

Run: `scripts/custom_compare.py /tmp/bench_baseline.mat.gz /tmp/bench_p1.mat.gz`
Expected: tolerance 5e-6, no differences

- [ ] **Step 5: Read profile report and compare vs baseline**

Read the report at `bench_reports/*-p1-i64-div.md`. Compare `write_matrix_value` CPU% in perf stat vs baseline. Fill in `Other hotspots observed`.

Expected: write_matrix_value drops from ~5.36% to ~2-3% CPU

- [ ] **Step 6: Run existing tests**

Run: `cargo test 2>&1`
Expected: All tests pass

- [ ] **Step 7: Commit**

```bash
git add src/io/writers/mod.rs bench_reports/
git commit -m "perf: use i64 fast path in write_matrix_value to avoid 128-bit division

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 3: Chromosome Name Interning (M2)

**Files:**
- Modify: `src/io/readers/bed.rs:80` (chrom field in BedRecord::parse)
- Modify: `src/io/readers/bed.rs:49-59` (BedRecord.chrom type)
- Modify: `src/pipeline/core/mod.rs` (all `record.chrom` comparisons — switch from `String` to `Arc<str>` equality)

- [ ] **Step 1: Add interner module to bed.rs**

At the top of `src/io/readers/bed.rs`, after the imports, add:

```rust
use std::cell::RefCell;
use std::collections::HashMap;

thread_local! {
    static CHROM_INTERNER: RefCell<HashMap<String, Arc<str>>> = RefCell::new(HashMap::new());
}

fn intern_chrom(s: String) -> Arc<str> {
    CHROM_INTERNER.with(|map| {
        let mut map = map.borrow_mut();
        if let Some(existing) = map.get(&s) {
            Arc::clone(existing)
        } else {
            let arc: Arc<str> = Arc::from(s.as_str());
            map.insert(s, arc.clone());
            arc
        }
    })
}
```

- [ ] **Step 2: Change BedRecord.chrom type**

Edit `BedRecord` struct at `src/io/readers/bed.rs:49-59`:

```rust
#[derive(Debug, Clone, PartialEq)]
pub struct BedRecord {
    pub chrom: Arc<str>,
    pub start: u32,
    pub end: u32,
    pub name: Option<String>,
    pub score: Option<f32>,
    pub score_raw: Option<String>,
    pub strand: Strand,
    pub strand_raw: Option<String>,
    pub extra_fields: Vec<String>,
}
```

- [ ] **Step 3: Apply interning in BedRecord::parse**

Edit `src/io/readers/bed.rs:80`, change:
```rust
let chrom = fields[0].to_string();
```
to:
```rust
let chrom = intern_chrom(fields[0].to_string());
```

- [ ] **Step 4: Update all sites that compare or clone chrom**

In `src/pipeline/core/mod.rs`, find all places that compare `record.chrom` and update:
- Line ~905: `item.record.chrom != current_chrom` — works with `Arc<str>` (PartialEq)
- Line ~917: `current_chrom = item.record.chrom.clone();` — `Arc<str>::clone()` is cheap ref-count increment
- Any `chrom.clone()` calls become cheap `Arc::clone()`

Run to find all chrom references:
```bash
grep -n '\.chrom' src/pipeline/core/mod.rs src/pipeline/reference_point.rs src/pipeline/scale_regions.rs
```

For each site, if it's a `.clone()` on `chrom`, it's now an `Arc::clone()` (cheap). If it's a comparison, `Arc<str>: PartialEq<str>` via Deref — check if explicit `.as_ref()` is needed.

Common pattern fix — `current_chrom` variable type:
```rust
// Before:
let mut current_chrom = String::new();
// After:
let mut current_chrom: Arc<str> = Arc::from("");
```

- [ ] **Step 5: Build and fix any type errors**

Run: `cargo build --release 2>&1`
Fix any compilation errors from `String` → `Arc<str>` mismatch. Key areas:
- `process_batch` in core/mod.rs (~line 986): `let chrom = &batch.items[0].2.chrom;` — `&Arc<str>` derefs to `&str` via Deref, should work
- `BigWigReader::values(chrom, ...)` takes `&str`, `Arc<str>` derefs to `str`

- [ ] **Step 6: Run profile_bench.sh**

```bash
./scripts/profile_bench.sh m2-chrom-intern \
  "interning chromosome names to share Arc<str> across records" \
  "load_groups -> BedRecord::parse -> chrom String allocation" \
  -- cargo run --release -- scale-regions -b 10 -p 4 --unscaled5prime 0 --unscaled3prime 0 --region test_data/encode_k562_atac.bed --scoreFileName test_data/sample1.bw test_data/sample2.bw test_data/sample3.bw -o /tmp/bench_m2.mat.gz
```

- [ ] **Step 7: Verify correct output and run tests**

Run:
```bash
scripts/custom_compare.py /tmp/bench_baseline.mat.gz /tmp/bench_m2.mat.gz
cargo test 2>&1
```

Expected: output matches, all tests pass, peak heap reduced ~20 MB

- [ ] **Step 8: Read profile report**

Read `bench_reports/*-m2-chrom-intern.md`. Compare peak heap vs baseline. Fill `Other hotspots observed`.

- [ ] **Step 9: Commit**

```bash
git add src/io/readers/bed.rs src/pipeline/core/mod.rs src/pipeline/reference_point.rs src/pipeline/scale_regions.rs bench_reports/
git commit -m "perf: intern chromosome names to share Arc<str> across BED records

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 4: block_cache Arc + CIR Node Cache Limit (A2 + M3)

**Files:**
- Modify: `src/io/readers/bwig.rs:291-298` (BigWigReader struct)
- Modify: `src/io/readers/bwig.rs:401-423` (get_or_cache_block)

- [ ] **Step 1: Change block_cache to store Arc<[u8]>**

Edit `BigWigReader` struct at `src/io/readers/bwig.rs:294-298`:

```rust
pub struct BigWigReader {
    shared: Arc<SharedBigWigReader>,
    cir_node_cache: HashMap<u64, Arc<CachedCirNode>>,
    block_cache: HashMap<(u64, u64), Arc<[u8]>>,
}
```

- [ ] **Step 2: Add CIR cache size limit constant**

At the top of `src/io/readers/bwig.rs`, near the other constants (search for `MAX_BLOCK_CACHE_ENTRIES`):

```rust
const MAX_CIR_CACHE_ENTRIES: usize = 1000;
```

- [ ] **Step 3: Update get_or_cache_block to store and return Arc<[u8]>**

Replace `get_or_cache_block` at `src/io/readers/bwig.rs:401-423`:

```rust
fn get_or_cache_block(
    &mut self,
    offset: u64,
    size: u64,
    work_buf: &mut Vec<u8>,
) -> io::Result<Arc<[u8]>> {
    let key = (offset, size);
    if let Some(data) = self.block_cache.get(&key) {
        return Ok(Arc::clone(data));
    }

    let raw = read_and_decompress(&self.shared.file, offset, size, work_buf)?;
    let data: Arc<[u8]> = Arc::from(raw.to_vec().into_boxed_slice());

    if !data.is_empty() {
        if self.block_cache.len() >= MAX_BLOCK_CACHE_ENTRIES {
            self.block_cache.clear();
        }
        self.block_cache.insert(key, Arc::clone(&data));
    }

    Ok(data)
}
```

- [ ] **Step 4: Add CIR cache size check in search_cir_tree**

Find the CIR cache insertion in `search_cir_tree` (around line 360-380 in bwig.rs). Find where `cir_node_cache.insert` is called and add a size check before it:

```rust
if self.cir_node_cache.len() >= MAX_CIR_CACHE_ENTRIES {
    self.cir_node_cache.clear();
}
self.cir_node_cache.insert(offset, Arc::new(node));
```

- [ ] **Step 5: Update callers of get_or_cache_block**

In `values()` at line ~340, the call is:
```rust
let data = self.get_or_cache_block(block.offset, block.size, &mut work_buf)?;
```
Now `data` is `Arc<[u8]>` instead of `Vec<u8>`. The `parse_block_values` function takes `&[u8]`, so `&data` (which derefs `Arc<[u8]>` to `&[u8]`) works. Later the `data.is_empty()` check also works via Deref.

No changes needed to callers — `Arc<[u8]>` derefs to `&[u8]`.

- [ ] **Step 6: Build**

Run: `cargo build --release 2>&1`
Expected: Compiles successfully

- [ ] **Step 7: Run profile_bench.sh**

```bash
./scripts/profile_bench.sh a2-arc-cache \
  "block_cache uses Arc<[u8]> to avoid clone; CIR cache has 1000-entry limit" \
  "BigWigReader::values -> get_or_cache_block -> data.clone()" \
  -- cargo run --release -- scale-regions -b 10 -p 4 --unscaled5prime 0 --unscaled3prime 0 --region test_data/encode_k562_atac.bed --scoreFileName test_data/sample1.bw test_data/sample2.bw test_data/sample3.bw -o /tmp/bench_a2.mat.gz
```

- [ ] **Step 8: Verify correct output and run tests**

Run:
```bash
scripts/custom_compare.py /tmp/bench_baseline.mat.gz /tmp/bench_a2.mat.gz
cargo test 2>&1
```

- [ ] **Step 9: Read profile report and commit**

```bash
git add src/io/readers/bwig.rs bench_reports/
git commit -m "perf: use Arc<[u8]> for block cache and add CIR cache size limit

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 5: Replace zune-inflate with zlib-rs Deflate (A1)

**Files:**
- Modify: `src/io/readers/bwig.rs:1-9` (imports)
- Modify: `src/io/readers/bwig.rs:449-458` (deflate call in read_and_decompress)

- [ ] **Step 1: Replace import**

Edit `src/io/readers/bwig.rs:9`:
```rust
// Remove:
use zune_inflate::DeflateDecoder;
// Add:
use flate2::Decompress;
```

The `flate2` crate is already in Cargo.toml with `zlib-rs` feature enabled.

- [ ] **Step 2: Replace the deflate call in read_and_decompress**

Replace lines 449-458 in `src/io/readers/bwig.rs`:

```rust
    if block[0] == 0x78 {
        // zlib compressed — use zlib-rs via flate2
        let mut decoder = Decompress::new(true); // true = zlib wrapper
        let mut decoded = Vec::with_capacity(buf_len * 4); // typical compression ratio
        decoder
            .decompress_vec(block, &mut decoded, flate2::FlushDecompress::Finish)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;

        let len = decoded.len();
        if len > work_buf.len() {
            work_buf.resize(len, 0);
        }
        work_buf[..len].copy_from_slice(&decoded);
        Ok(&work_buf[..len])
    } else {
        // Uncompressed — already in work_buf
        Ok(&work_buf[..buf_len])
    }
```

- [ ] **Step 3: Build**

Run: `cargo build --release 2>&1`
Expected: Compiles. If `decompress_vec` API differs, adjust per flate2 docs.

- [ ] **Step 4: Run profile_bench.sh — zlib-rs variant**

```bash
./scripts/profile_bench.sh a1-zlib-rs \
  "replace zune-inflate with zlib-rs deflate for decompression" \
  "BigWigReader::values -> read_and_decompress -> deflate" \
  -- cargo run --release -- scale-regions -b 10 -p 4 --unscaled5prime 0 --unscaled3prime 0 --region test_data/encode_k562_atac.bed --scoreFileName test_data/sample1.bw test_data/sample2.bw test_data/sample3.bw -o /tmp/bench_a1_zlib.mat.gz
```

- [ ] **Step 5: Verify correct output**

Run: `scripts/custom_compare.py /tmp/bench_baseline.mat.gz /tmp/bench_a1_zlib.mat.gz`
Expected: tolerance 5e-6, no differences

- [ ] **Step 6: Remove zune-inflate from Cargo.toml**

Edit `Cargo.toml`, remove:
```diff
- zune-inflate = "0.2"
```

- [ ] **Step 7: Build again and run tests**

Run:
```bash
cargo build --release 2>&1
cargo test 2>&1
```

- [ ] **Step 8: Read profile report and commit**

Check that `DeflateDecoder::new` / `zune_inflate` CPU% dropped from ~3% in perf stat.

```bash
git add src/io/readers/bwig.rs Cargo.toml bench_reports/
git commit -m "perf: replace zune-inflate with zlib-rs via flate2 for bigWig deflate

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 6: aggregate_slice f32 + Small-Slice Optimization (P2)

**Files:**
- Modify: `src/pipeline/core/mod.rs:449-531`

- [ ] **Step 1: Replace aggregate_slice with f32 accumulation and small-slice fast paths**

Replace the entire `aggregate_slice` function at `src/pipeline/core/mod.rs:449-531`:

```rust
fn aggregate_slice(slice: &[f32], average_type: AverageTypeBins) -> Option<f32> {
    let len = slice.len();
    if len == 0 {
        return None;
    }

    match average_type {
        AverageTypeBins::Mean => {
            let mut sum = 0.0f32;
            let mut count = 0u32;
            for &value in slice {
                if !value.is_nan() {
                    sum += value;
                    count += 1;
                }
            }
            if count == 0 { None } else { Some(sum / count as f32) }
        }
        AverageTypeBins::Sum => {
            let mut sum = 0.0f32;
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
            let mut min = f32::INFINITY;
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
            let mut max = f32::NEG_INFINITY;
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
            let mut sum = 0.0f32;
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
            let mean = sum / count as f32;
            let mut variance_sum = 0.0f64;
            for &value in slice {
                if !value.is_nan() {
                    let delta = value as f64 - mean as f64;
                    variance_sum += delta * delta;
                }
            }
            Some((variance_sum / count as f64).sqrt() as f32)
        }
        AverageTypeBins::Median => {
            let mut values: Vec<f32> = slice
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

Key changes from original:
- Mean/Sum/Min/Max use `f32` accumulation instead of `f64` — sufficient precision for genomics data
- Std keeps `f64` for variance to avoid catastrophic cancellation
- Median keeps its allocation (small N per call)
- Removed the pre-check `slice.len() <= 16` unrolling — the compiler auto-vectorizes simple f32 loops

- [ ] **Step 2: Build**

Run: `cargo build --release 2>&1`

- [ ] **Step 3: Run profile_bench.sh**

```bash
./scripts/profile_bench.sh p2-aggregate-f32 \
  "f32 accumulation in aggregate_slice (mean/sum/min/max)" \
  "aggregate_slice -> f64 cast + f64 accumulation" \
  -- cargo run --release -- scale-regions -b 10 -p 4 --unscaled5prime 0 --unscaled3prime 0 --region test_data/encode_k562_atac.bed --scoreFileName test_data/sample1.bw test_data/sample2.bw test_data/sample3.bw -o /tmp/bench_p2.mat.gz
```

- [ ] **Step 4: Verify correct output**

Run: `scripts/custom_compare.py /tmp/bench_baseline.mat.gz /tmp/bench_p2.mat.gz`
Expected: tolerance 5e-6. f32 accumulation may introduce very minor (< 1e-6) differences.

- [ ] **Step 5: Run existing tests**

Run: `cargo test 2>&1`
Expected: All tests pass, especially `aggregate_slice_ignores_nans`

- [ ] **Step 6: Commit**

```bash
git add src/pipeline/core/mod.rs bench_reports/
git commit -m "perf: use f32 accumulation in aggregate_slice for mean/sum/min/max

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 7: work_buf Reuse Across values() Calls (A3)

**Files:**
- Modify: `src/io/readers/bwig.rs:291-298` (BigWigReader struct)
- Modify: `src/io/readers/bwig.rs:322-348` (values method)
- Modify: `src/io/readers/bwig.rs:401-423` (get_or_cache_block signature)

- [ ] **Step 1: Add work_buf to BigWigReader struct**

Edit `src/io/readers/bwig.rs:294-298`:

```rust
pub struct BigWigReader {
    shared: Arc<SharedBigWigReader>,
    cir_node_cache: HashMap<u64, Arc<CachedCirNode>>,
    block_cache: HashMap<(u64, u64), Arc<[u8]>>,
    work_buf: Vec<u8>,
}
```

- [ ] **Step 2: Initialize work_buf in constructors**

Edit `from_shared` at line 310-316:

```rust
pub fn from_shared(shared: Arc<SharedBigWigReader>) -> Self {
    let uncompress_buf_size = shared.uncompress_buf_size;
    Self {
        shared,
        cir_node_cache: HashMap::new(),
        block_cache: HashMap::new(),
        work_buf: Vec::with_capacity(uncompress_buf_size),
    }
}
```

- [ ] **Step 3: Use struct field in values()**

In `values()` at line ~337, replace:
```rust
let mut work_buf = vec![0u8; self.shared.uncompress_buf_size];
```
with nothing (remove the local allocation). Then update the `get_or_cache_block` call at line ~340 to pass `&mut self.work_buf`:

```rust
let data = self.get_or_cache_block(block.offset, block.size)?;
```

- [ ] **Step 4: Update get_or_cache_block to use self.work_buf**

Change the signature at line 401-406:
```rust
fn get_or_cache_block(
    &mut self,
    offset: u64,
    size: u64,
) -> io::Result<Arc<[u8]>> {
```

And update the `read_and_decompress` call inside:
```rust
let raw = read_and_decompress(&self.shared.file, offset, size, &mut self.work_buf)?;
```

- [ ] **Step 5: Build**

Run: `cargo build --release 2>&1`
Expected: Compiles. If any other caller of `get_or_cache_block` passes `work_buf`, update it.

- [ ] **Step 6: Run profile_bench.sh**

```bash
./scripts/profile_bench.sh a3-work-buf-reuse \
  "reuse work_buf across values() calls instead of per-call allocation" \
  "BigWigReader::values -> vec![0u8; uncompress_buf_size] per call" \
  -- cargo run --release -- scale-regions -b 10 -p 4 --unscaled5prime 0 --unscaled3prime 0 --region test_data/encode_k562_atac.bed --scoreFileName test_data/sample1.bw test_data/sample2.bw test_data/sample3.bw -o /tmp/bench_a3.mat.gz
```

- [ ] **Step 7: Verify correct output and run tests**

```bash
scripts/custom_compare.py /tmp/bench_baseline.mat.gz /tmp/bench_a3.mat.gz
cargo test 2>&1
```

- [ ] **Step 8: Commit**

```bash
git add src/io/readers/bwig.rs bench_reports/
git commit -m "perf: reuse work_buf across BigWigReader::values calls

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 8: Arc<BedRecord> in RegionTask (M1)

**Files:**
- Modify: `src/pipeline/core/mod.rs:168-174` (RegionTask struct)
- Modify: `src/pipeline/core/mod.rs:828-834` (WorkItem struct — record type)
- Modify: `src/pipeline/core/mod.rs:884-891` (CoalescedBatch.items type)
- Modify: `src/pipeline/core/mod.rs` (process_batch, create_batches — BedRecord → Arc<BedRecord>)
- Modify: `src/pipeline/reference_point.rs` (RegionTask construction)
- Modify: `src/pipeline/scale_regions.rs` (RegionTask construction)

- [ ] **Step 1: Change RegionTask to use Arc<BedRecord>**

Edit `RegionTask` at `src/pipeline/core/mod.rs:168-174`:

```rust
#[derive(Clone)]
pub struct RegionTask {
    pub index: usize,
    pub group_index: usize,
    pub record: Arc<BedRecord>,
}
```

- [ ] **Step 2: Change WorkItem to use Arc<BedRecord>**

Edit `WorkItem` at `src/pipeline/core/mod.rs:828-834`:

```rust
struct WorkItem {
    orig_idx: usize,
    group_index: usize,
    record: Arc<BedRecord>,
    query_start: i64,
    query_end: i64,
}
```

This is critical — `WorkItem` is populated from `RegionTask.record` in `execute_mode` Phase 1. After T8, both hold `Arc<BedRecord>`, so the move is cheap (ref-count bump, no clone needed).

- [ ] **Step 3: Update RegionTask construction in scale_regions.rs**

Find where `RegionTask` is constructed in `scale_regions.rs` (search for `RegionTask {`). Wrap `record` in `Arc::new(...)`:

```rust
RegionTask {
    index,
    group_index,
    record: Arc::new(record),
}
```

- [ ] **Step 4: Update RegionTask construction in reference_point.rs**

Same change — find `RegionTask {` and wrap `record: Arc::new(record)`.

- [ ] **Step 5: Update CoalescedBatch to use Arc<BedRecord>**

At `src/pipeline/core/mod.rs:884-891`:

```rust
struct CoalescedBatch {
    /// Items in original sorted order: (orig_idx, group_index, record).
    items: Vec<(usize, usize, Arc<BedRecord>)>,
    /// Start of the merged query window (minimum of all item windows).
    query_start: i64,
    /// End of the merged query window (maximum of all item windows).
    query_end: i64,
}
```

- [ ] **Step 6: Update create_batches to work with Arc<BedRecord>**

In `create_batches` at line ~893, update the Vec types:
```rust
let mut current_items: Vec<(usize, usize, Arc<BedRecord>)> = Vec::new();
```

The `item.record.chrom.clone()` becomes `Arc::clone(&item.record.chrom)` (or since `chrom` is `Arc<str>`, just `Arc::clone(&item.record.chrom)`).

- [ ] **Step 7: Update process_batch to use Arc<BedRecord>**

In `process_batch`, the pattern `for (orig_idx, group_index, record) in batch.items` now iterates `Arc<BedRecord>`. Calls like `mode.plan_for(&record, metadata)` work since `Arc<BedRecord>` derefs to `BedRecord`.

The `mode.postprocess_row(record, ...)` at line ~972 (which takes ownership of `record`) becomes:
```rust
mode.postprocess_row(Arc::unwrap_or_clone(record), flat, sc, bc, metadata)
```

- [ ] **Step 8: Update record.move into batch items**

In `core/mod.rs` around line ~920 where items are pushed:
```rust
current_items.push((item.orig_idx, item.group_index, item.record));
```
`item.record` is already `Arc<BedRecord>` (from step 2-3).

- [ ] **Step 9: Build and fix any type errors**

Run: `cargo build --release 2>&1`
Fix compilation errors. Main areas: any place that takes ownership of `BedRecord` needs `Arc::unwrap_or_clone`.

- [ ] **Step 10: Run profile_bench.sh**

```bash
./scripts/profile_bench.sh m1-arc-bedrecord \
  "Arc<BedRecord> to avoid cloning BedRecord into each RegionTask" \
  "load_groups -> RegionTask.record.clone() -> String allocations" \
  -- cargo run --release -- scale-regions -b 10 -p 4 --unscaled5prime 0 --unscaled3prime 0 --region test_data/encode_k562_atac.bed --scoreFileName test_data/sample1.bw test_data/sample2.bw test_data/sample3.bw -o /tmp/bench_m1.mat.gz
```

- [ ] **Step 11: Verify correct output and run tests**

```bash
scripts/custom_compare.py /tmp/bench_baseline.mat.gz /tmp/bench_m1.mat.gz
cargo test 2>&1
```
Expected: Peak heap reduced 50-100 MB vs baseline

- [ ] **Step 12: Commit**

```bash
git add src/pipeline/core/mod.rs src/pipeline/reference_point.rs src/pipeline/scale_regions.rs bench_reports/
git commit -m "perf: share BedRecord via Arc to avoid clone per RegionTask

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 9: Sparse Detection + CoalesceStrategy (P3)

**Files:**
- Modify: `src/pipeline/core/mod.rs:838-936` (estimate_coalesce_gap, create_batches, CoalescedBatch)
- Modify: `src/pipeline/core/mod.rs:1174-1186` (coalescing call site in execute_mode)
- Modify: `src/pipeline/scale_regions.rs:212-218` (streaming + sort_regions check)
- Modify: `src/pipeline/reference_point.rs:186-194` (streaming + sort_regions check)

- [ ] **Step 1: Define CoalesceStrategy enum**

Add at `src/pipeline/core/mod.rs` near line 838 (before `estimate_coalesce_gap`):

```rust
const COALESCE_CLAMP_MAX: i64 = 2000;

/// Strategy chosen based on estimated coalescing gap.
enum CoalesceStrategy {
    /// Normal coalescing with the given gap threshold.
    Coalesce(i64),
    /// Sparse data: skip coalescing, each work item forms its own 1-item batch.
    NoCoalesce,
}
```

- [ ] **Step 2: Modify create_batches to accept CoalesceStrategy**

Replace `create_batches` at line 896:

```rust
fn create_batches(work_items: Vec<WorkItem>, strategy: &CoalesceStrategy) -> Vec<CoalescedBatch> {
    match strategy {
        CoalesceStrategy::Coalesce(coalesce_gap) => {
            create_coalesced_batches(work_items, *coalesce_gap)
        }
        CoalesceStrategy::NoCoalesce => {
            create_per_item_batches(work_items)
        }
    }
}

fn create_per_item_batches(work_items: Vec<WorkItem>) -> Vec<CoalescedBatch> {
    work_items
        .into_iter()
        .map(|item| CoalescedBatch {
            query_start: item.query_start,
            query_end: item.query_end,
            items: vec![(item.orig_idx, item.group_index, item.record)],
        })
        .collect()
}

fn create_coalesced_batches(work_items: Vec<WorkItem>, coalesce_gap: i64) -> Vec<CoalescedBatch> {
    let mut batches = Vec::new();
    let mut current_chrom: Arc<str> = Arc::from("");
    let mut current_items: Vec<(usize, usize, Arc<BedRecord>)> = Vec::new();
    let mut batch_start: i64 = 0;
    let mut batch_end: i64 = 0;

    for item in work_items {
        if current_items.is_empty() {
            current_chrom = Arc::clone(&item.record.chrom);
            batch_start = item.query_start;
            batch_end = item.query_end;
            current_items.push((item.orig_idx, item.group_index, item.record));
        } else if item.record.chrom != current_chrom
            || item.query_start > batch_end.saturating_add(coalesce_gap)
        {
            batches.push(CoalescedBatch {
                items: std::mem::take(&mut current_items),
                query_start: batch_start,
                query_end: batch_end,
            });
            current_chrom = Arc::clone(&item.record.chrom);
            batch_start = item.query_start;
            batch_end = item.query_end;
            current_items.push((item.orig_idx, item.group_index, item.record));
        } else {
            batch_end = batch_end.max(item.query_end);
            current_items.push((item.orig_idx, item.group_index, item.record));
        }
    }

    if !current_items.is_empty() {
        batches.push(CoalescedBatch {
            items: current_items,
            query_start: batch_start,
            query_end: batch_end,
        });
    }

    batches
}
```

- [ ] **Step 3: Determine CoalesceStrategy in execute_mode**

Replace lines 1179-1186 in `src/pipeline/core/mod.rs`:

```rust
    // ── Phase 3.5: Determine coalescing strategy ────────────────────────
    let coalesce_gap = estimate_coalesce_gap(&work_items);
    let strategy = if coalesce_gap >= COALESCE_CLAMP_MAX {
        // Sparse dataset: gap at clamp ceiling means regions are far apart.
        // Skipping coalescing avoids large-window reads that pull in
        // irrelevant data between isolated regions.
        CoalesceStrategy::NoCoalesce
    } else {
        CoalesceStrategy::Coalesce(coalesce_gap)
    };
    let batches = create_batches(work_items, &strategy);
    eprintln!(
        "[coalesce-gap] strategy={:?} batches={} items={} ratio={:.2}",
        match &strategy { CoalesceStrategy::Coalesce(g) => format!("coalesce({g})"), CoalesceStrategy::NoCoalesce => "no-coalesce".into() },
        batches.len(),
        task_count,
        batches.len() as f64 / task_count as f64
    );
```

- [ ] **Step 4: Fix streaming check to allow SortRegions::No**

In `src/io/writers/mod.rs:53-77`, update `should_use_streaming_for_plan`:

```rust
pub fn should_use_streaming_for_plan(
    row_count: usize,
    sample_count: usize,
    bin_count: usize,
    sort_is_keep: bool,
    io: &IoOptions,
) -> bool {
    if io.matrix_values_output.is_some() || io.sorted_regions_output.is_some() {
        return false;
    }

    if !sort_is_keep {
        return false;
    }

    if row_count == 0 {
        return false;
    }

    let cell_count = row_count
        .saturating_mul(sample_count)
        .saturating_mul(bin_count);

    cell_count >= STREAMING_CELL_THRESHOLD
}
```

Now update callers in `scale_regions.rs:212-218` and `reference_point.rs:186-194` to pass `matches!(general.sort_regions, SortRegions::Keep | SortRegions::No)` instead of just `SortRegions::Keep`:

In `scale_regions.rs` line ~216:
```rust
matches!(general.sort_regions, SortRegions::Keep | SortRegions::No),
```

In `reference_point.rs` line ~190:
```rust
matches!(general.sort_regions, SortRegions::Keep | SortRegions::No),
```

- [ ] **Step 5: Build**

Run: `cargo build --release 2>&1`

- [ ] **Step 6: Run profile_bench.sh with sparse data**

```bash
./scripts/profile_bench.sh p3-sparse-coalesce \
  "sparse detection: skip coalescing when gap >= 2000, fix SortRegions::No streaming" \
  "create_batches -> estimate_coalesce_gap -> CoalesceStrategy decision" \
  -- cargo run --release -- scale-regions -b 10 -p 4 --unscaled5prime 0 --unscaled3prime 0 --region test_data/encode_k562_atac.bed --scoreFileName test_data/sample1.bw test_data/sample2.bw test_data/sample3.bw -o /tmp/bench_p3.mat.gz
```

- [ ] **Step 7: Also test with a dense dataset to ensure no regression**

If available, run with a dense BED file (e.g., tiling windows at 100bp). Verify coalescing still activates for dense data.

- [ ] **Step 8: Verify correct output and run tests**

```bash
scripts/custom_compare.py /tmp/bench_baseline.mat.gz /tmp/bench_p3.mat.gz
cargo test 2>&1
```

- [ ] **Step 9: Commit**

```bash
git add src/pipeline/core/mod.rs src/pipeline/scale_regions.rs src/pipeline/reference_point.rs src/io/writers/mod.rs bench_reports/
git commit -m "perf: skip coalescing for sparse datasets and fix SortRegions::No streaming

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 10: Channel Pipeline for Compute/IO Separation (AR2)

**Files:**
- Modify: `src/pipeline/core/mod.rs:1107-1250` (execute_mode — Phase 5-7 replacement)
- Modify: `src/pipeline/mod.rs:51-62` (remove spawn_writer_thread or keep for in-memory fallback)
- Modify: `src/pipeline/scale_regions.rs:222-267` (streaming + in-memory paths)
- Modify: `src/pipeline/reference_point.rs:196-280` (streaming + in-memory paths)

**WARNING: This is the largest change. It touches the execute_mode function, scale_regions, and reference_point. Builds on T9's CoalesceStrategy.**

- [ ] **Step 1: Add import for mpsc channel**

In `src/pipeline/core/mod.rs` at line 1, add:
```rust
use std::sync::mpsc;
```

- [ ] **Step 2: Define channel message type**

Add near the `WorkItem` struct at line ~827:

```rust
/// Channel message: a single result row from a processed batch, tagged with
/// its original index for reordering at the writer end.
type BatchResult = (usize, usize, Option<MatrixRow>);  // (orig_idx, group_index, row)
```

- [ ] **Step 3: Rewrite execute_mode signature and body (Phase 4-7)**

Replace `execute_mode` from line 1107 through line 1250:

```rust
pub fn execute_mode<M, C, F>(
    tasks: Vec<RegionTask>,
    general: &GeneralOptions,
    sample_paths: Arc<Vec<PathBuf>>,
    collector: C,
    thread_count: usize,
    mode: &M,
    metadata: Arc<M::Metadata>,
    header_builder: F,
    group_count: usize,
) -> Result<C::Output>
where
    M: PipelineMode,
    C: RowCollector + Send + 'static,
    F: FnOnce(Vec<usize>) -> Result<MatrixHeader> + Send + 'static,
{
    let task_count = tasks.len();

    // ── Phase 1: Pre-compute sort keys for I/O locality ──────────────────
    let mut work_items: Vec<WorkItem> = tasks
        .into_iter()
        .map(|task| {
            let plan = mode.plan_for(&task.record, metadata.as_ref());
            WorkItem {
                orig_idx: task.index,
                group_index: task.group_index,
                record: task.record,
                query_start: plan.window_start(),
                query_end: plan.window_end(),
            }
        })
        .collect();

    // ── Phase 2: Sort by (chrom, window_start, window_end) ────────────────
    work_items.sort_by(|a, b| {
        a.record
            .chrom
            .cmp(&b.record.chrom)
            .then(a.query_start.cmp(&b.query_start))
            .then(a.query_end.cmp(&b.query_end))
    });

    if task_count == 0 {
        let header = header_builder(vec![0; group_count])?;
        return collector.finalize(header);
    }

    // ── Phase 3: Open shared bigWig readers once ────────────────────────
    let shared_readers = Arc::new(
        sample_paths
            .iter()
            .map(|path| {
                SharedBigWigReader::open(path)
                    .map(Arc::new)
                    .with_context(|| {
                        format!("Failed to open bigWig file '{}'", path.display())
                    })
            })
            .collect::<Result<Vec<_>>>()?,
    );

    // ── Phase 3.5: Determine coalescing strategy ────────────────────────
    let coalesce_gap = estimate_coalesce_gap(&work_items);
    let strategy = if coalesce_gap >= COALESCE_CLAMP_MAX {
        CoalesceStrategy::NoCoalesce
    } else {
        CoalesceStrategy::Coalesce(coalesce_gap)
    };
    let batches = create_batches(work_items, &strategy);
    eprintln!(
        "[coalesce-gap] strategy={:?} batches={} items={} ratio={:.2}",
        match &strategy {
            CoalesceStrategy::Coalesce(g) => format!("coalesce({g})"),
            CoalesceStrategy::NoCoalesce => "no-coalesce".into(),
        },
        batches.len(),
        task_count,
        batches.len() as f64 / task_count as f64
    );

    // ── Phase 4: Spawn writer thread, create channel ────────────────────
    let (tx, rx) = mpsc::sync_channel::<BatchResult>(256);

    let writer_handle = std::thread::Builder::new()
        .name("matrix-writer".into())
        .spawn(move || {
            let mut next_idx: usize = 0;
            let mut pending: std::collections::BTreeMap<usize, (usize, Option<MatrixRow>)> =
                std::collections::BTreeMap::new();
            let mut collector = collector;
            let mut group_counts = vec![0usize; group_count];

            for (orig_idx, group_index, row) in rx {
                pending.insert(orig_idx, (group_index, row));
                while let Some(entry) = pending.remove(&next_idx) {
                    let (grp, row_opt) = entry;
                    if let Some(row) = row_opt {
                        collector.on_row(row)?;
                    }
                    group_counts[grp] += 1;
                    next_idx += 1;
                }
            }

            // All senders dropped — channel closed
            let header = header_builder(group_counts)?;
            collector.finalize(header)
        })
        .context("Failed to spawn writer thread")?;

    // ── Phase 5: Parallel processing over batches ──────────────────────
    let pool = ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build()
        .context("Failed to initialise rayon thread pool")?;

    let sample_paths_for_workers = Arc::clone(&sample_paths);
    let shared_for_workers = Arc::clone(&shared_readers);
    let metadata_ref = metadata.as_ref();

    pool.install(|| {
        batches
            .into_par_iter()
            .map_init(
                move || {
                    WorkerSamples::from_shared(
                        Arc::clone(&sample_paths_for_workers),
                        Arc::clone(&shared_for_workers),
                    )
                },
                |worker_samples, batch| {
                    let samples = worker_samples.samples()?;
                    let results = process_batch(
                        samples.as_mut_slice(),
                        batch,
                        mode,
                        general,
                        metadata_ref,
                    )?;
                    for result in results {
                        if tx.send(result).is_err() {
                            // Writer thread terminated — stop processing
                            break;
                        }
                    }
                    Ok(())
                },
            )
            .collect::<Vec<Result<()>>>()
    });

    // Drop the last sender so the writer thread can exit
    drop(tx);
    drop(shared_readers);

    // ── Phase 6: Collect any errors from workers ────────────────────────
    // (Errors from process_batch are propagated via Result<()> in the collect above;
    //  the for loop below would have bailed early on first error.)

    // Wait for writer thread and get its result
    match writer_handle.join() {
        Ok(Ok(output)) => Ok(output),
        Ok(Err(e)) => Err(e),
        Err(panic) => std::panic::resume_unwind(panic),
    }
}
```

- [ ] **Step 4: Update scale_regions.rs callers**

In `scale_regions.rs`, the streaming path now works differently. The `execute_mode` still takes `collector` and `header_builder`:

Streaming path (line ~234-266): Unchanged — passes `FileCollector` + `header_builder` to `execute_mode`.

In-memory path (line ~269-298): Unchanged — passes `InMemoryCollector` + `header_builder` to `execute_mode`.

But wait: `execute_mode` now spawns a writer thread internally. The `collector` is moved into the writer thread. So:
- Streaming: `FileCollector` writes to gzip inside writer thread — compute/IO overlap achieved
- In-memory: `InMemoryCollector` collects rows, `finalize` returns `MatrixData`. Then `sort_groups` is called on the returned MatrixData

The in-memory caller needs adjustment — `sort_groups` call moves after `execute_mode` returns:

```rust
let mut matrix = core::execute_mode(
    tasks, general, sample_paths, collector, thread_count,
    &mode, metadata, header_builder, group_labels.len(),
)?;
// sort_groups runs on main thread after writer thread collected all rows
let sort_sample_indices =
    core::normalize_sort_sample_indices(general.sort_using_samples.as_ref(), sample_count)?;
matrix.sort_groups(
    general.sort_regions,
    general.sort_using,
    sort_sample_indices.as_deref(),
)?;
```

This matches the existing code — no change needed.

- [ ] **Step 5: Same updates for reference_point.rs**

Apply the same pattern: `execute_mode` call stays the same, `sort_groups` call after it stays the same. The internal channel pipeline is transparent to callers.

- [ ] **Step 6: Clean up pipeline/mod.rs**

The `spawn_writer_thread` function at `src/pipeline/mod.rs:51-62` is now only used for the in-memory path's final `write_outputs`. It can stay as-is or be simplified. For now, leave it — the in-memory path still calls `spawn_writer_thread(matrix, io)` after `sort_groups`.

- [ ] **Step 7: Build**

Run: `cargo build --release 2>&1`
Expected: Compiles. Expect type errors around `CoalesceStrategy` import and `CoalescedBatch.items` type (Arc<BedRecord> from T8).

- [ ] **Step 8: Run profile_bench.sh**

```bash
./scripts/profile_bench.sh ar2-channel-pipeline \
  "channel-based compute/IO pipeline with BTreeMap reordering" \
  "execute_mode Phase 5-6: result_slots + sequential write -> sync_channel + writer thread" \
  -- cargo run --release -- scale-regions -b 10 -p 4 --unscaled5prime 0 --unscaled3prime 0 --region test_data/encode_k562_atac.bed --scoreFileName test_data/sample1.bw test_data/sample2.bw test_data/sample3.bw -o /tmp/bench_ar2.mat.gz
```

- [ ] **Step 9: Verify correct output and run tests**

```bash
scripts/custom_compare.py /tmp/bench_baseline.mat.gz /tmp/bench_ar2.mat.gz
cargo test 2>&1
```

Expected: Output identical, streaming path shows reduced peak memory (no result_slots).

- [ ] **Step 10: Commit**

```bash
git add src/pipeline/core/mod.rs src/pipeline/mod.rs src/pipeline/scale_regions.rs src/pipeline/reference_point.rs bench_reports/
git commit -m "perf: replace result_slots with channel-based compute/IO pipeline

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 11: thread_local Coverage Buffer Pool (M4)

**Files:**
- Modify: `src/pipeline/core/mod.rs:988-1033` (sample_coverages allocation in process_batch)

- [ ] **Step 1: Add thread_local buffer pool**

At the top of `process_batch` function or as a module-level item near line ~980:

```rust
use std::cell::RefCell;

thread_local! {
    /// Reusable coverage buffers per worker thread. Each entry is a
    /// sample's coverage Vec; the outer Vec is resized to sample_count
    /// and each inner Vec is resized to window_len per batch.
    static COVERAGE_POOL: RefCell<Vec<Vec<f32>>> = RefCell::new(Vec::new());
}

fn take_coverage_buffers(sample_count: usize, window_len: usize, default_fill: f32) -> Vec<Vec<f32>> {
    COVERAGE_POOL.with(|pool| {
        let mut bufs = pool.borrow_mut();
        bufs.resize_with(sample_count, Vec::new);
        for buf in bufs.iter_mut() {
            buf.clear();
            buf.resize(window_len, default_fill);
        }
        std::mem::take(&mut *bufs)
    })
}

fn return_coverage_buffers(bufs: Vec<Vec<f32>>) {
    COVERAGE_POOL.with(|pool| {
        *pool.borrow_mut() = bufs;
    });
}
```

- [ ] **Step 2: Replace sample_coverages allocation in process_batch**

Find lines 989-1033 in `src/pipeline/core/mod.rs`. Replace the `let mut sample_coverages: Vec<Vec<f32>> = Vec::with_capacity(sample_count);` and the subsequent loop that allocates per-sample `vec![default_fill; window_len]`.

Replace with:

```rust
    // ── ONE bigWig read per sample for the entire merged window ────────
    let mut sample_coverages = take_coverage_buffers(sample_count, window_len, default_fill);

    for (sample_idx, sample) in samples.iter_mut().enumerate() {
        let chrom_length = match sample.chrom_length(chrom) {
            Some(l) => l,
            None => continue,  // buffer already filled with default_fill
        };

        let fetch_start = clamp_coordinate(batch.query_start, chrom_length);
        let fetch_end = clamp_coordinate(batch.query_end, chrom_length);

        if fetch_start >= fetch_end {
            continue;  // buffer already filled with default_fill
        }

        let intervals = sample
            .reader_mut()
            .values(chrom, fetch_start, fetch_end)
            .map_err(anyhow::Error::new)
            .with_context(|| {
                format!(
                    "Failed to read bigWig intervals for '{}' in '{}'",
                    chrom,
                    sample.path().display()
                )
            })?;

        let cov = &mut sample_coverages[sample_idx];
        for v in intervals {
            let rs = i64::from(v.start)
                .saturating_sub(batch.query_start)
                .max(0);
            let re = i64::from(v.end)
                .saturating_sub(batch.query_start)
                .min(window_span)
                .max(0);
            if rs < re {
                cov[rs as usize..re as usize].fill(v.value);
            }
        }
    }

    // ── Extract per-region bins from the pre-read coverage buffers ─────
    // (existing code, unchanged)
    let nan_after_end = mode.nan_after_end(metadata);
    let mut results = Vec::with_capacity(item_count);
    // ... existing per-item extraction loop ...

    // After extraction, return buffers to pool
    return_coverage_buffers(sample_coverages);

    Ok(results)
```

Important: The extraction loop at the end uses `sample_coverages` — make sure the buffers are returned to pool AFTER the extraction loop finishes.

- [ ] **Step 3: Build**

Run: `cargo build --release 2>&1`

- [ ] **Step 4: Run profile_bench.sh**

```bash
./scripts/profile_bench.sh m4-coverage-pool \
  "thread_local coverage buffer pool to reuse Vec<Vec<f32>> across batches" \
  "process_batch -> vec![default_fill; window_len] per sample per batch" \
  -- cargo run --release -- scale-regions -b 10 -p 4 --unscaled5prime 0 --unscaled3prime 0 --region test_data/encode_k562_atac.bed --scoreFileName test_data/sample1.bw test_data/sample2.bw test_data/sample3.bw -o /tmp/bench_m4.mat.gz
```

- [ ] **Step 5: Verify correct output and run tests**

```bash
scripts/custom_compare.py /tmp/bench_baseline.mat.gz /tmp/bench_m4.mat.gz
cargo test 2>&1
```

- [ ] **Step 6: Commit**

```bash
git add src/pipeline/core/mod.rs bench_reports/
git commit -m "perf: reuse coverage buffers via thread_local pool across batches

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Orchestrator Dispatch Order

Tasks within each phase can be dispatched in parallel unless they share files:

| Phase | Parallel group | Tasks |
|-------|---------------|-------|
| 0 | — | T0 (sequential, infrastructure) |
| 1 | Group A | T1, T3 (different files: Cargo.toml/bigwig.rs vs bed.rs) |
| 1 | After A | T2 (writers/mod.rs), T4 (bwig.rs) — independent of each other |
| 2 | Group B | T5, T6, T8 — different primary files |
| 2 | After B | T7 (bwig.rs — waits for T5 which also touches bwig.rs) |
| 3 | Sequential | T9 → T10 → T11 (all touch core/mod.rs, T10 depends on T9's CoalesceStrategy) |

**File conflict resolution:**
- `bwig.rs`: T4 → T5 → T7 (sequential dispatch)
- `core/mod.rs`: T3 → T6 → T8 → T9 → T10 → T11 (T3 and T8 touch BedRecord types; later tasks need those types; T9 needs CoalesceStrategy from T9; T10 needs T9's types; T11 touches process_batch which is also modified by T10)

For simplicity, the orchestrator should dispatch Phase 2 as two sequential waves and Phase 3 fully sequentially.
