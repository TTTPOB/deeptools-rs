# Refactor Large Files & Shared Block Cache

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Improve code organization by splitting oversized modules and eliminate per-worker block cache duplication with a shared DashMap cache (500 entry cap).

**Architecture:** Pure refactoring of module boundaries (no behavioral changes) followed by replacing per-worker `HashMap<(u64,u64), Arc<[u8]>>` block caches with a single `Arc<DashMap>` shared across all workers. The shared cache uses a simple "skip insert when full" eviction policy.

**Tech Stack:** Rust, `dashmap` crate (new dependency), existing `flate2`/`rayon`/`anyhow` stack.

---

## File Structure

### Files to create (refactoring splits):
- `src/pipeline/core/traits.rs` — SignalBin, RegionPlan, PipelineMode trait definitions + ModeTag
- `src/pipeline/core/regions.rs` — Group, RegionTask, load_groups, derive_sample_labels, normalize_sort_sample_indices, BED/GTF parsing helpers
- `src/pipeline/core/executor.rs` — execute_mode, OutputStrategy, WorkItem, into_chunks, input_order_is_compute_sorted
- `src/pipeline/core/samples.rs` — Sample, WorkerSamples (bigwig handle wrappers)
- `src/pipeline/zones/mod.rs` — re-exports + non-metagene zone logic (ReferenceBin, ScaleBin, plans, helpers)
- `src/pipeline/zones/metagene.rs` — metagene module extracted verbatim
- `src/io/writers/mod.rs` — re-exports only
- `src/io/writers/matrix_gz.rs` — write_matrix_gz, write_matrix_gz_streaming, header helpers, StreamingMatrixWriter
- `src/io/writers/auxiliary.rs` — write_matrix_values, write_sorted_regions, helpers
- `src/io/writers/formatting.rs` — write_matrix_row, write_matrix_value, write_scaled_i64/i128, ROW_BUFFER
- `src/pipeline/run.rs` — shared generic run_pipeline function

### Files to modify:
- `src/pipeline/core/mod.rs` — becomes thin re-export hub
- `src/pipeline/zones.rs` → replaced by `src/pipeline/zones/mod.rs` + `metagene.rs`
- `src/pipeline/reference_point.rs` — delegates to run_pipeline
- `src/pipeline/scale_regions.rs` — delegates to run_pipeline
- `src/io/readers/bwig.rs` — remove per-worker block_cache, accept shared cache
- `src/pipeline/mod.rs` — add `mod run;`
- `Cargo.toml` — add `dashmap` dependency

### Files to create (shared cache):
- `src/io/readers/block_cache.rs` — SharedBlockCache type wrapping DashMap

---

## Task 1: Split `pipeline/core/mod.rs` into submodules

**Files:**
- Create: `src/pipeline/core/traits.rs`
- Create: `src/pipeline/core/regions.rs`
- Create: `src/pipeline/core/executor.rs`
- Create: `src/pipeline/core/samples.rs`
- Modify: `src/pipeline/core/mod.rs`

- [ ] **Step 1: Create `traits.rs` with trait definitions**

Move lines 16-65 and 415-441 from `core/mod.rs` into `traits.rs`:

```rust
// src/pipeline/core/traits.rs
use std::fmt;
use std::sync::Arc;

use anyhow::{Result, bail};

use crate::config::GeneralOptions;
use crate::io::BedRecord;
use crate::pipeline::matrix::{MatrixHeader, MatrixRow};

pub trait SignalBin {
    fn start(&self) -> i64;
    fn end(&self) -> i64;
    fn beyond_region(&self) -> bool;
}

#[derive(Debug, Clone, Copy)]
pub enum ModeTag {
    ReferencePoint,
    ScaleRegions,
}

impl fmt::Display for ModeTag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let label = match self {
            ModeTag::ReferencePoint => "reference-point",
            ModeTag::ScaleRegions => "scale-regions",
        };
        f.write_str(label)
    }
}

pub fn ensure_positive(value: u32, flag: &str, mode: ModeTag) -> Result<()> {
    if value == 0 {
        bail!("[{mode}] {flag} must be a positive integer");
    }
    Ok(())
}

pub fn ensure_multiple(bin_size: u32, distance: u32, flag: &str, mode: ModeTag) -> Result<()> {
    if distance % bin_size != 0 {
        bail!("[{mode}] {flag} ({distance}) must be a multiple of the bin size ({bin_size})");
    }
    Ok(())
}

pub trait RegionPlan {
    type Bin: SignalBin;

    fn window_start(&self) -> i64;
    fn window_end(&self) -> i64;
    fn bins(&self) -> &[Self::Bin];

    fn included_intervals(&self) -> Option<&[(i64, i64)]> {
        None
    }
}

pub trait PipelineMode: Sync {
    type Plan: RegionPlan;
    type Metadata: Send + Sync;

    fn validate(&self, general: &GeneralOptions) -> Result<Self::Metadata>;
    fn total_bins(&self, metadata: &Self::Metadata) -> usize;
    fn plan_for(&self, record: &BedRecord, metadata: &Self::Metadata) -> Self::Plan;
    fn nan_after_end(&self, metadata: &Self::Metadata) -> bool;
    fn postprocess_row(
        &self,
        record: BedRecord,
        values: Vec<f32>,
        sample_count: usize,
        bin_count: usize,
        metadata: &Self::Metadata,
    ) -> MatrixRow;
    fn build_header(
        &self,
        general: &GeneralOptions,
        metadata: &Self::Metadata,
        sample_labels: &[String],
        group_labels: &[String],
        group_counts: &[usize],
        thread_count: usize,
        sample_count: usize,
    ) -> MatrixHeader;
}
```

- [ ] **Step 2: Create `samples.rs` with Sample/WorkerSamples**

Move lines 67-154 from `core/mod.rs`:

```rust
// src/pipeline/core/samples.rs
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context, Result, anyhow};

use crate::io::{BigWigReader, SharedBigWigReader};

pub struct Sample {
    path: PathBuf,
    reader: BigWigReader,
    chrom_lengths: HashMap<String, u32>,
}

impl Sample {
    pub fn open(path: &Path) -> Result<Self> {
        let reader = BigWigReader::open(path)
            .with_context(|| format!("Failed to open bigWig file '{}'", path.display()))?;
        let chrom_lengths = reader
            .chroms()
            .iter()
            .map(|chrom| (chrom.name.clone(), chrom.length))
            .collect::<HashMap<_, _>>();

        Ok(Self {
            path: path.to_path_buf(),
            reader,
            chrom_lengths,
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
        let chrom_lengths = shared
            .chroms()
            .iter()
            .map(|chrom| (chrom.name.clone(), chrom.length))
            .collect::<HashMap<_, _>>();
        Self {
            path,
            reader: BigWigReader::from_shared(shared),
            chrom_lengths,
        }
    }

    pub fn chrom_length(&self, chrom: &str) -> Option<u32> {
        self.chrom_lengths.get(chrom).copied()
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

- [ ] **Step 3: Create `regions.rs` with group loading logic**

Move lines 156-413 from `core/mod.rs`:

```rust
// src/pipeline/core/regions.rs
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context, Result, bail};

use crate::config::{GeneralOptions, GtfOptions};
use crate::io::{BedReadError, BedRecord, load_gtf_records};

pub struct Group {
    pub label: String,
    pub records: Vec<BedRecord>,
}

#[derive(Clone)]
pub struct RegionTask {
    pub index: usize,
    pub group_index: usize,
    pub record: Arc<BedRecord>,
}

pub fn load_groups(paths: &[PathBuf], gtf: &GtfOptions) -> Result<Vec<Group>> {
    // ... (exact content from lines 174-201 of current core/mod.rs)
    todo!("move verbatim from core/mod.rs")
}

pub fn derive_sample_labels(paths: &[PathBuf], general: &GeneralOptions) -> Result<Vec<String>> {
    // ... (exact content from lines 204-220)
    todo!("move verbatim from core/mod.rs")
}

pub fn normalize_sort_sample_indices(
    raw: Option<&Vec<usize>>,
    sample_count: usize,
) -> Result<Option<Vec<usize>>> {
    // ... (exact content from lines 222-247)
    todo!("move verbatim from core/mod.rs")
}

// All private helpers: parse_grouped_bed, parse_grouped_gtf, finalize_group,
// next_unique_label, label_from_path, bed_file_label, infer_region_format, RegionFormat
// ... (lines 258-413, move verbatim)
```

- [ ] **Step 4: Create `executor.rs` with execute_mode**

Move lines 456-772 from `core/mod.rs`:

```rust
// src/pipeline/core/executor.rs
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{Context, Result};
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

use crate::config::{GeneralOptions, SortRegions};
use crate::io::{BedRecord, SharedBigWigReader};
use crate::pipeline::matrix::{MatrixHeader, MatrixRow};

use super::coalesce::{CoalesceStrategy, COALESCE_CLAMP_MAX, create_batches, estimate_coalesce_gap};
use super::collector::{GroupBucketCollector, RowCollector};
use super::regions::RegionTask;
use super::samples::WorkerSamples;
use super::traits::PipelineMode;
use super::worker::process_batch;

enum OutputStrategy {
    StreamOrdered,
    InMemoryKeep,
    InMemoryGroupBucket,
}

fn into_chunks<T>(items: Vec<T>, chunk_size: usize) -> Vec<Vec<T>> {
    // ... (verbatim from lines 468-479)
    todo!("move verbatim")
}

type BatchResult = (usize, usize, Option<MatrixRow>);

pub(crate) struct WorkItem {
    pub(crate) orig_idx: usize,
    pub(crate) group_index: usize,
    pub(crate) record: Arc<BedRecord>,
    pub(crate) query_start: i64,
    pub(crate) query_end: i64,
}

fn input_order_is_compute_sorted(items: &[WorkItem]) -> bool {
    // ... (verbatim from lines 492-499)
    todo!("move verbatim")
}

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
    // ... (verbatim from lines 501-772)
    todo!("move verbatim")
}
```

- [ ] **Step 5: Rewrite `core/mod.rs` as re-export hub**

```rust
// src/pipeline/core/mod.rs
pub mod traits;
pub mod samples;
pub mod regions;
mod executor;
mod collector;
mod coalesce;
mod worker;

pub use traits::{SignalBin, RegionPlan, PipelineMode, ModeTag, ensure_positive, ensure_multiple};
pub use samples::{Sample, WorkerSamples};
pub use regions::{Group, RegionTask, load_groups, derive_sample_labels, normalize_sort_sample_indices};
pub use executor::execute_mode;
pub use collector::{RowCollector, InMemoryCollector, FileCollector, GroupBucketCollector};
pub use worker::compute_row;

#[cfg(test)]
mod tests {
    // move existing tests here — they test cross-module integration
    // ... (verbatim from lines 774-897)
}
```

- [ ] **Step 6: Run tests to verify**

Run: `cargo test`
Expected: All existing tests pass with no behavioral change.

- [ ] **Step 7: Commit**

```bash
git add src/pipeline/core/
git commit -m "refactor: split pipeline/core/mod.rs into traits, samples, regions, executor submodules"
```

---

## Task 2: Split `pipeline/zones.rs` into `zones/mod.rs` + `zones/metagene.rs`

**Files:**
- Create: `src/pipeline/zones/mod.rs`
- Create: `src/pipeline/zones/metagene.rs`
- Delete: `src/pipeline/zones.rs`

- [ ] **Step 1: Create directory and move metagene module**

Create `src/pipeline/zones/metagene.rs` with the content of the current `mod metagene { ... }` block (lines 441-1119 of zones.rs), removing the `mod metagene {` wrapper and making it a standalone module:

```rust
// src/pipeline/zones/metagene.rs
use super::{
    ReferenceBin, ScaleBin, ScaleRegionsPlan, collect_window_bounds,
    intervals_to_bins, intervals_total_length,
};
use crate::config::ScaleRegionsOptions;
use crate::io::{BedRecord, Strand};
use crate::config::ReferencePoint;

pub fn reference_bins(
    // ... (exact content from metagene::reference_bins, lines 449-503)
) -> Option<Vec<ReferenceBin>> {
    todo!("move verbatim")
}

pub fn scale_bins(
    // ... (exact content from metagene::scale_bins, lines 506-573)
) -> Option<ScaleRegionsPlan> {
    todo!("move verbatim")
}

// All private helpers: build_tss, build_tes, build_center,
// append_reference_bins, build_scale_positive, build_scale_negative,
// append_scale_bins, take_from_start, take_from_end, chop_regions,
// chop_regions_from_middle
// ... (move verbatim)
```

- [ ] **Step 2: Create `zones/mod.rs` with remaining logic**

```rust
// src/pipeline/zones/mod.rs
pub(crate) mod metagene;

// Everything from zones.rs EXCEPT the `mod metagene { ... }` block:
// - use statements (lines 1-8)
// - ReferenceBin, ReferencePointPlan, impls (lines 10-93)
// - ScaleBin, ScaleRegionsPlan, impls (lines 96-234)
// - build_bins, reference_coordinate, bin_boundaries, bin_beyond_region, append_bins (lines 236-353)
// - coordinate_from_offset, intervals_total_length, intervals_to_bins, collect_window_bounds (lines 355-439)
// - tests module (lines 1121-1306)

// Make intervals_total_length and related helpers pub(crate) so metagene.rs can use them
pub(crate) fn intervals_total_length(intervals: &[(i64, i64)]) -> i64 { ... }
pub(crate) fn intervals_to_bins(intervals: &[(i64, i64)], bin_count: usize) -> Vec<(i64, i64, bool)> { ... }
pub(crate) fn collect_window_bounds(intervals: &[(i64, i64)], window: (i64, i64)) -> (i64, i64) { ... }
```

- [ ] **Step 3: Update `pipeline/mod.rs` to use directory module**

No change needed — `pub mod zones;` resolves to `zones/mod.rs` automatically when the directory exists.

- [ ] **Step 4: Delete old `src/pipeline/zones.rs`**

```bash
rm src/pipeline/zones.rs
```

- [ ] **Step 5: Run tests to verify**

Run: `cargo test`
Expected: All tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/pipeline/zones/ && git rm src/pipeline/zones.rs
git commit -m "refactor: split pipeline/zones into zones/mod.rs and zones/metagene.rs"
```

---

## Task 3: Split `io/writers/mod.rs` into submodules

**Files:**
- Create: `src/io/writers/formatting.rs`
- Create: `src/io/writers/matrix_gz.rs`
- Create: `src/io/writers/auxiliary.rs`
- Modify: `src/io/writers/mod.rs`

- [ ] **Step 1: Create `formatting.rs`**

Move the row serialization logic (lines 317-598 of current writers/mod.rs):

```rust
// src/io/writers/formatting.rs
use std::cell::RefCell;
use std::io::{self, Write};

use anyhow::{Context, Result};
use itoa::Buffer;

use crate::pipeline::matrix::MatrixRow;

thread_local! {
    static ROW_BUFFER: RefCell<Vec<u8>> = RefCell::new(Vec::with_capacity(32768));
}

pub fn write_matrix_row<W: Write>(writer: &mut W, row: &MatrixRow) -> Result<()> {
    // ... (verbatim from lines 317-383)
    todo!("move verbatim")
}

pub fn write_matrix_value<W: Write>(writer: &mut W, value: f32) -> io::Result<()> {
    // ... (verbatim from lines 503-531)
    todo!("move verbatim")
}

#[inline]
fn write_scaled_i64<W: Write>(writer: &mut W, scaled: i64) -> io::Result<()> {
    // ... (verbatim)
    todo!("move verbatim")
}

#[inline]
fn write_scaled_i128<W: Write>(writer: &mut W, scaled: i128) -> io::Result<()> {
    // ... (verbatim)
    todo!("move verbatim")
}

fn format_plain_value(value: f32) -> String {
    // ... (verbatim)
    todo!("move verbatim")
}

pub fn write_plain_row<W: Write>(writer: &mut W, row: &MatrixRow) -> Result<()> {
    // ... (verbatim from lines 434-446)
    todo!("move verbatim")
}

#[cfg(test)]
mod tests {
    // ... (verbatim from lines 600-669, the formatting tests)
}
```

- [ ] **Step 2: Create `matrix_gz.rs`**

Move streaming/one-shot gz writers (lines 15-298):

```rust
// src/io/writers/matrix_gz.rs
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Seek, SeekFrom, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};
use flate2::write::GzEncoder;
use flate2::{Compression, GzBuilder};
use tempfile::{NamedTempFile, TempPath};

use crate::pipeline::matrix::{MatrixData, MatrixHeader, MatrixRow};

use super::formatting::write_matrix_row;

pub const STREAMING_CELL_THRESHOLD: usize = 100_000;
const RESERVED_HEADER_COMPRESSED: usize = 8192;
const RESERVED_HEADER_PAYLOAD: usize = RESERVED_HEADER_COMPRESSED - 23;

pub fn write_matrix_gz(path: &Path, matrix: &MatrixData) -> Result<()> { ... }
pub fn write_matrix_gz_streaming(path: &Path, matrix: &mut MatrixData) -> Result<()> { ... }
pub fn build_padded_header_payload(header: &MatrixHeader) -> Result<Vec<u8>> { ... }
pub fn ensure_streaming_header_capacity(header: &MatrixHeader) -> Result<()> { ... }

pub struct StreamingMatrixWriter { ... }
impl StreamingMatrixWriter {
    pub fn start(path: &Path) -> Result<Self> { ... }
    pub fn write_row(&mut self, row: &MatrixRow) -> Result<()> { ... }
    pub fn finish(self, header: &MatrixHeader) -> Result<()> { ... }
    pub fn abort(self) { ... }
}

// Private helpers: build_header_line_from_header, pad_header_payload,
// placeholder_header_payload, write_header_member, rewrite_header_member,
// write_header_line, spool_rows
```

- [ ] **Step 3: Create `auxiliary.rs`**

Move matrix_values and sorted_regions writers (lines 385-500):

```rust
// src/io/writers/auxiliary.rs
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};

use crate::pipeline::matrix::{MatrixData, MatrixRow};

use super::formatting::write_plain_row;

pub fn write_matrix_values(path: &Path, matrix: &MatrixData) -> Result<()> { ... }
pub fn write_sorted_regions(path: &Path, matrix: &MatrixData) -> Result<()> { ... }

// Private helpers: write_matrix_values_header, group_label_for_index, diff
```

- [ ] **Step 4: Rewrite `writers/mod.rs` as re-export hub**

```rust
// src/io/writers/mod.rs
mod formatting;
mod matrix_gz;
mod auxiliary;

use anyhow::Result;

use crate::config::{IoOptions, SortRegions};
use crate::pipeline::matrix::MatrixData;

pub use matrix_gz::{
    StreamingMatrixWriter, build_padded_header_payload, ensure_streaming_header_capacity,
    STREAMING_CELL_THRESHOLD,
};

pub fn should_use_streaming_for_plan(
    row_count: usize,
    sample_count: usize,
    bin_count: usize,
    sort_regions: SortRegions,
    io: &IoOptions,
) -> bool {
    if io.matrix_values_output.is_some() || io.sorted_regions_output.is_some() {
        return false;
    }
    if !matches!(sort_regions, SortRegions::Keep | SortRegions::No) {
        return false;
    }
    if row_count == 0 {
        return false;
    }
    let cell_count = row_count.saturating_mul(sample_count).saturating_mul(bin_count);
    cell_count >= STREAMING_CELL_THRESHOLD
}

pub fn write_outputs(mut matrix: MatrixData, io: &IoOptions) -> Result<()> {
    let sort_regions_str = &matrix.header.sort_regions;
    let sort_ok = sort_regions_str == "keep" || sort_regions_str == "no";
    let should_stream = should_use_streaming_for_plan(
        matrix.rows.len(),
        matrix.sample_count,
        matrix.bin_count,
        if sort_ok { SortRegions::Keep } else { SortRegions::Descend },
        io,
    );

    if should_stream {
        matrix_gz::write_matrix_gz_streaming(&io.matrix_output, &mut matrix)?;
        return Ok(());
    }

    matrix_gz::write_matrix_gz(&io.matrix_output, &matrix)?;

    if let Some(path) = &io.matrix_values_output {
        auxiliary::write_matrix_values(path, &matrix)?;
    }
    if let Some(path) = &io.sorted_regions_output {
        auxiliary::write_sorted_regions(path, &matrix)?;
    }

    Ok(())
}
```

- [ ] **Step 5: Run tests to verify**

Run: `cargo test`
Expected: All tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/io/writers/
git commit -m "refactor: split io/writers into formatting, matrix_gz, auxiliary submodules"
```

---

## Task 4: Extract shared `run_pipeline` from reference_point/scale_regions

**Files:**
- Create: `src/pipeline/run.rs`
- Modify: `src/pipeline/mod.rs`
- Modify: `src/pipeline/reference_point.rs`
- Modify: `src/pipeline/scale_regions.rs`

- [ ] **Step 1: Create `src/pipeline/run.rs` with generic run function**

Both `reference_point::run` and `scale_regions::run` follow the same pattern. Extract the shared skeleton:

```rust
// src/pipeline/run.rs
use std::sync::Arc;

use anyhow::Result;

use crate::config::{GeneralOptions, GtfOptions, IoOptions};
use crate::io::writers;
use crate::pipeline::core::{
    self, FileCollector, InMemoryCollector, PipelineMode, RegionTask,
};
use crate::pipeline::matrix::MatrixHeader;

use super::RunOutcome;

pub fn run_pipeline<M>(
    mode: M,
    general: &GeneralOptions,
    io: &IoOptions,
    gtf: &GtfOptions,
) -> Result<RunOutcome>
where
    M: PipelineMode + Clone + Send + 'static,
    M::Metadata: Clone + 'static,
{
    let metadata = Arc::new(mode.validate(general)?);

    let sample_labels = core::derive_sample_labels(&io.scores, general)?;
    let sample_count = sample_labels.len();

    let groups = core::load_groups(&io.regions, gtf)?;
    let group_labels: Vec<String> = groups.iter().map(|g| g.label.clone()).collect();
    let group_capacity: Vec<usize> = groups.iter().map(|g| g.records.len()).collect();

    let mut tasks = Vec::new();
    for (group_index, group) in groups.into_iter().enumerate() {
        for record in group.records {
            let index = tasks.len();
            tasks.push(RegionTask {
                index,
                group_index,
                record: Arc::new(record),
            });
        }
    }

    let thread_count = std::cmp::max(1, general.number_of_processors.resolve() as usize);
    let total_bins = mode.total_bins(metadata.as_ref());
    let row_count = tasks.len();
    let should_stream = writers::should_use_streaming_for_plan(
        row_count,
        sample_count,
        total_bins,
        general.sort_regions,
        io,
    );

    let sample_paths = Arc::new(io.scores.clone());

    if should_stream {
        let header_estimate = mode.build_header(
            general,
            metadata.as_ref(),
            &sample_labels,
            &group_labels,
            &group_capacity,
            thread_count,
            sample_count,
        );
        writers::ensure_streaming_header_capacity(&header_estimate)?;

        let writer = writers::StreamingMatrixWriter::start(&io.matrix_output)?;
        let collector = FileCollector::new(writer);
        let header_builder = {
            let general = general.clone();
            let sample_labels = sample_labels.clone();
            let group_labels = group_labels.clone();
            let metadata = Arc::clone(&metadata);
            let mode = mode.clone();
            move |group_counts: Vec<usize>| -> Result<MatrixHeader> {
                Ok(mode.build_header(
                    &general,
                    metadata.as_ref(),
                    &sample_labels,
                    &group_labels,
                    &group_counts,
                    thread_count,
                    sample_count,
                ))
            }
        };

        core::execute_mode(
            tasks,
            general,
            Arc::clone(&sample_paths),
            collector,
            thread_count,
            &mode,
            Arc::clone(&metadata),
            header_builder,
            group_labels.len(),
        )?;
        return Ok(RunOutcome::Streamed);
    }

    let collector = InMemoryCollector::with_capacity(row_count, sample_count, total_bins);
    let header_builder = {
        let general = general.clone();
        let sample_labels = sample_labels.clone();
        let group_labels = group_labels.clone();
        let metadata = Arc::clone(&metadata);
        let mode = mode.clone();
        move |group_counts: Vec<usize>| -> Result<MatrixHeader> {
            Ok(mode.build_header(
                &general,
                metadata.as_ref(),
                &sample_labels,
                &group_labels,
                &group_counts,
                thread_count,
                sample_count,
            ))
        }
    };

    let mut matrix = core::execute_mode(
        tasks,
        general,
        sample_paths,
        collector,
        thread_count,
        &mode,
        metadata,
        header_builder,
        group_labels.len(),
    )?;

    let sort_sample_indices =
        core::normalize_sort_sample_indices(general.sort_using_samples.as_ref(), sample_count)?;

    matrix.sort_groups(
        general.sort_regions,
        general.sort_using,
        sort_sample_indices.as_deref(),
    )?;

    matrix.prune_zero_rows();

    Ok(RunOutcome::Matrix(matrix))
}
```

- [ ] **Step 2: Update `pipeline/mod.rs`**

Add the new module:

```rust
// Add after existing mod declarations:
mod run;
pub(crate) use run::run_pipeline;
```

- [ ] **Step 3: Simplify `reference_point.rs`**

Replace the `run` function body with a call to `run_pipeline`:

```rust
// src/pipeline/reference_point.rs
use anyhow::Result;

use crate::config::{GeneralOptions, GtfOptions, IoOptions, ReferencePointOptions};

use super::RunOutcome;
use super::run_pipeline;

// Keep ReferencePointMode struct and PipelineMode impl as-is (lines 18-152)
// ...

pub fn run(
    general: &GeneralOptions,
    io: &IoOptions,
    gtf: &GtfOptions,
    options: &ReferencePointOptions,
) -> Result<RunOutcome> {
    let mode = ReferencePointMode::new(options.clone(), gtf.keep_exons);
    run_pipeline(mode, general, io, gtf)
}
```

- [ ] **Step 4: Simplify `scale_regions.rs`**

Same treatment:

```rust
pub fn run(
    general: &GeneralOptions,
    io: &IoOptions,
    gtf: &GtfOptions,
    options: &ScaleRegionsOptions,
) -> Result<RunOutcome> {
    let mode = ScaleRegionsMode::new(options.clone(), gtf.keep_exons);
    run_pipeline(mode, general, io, gtf)
}
```

- [ ] **Step 5: Run tests to verify**

Run: `cargo test`
Expected: All tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/pipeline/run.rs src/pipeline/mod.rs src/pipeline/reference_point.rs src/pipeline/scale_regions.rs
git commit -m "refactor: extract shared run_pipeline to eliminate duplicate run() logic"
```

---

## Task 5: Implement shared block cache with DashMap

**Files:**
- Create: `src/io/readers/block_cache.rs`
- Modify: `src/io/readers/bwig.rs`
- Modify: `src/io/readers/mod.rs`
- Modify: `Cargo.toml`

- [ ] **Step 1: Add `dashmap` dependency**

```toml
# In Cargo.toml [dependencies], add:
dashmap = "6"
```

- [ ] **Step 2: Create `block_cache.rs`**

```rust
// src/io/readers/block_cache.rs
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use dashmap::DashMap;

const MAX_SHARED_BLOCK_CACHE_ENTRIES: usize = 500;

pub struct SharedBlockCache {
    map: DashMap<(u64, u64), Arc<[u8]>>,
    len: AtomicUsize,
}

impl SharedBlockCache {
    pub fn new() -> Self {
        Self {
            map: DashMap::with_capacity(MAX_SHARED_BLOCK_CACHE_ENTRIES),
            len: AtomicUsize::new(0),
        }
    }

    pub fn get(&self, key: &(u64, u64)) -> Option<Arc<[u8]>> {
        self.map.get(key).map(|entry| Arc::clone(entry.value()))
    }

    pub fn insert(&self, key: (u64, u64), value: Arc<[u8]>) {
        if self.len.load(Ordering::Relaxed) >= MAX_SHARED_BLOCK_CACHE_ENTRIES {
            return;
        }
        if self.map.contains_key(&key) {
            return;
        }
        self.map.insert(key, value);
        self.len.fetch_add(1, Ordering::Relaxed);
    }
}
```

- [ ] **Step 3: Update `readers/mod.rs` to export SharedBlockCache**

```rust
// src/io/readers/mod.rs
pub mod bed;
pub mod bwig;
pub mod gtf;
pub mod block_cache;

pub use bed::{BedReadError, BedReader, BedRecord, Strand};
pub use bwig::{BigWigReadError, BigWigReader, BigWigValue, ChromInfo, SharedBigWigReader};
pub use block_cache::SharedBlockCache;
pub use gtf::load_gtf_records;
```

- [ ] **Step 4: Modify `SharedBigWigReader` to hold shared cache**

In `src/io/readers/bwig.rs`, add `block_cache: Arc<SharedBlockCache>` to `SharedBigWigReader`:

```rust
use super::block_cache::SharedBlockCache;

pub struct SharedBigWigReader {
    file: File,
    uncompress_buf_size: usize,
    chroms: Vec<ChromInfo>,
    chrom_id_by_name: Vec<(String, u32)>,
    cir_tree_root: u64,
    block_cache: Arc<SharedBlockCache>,
}

impl SharedBigWigReader {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, BigWigReadError> {
        Self::open_with_cache(path, Arc::new(SharedBlockCache::new()))
    }

    pub fn open_with_cache(
        path: impl AsRef<Path>,
        block_cache: Arc<SharedBlockCache>,
    ) -> Result<Self, BigWigReadError> {
        // ... existing parsing logic ...
        Ok(Self {
            file,
            uncompress_buf_size,
            chroms,
            chrom_id_by_name,
            cir_tree_root,
            block_cache,
        })
    }

    pub fn block_cache(&self) -> &Arc<SharedBlockCache> {
        &self.block_cache
    }

    // ... rest unchanged
}
```

- [ ] **Step 5: Remove per-worker block_cache from BigWigReader**

Replace the per-worker `block_cache: HashMap<(u64, u64), Arc<[u8]>>` with a reference to the shared cache:

```rust
pub struct BigWigReader {
    shared: Arc<SharedBigWigReader>,
    cir_node_cache: HashMap<u64, Arc<CachedCirNode>>,
    // REMOVED: block_cache: HashMap<(u64, u64), Arc<[u8]>>,
    work_buf: Vec<u8>,
    decode_buf: Vec<u8>,
    values_buf: Vec<BigWigValue>,
    blocks_buf: Vec<Block>,
    remaining_buf: VecDeque<u64>,
    pub stats: BigWigReaderStats,
}

impl BigWigReader {
    pub fn from_shared(shared: Arc<SharedBigWigReader>) -> Self {
        let uncompress_buf_size = shared.uncompress_buf_size;
        Self {
            shared,
            cir_node_cache: HashMap::new(),
            // REMOVED: block_cache
            work_buf: Vec::with_capacity(uncompress_buf_size),
            decode_buf: Vec::new(),
            values_buf: Vec::new(),
            blocks_buf: Vec::new(),
            remaining_buf: VecDeque::new(),
            stats: BigWigReaderStats::default(),
        }
    }

    fn get_or_cache_block(
        &mut self,
        offset: u64,
        size: u64,
    ) -> io::Result<Arc<[u8]>> {
        let key = (offset, size);
        let cache = &self.shared.block_cache;

        if let Some(data) = cache.get(&key) {
            self.stats.block_cache_hits += 1;
            return Ok(data);
        }

        self.stats.block_cache_misses += 1;
        let raw = read_and_decompress(
            &self.shared.file,
            offset,
            size,
            &mut self.work_buf,
            &mut self.decode_buf,
        )?;
        self.stats.decoded_bytes += raw.len() as u64;
        let data: Arc<[u8]> = Arc::from(raw);

        if !data.is_empty() {
            cache.insert(key, Arc::clone(&data));
        }

        Ok(data)
    }
}
```

- [ ] **Step 6: Update `execute_mode` to pass a shared cache**

In `executor.rs` (or wherever `shared_readers` are opened), create one `SharedBlockCache` and pass it to all readers:

```rust
// In execute_mode, Phase 3:
let block_cache = Arc::new(SharedBlockCache::new());
let shared_readers = Arc::new(
    sample_paths
        .iter()
        .map(|path| {
            SharedBigWigReader::open_with_cache(path, Arc::clone(&block_cache))
                .map(Arc::new)
                .with_context(|| {
                    format!("Failed to open bigWig file '{}'", path.display())
                })
        })
        .collect::<Result<Vec<_>>>()?,
);
```

- [ ] **Step 7: Remove `block_cache_clears` from stats**

Update `BigWigReaderStats` to remove `block_cache_clears` (no longer applicable):

```rust
#[derive(Debug, Default)]
pub struct BigWigReaderStats {
    pub values_calls: u64,
    pub values_returned: u64,
    pub blocks_per_query_total: u64,
    pub cir_cache_hits: u64,
    pub cir_cache_misses: u64,
    pub cir_cache_clears: u64,
    pub block_cache_hits: u64,
    pub block_cache_misses: u64,
    // REMOVED: block_cache_clears
    pub decoded_bytes: u64,
}
```

- [ ] **Step 8: Run tests to verify**

Run: `cargo test`
Expected: All tests pass.

- [ ] **Step 9: Run `cargo clippy` to check for warnings**

Run: `cargo clippy -- -D warnings`
Expected: No errors.

- [ ] **Step 10: Commit**

```bash
git add Cargo.toml src/io/readers/
git commit -m "feat: replace per-worker block cache with shared DashMap cache (500 entry cap)"
```

---

## Task 6: Update `io/mod.rs` exports and final cleanup

**Files:**
- Modify: `src/io/mod.rs`

- [ ] **Step 1: Add SharedBlockCache to io re-exports if needed externally**

```rust
// src/io/mod.rs
pub mod readers;
pub mod writers;

pub use readers::bed::{BedReadError, BedReader, BedRecord, Strand};
pub use readers::bwig::{BigWigReadError, BigWigReader, BigWigValue, ChromInfo, SharedBigWigReader};
pub use readers::block_cache::SharedBlockCache;
pub use readers::gtf::load_gtf_records;
```

- [ ] **Step 2: Full test + build check**

Run: `cargo test && cargo build --release`
Expected: All pass, release build succeeds.

- [ ] **Step 3: Commit any final fixups**

```bash
git add -A
git commit -m "refactor: final cleanup after module restructuring"
```
