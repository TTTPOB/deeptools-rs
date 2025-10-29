# Rust `computeMatrix` Reimplementation Plan

## Objectives
- Maintain CLI parity with DeepTools `computeMatrix` (`reference-point` and `scale-regions`) including flag spelling, defaults, aliases, mutually exclusive groups, and help semantics. Any unavoidable deviation must be documented alongside user-visible implications.
- Reproduce DeepTools matrix outputs byte-for-byte (JSON header + BED-like rows) for supported options so that downstream tools (`plotHeatmap`, `plotProfile`, `computeMatrixOperations`) remain interoperable.
- Mirror Python behaviour across edge cases: strand-aware windows, NaN vs zero padding, unscaled zones, reference-point handling, sorting/threshold logic, and headers for legacy compatibility.
- Build an idiomatic Rust architecture that leverages existing crates (`bigtools` for bigWig IO, `rayon` for parallelism) while keeping the codebase modular and testable.

## Reference Python Implementation (Ground Truth)

### CLI & Parameter Normalisation (`deeptools/deeptools/computeMatrix.py`)
- Uses `argparse` with two subcommands sharing required (`-R/-S`), output, optional, and GTF option groups. `process_args()` massages negative upstream/downstream values, enforces reference-point preconditions, and expands flags into a `parameters` dict consumed later.
- CLI flags influence both runtime behaviour (e.g. `--missingDataAsZero`, `--skipZeros`, `--nanAfterEnd`) and metadata persisted in the output header. `--sortRegions keep` triggers a later BED reordering via `computeMatrixOperations.sortMatrix`.
- After parsing, `main()` creates a `parameters` dict with canonical key names (e.g. `'bin size'`, `'ref point'`, `'unscaled 5 prime'`) and instantiates `heatmapper.heatmapper`.

### Heatmapper Orchestration (`heatmapper.heatmapper.computeMatrix`)
- Performs validation that body/upstream/downstream lengths are multiples of `bin size` and that unscaled regions are only used with `scale-regions`.
- Handles GTF-specific flags (`transcriptID`, `exonID`, `keepExons`, `quiet`) and loads chromosome sizes with `getScorePerBigWigBin.getChromSizes`.
- Delegates parallel computation to `mapReduce.mapReduce`, which:
  - Splits chromosomes into ~100 kb chunks (configurable) and intersects with region batches from BED/GTF (handled by `deeptoolsintervals.GTF` to supply exon lists, group IDs, transcript labels, and strand).
  - Sends each batch to `heatmapper.compute_sub_matrix_worker`, passing the static `parameters` and full score file list.
- Collates worker results (each returns `sub_matrix`, matching region metadata, and a count of “no score” regions), sorts regions within groups, masks invalid floats via `np.ma`, assembles sample/group boundaries, and builds the `_matrix` helper. `skip zeros` optionally removes all-zero rows post assembly.

### Worker-Level Binning (`heatmapper.compute_sub_matrix_worker`)
- Opens each bigWig file for the current task and pre-computes matrix column count (`n_samples * total_bins`). Initializes `nan`-filled arrays unless `missing data as zero` is requested.
- For each transcript:
  - Derives strand-specific upstream/downstream intervals and invokes `chopRegions` or `chopRegionsFromMiddle` to carve exonic coordinates into five logical zones (`upstream`, `unscaled5'`, `body`, `unscaled3'`, `downstream`) or two zones for reference-point mode.
  - Tracks padding shortfalls (`padLeft`, `padRight`) to optionally fill with NaN or extend genomic intervals depending on `nan after end`.
  - Skips regions shorter than the bin size in scale mode (unless treated as missing data).
  - Builds a `zones` list pairing ordered genomic intervals with expected bin counts, adjusts for off-chromosome overflow (`trimZones`), and records per-zone bin totals (`a`–`e`).
- Delegates to `heatmapper.coverage_from_bigwig` / `coverage_from_array` to fetch signal values for each interval, convert to a contiguous array (NaN padded when necessary), and collapse into per-bin aggregates via `my_average` (mean/median/etc selectable). Applies min/max thresholds and early exits for all-zero rows (`skip zeros`) or regions outside thresholds.

### Coverage Extraction & Averaging
- `coverage_from_bigwig` marshals contiguous slices from `pyBigWig`, performing bounds checking, chromosome name normalisation (`change_chrom_names`), NaN padding when requests extend outside chromosome bounds, and zero-fills when `missingDataAsZero` is set.
- `coverage_from_array` bins the flattened coverage array into requested bin counts using either vectorized reshape (`reshapeZones`) or an accumulation loop to ensure identical bin counts despite shortfall padding.
- `my_average` applies numpy masked operations for the chosen statistic (mean/median/min/max/sum/std), propagating NaNs if an entire bin lacks valid data.

### Matrix Assembly & Post-processing
- `_matrix` stores the master array, region metadata (including group IDs and transcript names), and lists of group/sample boundaries; provides region-level sorting, cluster-based regrouping, silhouette computation, and submatrix retrieval.
- Sorting honours `sort using` (mean/median/max/min/sum/region_length) and `sortRegions` directives, optionally restricted to samples specified by `--sortUsingSamples`.
- `skip zeros` removes rows where all samples sum to zero, while threshold filters occur in the worker before rows are ever appended.

### Output & Auxiliary Operations
- `save_matrix` writes a gzipped file with a JSON header (prefixed by `@`) and BED-ish rows. Legacy “special_params” are homogenized to per-sample lists for backward compatibility.
- `save_matrix_values`/`save_tabulated_values` compute per-group/tabular outputs with tick labels derived from upstream/downstream/body lengths.
- `save_BED` writes sorted regions after filtering; `computeMatrixOperations.sortMatrix` can resort existing outputs when `--sortRegions keep` is used.
- Auxiliary modules (`mapReduce`, `getScorePerBigWigBin`, `heatmapper_utilities`) supply genome chunking, chromosome size lookup, tick calculation, and label management.

## Rust Architecture Plan

### Module Layout & Responsibilities
- `cli`: Clap-based argument parsing mirroring Python’s flags, aliases, and help text. Should emit a fully typed `Config` (see `src/config.rs`) with validation applied during parsing to fail fast.
- `config`: Translate CLI input into mode-specific structs (`ReferencePointOptions`, `ScaleRegionsOptions`) plus shared `GeneralOptions`/`IoOptions`. Normalise defaults (smart labels, processor counts, bin size validation) and expose derived values (bin counts, label sets).
- `io`: Provide submodules for BED/GTF parsing (`bed`, `gtf`), bigWig handling (`bigwig`), and output serialization. Implement `BedRecord`/`GtfTranscript` with group tracking, optional exon retention, and strand info. Wrap `bigtools` handles with buffered, thread-safe readers capable of slicing coverage segments.
- `pipeline`: House orchestrators for each mode (`reference_point`, `scale_regions`) plus shared helpers (`zones`, `binning`, `workers`, `matrix`). Provide an abstraction similar to Python’s `_matrix` to manage boundaries, sorting, and metadata.
- `pipeline::matrix`: Already contains `MatrixHeader`/`MatrixRow`/`MatrixData`; extend to include sorting, skip-zero filtering, and grouping utilities to mimic Python behaviour.

### Configuration Flow
1. `main` calls `cli::Cli::parse()`, then `Config::try_from()` to build typed structs.
2. `Config::execute()` dispatches to `pipeline::reference_point::run` or `pipeline::scale_regions::run`.
3. Each pipeline performs preflight validation equivalent to Python’s early exits (bin multiples, unscaled sanity, nan/threshold combinations) before scheduling work.

### Region Ingestion & Grouping
- Implement a `RegionLoader` that reads multiple BED/GTF files, handling `#` group delimiters, default group names (file stem) and unique label enforcement. For GTF mode, reuse or port logic from `deeptoolsintervals.GTF` to extract transcripts with exon lists, optional filtering by `--transcriptID`, `--exonID`, and `--keepExons`.
- Produce a `Vec<Group>` where each group holds ordered `Region` structs (`chrom`, `Vec<Interval>`, `name`, `strand`, `score`, `group_idx`), plus a global label list.

### BigWig IO Layer
- `BigWigReader` should wrap `bigtools::bigwig::BigWigRead`. Cache chromosome lengths and expose:
  - `values(chrom: &str, start: u32, end: u32) -> Result<Vec<f32>>` returning contiguous signal with NaN padding if outside bounds and optional zero fill when `missing_data_as_zero` is true.
  - `summaries(chrom, bins, aggregator)` for bin-level summaries, leveraging `bigtools::BBIRead::summarize` where possible but falling back to manual slicing to match Python semantics (especially across split exons and NaN padding).
- Manage a per-thread reader pool using `rayon::ThreadPool` thread-local storage to avoid reopening files for every region. Each worker receives a `SampleCache` containing the reader and `chrom_lengths`.

### Task Scheduling & Parallelism
- Emulate `mapReduce` by chunking transcripts by chromosome and optionally chunk size (e.g., 100 kb windows) to distribute work. Implement a `Task` struct containing `chrom`, interval span, and the subset of regions overlapping the chunk (trimmed to chunk bounds).
- Use `rayon::scope` or `par_iter` over tasks, ensuring deterministic ordering by collecting results and sorting by original group indices after parallel execution (matching Python’s `sorted(zip(groups, ...))` logic).
- Honour `GeneralOptions.processor_request` to set rayon thread pool size. Provide a sequential fallback for `proc_number == 1`.

### Zone Construction & Bin Layout
- Implement helpers mirroring `chopRegions`, `chopRegionsFromMiddle`, and `trimZones` in Rust:
  - Accept exonic interval lists and return `ZonePlan { zones: Vec<Zone>, bin_counts: Vec<usize>, pad_left/right: usize }`.
  - Support scale-regions (five zones) and reference-point variants (`TSS`, `TES`, `center`), respecting strand orientation, unscaled flanks, and `nan_after_end` semantics.
- Provide deterministic bin count calculations so that `total_bins = upstream + unscaled5 + body + unscaled3 + downstream` divided by `bin_size` matches Python even when padding occurs.

### Signal Sampling Pipeline
- For each region:
  - Build contiguous slices per zone, requesting coverage via `BigWigReader`.
  - Stitch slices into a single `Vec<f32>` (or `Vec<f64>` if precision issues appear) while inserting NaNs for missing bases and respecting padding instructions.
  - Aggregate contiguous coverage into bins using the configured statistic. Implement `AverageType` enum to mirror Python’s options, relying on manual numerics to match masked operations (NaNs ignored unless `missing_data_as_zero`).
  - Apply `min_threshold` / `max_threshold` tests across the full flattened row; skip rows with all zeros if `skip_zeros` is set. Record `no_score` counts to log warnings comparable to Python.
- Return per-sample vectors (`Vec<Vec<f32>>`) so the matrix assembly can flatten later.

### Matrix Model & Sorting
- Extend `MatrixData` with methods to:
  - Compute `group_boundaries`, `sample_boundaries`, and `sample_labels` (respecting `--samplesLabel` / `--smartLabels`).
  - Remove rows with zero sums when `skip_zeros` is true.
  - Sort groups using `sort_regions`/`sort_using`, optionally restricted to `sort_using_samples`, mimicking `_matrix.sort_groups` ordering and tie behaviour.
  - Produce derived artefacts: sorted BED rows, per-group statistics, cluster placeholders (`hmcluster`) for future parity.

### Output Serialisation & Interop
- Implement `writer::write_matrix_gz` that:
  - Serialises `MatrixHeader` as JSON prefixed with `'@'` and writes gzipped output using `flate2::write::GzEncoder`.
  - Emits BED-like lines with comma-joined start/end coordinates per exon, consistent score/strand fields, and formatted bin values using `"{:.6}"` to match Python’s `%f`.
- Support optional `outFileNameMatrix` by writing tabulated values (`averageType` per column) and `outFileSortedRegions` by emitting filtered BED entries in processing order.
- Ensure header `special_params` (`unscaled 5 prime`, `body`, `downstream`, `upstream`, `ref point`, `bin size`) are emitted as per-sample lists, even when scalars, to remain backward compatible.

### Diagnostics, Logging & Error Handling
- Mirror Python’s warning thresholds (e.g., >75% regions lacking score) using `stderr` warnings.
- Provide explicit errors for invalid input (bin multiples, mismatched chromosome names unless resolvable by `change_chrom_names` equivalent).
- Honour `--quiet` by suppressing per-region messages while still emitting fatal errors.

## Validation Strategy
- Use `pixi` to provision the DeepTools reference environment and fixtures (`deeptools/deeptools/test/test_data/`). Add reproducible scripts similar to `scripts/reference_point_regression.py` for both modes.
- Generate regression artefacts via `pixi run computeMatrix ...` and compare against Rust output (gzip contents, JSON header equality, per-value diff within tolerance). Integrate into `cargo test` as ignored tests gated by an environment variable to avoid requiring Python on every CI run.
- Add unit tests for zone splitting (`chop_regions`, `trim_zones`), coverage padding, threshold filtering, and matrix sorting to ensure deterministic behaviour.

## Performance & Future Work
- After achieving functional parity, benchmark against Python on large datasets to validate throughput. Profiling targets include bigWig fetch batching, rayon scheduling overhead, and potential columnar accumulation strategies.
- Future enhancements: implement clustering (`hmcluster`), silhouette scores, memory pooling for coverage buffers, and optional caching of bigWig blocks.

## Immediate Next Steps
- Document full CLI flag mapping (confirmed defaults, alias coverage) and ensure `config` layer encodes them.
- Finalise region ingestion abstractions (BED groups + upcoming GTF support) and map them to existing `BedRecord`.
- Flesh out zone-building helpers in Rust and cover them with unit tests using fixtures derived from Python’s `chopRegions` behaviour.
- Port the worker pipeline for `reference-point` mode end-to-end, validating against `scripts/reference_point_regression.py`.
- Establish regression harness via `pixi` to automate reference matrix generation for CI/local development.

## Priority
- Focus on BED-based `reference-point` parity first; defer full GTF handling, clustering, and `scale-regions` specific features until after the regression suite passes on BED fixtures.

## Implementation Order (Macro Level)
- [ ] `cli` + `config`: initial structs exist, but CLI parity, alias coverage, and validation gaps still need audit/patching.
- [ ] `io::bed` + grouping utilities: baseline BED parser present, yet group handling and GTF integration require completion.
- [ ] `io::bigwig`: thin wrapper in place; needs NaN/zero semantics, caching strategy, and error parity with DeepTools (currently implemented ad-hoc in the pipeline).
- [x] `pipeline::zones`: reference-point bin layout helper landed in `pipeline::zones`.
- [ ] `pipeline::reference_point`: worker now honours zone plan and average-type settings but still lacks concurrency and diagnostics.
- [ ] `io::writers`: gzip/tab outputs implemented; metadata normalisation (e.g., special param expansion) still pending despite `@` header prefix fix.
- [ ] `pipeline::matrix`: data structs exist; sorting, boundary helpers, and skip-zero removal logic unimplemented.
- [ ] Regression harness (`scripts/` + integration tests): Python comparison script exists but not wired into tests/CI.
- [ ] `pipeline::scale_regions`: not implemented.
- [ ] Advanced features (clustering, diagnostics polish) and performance tuning: pending.

## Implementation Order (Detailed within Reference-Point Milestone)
- [x] Zone helper port (`chop_regions`, `chop_regions_from_middle`, `trim_zones`): coverage windows now generated via `ReferencePointPlan`.
- [x] BigWig sampling layer: dense window reconstruction handles NaN/zero padding prior to bin aggregation.
- [x] Per-region worker: integrates zone plan, average-type aggregations, and scale factor application for each sample.
- [x] Task scheduler: rayon-based chunking not implemented; processing remains single-threaded.
- [x] Matrix assembly: needs boundary computation, skip-zero pruning, and sort hooks beyond basic struct fill.
- [x] Output serialization: header prefix now matches DeepTools; still need legacy list normalisation and value formatting review.
- [ ] Regression GLUE: script available, but cargo test/integration gating and automated diffing are outstanding.

### Task Scheduler Plan (Rayon)
- Partition the region list into deterministic chunks (target ~64 regions each) while preserving group boundaries so downstream sorting remains stable; expose chunk sizing via `Config.scheduler.chunk_size`.
- Convert the sequential region iterator in `pipeline::reference_point::execute` into a `rayon::ThreadPool`-backed `par_bridge`, yielding `TaskPayload` structs that bundle zones, sample handles, and provenance metadata.
- Provide a scoped resource manager that hands each worker thread a `Vec<BigWigCache>` built with `rayon::ThreadPool::install` to avoid `Send` conflicts and to reuse file handles across tasks.
- Ensure worker results implement `Send`/`Sync` by moving owned buffers into a `TaskResult` struct; aggregate with `rayon::iter::ParallelIterator::reduce` so the hot path stays lock-free except for a bounded `crossbeam` channel used to stream progress updates.
- Propagate errors via `Result<TaskResult>` and surface them with `rayon::join_context` so early exits cancel sibling jobs; wrap in a thin `scheduler::execute_parallel` helper to keep orchestration code testable with a single-threaded fallback.
