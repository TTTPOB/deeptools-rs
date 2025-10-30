# Streaming `matrix.gz` Investigation

## Current Behaviour Snapshot
- `reference_point::run` and `scale_regions::run` both return a `MatrixData` that owns the full `Vec<MatrixRow>` plus header metadata (`src/pipeline/reference_point.rs:160`, `src/pipeline/scale_regions.rs:128`).
- `write_outputs` currently serialises on the caller thread, feeding a `flate2::write::GzEncoder` row-by-row after sorting and zero-row pruning (`src/io/writers/mod.rs:13`–`55`).
- Because the header depends on final row counts and ordering, we keep every row in memory until all computation completes.

Peak memory therefore scales with `rows.len() * samples * bins`, matching Python’s behaviour and blocking very large inputs.

## Constraints That Block Naïve Streaming
- **Header ordering** – the JSON header embeds `group_boundaries` and `sample_boundaries`, which are only known after skip-zero/threshold pruning (and after optional sorting).
- **Sorted outputs** – `sortRegions` can reorder rows within each group; we must not stream rows until the sort keys are final, or we break byte-for-byte parity.
- **Shared consumers** – auxiliary sinks (`matrixValues`, `sortedRegions`) need the same ordered data.

## Streaming Implementation (stored-header approach)
1. Spill rows for large, unsorted workloads into a temporary plain-text file while keeping header metadata in memory. Dropping the `Vec<MatrixRow>` after spooling frees the dominant heap allocation even before gzip emission.
2. Write the gzip output as **two members** whenever the destination handle is seekable:
   - Member #1: a fixed-width header block built with `Compression::none()` so the encoder emits DEFLATE *stored* blocks. The payload is padded with spaces (with `@...\n` retained) to a reserved size of 4 096 compressed bytes (`payload = 4 096 - 23`). Because stored blocks have deterministic overhead (10-byte gzip header + 5-byte block + payload + 8-byte trailer), the member size is stable.
   - Member #2: the streamed matrix body compressed with the usual default level. Rows are replayed from the temporary file via `io::copy` into a fresh `GzBuilder::new().mtime(0)` encoder to keep timestamp-stable output.
3. After the data member finishes, seek back to byte 0 and rewrite member #1 with the *real* header (same padded length). Re-encoding the stored block updates CRC/ISIZE without disturbing the data member since the byte count is unchanged.
4. For small matrices or when auxiliary sinks are requested, fall back to the in-memory writer so we avoid the temp-file overhead.

This workflow stays fully standards-compliant—no fake padding members—while still allowing a single-pass stream for the heavy data portion.

## Header Size Guardrails
- The reserved payload (8 192 − 23 = 8 169 bytes) comfortably exceeds current header sizes; if `@header\n` ever outgrows that, we abandon the streaming path and render everything in-memory instead of risking corruption.
- JSON tolerates trailing spaces, so padding before the newline keeps downstream parsers happy.

## Implementation Notes
- A quick `seek(SeekFrom::Current(0))` probe guards against pipes/FIFOs; if the check fails we keep the legacy path.
- The temporary matrix stays as a simple newline-delimited text file so auxiliary writers can reuse it later if we add streaming support for `matrixValues`/`sortedRegions`.
- `mtime(0)` and no filename/comment yield deterministic gzip members for reproducibility.

## Compatibility Notes
- Consumers that only understand single-member gzips continue to work—the stored header is just a regular gzip member containing ASCII text.
- Split outputs (`--matrixFile`, `--matrixValues`, `--sortedRegions`) still match Python byte-for-byte because the payload copied into member #2 is identical to the in-memory writer.
