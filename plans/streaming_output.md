# Streaming `matrix.gz` Investigation

## Current Behaviour Snapshot
- `reference_point::run` and `scale_regions::run` now return a `RunOutcome`, choosing between a fully materialised `MatrixData` and a streaming write. When `RunOutcome::Streamed` is selected, rows are flushed to disk while computation is still in flight instead of waiting for the matrix to finish (`src/pipeline/reference_point.rs`, `src/pipeline/scale_regions.rs`).
- Streaming is gated up front on three conditions: `sortRegions keep`, no auxiliary outputs (`matrixValues`, `sortedRegions`), and an estimated cell count of at least 100 000 (`writers::should_use_streaming_for_plan`). Other workloads continue to use the in-memory writer.
- For streaming runs, the output `matrix.gz` grows immediately after rayon workers start. Ordering is preserved via per-row indices, and group counts are tracked on the fly for the final header.

## Streaming Architecture
1. The mode runner builds a task list (one per BED record) and computes per-group capacities. A tentative header using those capacities is assembled and passed through `writers::ensure_streaming_header_capacity` to confirm it fits the 8 169-byte reserved payload.
2. `StreamingMatrixWriter::start` opens the destination, probes for `Seek` support, writes a placeholder stored-member header, and exposes a gzip encoder for the data member.
3. Rayon workers emit `RegionResult` structs as soon as each row finishes. Results flow through an `mpsc` channel to a dedicated “matrix-streamer” thread. A `BTreeMap` reorders out-of-order results so only contiguous indices are written, guaranteeing deterministic row order even with parallel execution.
4. As rows arrive, the streamer writes them directly into the data member, updates per-group counts, and discards the matrix payload—no temporary spooling. Once the channel closes, remaining buffered rows flush, and a final `MatrixHeader` is constructed from the accumulated counts.
5. The streamer rewrites the stored header member in place using the padded payload produced by `build_padded_header_payload`, finalises the gzip encoder, and surfaces any I/O errors back to the caller. In-memory mode still drives `write_outputs`, preserving the auxiliary artifacts and sorting capabilities that depend on a full `MatrixData`.

## Header Guardrails
- The reserved payload (8 192 − 23 bytes) matches the earlier design. If a tentative header exceeds the budget, streaming is bypassed so the legacy writer can handle the oversized JSON without risking corruption.
- The final rewrite uses the same deterministic gzip layout (`mtime(0)`, no filename/comment), keeping byte-for-byte parity with the previous implementation.

## Compatibility Notes
- Sorting modes other than `keep`, or requests for additional sinks, transparently fall back to the existing aggregate-and-write path.
- Skip-zero behaviour, threshold filtering, and negative-strand reversal all occur before a row is sent to the streamer, so downstream consumers see identical data regardless of the code path.
- Existing tests and regression harnesses continue to operate on the in-memory path by default; streaming is currently targeted at large unsorted jobs to reduce peak latency and memory pressure.
