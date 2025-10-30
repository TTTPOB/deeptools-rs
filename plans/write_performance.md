# Write Path Performance – 2025-10-30

## Observations
- `write_matrix_row` in `src/io/writers/mod.rs:279` converts every float to a heap allocated `String` via `format!("{value:.6}")`, and the outer `write!` macro performs an additional formatting pass. The flame graph’s tall stack of `core::fmt::write`, `std::io::Write::write_all`, and `std::io::default_write_fmt` blocks confirms this path dominates CPU time when emitting rows.
- The region metadata fields (`name`, `score`) follow the same pattern (`format!("{score:.6}")`), so each row pays multiple small allocations before any bytes reach the encoder.
- Streaming output feeds a `flate2::write::GzEncoder<File>` constructed with `Compression::default()` (libz level 6). The profile shows a wide band under `flate2::deflate::compress_inner`, indicating the compressor competes with formatting work for total runtime.
- When the streaming fast-path is exercised, the gzip encoder is not wrapped in a buffered writer, so we rely solely on the flate2 internal buffer. Small formatted writes from the per-value `write!` calls can therefore churn the encoder’s deflate state and amplify compression cost.

## Progress
- 2025-10-30: `write_matrix_row` now emits matrix values via a stack-buffered formatter (manual fixed-point conversion + `itoa`) to remove the per-cell `String` allocation path.

## Priorities
1. **Eliminate per-value `String` allocations** by switching to stack-based formatters (e.g. `ryu::Buffer` for floats, `itoa` for integers) and writing the resulting byte slice directly. This tackles the `core::fmt` hotspot and should reduce allocator pressure significantly.
2. **Batch row serialization** so an entire line is assembled in a reusable `Vec<u8>` (or `SmallVec`) before a single `write_all`. This reduces `write!` overhead, improves cache locality, and feeds larger chunks to gzip.
3. **Expose a configurable gzip level** (with `Compression::new(1)` as the default for performance runs) or detect when output is not size-sensitive. Dropping from level 6 to 1 can roughly halve `compress_inner` time while still generating compressed output.
4. **Consider buffering the encoder** (`BufWriter<GzEncoder<File>>` or `GzEncoder<BufWriter<File>>`) once formatting costs drop, to cut down on syscalls and give deflate larger blocks.
5. **Profile the plain (uncompressed) writer path** after the formatting improvements to ensure no hidden hotspots remain before re-tuning compression.
