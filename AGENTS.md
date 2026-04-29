## ref impl test plan
- use pixi to create python test environment
- include bioconda channel, deeptools 3.5.6 package
- find relavant test data, generate an reference output using `pixi run computeMatrix ...`
- at the beginning of each session, read `plans/overall_status.md` to understand current progress, after each code change, update `plans/overall_status.md` to reflect current status

## architecture (as of v0.3.0)
- **two execution paths**: StreamOrdered (sort=no / keep+already_sorted) and HybridBucket (keep+coalesced / ascend / descend)
- **all output through FileCollector**: main gzip + optional auxiliary writers (outFileNameMatrix, outFileSortedRegions)
- **file spilling**: HybridBucketCollector accumulates rows per group, spills to temp files when memory > 1GB (injectable threshold for testing), finalize via mmap readback
- **compare_matrix binary**: dev tool for comparing matrix files (plain/gzip/multi-member), subcommands: header, values, diff
- **Rust integration tests**: `tests/python_compatibility.rs` — 10 scenarios using `compare_matrix diff`
- **removed**: InMemoryCollector, GroupBucketCollector, MatrixData, write_outputs(), RunOutcome::Matrix, spawn_writer_thread

## workspace
for openai codex: you are configured to run code in sandbox so that network access may not be reliable. when you try to run `cargo update` or `cargo metadata`, you may see errors. for `cargo update`, you don't need to run it. `rust-analyzer` will handle by automatically lock. `cargo metadata`, you may want to use it to get doc or codes, try to use rg to find them under local directory, as `rust-analyzer` may have already indexed them. 