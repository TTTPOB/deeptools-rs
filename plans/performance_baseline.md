# Performance Baseline

## Environment
- CPU: 12 cores
- RAM: ~8 GB available
- OS: Debian Trixie (WSL2)
- Rust: release profile (opt-level=3)

## Parity Tests: 10/10 PASS
All 10 python-compatibility scenarios pass with tolerance ≤ 5e-06.

## ENCODE Benchmark (K562 ATAC-seq, 269,800 regions, 3 bigWig samples)

### Reference-point mode (center, +/-100bp, bin=10, 20 bins/sample)
- Wall time: 8.59s
- User time: 24.00s
- System time: 5.22s
- Max RSS: 425,784 KB (~416 MB)
- Speedup vs Python: 27.95x

### Scale-regions mode (body=200, +/-100bp, unscaled=50/50, bin=10, 50 bins/sample)
- Wall time: 13.70s
- User time: 39.99s
- System time: 5.18s
- Max RSS: 425,060 KB (~415 MB)
- Speedup vs Python: 36.74x
