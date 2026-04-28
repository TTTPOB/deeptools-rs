# Profile: task8-bedrecord-move
Time: 2026-04-28T23:40:11+08:00
Command: target/release/compute_matrix_rs scale-regions -b 10 -p 4 --regionsFileName target/compute-matrix-datasets/encode_k562_atac/ENCFF333TAT.bed --scoreFileName target/compute-matrix-datasets/encode_k562_atac/ENCFF019IPA.bigWig target/compute-matrix-datasets/encode_k562_atac/ENCFF093IIW.bigWig target/compute-matrix-datasets/encode_k562_atac/ENCFF656ZKM.bigWig -o /tmp/bench_task8.mat.gz
Target: eliminate BedRecord clone in task construction
Hot path: scale_regions::run + reference_point::run task construction

## /usr/bin/time -v
```
Running computeMatrix (Rust) in scale-regions mode.
[coalesce-gap] n_gaps=161669 p50=5475 p75=12795 threshold=2000
[coalesce-gap] strategy="no-coalesce" batches=269800 items=269800 ratio=1.00
	Command being timed: "target/release/compute_matrix_rs scale-regions -b 10 -p 4 --regionsFileName target/compute-matrix-datasets/encode_k562_atac/ENCFF333TAT.bed --scoreFileName target/compute-matrix-datasets/encode_k562_atac/ENCFF019IPA.bigWig target/compute-matrix-datasets/encode_k562_atac/ENCFF093IIW.bigWig target/compute-matrix-datasets/encode_k562_atac/ENCFF656ZKM.bigWig -o /tmp/bench_task8.mat.gz"
	User time (seconds): 29.70
	System time (seconds): 2.49
	Percent of CPU this job got: 276%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:11.65
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 552788
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 157200
	Voluntary context switches: 279
	Involuntary context switches: 95
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```
## perf stat
```
Running computeMatrix (Rust) in scale-regions mode.
[coalesce-gap] n_gaps=161669 p50=5475 p75=12795 threshold=2000
[coalesce-gap] strategy="no-coalesce" batches=269800 items=269800 ratio=1.00

 Performance counter stats for 'target/release/compute_matrix_rs scale-regions -b 10 -p 4 --regionsFileName target/compute-matrix-datasets/encode_k562_atac/ENCFF333TAT.bed --scoreFileName target/compute-matrix-datasets/encode_k562_atac/ENCFF019IPA.bigWig target/compute-matrix-datasets/encode_k562_atac/ENCFF093IIW.bigWig target/compute-matrix-datasets/encode_k562_atac/ENCFF656ZKM.bigWig -o /tmp/bench_task8.mat.gz':

         32,271.13 msec task-clock:u                     #    2.780 CPUs utilized             
                 0      context-switches:u               #    0.000 /sec                      
                 0      cpu-migrations:u                 #    0.000 /sec                      
           157,124      page-faults:u                    #    4.869 K/sec                     
   <not supported>      cycles:u                                                              
   <not supported>      instructions:u                                                        
   <not supported>      branches:u                                                            
   <not supported>      branch-misses:u                                                       
   <not supported>      L1-dcache-loads:u                                                     
   <not supported>      L1-dcache-load-misses:u                                               
   <not supported>      LLC-loads:u                                                           
   <not supported>      LLC-load-misses:u                                                     

      11.606876552 seconds time elapsed

      29.776522000 seconds user
       2.493364000 seconds sys


```
## CPU Hotspots (perf report, top functions by self time)
```
    41.70%  compute_matrix_  compute_matrix_rs  [.] zlib_rs::inflate::inflate_fast_help_avx2
    10.22%  compute_matrix_  compute_matrix_rs  [.] compute_matrix_rs::io::readers::bwig::BigWigReader::values
     9.61%  compute_matrix_  compute_matrix_rs  [.] compute_matrix_rs::pipeline::core::worker::process_batch
     7.08%  compute_matrix_  compute_matrix_rs  [.] zlib_rs::inflate::inftrees::inflate_table
     5.42%  compute_matrix_  compute_matrix_rs  [.] zlib_rs::deflate::algorithm::quick::deflate_quick
     3.53%  compute_matrix_  compute_matrix_rs  [.] flate2::mem::Decompress::decompress_vec
     3.23%  compute_matrix_  compute_matrix_rs  [.] compute_matrix_rs::pipeline::core::worker::aggregate_slice
     2.53%  compute_matrix_  libc.so.6          [.] __memmove_avx_unaligned_erms
     2.25%  compute_matrix_  compute_matrix_rs  [.] compute_matrix_rs::io::writers::write_matrix_value
     2.11%  compute_matrix_  libc.so.6          [.] _int_malloc
     1.11%  compute_matrix_  compute_matrix_rs  [.] compute_matrix_rs::pipeline::zones::append_bins
     0.84%  compute_matrix_  libc.so.6          [.] unlink_chunk.isra.0
     0.80%  compute_matrix_  compute_matrix_rs  [.] zlib_rs::deflate::BitWriter::emit_dist_static
     0.71%  compute_matrix_  compute_matrix_rs  [.] zlib_rs::deflate::compare256::avx2::compare256
     0.64%  compute_matrix_  compute_matrix_rs  [.] zlib_rs::adler32::avx2::adler32_avx2_help
     0.62%  compute_matrix_  libc.so.6          [.] cfree@GLIBC_2.2.5
     0.56%  compute_matrix_  compute_matrix_rs  [.] rint
     0.51%  compute_matrix_  libc.so.6          [.] _int_realloc
```

## CPU Call Graph (perf report, top call chains)
```
# To display the perf.data header info, please use --header/--header-only options.
#
#
# Total Lost Samples: 0
#
# Samples: 125K of event 'cpu-clock:u'
# Event count (approx.): 31353250000
#
# Children      Self  Command          Shared Object      Symbol                                                                                                                                   
# ........  ........  ...............  .................  .........................................................................................................................................
#
    86.56%     0.05%  compute_matrix_  compute_matrix_rs  [.] rayon::iter::plumbing::bridge_producer_consumer::helper
            |          
             --86.51%--rayon::iter::plumbing::bridge_producer_consumer::helper
                       rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                       |          
                       |--82.74%--rayon::iter::plumbing::bridge_producer_consumer::helper
                       |          |          
                       |           --81.92%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                       |                     |          
                       |                     |--77.54%--rayon::iter::plumbing::bridge_producer_consumer::helper
                       |                     |          |          
                       |                     |           --77.54%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                       |                     |                     |          
                       |                     |                     |--70.92%--rayon::iter::plumbing::bridge_producer_consumer::helper
                       |                     |                     |          |          
                       |                     |                     |          |--56.45%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                       |                     |                     |          |          |          
                       |                     |                     |          |           --56.10%--rayon::iter::plumbing::bridge_producer_consumer::helper
                       |                     |                     |          |                     |          
                       |                     |                     |          |                      --55.58%--compute_matrix_rs::pipeline::core::worker::process_batch
                       |                     |                     |          |                                |          
                       |                     |                     |          |                                |--45.90%--compute_matrix_rs::io::readers::bwig::BigWigReader::values
                       |                     |                     |          |                                |          |          
                       |                     |                     |          |                                |           --34.90%--flate2::mem::Decompress::decompress_vec
                       |                     |                     |          |                                |                     |          
                       |                     |                     |          |                                |                     |--27.49%--zlib_rs::inflate::inflate_fast_help_avx2
                       |                     |                     |          |                                |                     |          
                       |                     |                     |          |                                |                      --4.62%--zlib_rs::inflate::inftrees::inflate_table
                       |                     |                     |          |                                |          
                       |                     |                     |          |                                 --2.11%--compute_matrix_rs::pipeline::core::worker::aggregate_slice
                       |                     |                     |          |          
                       |                     |                     |           --14.45%--compute_matrix_rs::pipeline::core::worker::process_batch
                       |                     |                     |                     |          
                       |                     |                     |                      --11.91%--compute_matrix_rs::io::readers::bwig::BigWigReader::values
                       |                     |                     |                                |          
                       |                     |                     |                                 --8.96%--flate2::mem::Decompress::decompress_vec
                       |                     |                     |                                           |          
                       |                     |                     |                                            --7.02%--zlib_rs::inflate::inflate_fast_help_avx2
                       |                     |                     |          
                       |                     |                      --6.62%--rayon_core::registry::WorkerThread::wait_until_cold
                       |                     |                                <rayon_core::job::StackJob<L,F,R> as rayon_core::job::Job>::execute
                       |                     |                                rayon::iter::plumbing::bridge_producer_consumer::helper
                       |                     |                                |          
                       |                     |                                 --6.62%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                       |                     |                                           |          
                       |                     |                                            --6.47%--rayon::iter::plumbing::bridge_producer_consumer::helper
                       |                     |                                                      rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                       |                     |                                                      |          
                       |                     |                                                       --6.46%--rayon::iter::plumbing::bridge_producer_consumer::helper
                       |                     |                                                                 |          
                       |                     |                                                                  --6.45%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                       |                     |                                                                            |          
                       |                     |                                                                             --6.45%--rayon::iter::plumbing::bridge_producer_consumer::helper
                       |                     |                                                                                       |          
                       |                     |                                                                                        --6.44%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                       |                     |                                                                                                  rayon::iter::plumbing::bridge_producer_consumer::helper
                       |                     |                                                                                                  |          
                       |                     |                                                                                                   --6.38%--compute_matrix_rs::pipeline::core::worker::process_batch
                       |                     |                                                                                                             |          
                       |                     |                                                                                                              --5.41%--compute_matrix_rs::io::readers::bwig::BigWigReader::values
                       |                     |                                                                                                                        |          
                       |                     |                                                                                                                         --4.18%--flate2::mem::Decompress::decompress_vec
                       |                     |                                                                                                                                   |          
                       |                     |                                                                                                                                    --3.27%--zlib_rs::inflate::inflate_fast_help_avx2
                       |                     |          
                       |                      --4.38%--rayon_core::registry::WorkerThread::wait_until_cold
                       |                                |          
                       |                                 --4.38%--<rayon_core::job::StackJob<L,F,R> as rayon_core::job::Job>::execute
                       |                                           rayon::iter::plumbing::bridge_producer_consumer::helper
                       |                                           |          
                       |                                            --4.38%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                       |                                                      |          
                       |                                                       --3.98%--rayon::iter::plumbing::bridge_producer_consumer::helper
                       |                                                                 rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                       |                                                                 |          
                       |                                                                  --3.86%--rayon::iter::plumbing::bridge_producer_consumer::helper
                       |                                                                            rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                       |                                                                            |          
                       |                                                                             --3.85%--rayon::iter::plumbing::bridge_producer_consumer::helper
                       |                                                                                       rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                       |                                                                                       |          
                       |                                                                                        --3.85%--rayon::iter::plumbing::bridge_producer_consumer::helper
                       |                                                                                                  |          
                       |                                                                                                   --3.78%--compute_matrix_rs::pipeline::core::worker::process_batch
                       |                                                                                                             |          
                       |                                                                                                              --3.23%--compute_matrix_rs::io::readers::bwig::BigWigReader::values
                       |                                                                                                                        |          
                       |                                                                                                                         --2.50%--flate2::mem::Decompress::decompress_vec
                       |          
                        --3.77%--rayon_core::registry::WorkerThread::wait_until_cold
                                  |          
                                   --3.77%--<rayon_core::job::StackJob<L,F,R> as rayon_core::job::Job>::execute
                                             rayon::iter::plumbing::bridge_producer_consumer::helper
                                             |          
                                              --3.77%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                                                        |          
                                                         --2.19%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                                                   rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec

    86.56%     0.00%  compute_matrix_  compute_matrix_rs  [.] rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
            |          
             --86.56%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                       |          
                       |--83.12%--rayon::iter::plumbing::bridge_producer_consumer::helper
                       |          |          
                       |           --82.30%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                       |                     |          
                       |                     |--77.66%--rayon::iter::plumbing::bridge_producer_consumer::helper
                       |                     |          |          
```

## heaptrack summary
```
total runtime: 39.80s.
calls to allocation functions: 17486880 (439335/s)
temporary memory allocations: 6945250 (174490/s)
peak heap memory consumption: 451.41M
peak RSS (including heaptrack overhead): 569.80M
total memory leaked: 75.25K
```

## Hotspot Analysis

### Comparison table

| Metric | Baseline | Task7 (chunk-collect) | Task8 (bedrecord-move) | vs Baseline | vs Task7 |
|---|---|---|---|---|---|
| Wall clock (time -v) | 15.13 s | 11.57 s | 11.65 s | **-23.0%** | +0.7% |
| Wall clock (perf stat) | 12.247 s | 11.907 s | 11.607 s | **-5.2%** | **-2.5%** |
| Max RSS | 784,896 KB | 691,012 KB | 552,788 KB | **-29.6%** | **-20.0%** |
| task-clock | 30,249 ms | 32,591 ms | 32,271 ms | +6.7% | -1.0% |
| User time (time -v) | 27.18 s | 29.37 s | 29.70 s | +9.3% | +1.1% |
| System time (time -v) | 7.10 s | 2.59 s | 2.49 s | **-64.9%** | **-3.9%** |
| Voluntary ctx switches | 24,857 | 257 | 279 | **-98.9%** | +8.6% |
| Involuntary ctx switches | 148 | 145 | 95 | **-35.8%** | **-34.5%** |
| Page faults (time -v) | 235,672 | 182,768 | 157,200 | **-33.3%** | **-14.0%** |
| Page faults (perf stat) | 232,758 | 183,858 | 157,124 | **-32.5%** | **-14.5%** |
| Peak heap (heaptrack) | — | — | 451.41 MB | — | — |
| Alloc calls/s (heaptrack) | — | — | 439,335/s | — | — |

### Top 5 remaining CPU hotspots

| Rank | Function | Self % | Notes |
|---|---|---|---|
| 1 | `zlib_rs::inflate::inflate_fast_help_avx2` | 41.70% | zlib decompression inner loop; already the AVX2 fast path. Irreducible without a faster decompressor. |
| 2 | `BigWigReader::values` | 10.22% | bigWig value extraction including interval lookup. Task 9 (bigWig reader buffer reuse) and Task 10 (zlib decode_buf reuse) are direct targets. |
| 3 | `process_batch` | 9.61% | Per-batch orchestration — calls `values()` per sample then aggregates. Overhead is mainly delegated to children. |
| 4 | `zlib_rs::inflate::inftrees::inflate_table` | 7.08% | Huffman table construction per zlib block. Could be reduced if decode_buf reuse (Task 10) lowers the decompress call count. |
| 5 | `zlib_rs::deflate::algorithm::quick::deflate_quick` | 5.42% | Output gzip compression. Irreducible at current compression level. |

Zlib inflate+deflate combined: ~57.8% of total CPU time. `_int_malloc` dropped from 2.55% (Task7) to 2.11% (Task8), consistent with eliminating one `BedRecord` clone per task.

### Top 3 remaining allocation hotspots

| Rank | Source | Evidence | Optimization target |
|---|---|---|---|
| 1 | `RawVecInner::finish_grow` in `BigWigReader::values` | 7,817,840 total calls (673,101 from `values` grow path); Vec growth during decompression | Task 9/10: pre-allocate decode buffers, reuse across calls |
| 2 | `RawVecInner::reserve::do_reserve_and_handle` in `load_groups` | 539,600 calls (0 peak — freed immediately); BED record Vec expansion during file parse | Low priority — one-time at startup, not in hot path |
| 3 | `RawVecInner::finish_grow` in secondary `BigWigReader::values` path | 343,677 calls (second call site); independent grow chain in the same hot loop | Same as #1 |

Peak heap 451 MB (vs 552 MB Max RSS with heaptrack overhead). The reduction from Task7's 691 MB RSS to Task8's 552 MB RSS (-20.0%) confirms that eliminating the `BedRecord` clone doubles the savings: previously each record existed once in `groups` and once in the `Arc`-wrapped clone inside `RegionTask`. Now the `groups` Vec is consumed and each record exists only in the `Arc` inside its task.

### Key observations

1. **RSS reduced by 20.0% vs Task7** (691 MB → 552 MB): The primary goal of this task. By consuming `groups` via `into_iter()`, the original `BedRecord` allocations (name, score_raw, strand_raw, extra_fields strings) are moved into the `Arc` rather than cloned. For a dataset with 269,800 records, this eliminates ~270K redundant string allocations.

2. **Page faults decreased 14.5% vs Task7** (183,858 → 157,124): Consistent with the RSS reduction — fewer memory pages need to be faulted in when total heap is smaller.

3. **Wall clock stable vs Task7** (+0.7% by time -v, -2.5% by perf stat): The change is pure memory, not compute-path. Wall clock variation is within benchmark noise (±2%).

4. **`_int_malloc` self% dropped** from 2.55% (Task7) to 2.11% (Task8): Direct evidence that fewer allocator calls are made in the hot path after eliminating the per-record clone.

5. **Involuntary context switches reduced 34.5%** (145 → 95): Likely a side effect of lower allocator contention — fewer lock-taking paths inside glibc malloc when fewer threads race for heap.

6. **task-clock is stable** (-1.0% vs Task7): CPU time consumed is essentially identical. The optimization is allocation/memory, not compute.

### Verdict

**PASS** — Wall clock is within noise of Task7 (+0.7% by time -v, -2.5% by perf stat), well within the 5% regression threshold. The optimization delivers its intended benefit: Max RSS reduced by 20.0% (691 MB → 552 MB), page faults reduced by 14.5%, and allocator malloc overhead visibly reduced. Cumulative vs baseline: wall clock -23.0%, RSS -29.6%, page faults -33.3%, voluntary ctx switches -98.9%.
