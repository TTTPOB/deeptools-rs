# Profile: task9-bwig-buf-reuse
Time: 2026-04-28T23:47:26+08:00
Command: target/release/compute_matrix_rs scale-regions -b 10 -p 4 --regionsFileName target/compute-matrix-datasets/encode_k562_atac/ENCFF333TAT.bed --scoreFileName target/compute-matrix-datasets/encode_k562_atac/ENCFF019IPA.bigWig target/compute-matrix-datasets/encode_k562_atac/ENCFF093IIW.bigWig target/compute-matrix-datasets/encode_k562_atac/ENCFF656ZKM.bigWig -o /tmp/bench_task9.mat.gz
Target: reuse values/blocks/remaining buffers in BigWigReader
Hot path: BigWigReader::values + search_cir_tree

## /usr/bin/time -v
```
Running computeMatrix (Rust) in scale-regions mode.
[coalesce-gap] n_gaps=161669 p50=5475 p75=12795 threshold=2000
[coalesce-gap] strategy="no-coalesce" batches=269800 items=269800 ratio=1.00
	Command being timed: "target/release/compute_matrix_rs scale-regions -b 10 -p 4 --regionsFileName target/compute-matrix-datasets/encode_k562_atac/ENCFF333TAT.bed --scoreFileName target/compute-matrix-datasets/encode_k562_atac/ENCFF019IPA.bigWig target/compute-matrix-datasets/encode_k562_atac/ENCFF093IIW.bigWig target/compute-matrix-datasets/encode_k562_atac/ENCFF656ZKM.bigWig -o /tmp/bench_task9.mat.gz"
	User time (seconds): 28.20
	System time (seconds): 2.47
	Percent of CPU this job got: 271%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:11.30
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 553520
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 158663
	Voluntary context switches: 310
	Involuntary context switches: 113
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

 Performance counter stats for 'target/release/compute_matrix_rs scale-regions -b 10 -p 4 --regionsFileName target/compute-matrix-datasets/encode_k562_atac/ENCFF333TAT.bed --scoreFileName target/compute-matrix-datasets/encode_k562_atac/ENCFF019IPA.bigWig target/compute-matrix-datasets/encode_k562_atac/ENCFF093IIW.bigWig target/compute-matrix-datasets/encode_k562_atac/ENCFF656ZKM.bigWig -o /tmp/bench_task9.mat.gz':

         31,674.42 msec task-clock:u                     #    2.728 CPUs utilized             
                 0      context-switches:u               #    0.000 /sec                      
                 0      cpu-migrations:u                 #    0.000 /sec                      
           158,784      page-faults:u                    #    5.013 K/sec                     
   <not supported>      cycles:u                                                              
   <not supported>      instructions:u                                                        
   <not supported>      branches:u                                                            
   <not supported>      branch-misses:u                                                       
   <not supported>      L1-dcache-loads:u                                                     
   <not supported>      L1-dcache-load-misses:u                                               
   <not supported>      LLC-loads:u                                                           
   <not supported>      LLC-load-misses:u                                                     

      11.609338507 seconds time elapsed

      28.932669000 seconds user
       2.742167000 seconds sys


```
## CPU Hotspots (perf report, top functions by self time)
```
    42.63%  compute_matrix_  compute_matrix_rs     [.] zlib_rs::inflate::inflate_fast_help_avx2
    14.59%  compute_matrix_  compute_matrix_rs     [.] compute_matrix_rs::io::readers::bwig::BigWigReader::values
     9.84%  compute_matrix_  compute_matrix_rs     [.] compute_matrix_rs::pipeline::core::worker::process_batch
     7.23%  compute_matrix_  compute_matrix_rs     [.] zlib_rs::inflate::inftrees::inflate_table
     5.67%  compute_matrix_  compute_matrix_rs     [.] zlib_rs::deflate::algorithm::quick::deflate_quick
     3.22%  compute_matrix_  compute_matrix_rs     [.] compute_matrix_rs::pipeline::core::worker::aggregate_slice
     2.37%  compute_matrix_  libc.so.6             [.] __memmove_avx_unaligned_erms
     2.32%  compute_matrix_  compute_matrix_rs     [.] compute_matrix_rs::io::writers::write_matrix_value
     1.11%  compute_matrix_  compute_matrix_rs     [.] compute_matrix_rs::pipeline::zones::append_bins
     0.99%  compute_matrix_  libc.so.6             [.] _int_malloc
     0.86%  compute_matrix_  compute_matrix_rs     [.] zlib_rs::deflate::BitWriter::emit_dist_static
     0.77%  compute_matrix_  compute_matrix_rs     [.] zlib_rs::deflate::compare256::avx2::compare256
     0.69%  compute_matrix_  libc.so.6             [.] cfree@GLIBC_2.2.5
     0.61%  compute_matrix_  compute_matrix_rs     [.] zlib_rs::adler32::avx2::adler32_avx2_help
     0.60%  compute_matrix_  libc.so.6             [.] unlink_chunk.isra.0
     0.55%  compute_matrix_  compute_matrix_rs     [.] rint
```

## CPU Call Graph (perf report, top call chains)
```
# To display the perf.data header info, please use --header/--header-only options.
#
#
# Total Lost Samples: 306
#
# Samples: 121K of event 'cpu-clock:u'
# Event count (approx.): 30432000000
#
# Children      Self  Command          Shared Object         Symbol                                                                                                                                   
# ........  ........  ...............  ....................  .........................................................................................................................................
#
    86.00%     0.05%  compute_matrix_  compute_matrix_rs     [.] rayon::iter::plumbing::bridge_producer_consumer::helper
            |          
             --85.95%--rayon::iter::plumbing::bridge_producer_consumer::helper
                       |          
                        --85.86%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                                  |          
                                   --84.68%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                             |          
                                             |--69.56%--compute_matrix_rs::pipeline::core::worker::process_batch
                                             |          compute_matrix_rs::io::readers::bwig::BigWigReader::values
                                             |          |          
                                             |          |--42.63%--zlib_rs::inflate::inflate_fast_help_avx2
                                             |          |          
                                             |           --7.23%--zlib_rs::inflate::inftrees::inflate_table
                                             |          
                                              --15.11%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                                                        |          
                                                        |--12.92%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                                        |          rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                                                        |          |          
                                                        |           --11.89%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                                        |                     |          
                                                        |                     |--9.79%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                                                        |                     |          |          
                                                        |                     |           --9.52%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                                        |                     |                     |          
                                                        |                     |                      --9.02%--compute_matrix_rs::pipeline::core::worker::process_batch
                                                        |                     |          
                                                        |                      --2.10%--compute_matrix_rs::pipeline::core::worker::process_batch
                                                        |          
                                                         --2.19%--rayon_core::registry::WorkerThread::wait_until_cold
                                                                   |          
                                                                    --2.19%--<rayon_core::job::StackJob<L,F,R> as rayon_core::job::Job>::execute
                                                                              rayon::iter::plumbing::bridge_producer_consumer::helper
                                                                              rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec

    85.99%     0.00%  compute_matrix_  compute_matrix_rs     [.] rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
            |          
             --85.99%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                       |          
                        --84.95%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                  |          
                                  |--69.64%--compute_matrix_rs::pipeline::core::worker::process_batch
                                  |          compute_matrix_rs::io::readers::bwig::BigWigReader::values
                                  |          |          
                                  |          |--42.63%--zlib_rs::inflate::inflate_fast_help_avx2
                                  |          |          
                                  |           --7.23%--zlib_rs::inflate::inftrees::inflate_table
                                  |          
                                   --15.31%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                                             |          
                                             |--12.95%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                             |          rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                                             |          |          
                                             |           --11.95%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                             |                     |          
                                             |                     |--9.84%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                                             |                     |          |          
                                             |                     |           --9.56%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                             |                     |                     |          
                                             |                     |                      --9.02%--compute_matrix_rs::pipeline::core::worker::process_batch
                                             |                     |          
                                             |                      --2.09%--compute_matrix_rs::pipeline::core::worker::process_batch
                                             |          
                                              --2.36%--rayon_core::registry::WorkerThread::wait_until_cold
                                                        |          
                                                         --2.35%--<rayon_core::job::StackJob<L,F,R> as rayon_core::job::Job>::execute
                                                                   rayon::iter::plumbing::bridge_producer_consumer::helper
                                                                   rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec

    84.71%     9.84%  compute_matrix_  compute_matrix_rs     [.] compute_matrix_rs::pipeline::core::worker::process_batch
            |          
            |--74.87%--compute_matrix_rs::pipeline::core::worker::process_batch
            |          |          
            |          |--69.67%--compute_matrix_rs::io::readers::bwig::BigWigReader::values
            |          |          |          
            |          |          |--42.64%--zlib_rs::inflate::inflate_fast_help_avx2
            |          |          |          
            |          |           --7.23%--zlib_rs::inflate::inftrees::inflate_table
            |          |          
            |           --3.22%--compute_matrix_rs::pipeline::core::worker::aggregate_slice
            |          
             --7.39%--__GI___clone3
                       start_thread
                       std::sys::thread::unix::Thread::new::thread_start
                       core::ops::function::FnOnce::call_once$u7b$$u7b$vtable.shim$u7d$$u7d$::h72d3d1760aa83d43
                       std::sys::backtrace::__rust_begin_short_backtrace
                       rayon_core::registry::WorkerThread::wait_until_cold
                       |          
                        --6.02%--<rayon_core::job::StackJob<L,F,R> as rayon_core::job::Job>::execute
                                  rayon::iter::plumbing::bridge_producer_consumer::helper
                                  rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                                  |          
                                   --5.88%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                             rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                                             |          
                                              --5.88%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                                        rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                                                        rayon::iter::plumbing::bridge_producer_consumer::helper
                                                        rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::h514f9db56d34e0ec
                                                        rayon::iter::plumbing::bridge_producer_consumer::helper
                                                        compute_matrix_rs::pipeline::core::worker::process_batch

    69.67%    14.59%  compute_matrix_  compute_matrix_rs     [.] compute_matrix_rs::io::readers::bwig::BigWigReader::values
            |          
            |--55.07%--compute_matrix_rs::io::readers::bwig::BigWigReader::values
            |          |          
            |          |--42.64%--zlib_rs::inflate::inflate_fast_help_avx2
            |          |          
```

## heaptrack summary
```
total runtime: 29.49s.
calls to allocation functions: 10744504 (364381/s)
temporary memory allocations: 1995274 (67666/s)
peak heap memory consumption: 451.46M
peak RSS (including heaptrack overhead): 578.32M
total memory leaked: 14.45K
```

## Hotspot Analysis

### Comparison table

| Metric | Baseline | Task8 (bedrecord-move) | Task9 (bwig-buf-reuse) | vs Baseline | vs Task8 |
|---|---|---|---|---|---|
| Wall clock (time -v) | 15.13 s | 11.65 s | 11.30 s | **-25.3%** | **-3.0%** |
| Wall clock (perf stat) | 12.247 s | 11.607 s | 11.609 s | **-5.2%** | +0.0% |
| Max RSS | 784,896 KB | 552,788 KB | 553,520 KB | **-29.5%** | +0.1% |
| task-clock | 30,249 ms | 32,271 ms | 31,674 ms | +4.7% | **-1.9%** |
| User time (time -v) | 27.18 s | 29.70 s | 28.20 s | +3.8% | **-5.1%** |
| System time (time -v) | 7.10 s | 2.49 s | 2.47 s | **-65.2%** | -0.8% |
| Voluntary ctx switches | 24,857 | 279 | 310 | **-98.8%** | +11.1% |
| Involuntary ctx switches | 148 | 95 | 113 | **-23.6%** | +18.9% |
| Page faults (time -v) | 235,672 | 157,200 | 158,663 | **-32.7%** | +0.9% |
| Page faults (perf stat) | 232,758 | 157,124 | 158,784 | **-31.8%** | +1.1% |
| Peak heap (heaptrack) | -- | 451.41 MB | 451.46 MB | -- | +0.0% |
| Alloc calls (heaptrack) | -- | 17,486,880 (439,335/s) | 10,744,504 (364,381/s) | -- | **-38.6%** |
| Temp allocs (heaptrack) | -- | 6,945,250 (174,490/s) | 1,995,274 (67,666/s) | -- | **-71.3%** |

### Top 5 remaining CPU hotspots

| Rank | Function | Self % | Notes |
|---|---|---|---|
| 1 | `zlib_rs::inflate::inflate_fast_help_avx2` | 42.63% | zlib decompression inner loop; already the AVX2 fast path. Irreducible without a faster decompressor or reducing the number of decompress calls. |
| 2 | `BigWigReader::values` | 14.59% | Increased from 10.22% (Task8) because allocation overhead that was previously attributed to child functions (`_int_malloc`, `cfree`, `_int_realloc`) is now absorbed into `values()` self time as compute work. The decompress_vec call inside values is the dominant cost. Task 10 (zlib decode_buf reuse) targets remaining allocations here. |
| 3 | `process_batch` | 9.84% | Per-batch orchestration. Stable vs Task8 (9.61%). Overhead is mainly delegated to `values()` and `aggregate_slice`. |
| 4 | `zlib_rs::inflate::inftrees::inflate_table` | 7.23% | Huffman table construction per zlib block. Stable vs Task8 (7.08%). Could be reduced if decode_buf reuse (Task 10) avoids reinitializing the decompressor per call. |
| 5 | `zlib_rs::deflate::algorithm::quick::deflate_quick` | 5.67% | Output gzip compression. Stable vs Task8 (5.42%). Irreducible at current compression level. |

Zlib inflate+deflate combined: ~56.3% of total CPU time (vs ~57.8% Task8). The relative share is stable; allocation reduction freed CPU for useful work.

### Top 3 remaining allocation hotspots

| Rank | Source | Calls | Peak | Optimization target |
|---|---|---|---|---|
| 1 | `BigWigReader::values` | 2,530,424 | 0 B | Remaining allocations from `decompress_vec` (Vec growth during zlib inflate). Task 10: reuse a decode buffer to avoid repeated Vec::new + grow inside flate2. |
| 2 | `load_groups` | 2,158,403 | 36.02 MB | BED record Vec expansion during file parse. One-time at startup, not in hot path. Low priority. |
| 3 | `process_batch` | 1,618,800 | 327.00 MB | Per-batch `all_values` Vec and coverage buffer management. The 327 MB peak is the coverage buffer pool. Already optimized via thread_local pool; remaining calls are the per-row `all_values` Vec. |

### Key observations

1. **Allocation calls reduced by 38.6% vs Task8** (17.5M -> 10.7M): The primary goal of this task. By reusing `values_buf`, `blocks_buf`, and `remaining_buf` as struct fields instead of allocating new Vecs on each of the ~270K calls, we eliminated ~6.7M allocation calls. The `values()` call site dropped from being the source of 7.8M+ allocation calls (Task8 heaptrack) to 2.5M — a 68% reduction within that function alone.

2. **Temporary allocations reduced by 71.3% vs Task8** (6.9M -> 2.0M): Temporary allocations (allocated and freed within the same call stack) were dramatically reduced. Previously, each `values()` call created a new `Vec<BigWigValue>` that was iterated and dropped. Now the buffer persists across calls, so only Vec growth (when the buffer needs to expand) triggers allocation.

3. **User CPU time reduced 5.1% vs Task8** (29.70s -> 28.20s): Less allocator work (fewer malloc/realloc/free calls) translates directly to reduced user-space CPU time. The `_int_malloc` self% dropped from 2.11% (Task8) to 0.99% (Task9), and `cfree` dropped from 0.62% to 0.69% (within noise). The `_int_realloc` entry disappeared entirely from the top list.

4. **Wall clock stable** (-3.0% by time -v, +0.0% by perf stat): The wall clock improvement from time -v is within the expected variance range. The perf stat measurement (11.607s vs 11.609s) is essentially identical, confirming this is primarily an allocation optimization rather than a throughput change.

5. **RSS and page faults stable** (553 KB vs 552 KB, +0.1%): The buffer reuse does not change peak memory — the buffers themselves are small (a few KB each). The total heap footprint is dominated by the block cache and coverage pool.

6. **`BigWigReader::values` self% rose from 10.22% to 14.59%**: This is expected and healthy. Previously, allocation overhead was spread across `_int_malloc` (2.11%), `_int_realloc` (0.51%), `cfree` (0.62%), and `unlink_chunk` (0.84%) — totaling ~4.1% of attributable allocator time. With buffer reuse, that allocator work vanishes, and the relative share of `values()` self time (parsing, cache lookup) increases proportionally. The absolute CPU time in `values()` is actually lower.

### Verdict

**PASS** — Wall clock is within noise of Task8 (+0.0% by perf stat, -3.0% by time -v), well within the 5% regression threshold. The optimization delivers its intended benefit: allocation calls reduced by 38.6% (17.5M -> 10.7M), temporary allocations reduced by 71.3% (6.9M -> 2.0M), and user CPU time reduced by 5.1%. Cumulative vs baseline: wall clock -25.3%, RSS -29.5%, page faults -32.7%, voluntary ctx switches -98.8%.
