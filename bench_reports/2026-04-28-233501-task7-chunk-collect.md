# Profile: task7-chunk-collect
Time: 2026-04-28T23:35:01+08:00
Command: target/release/compute_matrix_rs scale-regions -b 10 -p 4 --regionsFileName target/compute-matrix-datasets/encode_k562_atac/ENCFF333TAT.bed --scoreFileName target/compute-matrix-datasets/encode_k562_atac/ENCFF019IPA.bigWig target/compute-matrix-datasets/encode_k562_atac/ENCFF093IIW.bigWig target/compute-matrix-datasets/encode_k562_atac/ENCFF656ZKM.bigWig -o /tmp/bench_task7.mat.gz
Target: replace channel+BTreeMap with chunk-collect dispatch
Hot path: execute_mode chunk processing + output dispatch

## /usr/bin/time -v
```
Running computeMatrix (Rust) in scale-regions mode.
[coalesce-gap] n_gaps=161669 p50=5475 p75=12795 threshold=2000
[coalesce-gap] strategy="no-coalesce" batches=269800 items=269800 ratio=1.00
	Command being timed: "target/release/compute_matrix_rs scale-regions -b 10 -p 4 --regionsFileName target/compute-matrix-datasets/encode_k562_atac/ENCFF333TAT.bed --scoreFileName target/compute-matrix-datasets/encode_k562_atac/ENCFF019IPA.bigWig target/compute-matrix-datasets/encode_k562_atac/ENCFF093IIW.bigWig target/compute-matrix-datasets/encode_k562_atac/ENCFF656ZKM.bigWig -o /tmp/bench_task7.mat.gz"
	User time (seconds): 29.37
	System time (seconds): 2.59
	Percent of CPU this job got: 276%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:11.57
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 691012
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 4
	Minor (reclaiming a frame) page faults: 182768
	Voluntary context switches: 257
	Involuntary context switches: 145
	Swaps: 0
	File system inputs: 528
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

 Performance counter stats for 'target/release/compute_matrix_rs scale-regions -b 10 -p 4 --regionsFileName target/compute-matrix-datasets/encode_k562_atac/ENCFF333TAT.bed --scoreFileName target/compute-matrix-datasets/encode_k562_atac/ENCFF019IPA.bigWig target/compute-matrix-datasets/encode_k562_atac/ENCFF093IIW.bigWig target/compute-matrix-datasets/encode_k562_atac/ENCFF656ZKM.bigWig -o /tmp/bench_task7.mat.gz':

         32,591.44 msec task-clock:u                     #    2.737 CPUs utilized             
                 0      context-switches:u               #    0.000 /sec                      
                 0      cpu-migrations:u                 #    0.000 /sec                      
           183,858      page-faults:u                    #    5.641 K/sec                     
   <not supported>      cycles:u                                                              
   <not supported>      instructions:u                                                        
   <not supported>      branches:u                                                            
   <not supported>      branch-misses:u                                                       
   <not supported>      L1-dcache-loads:u                                                     
   <not supported>      L1-dcache-load-misses:u                                               
   <not supported>      LLC-loads:u                                                           
   <not supported>      LLC-load-misses:u                                                     

      11.906601034 seconds time elapsed

      29.776949000 seconds user
       2.808466000 seconds sys


```
## CPU Hotspots (perf report, top functions by self time)
```
    40.74%  compute_matrix_  compute_matrix_rs     [.] zlib_rs::inflate::inflate_fast_help_avx2
     9.83%  compute_matrix_  compute_matrix_rs     [.] compute_matrix_rs::io::readers::bwig::BigWigReader::values
     9.59%  compute_matrix_  compute_matrix_rs     [.] compute_matrix_rs::pipeline::core::worker::process_batch
     6.93%  compute_matrix_  compute_matrix_rs     [.] zlib_rs::inflate::inftrees::inflate_table
     5.14%  compute_matrix_  compute_matrix_rs     [.] zlib_rs::deflate::algorithm::quick::deflate_quick
     3.66%  compute_matrix_  compute_matrix_rs     [.] flate2::mem::Decompress::decompress_vec
     3.23%  compute_matrix_  compute_matrix_rs     [.] compute_matrix_rs::pipeline::core::worker::aggregate_slice
     2.58%  compute_matrix_  libc.so.6             [.] __memmove_avx_unaligned_erms
     2.55%  compute_matrix_  libc.so.6             [.] _int_malloc
     2.11%  compute_matrix_  compute_matrix_rs     [.] compute_matrix_rs::io::writers::write_matrix_value
     1.15%  compute_matrix_  compute_matrix_rs     [.] compute_matrix_rs::pipeline::zones::append_bins
     0.99%  compute_matrix_  libc.so.6             [.] unlink_chunk.isra.0
     0.76%  compute_matrix_  compute_matrix_rs     [.] zlib_rs::deflate::BitWriter::emit_dist_static
     0.74%  compute_matrix_  libc.so.6             [.] cfree@GLIBC_2.2.5
     0.73%  compute_matrix_  compute_matrix_rs     [.] zlib_rs::deflate::compare256::avx2::compare256
     0.65%  compute_matrix_  compute_matrix_rs     [.] zlib_rs::adler32::avx2::adler32_avx2_help
     0.62%  compute_matrix_  libc.so.6             [.] _int_realloc
     0.55%  compute_matrix_  compute_matrix_rs     [.] rint
     0.51%  compute_matrix_  libc.so.6             [.] _int_free_create_chunk
```

## CPU Call Graph (perf report, top call chains)
```
# To display the perf.data header info, please use --header/--header-only options.
#
#
# Total Lost Samples: 0
#
# Samples: 131K of event 'cpu-clock:u'
# Event count (approx.): 32771250000
#
# Children      Self  Command          Shared Object         Symbol                                                                                                                                   
# ........  ........  ...............  ....................  .........................................................................................................................................
#
    86.46%     0.07%  compute_matrix_  compute_matrix_rs     [.] rayon::iter::plumbing::bridge_producer_consumer::helper
            |          
             --86.39%--rayon::iter::plumbing::bridge_producer_consumer::helper
                       |          
                        --86.39%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::hc4312420b96fdb6d
                                  |          
                                  |--80.25%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                  |          |          
                                  |           --79.32%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::hc4312420b96fdb6d
                                  |                     |          
                                  |                     |--72.57%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                  |                     |          |          
                                  |                     |           --72.57%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::hc4312420b96fdb6d
                                  |                     |                     |          
                                  |                     |                     |--65.53%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                  |                     |                     |          |          
                                  |                     |                     |          |--53.73%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::hc4312420b96fdb6d
                                  |                     |                     |          |          |          
                                  |                     |                     |          |           --53.12%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                  |                     |                     |          |                     |          
                                  |                     |                     |          |                      --52.57%--compute_matrix_rs::pipeline::core::worker::process_batch
                                  |                     |                     |          |                                |          
                                  |                     |                     |          |                                |--43.37%--compute_matrix_rs::io::readers::bwig::BigWigReader::values
                                  |                     |                     |          |                                |          |          
                                  |                     |                     |          |                                |          |--32.66%--flate2::mem::Decompress::decompress_vec
                                  |                     |                     |          |                                |          |          |          
                                  |                     |                     |          |                                |          |          |--25.49%--zlib_rs::inflate::inflate_fast_help_avx2
                                  |                     |                     |          |                                |          |          |          
                                  |                     |                     |          |                                |          |           --4.38%--zlib_rs::inflate::inftrees::inflate_table
                                  |                     |                     |          |                                |          |          
                                  |                     |                     |          |                                |           --2.10%--alloc::raw_vec::RawVec<T,A>::grow_one
                                  |                     |                     |          |                                |                     |          
                                  |                     |                     |          |                                |                      --2.06%--alloc::raw_vec::RawVecInner<A>::finish_grow
                                  |                     |                     |          |                                |          
                                  |                     |                     |          |                                 --2.00%--compute_matrix_rs::pipeline::core::worker::aggregate_slice
                                  |                     |                     |          |          
                                  |                     |                     |           --11.79%--compute_matrix_rs::pipeline::core::worker::process_batch
                                  |                     |                     |                     |          
                                  |                     |                     |                      --9.72%--compute_matrix_rs::io::readers::bwig::BigWigReader::values
                                  |                     |                     |                                |          
                                  |                     |                     |                                 --7.24%--flate2::mem::Decompress::decompress_vec
                                  |                     |                     |                                           |          
                                  |                     |                     |                                            --5.66%--zlib_rs::inflate::inflate_fast_help_avx2
                                  |                     |                     |          
                                  |                     |                      --7.04%--rayon_core::registry::WorkerThread::wait_until_cold
                                  |                     |                                |          
                                  |                     |                                 --7.04%--<rayon_core::job::StackJob<L,F,R> as rayon_core::job::Job>::execute
                                  |                     |                                           rayon::iter::plumbing::bridge_producer_consumer::helper
                                  |                     |                                           rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::hc4312420b96fdb6d
                                  |                     |                                           |          
                                  |                     |                                            --6.92%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                  |                     |                                                      |          
                                  |                     |                                                       --6.92%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::hc4312420b96fdb6d
                                  |                     |                                                                 |          
                                  |                     |                                                                  --6.89%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                  |                     |                                                                            rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::hc4312420b96fdb6d
                                  |                     |                                                                            rayon::iter::plumbing::bridge_producer_consumer::helper
                                  |                     |                                                                            |          
                                  |                     |                                                                             --6.88%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::hc4312420b96fdb6d
                                  |                     |                                                                                       |          
                                  |                     |                                                                                        --6.88%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                  |                     |                                                                                                  |          
                                  |                     |                                                                                                   --6.80%--compute_matrix_rs::pipeline::core::worker::process_batch
                                  |                     |                                                                                                             |          
                                  |                     |                                                                                                              --5.89%--compute_matrix_rs::io::readers::bwig::BigWigReader::values
                                  |                     |                                                                                                                        |          
                                  |                     |                                                                                                                         --4.55%--flate2::mem::Decompress::decompress_vec
                                  |                     |                                                                                                                                   |          
                                  |                     |                                                                                                                                    --3.57%--zlib_rs::inflate::inflate_fast_help_avx2
                                  |                     |          
                                  |                      --6.75%--rayon_core::registry::WorkerThread::wait_until_cold
                                  |                                <rayon_core::job::StackJob<L,F,R> as rayon_core::job::Job>::execute
                                  |                                rayon::iter::plumbing::bridge_producer_consumer::helper
                                  |                                rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::hc4312420b96fdb6d
                                  |                                |          
                                  |                                 --6.03%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                  |                                           rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::hc4312420b96fdb6d
                                  |                                           |          
                                  |                                            --5.90%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                  |                                                      rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::hc4312420b96fdb6d
                                  |                                                      |          
                                  |                                                       --5.88%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                  |                                                                 rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::hc4312420b96fdb6d
                                  |                                                                 |          
                                  |                                                                  --5.87%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                  |                                                                            |          
                                  |                                                                             --5.76%--compute_matrix_rs::pipeline::core::worker::process_batch
                                  |                                                                                       |          
                                  |                                                                                        --4.84%--compute_matrix_rs::io::readers::bwig::BigWigReader::values
                                  |                                                                                                  |          
                                  |                                                                                                   --3.69%--flate2::mem::Decompress::decompress_vec
                                  |                                                                                                             |          
                                  |                                                                                                              --2.92%--zlib_rs::inflate::inflate_fast_help_avx2
                                  |          
                                   --6.14%--rayon_core::registry::WorkerThread::wait_until_cold
                                             |          
                                              --6.14%--<rayon_core::job::StackJob<L,F,R> as rayon_core::job::Job>::execute
                                                        |          
                                                         --6.14%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                                                   |          
                                                                    --6.14%--rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::hc4312420b96fdb6d
                                                                              |          
                                                                               --4.27%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                                                                         rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::hc4312420b96fdb6d
                                                                                         |          
                                                                                          --3.94%--rayon::iter::plumbing::bridge_producer_consumer::helper
                                                                                                    rayon_core::join::join_context::_$u7b$$u7b$closure$u7d$$u7d$::hc4312420b96fdb6d
                                                                                                    |          
                                                                                                     --3.79%--rayon::iter::plumbing::bridge_producer_consumer::helper
```

## heaptrack summary
```
(no heaptrack data captured for this run)
```

## Hotspot Analysis

### Comparison table

| Metric | Baseline | Previous (m4-pool) | Task7 (chunk-collect) | vs Baseline | vs Previous |
|---|---|---|---|---|---|
| Wall clock (time -v) | 15.13 s | 11.51 s | 11.57 s | **-23.5%** | +0.5% |
| Wall clock (perf stat) | 12.247 s | 12.059 s | 11.907 s | **-2.8%** | -1.3% |
| Max RSS | 784,896 KB | 515,576 KB | 691,012 KB | **-12.0%** | **+34.0%** |
| task-clock | 30,249 ms | 34,670 ms | 32,591 ms | +7.7% | **-6.0%** |
| User time (time -v) | 27.18 s | 29.05 s | 29.37 s | +8.1% | +1.1% |
| System time (time -v) | 7.10 s | 5.87 s | 2.59 s | **-63.5%** | **-55.9%** |
| Voluntary ctx switches | 24,857 | 189,810 | 257 | **-99.0%** | **-99.9%** |
| Involuntary ctx switches | 148 | 105 | 145 | -2.0% | +38.1% |
| Page faults (time -v) | 235,672 | 142,224 | 182,768 | **-22.5%** | **+28.5%** |
| Page faults (perf stat) | 232,758 | 166,714 | 183,858 | **-21.0%** | **+10.3%** |

### Top 5 remaining CPU hotspots

| Rank | Function | Self % | Notes |
|---|---|---|---|
| 1 | `zlib_rs::inflate::inflate_fast_help_avx2` | 40.74% | zlib decompression inner loop; irreducible — this is the AVX2-optimized fast path already. Not a target unless switching decompressors. |
| 2 | `BigWigReader::values` | 9.83% | bigWig value extraction including interval lookup and buffer management. Task 9 (bigWig reader buffer reuse) and Task 10 (zlib decode_buf reuse) are direct targets for reducing overhead here. |
| 3 | `process_batch` | 9.59% | Per-batch orchestration: calls values() per sample, then aggregates. The chunk-collect rewrite itself is this function. Overhead is dispatch + per-zone iteration — mostly delegated to children. |
| 4 | `zlib_rs::inflate::inftrees::inflate_table` | 6.93% | Huffman table construction per zlib block. Could be reduced if decode_buf reuse (Task 10) reduces the number of decompress calls by reusing inflated buffers across nearby regions. |
| 5 | `zlib_rs::deflate::algorithm::quick::deflate_quick` | 5.14% | Output gzip compression. Irreducible for current compression level; could be eliminated by writing uncompressed output or using a faster compressor (e.g., zstd). Low priority. |

Zlib inflate+deflate combined: ~56.5% of total CPU time. The decompression path (`inflate_fast_help_avx2` + `inftrees::inflate_table` + `Decompress::decompress_vec` = ~51.3%) dominates. The best remaining optimization target is reducing redundant decompressions via decode_buf reuse (Task 10) and bigWig reader buffer reuse (Task 9).

Application-level hotspots (`process_batch` 9.59% + `aggregate_slice` 3.23% + `write_matrix_value` 2.11% + `append_bins` 1.15%) sum to ~16.1%, down from the channel pipeline overhead in previous runs.

### Top 3 remaining allocation hotspots

No heaptrack data was captured for this run. Based on the perf call graph, the top allocation sites are:

| Rank | Source | Evidence | Optimization target |
|---|---|---|---|
| 1 | `alloc::raw_vec::RawVec::grow_one` under `BigWigReader::values` | 2.10% children in call graph; Vec growth during decompression buffer expansion | Task 9/10: pre-allocate decode buffers and reuse across calls |
| 2 | `_int_malloc` (glibc) | 2.55% self time; general heap allocation from repeated Vec/String creation in the hot path | Task 9: bigWig reader buffer reuse would reduce malloc pressure |
| 3 | `_int_realloc` + `_int_free_create_chunk` + `cfree` (glibc) | 0.62% + 0.51% + 0.74% = 1.87% combined; realloc churn from growing buffers | Same as above — pre-sized buffers would eliminate realloc chains |

### Key observations

1. **Wall clock is stable** (+0.5% vs previous, -23.5% vs baseline): The chunk-collect rewrite did not regress wall clock performance. The perf stat wall clock shows a slight improvement (-1.3% vs previous) within noise margin.

2. **System time dramatically reduced** (-55.9% vs previous, -63.5% vs baseline): The chunk-collect dispatch replaces the sync_channel + BTreeMap writer thread with direct chunk collection. This eliminates the futex syscalls from channel send/recv and the kernel-level synchronization overhead that was the primary cost of the channel pipeline.

3. **Voluntary context switches collapsed** (189,810 -> 257, -99.9%): The writer thread's channel-blocking pattern is gone. The 257 remaining voluntary context switches are from standard rayon thread pool scheduling. This confirms the channel was the source of the massive context switch count in previous runs.

4. **Max RSS regressed** (+34.0% vs previous, 515 MB -> 691 MB): The chunk-collect approach collects results in memory before dispatching to the writer, unlike the streaming channel pipeline that wrote rows as soon as they were available. This is a trade-off: lower synchronization overhead at the cost of higher peak memory. Still 12.0% below baseline (784 MB).

5. **Page faults increased** (+28.5% vs previous): Consistent with the RSS increase — more memory touched means more minor page faults. The increase tracks proportionally with the RSS delta.

6. **task-clock decreased** (-6.0% vs previous): Less total CPU time consumed. The sync_channel overhead (~2.5s of user time estimated in the ar2 report) is gone, and the system time savings (3.28s) is the largest single improvement.

### Verdict

**PASS** — Wall clock is within noise of previous (+0.5% by time -v, -1.3% by perf stat), well within the 5% regression threshold. The optimization successfully eliminates channel synchronization overhead (system time -55.9%, voluntary ctx switches -99.9%, task-clock -6.0%). The RSS regression (+34.0% vs previous) is a known trade-off of chunk collection vs streaming, and the absolute value (691 MB) remains 12.0% below the baseline (784 MB). No wall clock regression detected.
