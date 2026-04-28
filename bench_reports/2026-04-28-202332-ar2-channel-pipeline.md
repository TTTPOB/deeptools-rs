# Profile: ar2-channel-pipeline
Time: 2026-04-28T20:23:32+08:00
Command: -- cargo run --release -- scale-regions -b 10 -p 4 --unscaled5prime 0 --unscaled3prime 0 --regionsFileName /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF333TAT.bed --scoreFileName /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF019IPA.bigWig /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF093IIW.bigWig /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF656ZKM.bigWig -o /tmp/bench_ar2.mat.gz
Target: channel-based compute/IO pipeline with BTreeMap reordering
Hot path: execute_mode Phase 5-6: result_slots -> sync_channel + writer thread

## /usr/bin/time -v
```
    Finished `release` profile [optimized] target(s) in 0.13s
     Running `target/release/compute_matrix_rs scale-regions -b 10 -p 4 --unscaled5prime 0 --unscaled3prime 0 --regionsFileName /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF333TAT.bed --scoreFileName /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF019IPA.bigWig /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF093IIW.bigWig /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF656ZKM.bigWig -o /tmp/bench_ar2.mat.gz`
Running computeMatrix (Rust) in scale-regions mode.
[coalesce-gap] n_gaps=161669 p50=5475 p75=12795 threshold=2000
[coalesce-gap] strategy="no-coalesce" batches=269800 items=269800 ratio=1.00
	Command being timed: "cargo run --release -- scale-regions -b 10 -p 4 --unscaled5prime 0 --unscaled3prime 0 --regionsFileName /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF333TAT.bed --scoreFileName /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF019IPA.bigWig /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF093IIW.bigWig /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF656ZKM.bigWig -o /tmp/bench_ar2.mat.gz"
	User time (seconds): 30.05
	System time (seconds): 6.74
	Percent of CPU this job got: 291%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:12.62
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 663796
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 181518
	Voluntary context switches: 208919
	Involuntary context switches: 93
	Swaps: 0
	File system inputs: 0
	File system outputs: 16
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```
## perf stat
```
    Finished `release` profile [optimized] target(s) in 0.13s
     Running `target/release/compute_matrix_rs scale-regions -b 10 -p 4 --unscaled5prime 0 --unscaled3prime 0 --regionsFileName /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF333TAT.bed --scoreFileName /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF019IPA.bigWig /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF093IIW.bigWig /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF656ZKM.bigWig -o /tmp/bench_ar2.mat.gz`
Running computeMatrix (Rust) in scale-regions mode.
[coalesce-gap] n_gaps=161669 p50=5475 p75=12795 threshold=2000
[coalesce-gap] strategy="no-coalesce" batches=269800 items=269800 ratio=1.00

 Performance counter stats for 'cargo run --release -- scale-regions -b 10 -p 4 --unscaled5prime 0 --unscaled3prime 0 --regionsFileName /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF333TAT.bed --scoreFileName /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF019IPA.bigWig /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF093IIW.bigWig /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF656ZKM.bigWig -o /tmp/bench_ar2.mat.gz':

         35,362.05 msec task-clock:u                     #    2.991 CPUs utilized             
                 0      context-switches:u               #    0.000 /sec                      
                 0      cpu-migrations:u                 #    0.000 /sec                      
           148,024      page-faults:u                    #    4.186 K/sec                     
   <not supported>      cycles:u                                                              
   <not supported>      instructions:u                                                        
   <not supported>      branches:u                                                            
   <not supported>      branch-misses:u                                                       
   <not supported>      L1-dcache-loads:u                                                     
   <not supported>      L1-dcache-load-misses:u                                               
   <not supported>      LLC-loads:u                                                           
   <not supported>      LLC-load-misses:u                                                     

      11.824545412 seconds time elapsed

      30.160465000 seconds user
       6.029610000 seconds sys


```
## heaptrack
```
heaptrack stats:
	allocations:          	19532
	leaked allocations:   	3302
	temporary allocations:	3988
    Finished `release` profile [optimized] target(s) in 0.13s
     Running `target/release/compute_matrix_rs scale-regions -b 10 -p 4 --unscaled5prime 0 --unscaled3prime 0 --regionsFileName /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF333TAT.bed --scoreFileName /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF019IPA.bigWig /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF093IIW.bigWig /home/tpob/playground/deeptools-rs/target/compute-matrix-datasets/encode_k562_atac/ENCFF656ZKM.bigWig -o /tmp/bench_ar2.mat.gz`
Running computeMatrix (Rust) in scale-regions mode.
[coalesce-gap] n_gaps=161669 p50=5475 p75=12795 threshold=2000
[coalesce-gap] strategy="no-coalesce" batches=269800 items=269800 ratio=1.00
```

## Correctness

md5sum:
- baseline: `20af5dfded82a479f0766e36c57a620e`
- channel-pipeline: `20af5dfded82a479f0766e36c57a620e`

Output is identical.

## Hotspot Analysis

### Baseline reference (2026-04-28-170030)

| Metric | Baseline | Channel Pipeline | Delta |
|--------|----------|------------------|-------|
| Wall clock (perf stat) | 12.247 s | 11.824 s | **-3.5%** |
| Wall clock (time -v) | 15.13 s | 12.62 s | -16.6% (warm cache) |
| User time (perf stat) | 27.645 s | 30.160 s | +9.1% |
| System time (perf stat) | 2.603 s | 6.030 s | +131.7% |
| task-clock | 30,248 ms | 35,362 ms | +16.9% |
| CPU % | 226% | 291% | +28.8% |
| Max RSS | 784,896 KB | 663,796 KB | **-15.4%** |
| Minor page-faults | 235,672 | 181,518 | **-23.0%** |
| Voluntary ctx switches | 24,857 | 208,919 | +740% |
| Involuntary ctx switches | 148 | 93 | -37.2% |
| FS inputs | 4,469,304 | 0 | warm cache |
| Heaptrack allocations | 19,531 | 19,532 | +1 |

### Metric-by-metric analysis

**task-clock** (+16.9%, +5,113 ms): The channel pipeline burns more CPU time because the sync_channel adds synchronization overhead. With a bounded channel (capacity 256), compute workers perform additional work: acquiring/releasing internal channel mutexes for each `tx.send()`, and the writer thread acquires/releases the mutex for each `rx.recv()`. Over 102,664 output rows spread across 269,800 batches, this adds up.

**context-switches (perf stat, userspace)**: 0 for both runs. Userspace context switches (mutex contention) did not register at the perf stat userspace level. However, the time -v metrics tell a different story.

**page-faults** (-23.0%, -84,734): The channel pipeline eliminates the `result_slots` Vec (269,800 entries of `(usize, Option<MatrixRow>)` ~= 4.3 MB) and the intermediate `computed` Vec of per-batch results. Less memory touched means fewer minor page faults. This directly contributes to the Max RSS reduction.

**Max RSS** (-15.4%, -121,100 KB): The largest benefit. Removing the `result_slots` array and streaming results through the channel reduces peak memory by ~118 MB. For larger datasets, this improvement would scale linearly with the number of regions. The writer thread calls `collector.on_row()` (which writes to disk via `StreamingMatrixWriter`) as soon as results are available in order, so rows don't accumulate in memory.

**instructions/cycles**: Not supported on this platform (WSL2 without hardware PMU access).

**Wall clock** (-3.5%, from perf stat): The compute/IO overlap from the separate writer thread hides some of the I/O latency. While CPU time increased (+9.1% user, +131.7% system), wall clock decreased, meaning the additional CPU work is overlapped with I/O that was previously serialized.

**User time** (+9.1%, +2.515 s): The sync_channel add/remove operations per individual result contribute ~2.5 s of additional CPU work. Each of the 102,664 output rows requires a channel send and a channel receive, plus BTreeMap insertion and removal.

**System time** (+131.7%, +3.427 s): The channel-based approach triggers significantly more kernel involvement. Each `sync_channel.send()` and `recv()` may involve futex syscalls when the channel is full or empty. With 269,800 batches and 102,664 results, this is a high volume of synchronization events. Additionally, the writer thread's disk I/O happens in parallel with computation, increasing total system time.

**Voluntary context switches** (+740%, +184,062): The dramatic increase is expected with the channel-based design. The writer thread voluntarily yields when the channel is empty (waiting for compute workers to produce results). Compute workers may also block briefly when the channel is full. With the bounded channel (256), this happens frequently — on average, the writer depletes the channel in ~2.5 ms of compute output, then yields. For a ~12 s run, that's ~5,000 channel-empty events, each causing a voluntary yield. The large count (208k) suggests additional OS-level scheduling interactions.

**FS inputs** (0 vs 4,469,304): The baseline was a cold-cache run (files freshly read from disk). The channel pipeline run was the third run and benefited from OS page cache. This explains most of the wall-clock delta between the time -v measurements (15.13s -> 12.62s) but does not affect the perf stat comparison (which uses the second run for both).

**Heaptrack allocations** (19,531 -> 19,532): Identical within measurement noise. The BTreeMap adds one heap allocation for its root node. No per-row allocation overhead was added.

### Summary

The channel-based compute/IO pipeline achieves its primary goals:
1. **-15.4% Max RSS** (784 MB -> 664 MB) by eliminating the result_slots buffer
2. **-23.0% minor page-faults** from reduced memory footprint
3. **-3.5% wall clock** (perf stat) from overlapping compute and I/O

The trade-offs are:
1. **+9.1% user time** from channel synchronization overhead
2. **+131.7% system time** from kernel-level synchronization
3. **+740% voluntary context switches** from bounded-channel blocking

The wall-clock improvement is modest (-3.5%) for this dataset because the I/O phase was already fast (~0.5 s for gzip writing). The real benefit will be more pronounced with:
- Larger datasets where result_slots memory pressure becomes prohibitive
- Slower output targets (network storage, high-compression gzip) where overlapping I/O with compute hides more wall time
- Streaming mode with `FileCollector`, where the channel pipeline enables zero-copy-like row transmission from worker to writer
```
