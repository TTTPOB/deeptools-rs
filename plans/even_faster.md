# compute_matrix_rs 架构/性能/内存 优化建议

## Context

基于 perf flamegraph + heaptrack 对 `scale-regions` 模式在 134,900 个区域、3 个 bigWig 文件、binSize=10 条件下的 profile 分析。测试使用 `-p 4` 和 ENCODE K562 ATAC-seq 真实数据。

**关键统计数据:**

| 指标 | 134,900 区域 (3样本) | 269,800 区域 (2样本) |
|------|---------------------|---------------------|
| 总运行时间 | 10.6s | ~24s |
| 峰值堆内存 | 460 MB | 704 MB |
| 峰值 RSS | 594 MB | 900 MB |
| 总分配次数 | 6,756,127 | 11,550,737 |
| 临时分配 | 2,159,984 (32%) | 3,652,587 (31.6%) |
| 泄漏 | 22 KB (可忽略) | 22 KB |

---

## CPU 热点分析 (perf report, 1,960 samples)

### Top 热点函数 (self time)

| 排名 | 占比 | 函数 |
|------|------|------|
| 1 | **49.64%** | `BigWigReader::values` — bigWig 数据读取+解压 |
| 2 | **9.13%** | `process_batch` — 批处理调度 |
| 3 | **7.91%** | `deflate_quick` (zlib-rs) — gzip 输出压缩 |
| 4 | **5.66%** | `aggregate_slice` — bin 聚合计算 |
| 5 | **5.36%** | `write_matrix_value` — 浮点数格式化 |
| 6 | **2.81%** | `zune_inflate::build_decode_table_inner` — zlib 解压表重建 |
| 7 | **1.17%** | `append_bins` — zone 构建 |

### BigWigReader::values 内部细分

- **zune_inflate 解压**: ~3.0% (包括 `build_decode_table_inner` 2.81%)
- **内存分配/重分配**: ~4.6% (`grow_one` + `finish_grow` + `realloc` + `memset`)
- **memcpy**: ~1.7%
- **malloc/calloc**: ~1.1%

### 输出阶段细分

- **gzip 压缩** (`deflate_quick`): 7.91% — zlib-rs 的 quick 算法
  - 注意: libdeflate 在这里测过，zlib-rs 稍微好点
- **浮点数格式化** (`write_matrix_value`): 5.36%
  - `__divti3` + `u128_div_rem`: 2.55% (128位整数除法)
  - `__fixdfti`: 0.87% (f64→i128 转换)
  - `rint`: 0.82% (四舍五入)

---

## 内存热点分析 (heaptrack)

### 主要分配来源

1. **`BigWigReader::values`** — 最大的分配来源
   - 每线程每次调用: `work_buf` (decompression buffer, ~128KB 分配)
   - 每区间: `Vec::push` → `grow_one` 触发 realloc (~346K 次)
   - `get_or_cache_block` 每次返回 `data.clone()` — 缓存命中时也 clone
   - `values` 返回的 `Vec<BigWigValue>` 每次都是新分配的

2. **`load_groups`** — 加载所有 BED 记录
   - 每个 BedRecord 包含: `chrom: String`(独立分配), `name: Option<String>`, `extra_fields: Vec<String>`
   - 每个 RegionTask 再 clone 一次完整的 BedRecord（数据重复存储）

3. **`process_batch`** — 批处理阶段
   - 每个 batch 每个样本: `Vec<f32>` 覆盖窗口 (window_span 个 float)
   - 每个区域: `all_values` (sample_count * bin_count 个 float)

### 临时分配分析 (32% 的分配是临时的)

主要来源:
- `BigWigReader::values` 内部的 Vec 扩容
- `work_buf` 每次分配后再 resize
- `get_or_cache_block` 的 `data.clone()` 返回值

### 内存占用组成估算 (134,900 区域)

| 组成部分 | 估算大小 |
|----------|---------|
| BedRecord 存储 (含 strings) | ~100 MB |
| WorkItem 列表 | ~30 MB |
| CoalescedBatch 覆盖缓冲 (per-worker) | ~20 MB |
| 矩阵值 (流式写入, 不常驻) | ~5 MB (缓冲) |
| BigWig reader caches (4 workers × CIR + block caches) | ~50 MB |
| Streaming 临时文件 spool | ~40 MB |
| 其他 (zlib buffers, 字符串重复等) | ~215 MB |
| **总计** | ~460 MB |

---

## 已确认的优化项 (用户确认)

### ✅ A1. zune-inflate → zlib-rs 替换

**问题:** `zune-inflate` 占用 ~3% CPU，每次 `DeflateDecoder::new()` 都重建查找表 (2.81%)。

**方案:** 用 zlib-rs (项目已有依赖) 的 deflate 解压 API 替换。依次实现 zlib-rs 和 C deflate 两个变体，每个变体跑 `profile_bench.sh` 对比 baseline，选较优者。

- 文件: `src/io/readers/bwig.rs:451-452`
- 验证: 每个变体跑 `profile_bench.sh`，对比 perf stat task-clock 和 heaptrack 解压路径的分配

### ✅ A2. block_cache 用 Arc<[u8]> 替代 Vec<u8> clone

**问题:** `get_or_cache_block` (bwig.rs:408-409) 缓存命中时 `data.clone()` 完整复制 buffer。

**方案:** `HashMap<(u64, u64), Arc<[u8]>>`，返回 `Arc<[u8]>` 的 clone（仅引用计数增加）。

- 文件: `src/io/readers/bwig.rs:297, 401-423`

### ✅ A3. work_buf 复用

**问题:** `values()` 每次调用都 `vec![0u8; uncompress_buf_size]`。52,976 batch × 3 样本 = ~159K 次 128KB 分配。

**方案:** 将 `work_buf` 放入 `BigWigReader` 结构体字段，跨调用复用。

- 文件: `src/io/readers/bwig.rs:322-348`

---

### ✅ P1. write_matrix_value: 消除 128 位整数除法

**问题:** `write_matrix_value` 消耗 5.36%，其中 2.55% 是 `__divti3`（128位除）。用 i128 存 `f64 * 1e6`，但基因组学数据的 f32 值绝对值几乎总 `< 1e7`。

**方案:** 快速路径用 i64:
```rust
if value.abs() < 1e7 && value.is_finite() {
    let scaled = (value as f64 * 1_000_000.0 + 0.5) as i64;
    // 整数和小数部分用 64-bit 除
}
```
超出范围 fallback 到 i128 或 format!()。

- 文件: `src/io/writers/mod.rs:495-543`
- 预期: 2-3% CPU 节省

### ✅ P2. aggregate_slice 优化

**问题:** 消耗 5.66% CPU，~81M 次调用。每次 f32→f64→f32 累加，对 binSize=10 的小切片来说 f64 精度过剩。

**方案:**
1. 用 f32 累加替代 f64（mean/sum/min/max 精度足够）
2. 小切片 (<16 元素) 直接展开循环
3. NaN 检测可以用更高效的方式

- 文件: `src/pipeline/core/mod.rs:449-531`
- 预期: 1-3% CPU 节省

### ✅ P3. 去除 coalescing 对稀疏 BED 的开销 — 用户新方案

**背景:**
- Coalescing 通过合并相邻 query window 为一个大范围 bigWig read 来减少 I/O 次数
- 当前 `estimate_coalesce_gap` 将阈值 clamp 在 [100, 2000]
- 对于**密集** BED（区域间距 < 2000bp），合并有效
- 对于**稀疏** BED（区域间距 > 2000bp），gap 打到 clamp 上界，合并不发生但**覆盖缓冲仍然被分配**，且每次合并的大窗口 read 解压更多无关数据

**用户的新方案:**
当 coalesce gap 达到 clamp 上界 (2000) 时，判定为稀疏数据集。此时:
1. **跳过 coalescing** — 避免大窗口解压，直接 per-item 读取（更省内存和 CPU）
2. **仍然使用 streaming gzip 输出** — 因为正好利用已有的 `StreamingMatrixWriter`（header placeholder + 流式写入 + 最后 patch header）

**实现方式 (B+C):**
- `estimate_coalesce_gap` 照常返回 gap 值（不改变其职责）
- `create_batches` 入口处判断 gap ≥ 2000，若是则设置 `CoalesceStrategy` enum:
  ```rust
  enum CoalesceStrategy {
      Coalesce(i64),    // gap threshold, normal path
      NoCoalesce,       // sparse: 1 item = 1 batch
  }
  ```
- `CoalesceStrategy` 作为 pipeline 属性，后续代码据此决策

**streaming 生效条件:**
只有当 `sort_regions` ∈ {Keep, No} 时 streaming 才可用。若 `sort_regions` 为 `Ascend`/`Descend`，则 `sort_groups` 需要对所有行按 signal 排序→随机访问→必须在内存中，streaming 不可用。

| sort_regions | streaming 可用？ | 原因 |
|---|---|---|
| `Keep`（默认） | 是 | sort_groups 跳过 |
| `No` | 是 | sort_groups 跳过（当前代码未允许，需修） |
| `Descend` | 否 | sort_groups 需要所有行在内存 |
| `Ascend` | 否 | 同上 |

最终决策:
```
coalesce_gap >= 2000 AND sort_regions ∈ {Keep, No} → CoalesceStrategy::NoCoalesce + streaming
否则 → CoalesceStrategy::Coalesce(gap) + 按原逻辑选 streaming 或 in-memory
```

- 文件: `src/pipeline/core/mod.rs:843-878`, `src/pipeline/core/mod.rs:1179-1180`
- 收益: 对稀疏 BED 数据集的 CPU 和内存都有改善

---

### ✅ M1. BedRecord 用 Arc 共享避免 clone

**问题:** 每个 BedRecord 被 clone 进入 RegionTask。269,800 条记录 = 所有 String 字段重复存储。

**方案:** `RegionTask.record: Arc<BedRecord>`，clone 只增加引用计数。

- 文件: `src/pipeline/core/mod.rs:168-174`, `src/pipeline/reference_point.rs:178`, `src/pipeline/scale_regions.rs:204`
- 预期: 节省 50-100 MB

### ✅ M2. Chromosome 名称 interning

**问题:** 269,800 条 BED 记录可能只有 ~25 个唯一染色体名，但每条记录独立分配 String。

**方案:** 小型 interner: `HashMap<String, Arc<str>>`，解析 BED 时将 chrom 字段 intern。

- 文件: `src/io/readers/bed.rs:80`, `src/pipeline/core/mod.rs:81`
- 预期: 节省 ~20 MB

### ✅ M3. CIR node cache 添加大小限制

**问题:** `cir_node_cache: HashMap<u64, Arc<CachedCirNode>>` 无淘汰策略。

**方案:** 添加 `MAX_CIR_CACHE_ENTRIES` 常量（如 1000），超限时 clear 或用 LRU。

- 文件: `src/io/readers/bwig.rs:296`

### ✅ M4. sample_coverages buffer 复用

**问题:** 每个 batch 都分配 `Vec<Vec<f32>>` (core/mod.rs:989)，大窗口 batch 的分配开销高。每个 sample 一个 `window_len` 大小的 f32 数组，rayon 并行时多个 worker 同时分配。

**方案:** `thread_local!` buffer pool — 每个 rayon worker 线程持有一份 `RefCell<Vec<Vec<f32>>>`，跨 batch 复用:
- capacity 够 → `clear()` + 重新填充
- capacity 不够 → resize
- 4 个 worker × 各独占一份 buffer，无锁
- 与 rayon `map_init` 的 per-thread 初始化模式自然契合

```rust
thread_local! {
    static COVERAGE_POOL: RefCell<Vec<Vec<f32>>> = RefCell::new(Vec::new());
}

fn get_coverage_buffers(sample_count: usize, window_len: usize, default_fill: f32) -> Vec<Vec<f32>> {
    COVERAGE_POOL.with(|pool| {
        let mut bufs = pool.borrow_mut();
        bufs.resize_with(sample_count, || Vec::new());
        for buf in bufs.iter_mut() {
            buf.clear();
            buf.resize(window_len, default_fill);
        }
        // return borrowed data... 
    })
}
```

注意: 返回值不能 borrow `RefCell`。方案是用 `mem::take` + process 完后放回 pool，或直接在 `map_init` 闭包内借用 pool 并就地消费。

- 文件: `src/pipeline/core/mod.rs:989-1033`
- 收益: 减少大批次场景下的分配开销

---

## 架构层面发现和建议

### AR1. 消除未使用的依赖和死代码

| 项目 | 状态 |
|------|------|
| `crossbeam-channel = "0.5"` | Cargo.toml 中声明，源码中零引用 |
| `bigtools = "0.5.6"` | 仅被未使用的 `io/readers/bigwig.rs` 引用 |
| `src/io/readers/bigwig.rs` | 不在模块树中 (`readers/mod.rs` 不声明它) |

**方案:** 删除 `src/io/readers/bigwig.rs`，从 Cargo.toml 移除 `bigtools` 和 `crossbeam-channel`（见下一条 AR2 的讨论）。

### AR2. 计算/IO 线程分离 — 统一 Channel 流水线

**当前状态:**

```
rayon workers (并行处理 batch)
  → 所有结果收集到 result_slots (Phase 5 scatter) ← 全量缓冲
  → 主线程顺序遍历 result_slots (Phase 6)
    → FileCollector::on_row (streaming) 或
    → InMemoryCollector::on_row (in-memory)
  → collector.finalize(header)
```

`execute_mode` 签名: `fn execute_mode(tasks, collector, header_builder, ...) -> Result<C::Output>`

- collector 和 header_builder 传入 execute_mode，内部消耗
- 主线程阻塞在 rayon pool.install() 直到所有 batch 处理完
- 写入完全串行在 Phase 6

**目标架构 — 统一 Channel 流水线:**

两条路径用同一套 channel 基础设施，collector 移入 writer thread。

```
scale_regions.rs:
  创建 collector (FileCollector / InMemoryCollector)
  创建 sync_channel(256)
  spawn writer thread:
    ← 从 rx 接收 (orig_idx, group_index, Option<MatrixRow>)
    ← BTreeMap 重排保证 on_row() 按 orig_idx 顺序
    ← 统计 group_counts → header_builder → collector.finalize()
  调用 execute_mode(tasks, tx, ...)
  join writer thread
```

`execute_mode` 签名变化:
```rust
// 之前
fn execute_mode(tasks, collector, header_builder, ...) -> Result<C::Output>

// 之后 — 不关心谁在消费行，只管安排 worker 处理并塞进 channel
fn execute_mode(tasks, tx: Sender<...>, ...) -> Result<()>
```

**Writer thread 核心逻辑:**

```rust
let mut next_idx = 0;
let mut pending: BTreeMap<usize, (usize, Option<MatrixRow>)> = BTreeMap::new();
let mut group_counts = vec![0; group_count];

for (orig_idx, group_index, row) in rx {
    pending.insert(orig_idx, (group_index, row));
    while let Some(entry) = pending.remove(&next_idx) {
        let (grp, row_opt) = entry;
        if let Some(row) = row_opt {
            collector.on_row(row)?;
        }
        group_counts[grp] += 1;
        next_idx += 1;
    }
}
let header = header_builder(group_counts)?;
collector.finalize(header)
```

**对两条路径的效果:**

| 路径 | 变化 |
|------|------|
| Streaming | rows 直接写入 gzip，gzip 压缩/IO 与 rayon compute 重叠。内存只保留 BTreeMap 乱序缓冲（远小于 result_slots） |
| In-memory | 代码路径统一，不再需要 result_slots。性能无变化（sort_groups 仍需全量数据），但结构简化 |

**与 P3 的联动:**
- 稀疏数据 → CoalesceStrategy::NoCoalesce + streaming → writer thread 直接写 gzip
- 密集数据 → CoalesceStrategy::Coalesce(gap) → 仍走 channel，但 batch 窗口大，result_slots 消除省内存

**channel capacity:** 初始设为 256，后续 profiling 有需要再调整。

- 文件: `src/pipeline/core/mod.rs:1179-1250`, `src/pipeline/mod.rs:51-62`
- 收益: streaming 路径 compute/IO 重叠；消除 result_slots 全量缓冲

### AR3. 为什么 sort_regions != "keep" 禁用 streaming

`sort_groups` (matrix.rs:275-338) 需要所有行在内存中按 group 排序。streaming 路径中行已序列化到 gzip 文件，无法重新排序。

streaming 实际可用的条件: `sort_regions` ∈ {Keep, No}，两者都跳过 `sort_groups`。当前代码 streaming 检查只允许 `Keep`，这是一个 bug — `No` 也应该允许 streaming。需要在 streaming 条件判断处修复。

---

## T0: Profile Harness — 标准化性能报告

每个 worker 任务完成后，通过统一脚本采集性能数据并生成结构化报告，agent 自行对比前后报告判断 PASS/FAIL。

### 脚本: `scripts/profile_bench.sh`

```bash
#!/bin/bash
# Usage: ./scripts/profile_bench.sh <name> <target> <hot-path> -- <command>
# Example: ./scripts/profile_bench.sh p1-i64-div \
#              "消除 i128 除法，降低 write_matrix_value CPU" \
#              "write_matrix_value → rint + __divti3" \
#              -- cargo run --release -- ...

TIMESTAMP=$(date +%Y-%m-%d-%H%M%S)
NAME="$1"
TARGET="$2"
HOT_PATH="$3"
shift 3

REPORT="bench_reports/${TIMESTAMP}-${NAME}.md"
mkdir -p bench_reports

cat > "$REPORT" << EOF
# Profile: ${NAME}
Time: $(date -Iseconds)
Command: $@
Target: ${TARGET}
Hot path: ${HOT_PATH}

EOF

# /usr/bin/time -v
echo "## /usr/bin/time -v" >> "$REPORT"
echo '```' >> "$REPORT"
/usr/bin/time -v "$@" 2>> "$REPORT"
echo '```' >> "$REPORT"

# perf stat
echo "## perf stat" >> "$REPORT"
echo '```' >> "$REPORT"
perf stat -d "$@" 2>> "$REPORT"
echo '```' >> "$REPORT"

# heaptrack
echo "## heaptrack" >> "$REPORT"
echo '```' >> "$REPORT"
heaptrack "$@" 2>> "$REPORT"
echo '```' >> "$REPORT"

# Cleanup raw profiler files
rm -f perf.data heaptrack.*.gz 2>/dev/null

echo "Report written to $REPORT"
```

### 报告目录结构

```
bench_reports/
  2026-04-28-143215-baseline.md
  2026-04-28-150830-p1-i64-div.md
  2026-04-28-154500-a2-arc-cache.md
  ...
```

### 报告模板（agent 首次运行后提炼并固定格式）

```markdown
# Profile: <name>
Time: <ISO timestamp>
Command: <cmd>
Target: <这个优化在解决什么问题>
Hot path: <涉及的函数/代码路径>

## /usr/bin/time -v
| Metric | Value |
|--------|-------|
| Elapsed (wall clock) | X.X s |
| Max RSS | X KB |

## perf stat
| Metric | Value |
|--------|-------|
| task-clock | X ms |
| cycles | X |
| instructions | X |
| IPC | X.XX |
| cache-misses | X |

## heaptrack
| Metric | Value |
|--------|-------|
| Peak heap | X MB |
| Total allocations | X |
| Temporary | X (X%) |

## Comparison vs baseline
| Metric | Before | After | Delta |
|--------|--------|-------|-------|
| ... | ... | ... | ...% |

## Other hotspots observed
<!-- 当前目标路径之外的 CPU/内存热点 (perf top-10 函数及占比) -->
| Rank | Function | % CPU |
|------|----------|-------|
| 1 | ... | ... |
| ...

## Verdict
PASS/FAIL — <具体说明>
```

### Agent 使用协议

1. 读基准报告 → 提取 baseline 指标
2. 读当前报告 → 提取 structured 表格 → 对比 delta
3. 聚焦 `Hot path` 判断目标瓶颈是否改善
4. 填写 `Other hotspots observed` — 记录目标路径之外观察到的热点，供最后一轮汇总判断是否还有可优化项
5. 写 `Verdict`: PASS (改善或退化 ≤5%) 或 FAIL (退化 >5%)，含具体对比数字
6. 报告提交到 git（原始 perf.data / heaptrack.*.gz 已在脚本中删除，不保留）

---

## 实施优先级

### Task 依赖总览

| Task | 优化项 | Phase | 依赖 | 文件 |
|------|--------|-------|------|------|
| T0 | Profile harness | 0 | 无 | 新文件 `scripts/profile_bench.sh` |
| T1 | AR1: 清除 dead deps | 1 | 无 | `Cargo.toml`, `bigwig.rs` |
| T2 | P1: i64 除法快速路径 | 1 | T0 | `writers/mod.rs:495-543` |
| T3 | M2: Chromosome interning | 1 | 无 | `bed.rs:80`, `core/mod.rs:81` |
| T4 | A2: block_cache Arc | 1 | 无 | `bwig.rs:297,401-423` |
| T5 | A1: zune-inflate → zlib-rs | 2 | T0 | `bwig.rs:451-452` |
| T6 | P2: aggregate_slice 优化 | 2 | T0 | `core/mod.rs:449-531` |
| T7 | A3: work_buf 复用 | 2 | 无 | `bwig.rs:322-348` |
| T8 | M1: Arc\<BedRecord\> | 2 | 无 | `core/mod.rs:168-174`, ref_point, scale_regions |
| T9 | P3: 稀疏检测 + CoalesceStrategy | 3 | T7, T8 | `core/mod.rs:843-878,1179-1180` |
| T10 | AR2: Channel 流水线 | 3 | T9 | `core/mod.rs:1179-1250`, `pipeline/mod.rs:51-62` |
| T11 | M4: sample_coverages 复用 | 3 | 无 | `core/mod.rs:989-1033` |

### Phase 0 — 基础设施

| 优化 | 预期收益 | 改动量 |
|------|---------|--------|
| T0: Profile harness | 标准化性能对比流程 | 新脚本 + bench_reports/ |

### Phase 1 — 快速获胜 (独立、低风险)

| 优化 | 预期收益 | 改动量 |
|------|---------|--------|
| P1: i64 除法 | 2-3% CPU | ~20 行 |
| M2: Chromosome interning | ~20 MB 内存 | ~50 行新模块 |
| A2: block_cache Arc | 1-2% CPU + less alloc | ~10 行 |
| AR1: 清除 dead deps | 编译时间缩短 | 删除 2 个文件 |

### Phase 2 — 核心优化

| 优化 | 预期收益 | 改动量 |
|------|---------|--------|
| A1: zune-inflate → zlib-rs | 2-3% CPU | ~30 行, 需 benchmark |
| P2: aggregate_slice 优化 | 1-3% CPU | ~50 行 |
| A3: work_buf 复用 | less alloc | ~15 行 |
| M1: Arc<BedRecord> | 50-100 MB | ~30 行 (类型替换) |

### Phase 3 — 架构改进

| 优化 | 说明 |
|------|------|
| P3: 稀疏 BED 跳过 coalescing + streaming | CoalesceStrategy enum + sort_regions 条件判断 |
| AR2: 计算/IO channel 分离 | 统一 sync_channel 流水线，消除 result_slots |
| M4: sample_coverages 复用 | thread_local buffer pool |

### 不做的

| 项目 | 原因 |
|------|------|
| libdeflate 替换 zlib-rs | 用户已测试，zlib-rs 稍好 |
| Coalesce gap 上限提高 | 用户已测试，对稀疏 BED 并不好，可能增加内存和 CPU |
| 流式 chunk 临时文件消除 | 当前 streaming 模式 (`StreamingMatrixWriter`) 已经实现了增量写入 gzip，不需要 spool。spool_rows 只用于非 streaming path 的 `should_use_streaming` fallback |

---

## 验证方法

1. `scripts/custom_compare.py` 回归测试确保输出数值一致性 (tolerance 5e-6)
2. `scripts/profile_bench.sh` 标准化性能报告，含 perf stat + heaptrack + /usr/bin/time -v
3. 多规模测试: 1万 / 10万 / 100万 区域
4. A1 (zune-inflate → zlib-rs) 需额外解压 benchmark 对比 zune-inflate vs zlib-rs vs bigtools C deflate
