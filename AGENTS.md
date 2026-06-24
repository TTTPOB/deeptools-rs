## 项目状态
- 目标：Rust 版 deepTools `computeMatrix`，覆盖 `reference-point` 和 `scale-regions`。
- 数值目标：输出矩阵和 deepTools 3.5.6 在支持用例内保持 ≤5e-6 绝对误差。
- 当前状态：主流程可用，metagene、GTF/BED12、blacklist、skipZeros、threshold、auxiliary outputs 都有回归覆盖。

## 当前分支状态
- 当前 harness 入口是 `scripts/harness.py`。
- 配置位于 `scripts/configs/`，按任务拆成多个小 JSON。
- CI 在 release build 后运行 `cargo test --test python_compatibility -- --test-threads=1`。
- 最近确认过的全量兼容性入口：`pixi run compat --quiet`，结果为 26/26 passed。

## 工具链
- Python/deepTools 环境走 pixi。`pixi.toml` 使用 `conda-forge` 和 `bioconda`，锁定 `deeptools ==3.5.6`。
- 不要直接调用 `python3 scripts/harness.py ...` 做常规验证；使用 pixi task。
- 不要运行 `cargo update`。需要依赖信息时，先用本地源码和 `rg` 查。
- 不要依赖网络一定可用。ENCODE 数据准备命令会下载大文件，运行前确认用户确实需要。

## 验证入口
- 兼容性全量验证：`pixi run compat`
- 单个兼容性用例：`pixi run compat --case <case_id>`
- Rust 集成测试：`cargo test --test python_compatibility -- --test-threads=1`
- 全部 Rust 测试：`cargo test`
- 重新生成 Python 参考产物：`pixi run regen-artifacts`
- 验证参考产物：`pixi run verify-artifacts`
- 性能烟测：`pixi run bench-smoke`
- ENCODE 数据准备：`pixi run prepare-data encode_k562_atac`
- ENCODE 性能测试：`pixi run encode`

## 配置文件
- 共享路径和比较参数：`scripts/configs/common.json`
- 兼容性用例：`scripts/configs/compat/core.json`、`blacklist.json`、`corner.json`、`metagene.json`
- 参考产物用例：`scripts/configs/artifacts.json`
- 性能和 ENCODE 用例：`scripts/configs/benchmarks.json`
- 数据集下载描述：`scripts/configs/datasets.json`
- 新增用例时，把它放进对应的小 JSON。不要重新合并成一个大文件。

## 架构
- 两条执行路径：
  - `StreamOrdered`：`sort=no` 或 `sort=keep` 且输入已经排序。
  - `HybridBucket`：`sort=keep` 需要重排，或 `ascend` / `descend`。
- 所有输出通过 `FileCollector`：主 gzip 矩阵、`outFileNameMatrix`、`outFileSortedRegions`。
- 大矩阵排序使用文件 spilling：`HybridBucketCollector` 按 group 收集行，超过阈值写临时文件，finalize 时 mmap 读回。
- `matrix_compare` 是共享库模块，负责解析 plain/gzip/multi-member matrix，比较 header、BED 字段和数值。
- `compare_matrix` binary 只作为开发 CLI。测试代码直接调用 `matrix_compare`。

## 兼容性事实
- Rust 集成测试从 `scripts/configs/common.json` 和 `scripts/configs/compat/*.json` 读取用例。
- 当前兼容性任务包含 26 个 manifest case。
- 默认比较忽略 header 里的 `proc number` 和 `scale`。
- `scale` header：Python 默认值可能写 int，Rust 总是写 float。下游工具不读回这个字段。
- `outFileNameMatrix` 第一行在 finalize 时原地重写 group counts；Rust 可能用 trailing spaces 填满预留宽度，解析兼容。
- `outFileSortedRegions` 的 BED12 `blockStart` 使用相对 `chromStart` 的标准 offset；Python/deeptoolsintervals 写自身输出格式里的坐标值。该差异按 intended 记录。
- deepTools 的多字符短 flag `-bs` 和 `-bl` 不受 clap 支持；Rust CLI 使用 `--bs/--binSize` 和 `--bl/--blackListFileName`。

## 性能基线
- 历史 ENCODE K562 ATAC，4 cores，bin 10：
  - reference-point center ±100 bp：Python 171.35s，Rust 17.90s，约 9.57x。
  - scale-regions body 200 ±100 bp unscaled 50/50：Python 346.56s，Rust 18.64s，约 18.59x。

## 编辑规则
- 代码注释写英文。
- 常规 harness 验证使用 pixi task。只有调试 `scripts/harness.py` 自身时才直接调用 Python。
- 改动后按影响面跑验证。涉及 harness/config 时至少跑 `pixi run compat --case reference_point_basic` 和相关 task。
- 提交要原子化，不要把配置迁移、逻辑变更、格式化、文档改动混成一个提交。
- 若改动会影响用例清单，更新对应 `scripts/configs/*.json`，再跑 pixi task。
