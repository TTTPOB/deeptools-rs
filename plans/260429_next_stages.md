1. remove num_cpus dep
2. ram usage optimization, when try store the matrix data in memory, try do file spilling, we may pre compute sort metrics and use that as an index. when we need to write the matrix, we can mmap the file, and use index to accelerate sorting. what is the common way of doing this sort related file spilling in rust/c?
3. 
❯ 我对当前的两个python管理的regression脚本不太满意，它们和现有的测试体系不在一起，融合很差，可否做下调查给个意见？

● 我来调查一下项目的测试体系和 Python 回归脚本的现状。

● Explore(Explore test and regression scripts)
  ⎿  Done (40 tool uses · 72.6k tokens · 2m 48s)
  (ctrl+o to expand)

  Read 2 files (ctrl+o to expand)

● 调查完成，以下是我的分析和建议。

  现状分析

  你有两套完全割裂的测试体系：

  ┌──────────┬───────────────────────────┬───────────────────────────────────────────┐
  │   维度   │  Rust 测试 (cargo test)   │        Python 回归脚本 (scripts/)         │
  ├──────────┼───────────────────────────┼───────────────────────────────────────────┤
  │ 框架     │ 标准 #[test]              │ 自定义（非 pytest/unittest）              │
  ├──────────┼───────────────────────────┼───────────────────────────────────────────┤
  │ 场景数   │ 分散在 6 个文件的单元测试 │ 10 个兼容性场景                           │
  ├──────────┼───────────────────────────┼───────────────────────────────────────────┤
  │ 运行方式 │ cargo test                │ pixi run python scripts/custom_compare.py │
  ├──────────┼───────────────────────────┼───────────────────────────────────────────┤
  │ CI 集成  │ 无（CI 只做 build）       │ 无                                        │
  ├──────────┼───────────────────────────┼───────────────────────────────────────────┤
  │ 代码量   │ 少量 inline 测试          │ ~1500 行 Python + YAML 配置               │
  ├──────────┼───────────────────────────┼───────────────────────────────────────────┤
  │ 依赖     │ 无额外依赖                │ pixi + deepTools 3.5.6 + pyyaml           │
  └──────────┴───────────────────────────┴───────────────────────────────────────────┘

  核心问题：Python 脚本做了大量基础设施工作（命令执行、缓存、矩阵解析、比较、报告），但完全游离于 cargo test 之外，开发者需要手动知道并运行它们。

  建议方案：将兼容性测试迁移为 Rust 集成测试

  关键发现：参考矩阵（10 个 .mat 文件）已经 vendor 在 deeptools/deeptools/test/test_heatmapper/ 目录下。这意味着核心的兼容性验证完全不需要 Python —— 只需要：

  1. 用相同参数运行 Rust 二进制
  2. 解析输出矩阵
  3. 与参考矩阵按 5e-6 容差比较

  具体做法：

  tests/
    python_compatibility.rs     # 10 个 #[test] 函数，每个对应一个场景

  每个测试函数：
  - 调用 Command::new(env!("CARGO_BIN_EXE_compute_matrix_rs")) 运行编译好的二进制
  - 输出到 tempfile
  - 用已有的矩阵解析逻辑（gzip 解压 + JSON header + 行比较）
  - 按 tolerance 断言

  优势：
  - 所有测试统一到 cargo test，一个命令跑完
  - 消除 ~1500 行 Python 基础设施代码
  - 可以轻松加入 CI（不需要 pixi/deepTools）
  - ENCODE 大数据集测试可用 #[ignore] 标记，按需运行

  Python 脚本的保留价值：
  - custom_compare.py 中的性能对比功能（Rust vs Python 计时）在开发阶段仍有用
  - 可以精简为一个纯粹的 benchmark/性能对比工具，而非测试入口

  权衡

  ┌────────────┬─────────────────────────────────────┬──────────────────────────────┐
  │            │        迁移到 Rust 集成测试         │ 保持 Python 但用 pytest 改写 │
  ├────────────┼─────────────────────────────────────┼──────────────────────────────┤
  │ 运行便捷性 │ cargo test 一条命令                 │ 仍需 pixi 环境               │
  ├────────────┼─────────────────────────────────────┼──────────────────────────────┤
  │ CI 可行性  │ 直接可用                            │ 需要配置 pixi                │
  ├────────────┼─────────────────────────────────────┼──────────────────────────────┤
  │ 维护成本   │ 矩阵解析需在 Rust 侧重写（~200 行） │ 改写为 pytest 约 ~300 行     │
  ├────────────┼─────────────────────────────────────┼──────────────────────────────┤
  │ 性能对比   │ 需单独工具                          │ 天然支持                     │
  └────────────┴─────────────────────────────────────┴──────────────────────────────┘

  我的建议是走 Rust 集成测试路线。矩阵解析逻辑不复杂（gzip 解压 + JSON header + tab 分隔值），你的项目已经有 flate2 和 serde_json 依赖，代码量会很少。Python
  脚本保留一个精简版做性能 benchmark 用就够了。

  要我出一个具体的实施计划吗？

我的想法是加一个新的只在测试或开发时编译的binary做matrix comparison.就是把当前的数值对比和header对比的东西搬进去。这个Module同时也作为测试的时候，对比两个矩阵是否相同的基础设施来使用
关于rust和Python comparison, 可以考虑用shell脚本，或者是换用另一种Python测试写法，例如你提到的pytest来做