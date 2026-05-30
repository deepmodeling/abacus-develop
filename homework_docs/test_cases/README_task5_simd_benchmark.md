# Task 5 SIMD Benchmark Workflow

这套脚本用于“题目 5：平面波 gather 操作 SIMD 向量化优化”的正式性能对比。它们的设计目标是：

- 同一套脚本可用于 baseline 原始分支和 SIMD 优化分支。
- 不依赖当前 SIMD 分支新增的源码函数、专有测试程序或专有 timer 名称。
- 如果日志中没有某个 timer 字段，汇总结果写为 `NA`，不因为字段缺失而退出。

## 文件说明

- `run_task5_simd_suite.sh`
  - 一键运行固定的题目 5 性能测试矩阵。
  - 作为 wrapper 调用 `run_task5_simd_benchmark.sh`，不会重复实现底层运行逻辑。
  - 为每个 `case x nproc x threads` 组合创建单独子目录，并生成 `suite_manifest.csv`、`suite_run_log.txt`。
- `run_task5_simd_benchmark.sh`
  - 运行指定 case 的 warmup 和正式 benchmark。
  - 保存每次运行的完整 `stdout/stderr`。
  - 记录环境信息、git 信息、OpenMP 环境变量和 MPI 版本。
  - 自动调用汇总脚本生成单次结果目录下的 CSV/Markdown。
- `collect_task5_simd_results.py`
  - 从一个或多个 benchmark 输出目录收集日志。
  - 生成原始结果 CSV、汇总结果 CSV 和报告可直接使用的 Markdown 表格。
  - 当同时输入 baseline 和 simd 两组目录，若配置匹配，会自动计算 `speedup_vs_baseline = baseline_median / current_median`。
- `collect_task5_simd_suite_results.py`
  - 扫描一个或多个 suite 输出目录下的全部子测试目录。
  - 生成 `suite_raw_results.csv`、`suite_summary.csv` 和 `suite_summary.md`。
  - 当同时输入 baseline suite 和 simd suite 时，会自动按同一机器、同一 case、同一 `nproc/threads` 计算加速比。

## 固定 suite 配置

`run_task5_simd_suite.sh` 在脚本开头固定了以下测试矩阵，便于你后续手动调整：

- Cases
  - `homework_docs/test_cases/gaas_small`
  - `homework_docs/test_cases/gaas_medium`
  - `homework_docs/test_cases/gaas_large`
- MPI 进程数
  - `1 2 4`
- OpenMP 线程数
  - `1 2 4 8`
- Warmup
  - `1`
- Repeat
  - `5`
- Timeout
  - `1800`

选择这两个 case 的原因是：

- `gaas_small` 比 `gaas_tiny` 更容易稳定观察到平面波 FFT / gather-scatter 开销，同时总运行时间仍然可控。
- `gaas_medium` 比 `gaas_small` 更能放大 gather/scatter 重排的影响，但又比再继续增大规模更适合作为整套矩阵的常规测试。
- `gaas_large` 在保持同一 GaAs 体系和相同物理参数的前提下，把 FFT 网格增大到 `64 x 64 x 64`，更适合作为题目 5 的高负载补充样例。
- `gaas_tiny` 太轻，不适合作为 SIMD 性能结论的主样例。
- `gaas_tiny_40Ry` 当前目录下输入与 `gaas_tiny` 没有形成更有代表性的额外负载，因此没有默认纳入 suite。

## 运行前准备

1. 手动切换到要测试的 git 分支。
2. 使用相同的 CMake 配置重新编译。
3. 确认 baseline 和 simd 两边使用相同的：
   - 编译器与编译选项；
   - MPI 实现；
   - OpenMP 线程设置；
   - 测试 case；
   - `nproc` / `threads` / `repeat`。

## 在 baseline 原始分支上运行

1. 切换分支。

```bash
git checkout <baseline-branch-or-commit>
```

2. 用和 SIMD 分支一致的 CMake 参数重新编译。推荐优先使用下面这套偏性能测试、且便于复现的配置：

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DENABLE_MPI=ON \
  -DUSE_OPENMP=ON \
  -DENABLE_NATIVE_OPTIMIZATION=OFF
cmake --build build -j4
```

如果实体机上没有装好 ELPA，通常这套配置也可以直接工作；如果你遇到 CMake 依赖检查失败，再补充：

```bash
-DUSE_ELPA=OFF
```

最重要的是 baseline 和 SIMD 两边必须使用完全相同的 CMake 参数。

3. 运行 suite 脚本，`--label` 设为 `baseline`。

```bash
bash homework_docs/test_cases/run_task5_simd_suite.sh \
  --repo-dir . \
  --build-dir build \
  --label baseline \
  --out-dir homework_docs/test_cases/task5_suite_runs/baseline_$(date +%Y%m%d_%H%M%S)
```

如果你已经知道可执行文件路径，也可以不用 `--build-dir`，改为：

```bash
bash homework_docs/test_cases/run_task5_simd_suite.sh \
  --repo-dir . \
  --abacus-bin ./build/abacus_basic_para \
  --label baseline \
  --out-dir homework_docs/test_cases/task5_suite_runs/baseline_$(date +%Y%m%d_%H%M%S)
```

## 在 SIMD 优化分支上运行

1. 切换到 SIMD 分支。

```bash
git checkout feat/SIMD
```

2. 用和 baseline 完全一致的 CMake 参数重新编译。推荐直接复用 baseline 分支已经验证通过的同一条命令：

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DENABLE_MPI=ON \
  -DUSE_OPENMP=ON \
  -DENABLE_NATIVE_OPTIMIZATION=OFF
cmake --build build -j4
```

3. 运行同一份 suite 脚本，唯一关键区别是 `--label simd`，以及输出目录改成另一处。

```bash
bash homework_docs/test_cases/run_task5_simd_suite.sh \
  --repo-dir . \
  --build-dir build \
  --label simd \
  --out-dir homework_docs/test_cases/task5_suite_runs/simd_$(date +%Y%m%d_%H%M%S)
```

## 输出目录内容

每次 suite 运行后的总输出目录中至少包含：

- `suite_manifest.csv`：每个子测试的索引，记录 case、`nproc`、`threads`、label、输出目录、exit code、开始/结束时间。
- `suite_run_log.txt`：suite 级别的运行日志。
- 每个子测试目录，如 `gaas_medium_np2_omp4/`
- `suite_raw_results.csv`：suite 级别原始结果汇总。
- `suite_summary.csv`：suite 级别汇总表。
- `suite_summary.md`：可直接复制到报告中的整套测试表格。

每个子测试目录中至少包含：

- `logs/*.stdout`：每次运行的完整标准输出日志。
- `logs/*.stderr`：每次运行的完整标准错误日志。
- `environment_info.txt`：主机、日期、git commit、分支名、编译器、MPI 版本、OpenMP 环境变量等。
- `run_manifest.csv`：每次运行的原始索引。
- `raw_results.csv`：按单次运行展开的结果表。
- `summary_results.csv`：按同一组配置汇总后的结果表。
- `summary_table.md`：可复制到报告中的 Markdown 表格。

默认行为是先 warmup 1 次，再做 `--repeat` 次正式运行。汇总只统计 `repeat` 运行，不把 warmup 算入最终性能数据。

## 合并 baseline 和 simd 结果

假设你已经得到两个 suite 输出目录：

- `homework_docs/test_cases/task5_suite_runs/baseline_20260530_210000`
- `homework_docs/test_cases/task5_suite_runs/simd_20260530_223000`

可以把它们一起交给 suite 汇总脚本：

```bash
python3 homework_docs/test_cases/collect_task5_simd_suite_results.py \
  --suite-dir homework_docs/test_cases/task5_suite_runs/baseline_20260530_210000 \
  --suite-dir homework_docs/test_cases/task5_suite_runs/simd_20260530_223000 \
  --out-dir homework_docs/test_cases/task5_suite_runs/combined_baseline_vs_simd
```

这样会生成：

- `combined_baseline_vs_simd/suite_raw_results.csv`
- `combined_baseline_vs_simd/suite_summary.csv`
- `combined_baseline_vs_simd/suite_summary.md`

## 如何生成报告中的性能表

1. 从 `suite_summary.csv` 中读取每个 `case x nproc x threads` 组合的 `median_time_s`、`mean_time_s`、`min_time_s`、`max_time_s`。
2. 若同时存在 `baseline` 和 `simd`，查看 `speedup_vs_baseline` 列。
3. 直接复制 `suite_summary.md` 到实验报告，再补充“测试方法说明”即可。

推荐在报告中明确写出：

- 测试机器名称与 CPU 信息；
- 使用的分支名和 commit；
- 编译器与编译参数；
- `nproc`、`OMP_NUM_THREADS`、warmup 次数、repeat 次数；
- 使用中位数作为主对比指标；
- `speedup = baseline_median / simd_median`。

## 注意事项

- WSL 上之前得到的旧 baseline 数据不应作为最终性能结论。
- 最终性能结论应来自同一台实体机、同一编译器、同一 MPI/OpenMP 环境、同一输入样例和同一套脚本。
- 脚本不会自动切换 git 分支，分支切换请手动完成。
- 脚本不会伪造缺失 timer；缺失时会在 CSV/Markdown 中写 `NA`，并在 `notes` 中提示。
- 为了保证可比性，baseline 和 simd 必须在同一台实体机上、同样的配置下分别执行。
