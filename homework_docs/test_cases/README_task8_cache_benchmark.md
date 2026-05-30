# Task 8 Cache Benchmark Workflow

这套脚本用于“题目 8：平面波预计算与缓存优化”的正式性能对比。设计目标是：

- 同一套脚本可用于 baseline 原始分支和题目 8 cache 优化分支。
- 不依赖题目 8 优化分支新增的源码函数、专有接口或专有测试程序。
- 不自动切换 git 分支，分支切换由你手动控制。
- 如果日志中没有某个细粒度 timer，汇总结果写为 `NA`，并在 `notes` 中说明 missing。

## 文件说明

- `run_task8_cache_suite.sh`
  - 一键运行固定的题目 8 性能测试矩阵。
  - 作为 wrapper 调用 `run_task8_cache_benchmark.sh`。
  - 为每个 `case x nproc x threads` 组合创建单独子目录，并生成 `suite_manifest.csv`、`suite_run_log.txt`。
- `run_task8_cache_benchmark.sh`
  - 运行指定 case 的 warmup 和正式 benchmark。
  - 保存每次运行的完整 `stdout/stderr`、`OUT.ABACUS` 归档副本、环境信息和 git 信息。
  - 自动调用汇总脚本生成单次结果目录下的 CSV/Markdown。
- `collect_task8_cache_results.py`
  - 从一个或多个 benchmark 输出目录收集日志。
  - 生成 `raw_results.csv`、`summary_results.csv` 和 `summary_table.md`。
  - 当同时输入 baseline 和 cache 两组目录时，若 `case_name + nproc + threads` 匹配，会自动计算 `speedup_vs_baseline = baseline_median / current_median`。
- `collect_task8_cache_suite_results.py`
  - 扫描一个或多个 suite 输出目录下的全部子测试目录。
  - 生成 `suite_raw_results.csv`、`suite_summary.csv` 和 `suite_summary.md`。
  - 支持把 baseline suite 和 cache suite 合并成同一张对比表。

## 默认 suite 配置

`run_task8_cache_suite.sh` 开头集中维护固定测试矩阵：

- Cases
  - `homework_docs/test_cases/gaas_small`
  - `homework_docs/test_cases/gaas_medium`
  - `gaas_large` 作为可选补充样例，默认先注释掉
- MPI 进程数
  - `1 2 4`
- OpenMP 线程数
  - `1 2 4`
- Warmup
  - `1`
- Repeat
  - `3`
- Timeout
  - `1800`

选择理由：

- `gaas_small` 和 `gaas_medium` 都会触发 PW 初始化、平面波相关数据构建、FFT 与 SCF 迭代，且运行时间相对可控。
- 它们比 `gaas_tiny` 更容易体现平面波重复计算减少与缓存复用的收益。
- `gaas_large` 更重，适合作为补充验证，但在部分机器上可能显著拉长整套测试时间，因此默认注释。

## 在 baseline 原始分支上运行

1. 手动切换到 baseline 分支或 commit。

```bash
git checkout <baseline-branch-or-commit>
```

2. 用和优化分支完全一致的 CMake 参数重新编译。推荐示例：

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

若环境需要，也可以补充 `-DUSE_ELPA=OFF`。关键点不是具体参数长什么样，而是 baseline 和题目 8 优化分支必须使用完全相同的 CMake 参数。

3. 运行 suite 脚本，`--label` 设为 `baseline`：

```bash
bash homework_docs/test_cases/run_task8_cache_suite.sh \
  --repo-dir . \
  --build-dir build \
  --label baseline \
  --out-dir homework_docs/test_cases/task8_cache_suite_runs/baseline_$(date +%Y%m%d_%H%M%S)
```

如果你更想直接传二进制路径，也可以用：

```bash
bash homework_docs/test_cases/run_task8_cache_suite.sh \
  --repo-dir . \
  --abacus-bin ./build/abacus_basic_para \
  --label baseline \
  --out-dir homework_docs/test_cases/task8_cache_suite_runs/baseline_$(date +%Y%m%d_%H%M%S)
```

## 在题目 8 优化分支上运行

1. 手动切换到题目 8 cache 优化分支。

```bash
git checkout <task8-cache-branch>
```

2. 使用与 baseline 完全一致的 CMake 参数重新编译。

3. 运行同一份 suite 脚本，`--label` 改为 `cache` 或 `optimized`：

```bash
bash homework_docs/test_cases/run_task8_cache_suite.sh \
  --repo-dir . \
  --build-dir build \
  --label cache \
  --out-dir homework_docs/test_cases/task8_cache_suite_runs/cache_$(date +%Y%m%d_%H%M%S)
```

## 输出目录内容

每次 suite 输出目录中至少包含：

- `suite_manifest.csv`
- `suite_run_log.txt`
- 每个子测试目录，例如 `gaas_medium_np2_omp4/`
- `suite_raw_results.csv`
- `suite_summary.csv`
- `suite_summary.md`

每个子测试目录中至少包含：

- `logs/*.stdout`
- `logs/*.stderr`
- `logs/*.OUT.ABACUS/` 或同名归档副本
- `environment_info.txt`
- `run_manifest.csv`
- `raw_results.csv`
- `summary_results.csv`
- `summary_table.md`

汇总只统计 `repeat` 运行，不把 warmup 计入最终性能数据。

## 合并 baseline 和优化结果

假设你已经得到两个 suite 输出目录：

- `homework_docs/test_cases/task8_cache_suite_runs/baseline_20260530_210000`
- `homework_docs/test_cases/task8_cache_suite_runs/cache_20260530_223000`

可用下面的命令生成最终对比表：

```bash
python3 homework_docs/test_cases/collect_task8_cache_suite_results.py \
  --suite-dir homework_docs/test_cases/task8_cache_suite_runs/baseline_20260530_210000 \
  --suite-dir homework_docs/test_cases/task8_cache_suite_runs/cache_20260530_223000 \
  --out-dir homework_docs/test_cases/task8_cache_suite_runs/combined_baseline_vs_cache
```

匹配键为：

- `case_name`
- `nproc`
- `threads`

加速比计算方式为：

```text
speedup_vs_baseline = baseline_median_time_s / optimized_median_time_s
```

这里使用的是 median time，而不是任意单次运行时间。

## 细粒度 timer 的解释

汇总表会尝试解析这些字段：

- `pw_init_s_median`
- `plane_wave_generation_s_median`
- `gcar_or_gk2_s_median`
- `real2recip_s_median`
- `recip2real_s_median`

这些列只在日志中存在对应通用 timer 时才会填数值。当前脚本使用的是保守候选匹配，例如：

- `setuptransform` 近似对应 `pw_init`
- `collect_local_pw` / `collect_uniqgg` / `distributeg` 近似对应 `plane_wave_generation`
- `setupIndGk` 近似对应 `gcar_or_gk2`

如果你的 ABACUS 日志没有这些 timer，结果会写成 `NA`，`notes` 中会记录 missing。这种情况下，报告中应明确说明：

- 当前版本日志未暴露更细粒度的平面波初始化或 `gk2/gcar` 计时；
- 因此结论主要依据总 wall-clock 时间与可解析的 FFT 往返计时；
- `NA` 表示“日志不可得”，不表示该阶段耗时为零。

## 注意事项

- baseline 和题目 8 优化分支必须使用相同的 CMake 参数、相同 case、相同 MPI/OpenMP 配置。
- WSL 上的性能数据不应作为最终结论，最终报告应以同一台实体机上的结果为准。
- 脚本不会自动切换 git 分支。
- 脚本不会伪造缺失 timer。
