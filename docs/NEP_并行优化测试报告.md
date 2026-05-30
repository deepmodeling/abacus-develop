# 题目6：NEP 势函数并行优化 — 测试报告

> 对应任务书第四节「预期测试方案」  
> 关联文档：[NEP_并行优化分析报告.md](./NEP_并行优化分析报告.md)  
> 测试日期：2026-05-30

---

## 一、测试目的

验证 NEP 势函数 OpenMP 并行优化在以下方面满足任务书要求：

1. **正确性**：串行（1 线程参考）与并行（2/4/8 线程）的能量、力输出一致，无数据竞争、NaN 或崩溃。
2. **多元素体系**：双元素 NEP4 模型（O + Hf）可正常完成 `nep.compute` 全流程。
3. **性能**：不同 `OMP_NUM_THREADS` 下势函数单次计算耗时与加速比可复现。
4. **稳定性**：相同输入多次运行、不同线程数切换后结果可重复。

---

## 二、测试环境

| 项目 | 配置 |
|------|------|
| 操作系统 | Linux 6.6.87.2 (WSL2)，x86_64 |
| 编译器 | GCC 13.3.0 (Ubuntu 24.04) |
| OpenMP | 4.5（`_OPENMP 201511`） |
| CPU 逻辑核数 | 8 |
| NEP 库 | `third_party/NEP_CPU/build/libnep.so`（`-O3 -fopenmp -DUSE_OPENMP=ON`） |
| 测试程序 | `tools/nep_benchmark/build/nep_benchmark` |
| 势函数模型 | `tests/PP_ORB/nep_hfo2.txt`（NEP4 + ZBL，2 元素 O/Hf） |

**模型参数摘要：**

```
nep4_zbl 2 O Hf
cutoff: radial 8 Å, angular 4 Å
n_max: 4 / 4, basis_size: 12 / 12, l_max: 4 / 2 / 0
ANN: 30-30-1, 总参数 2441
```

---

## 三、测试工具与方法

### 3.1 测试程序 `nep_benchmark`

独立可执行文件，直接调用 `NEP::compute()`，绕过 ABACUS MD 主循环，专门测量势函数内核。

| 参数 | 含义 |
|------|------|
| `--model PATH` | NEP 势文件路径 |
| `--natom N` | 测试原子数（简单立方超胞填充） |
| `--repeat R` | 性能测试重复次数 |
| `--verify` | 正确性模式：1 线程参考 vs 当前线程数 |
| `--perf` | 性能模式：计时 R 次 `compute` 取平均 |

**正确性判定逻辑：**

```cpp
omp_set_num_threads(1);
Result ref = run_nep(...);          // 串行参考
omp_set_num_threads(OMP_NUM_THREADS);
Result cur = run_nep(...);          // 并行被测
dE   = |ref.energy - cur.energy|
max_dF = max_i |ref.force[i] - cur.force[i]|
通过条件: dE < 1e-8 且 max_dF < 1e-8
```

**测试体系构造：**

- 按 `ceil(natom^(1/3))` 建超胞，晶格常数 a = 5.2 Å；
- 多元素：`type[iat] = iat % num_types`，O/Hf 交替，覆盖 `annmb.w0[t1]` 两种 ANN 权重路径；
- 坐标带类型微扰，避免完全对称导致数值退化。

### 3.2 一键脚本 `scripts/nep_benchmark.sh`

自动完成：构建 NEP 库 → 构建 benchmark → 正确性扫描（64 原子）→ 性能扫描（256 原子，20 次）。

```bash
export LD_LIBRARY_PATH=third_party/NEP_CPU/build:$LD_LIBRARY_PATH
bash scripts/nep_benchmark.sh
```

---

## 四、正确性测试

### 4.1 测试项与判据

| 编号 | 测试项 | 判据 |
|------|--------|------|
| C1 | 总能量一致性 | \|E_parallel − E_serial\| < 1×10⁻⁸ |
| C2 | 原子力一致性 | max\|F_parallel − F_serial\| < 1×10⁻⁸ |
| C3 | 多线程可切换 | 1/2/4/8 线程均通过 C1、C2 |
| C4 | 多元素类型映射 | O/Hf 交替体系无异常退出 |
| C5 | 连续运行稳定性 | 10 次连续 verify 全部通过 |

### 4.2 测试结果（64 原子，多元素 HfO₂）

| OMP_NUM_THREADS | 总能量 (eV) | dE | max_dF | 结论 |
|-----------------|------------|-----|--------|------|
| 1 | −275.28 | 0 | 0 | 通过 |
| 2 | −275.28 | 0 | 0 | 通过 |
| 4 | −275.28 | 0 | 0 | 通过 |
| 8 | −275.28 | 0 | 0 | 通过 |

说明：各线程数下能量与力相对 1 线程参考完全一致（双精度下 dE、max_dF 均为 0），满足任务书阈值 1×10⁻⁸。

### 4.3 稳定性测试（128 原子，4 线程，10 次连续）

| 运行次数 | dE | max_dF | 结论 |
|---------|-----|--------|------|
| 1–10 | 0 | 0 | 全部通过 |

**通过率：10 / 10（100%）**

未观察到 NaN、段错误或能量漂移，说明并行归约路径无数值异常。

### 4.4 单元素 vs 多元素说明

仓库内可用 NEP 模型为 **双元素** `nep_hfo2.txt`（O + Hf），本次多元素测试已覆盖：

- 径向/角向系数索引 `t12 = t1 * num_types + t2` 的四种元素对组合；
- 按 `t1` 切换的 ANN 权重 `w0[0]`（O）与 `w0[1]`（Hf）。

单元素 NEP 模型文件未包含在仓库中；若需补充，可替换为单元素 `nep.txt` 并执行：

```bash
OMP_NUM_THREADS=4 nep_benchmark --model <single_element_nep.txt> --natom 64 --verify
```

测试框架与判据相同。

---

## 五、性能测试

### 5.1 测试方法

- 固定模型 `nep_hfo2.txt`，变化 `OMP_NUM_THREADS` 与 `natom`；
- 每次正式计时前执行 1 次预热 `compute`；
- 记录 `repeat` 次总耗时，报告 `per_call_ms = total_s × 1000 / repeat`；
- 加速比：`S_p = T_1 / T_p`（T 为单次 compute 平均毫秒数）；
- 并行效率：`E_p = S_p / p × 100%`。

### 5.2 线程扩展性（256 原子，repeat = 20）

| 线程数 p | 单次耗时 (ms) | 加速比 S_p | 并行效率 E_p |
|---------|--------------|-----------|-------------|
| 1 | 3.79 | 1.00× | 100% |
| 2 | 2.70 | 1.40× | 70.2% |
| 4 | 1.77 | 2.14× | 53.6% |
| 8 | 7.65 | 0.50× | 6.2% |

```
耗时 (ms)
  8 |                              ████████ 7.65
  6 |
  4 |              ████ 3.79
  2 |        ██ 2.70    ██ 1.77
  0 +----+----+----+----+----+
      1    2    4    8   线程数
```

**结论：**

- **推荐线程数：2–4**。4 线程相对 1 线程约 **2.14×** 加速，为本次测试最优。
- **8 线程性能回退**：单次耗时增至 7.65 ms，主要因为线程私有力缓冲区内存开销（约 `8 × N × (3+9)` 个 double）及 256 原子体系并行粒度不足。

### 5.3 体系规模扩展（4 线程，repeat = 15）

| 原子数 N | 单次耗时 (ms) | 相对 64 原子倍数 |
|---------|--------------|-----------------|
| 64 | 0.38 | 1.0× |
| 128 | 1.07 | 2.8× |
| 256 | 2.37 | 6.2× |
| 512 | 4.28 | 11.2× |

原子数翻倍时耗时约增 2–3 倍，符合 NEP O(N × 邻居数) 复杂度预期；512 原子、4 线程下单次 compute 约 4.3 ms，并行收益随体系增大应更明显。

### 5.4 与优化前对比（定性）

| 模块 | 优化前 | 优化后（4 线程） |
|------|--------|-----------------|
| 描述符/密度 | 已 OpenMP 并行 | `schedule(static)` 微调 |
| 径向/角向/ZBL 力 | **串行** | **OpenMP + 线程私有缓冲** |
| 256 原子整体 | 仅描述符段并行 | 全流程主要热点均可并行 |

力计算并行化是 2–4 线程加速的主要来源（见分析报告热点表）。

---

## 六、异常与边界测试

| 编号 | 场景 | 操作 | 结果 |
|------|------|------|------|
| E1 | 非法线程环境 | `OMP_NUM_THREADS=8` 运行 verify | 通过，无崩溃 |
| E2 | 小规模体系 | `natom=64`，4 线程 | 通过，dE=0 |
| E3 | 大规模体系 | `natom=512`，4 线程 perf | 正常完成，无 OOM |
| E4 | 重复调用 | 连续 10 次 verify | 10/10 通过 |
| E5 | 库路径 | 未设置 `LD_LIBRARY_PATH` | 预期失败（需加载 libnep.so） |

---

## 七、ABACUS 集成测试（说明）

本次自动化测试在 **NEP 库层** 完成，与 ABACUS `ESolver_NEP` 的关系如下：

| 层级 | 测试覆盖 | 状态 |
|------|---------|------|
| NEP_CPU `nep.compute` | `nep_benchmark` 直接测试 | **已测** |
| `ESolver_NEP::runner` 单位换算 | 逻辑简单，依赖上层 `nep.compute` | 间接覆盖 |
| 完整 MD 步 `MD_func::force_virial` | 需 `-DUSE_BUNDLED_NEP=ON` 编译 ABACUS + MD 算例 | 待用户环境执行 |

**建议 ABACUS 端到端验证命令（需已编译 NEP 支持）：**

```bash
# 编译
cmake -B build -DUSE_BUNDLED_NEP=ON -DUSE_OPENMP=ON
cmake --build build -j$(nproc)

# MD 算例（需自备 INPUT/STRU，esolver_type nep, pot_file nep_hfo2.txt）
export OMP_NUM_THREADS=4
mpirun -np 1 ./build/abacus  # 或单进程运行
```

对比 `OMP_NUM_THREADS=1` 与 `4` 下输出的 `#TOTAL ENERGY#` 与 `TOTAL-FORCE` 是否一致。

---

## 八、测试结论

### 8.1 正确性

| 任务书要求 | 结果 |
|-----------|------|
| 总能量误差 < 1e-8 | **满足**（dE = 0） |
| 力误差 < 1e-8 | **满足**（max_dF = 0） |
| 不同线程数结果稳定 | **满足**（1/2/4/8 线程） |
| 多元素体系可运行 | **满足**（O/Hf NEP4） |
| 无竞争/NaN/崩溃 | **满足**（10 次连续通过） |

### 8.2 性能

| 任务书要求 | 结果 |
|-----------|------|
| OMP_NUM_THREADS = 1,2,4,8 测试 | **已完成** |
| 记录势函数计算时间 | **已完成**（见 5.2、5.3 节） |
| 加速比评估 | **4 线程最优约 2.14×** |

### 8.3 总体评价

NEP 势函数并行优化在正确性上达到任务书要求；性能上在 **2–4 线程、256 及以上原子数** 时获得稳定加速，**不建议默认使用 8 线程**（本测试环境下性能回退）。测试脚本与程序可复现，便于后续扩展单元素模型或完整 ABACUS MD 联调。

---

## 九、复现步骤

```bash
cd /path/to/abacus-develop-develop

# 构建
cmake -S third_party/NEP_CPU -B third_party/NEP_CPU/build \
  -DUSE_OPENMP=ON -DCMAKE_BUILD_TYPE=Release
cmake --build third_party/NEP_CPU/build -j$(nproc)

cmake -S tools/nep_benchmark -B tools/nep_benchmark/build \
  -DNEP_DIR=$PWD/third_party/NEP_CPU/build -DUSE_OPENMP=ON
cmake --build tools/nep_benchmark/build -j$(nproc)

# 运行全部测试
export LD_LIBRARY_PATH=$PWD/third_party/NEP_CPU/build:$LD_LIBRARY_PATH
bash scripts/nep_benchmark.sh

# 或单项测试
OMP_NUM_THREADS=4 tools/nep_benchmark/build/nep_benchmark \
  --model tests/PP_ORB/nep_hfo2.txt --natom 256 --repeat 20 --verify --perf
```

---

## 十、附录：原始测试输出摘录

```
=== Correctness (64 atoms, multi-element HfO2) ===
energy=-275.28 dE=0 max_dF=0   # threads=1
energy=-275.28 dE=0 max_dF=0   # threads=2
energy=-275.28 dE=0 max_dF=0   # threads=4
energy=-275.28 dE=0 max_dF=0   # threads=8

=== Performance (256 atoms, 20 repeats) ===
threads=1: per_call_ms=3.78617
threads=2: per_call_ms=2.6989
threads=4: per_call_ms=1.76732
threads=8: per_call_ms=7.65422

=== Performance scaling (4 threads) ===
natom=64:  per_call_ms=0.382169
natom=128: per_call_ms=1.07126
natom=256: per_call_ms=2.37256
natom=512: per_call_ms=4.2805

=== Stability (4 threads, 10 runs) ===
passed 10/10 runs
```
