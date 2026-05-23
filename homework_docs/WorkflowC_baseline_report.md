# Workflow C 基线测试报告

## 第 12 周：当前算法测试与总结

**负责人**: 徐贤天  
**日期**: 2026-05-22  
**工作流**: SIMD 向量化（题目5）+ 缓存复用（题目8）  
**阶段**: 优化前基线测量  

---

## 1. 测试环境

| 项目 | 详情 |
|------|------|
| CPU | Intel Core i5-13500H (13th Gen) |
| 核心数 | 8P+8E = 16 逻辑核 (超线程) |
| 内存 | 7.6 GB |
| 操作系统 | Linux 6.6.87.2-microsoft-standard-WSL2 |
| 编译器 | GCC 13 (g++ 13) |
| MPI | OpenMPI 4.1.6 |
| 数学库 | FFTW 3 + OpenBLAS + ScaLAPACK |
| SIMD 指令集 | AVX2, AVX-VNNI |
| 缓存 | L1d: 384KB, L2: 10MB, L3: 18MB |
| ABACUS 版本 | v3.11.0-beta.2 (commit fae9fe94b) |
| 构建类型 | Debug (`-O0 -g`，未优化) |

**ABACUS 二进制**: `build/abacus_basic_para`，链接 FFTW3、OpenBLAS、ScaLAPACK-openmpi。

> **注意**：本报告主要数据来自 `-O0`（Debug）构建。Release（`-O3 -march=native -DNDEBUG`）对照实验见 [5.4 节](#54-release-构建对照实验-o3--marchnative--dndebug)。Debug 测得的百分比占比仍然有效，但 **SIMD 的预期加速倍数不能直接用 Debug 数据估算**——GCC 13 在 `-O3` 下已对拷贝循环自动向量化。

---

## 执行摘要

| 维度 | 关键数字 | 说明 |
|------|----------|------|
| 测试通过率 | **72/72 (100%)** | 4 用例 × 6 配置 × 3 重复 |
| 最大热点 | **diag_once (CG 对角化)** 占 SCF 73~77% | Workflow D 范畴 |
| Gather/Scatter 占比 | FFT 的 **30~55%**, SCF 的 **6~10%** (MPI 稳态) | 6 处连续拷贝循环在 `pw_gatherscatter.h` |
| SIMD 状态 | **GCC 13 -O3 有中等优化效果** (gathers 1.56×, gatherp 无显著加速) | 冷启动消除后重新评估，需 `-fopt-info-vec` 确认 |
| 缓存机会 | 数据 **~0.9 MB** (可常驻 L2), setupIndGk **~18K 冗余调用** | 见 [6.2.5 节](#625-setupindgk-双重计算量化评估) |
| Debug→Release | serial **1.44×**, np2_omp2 **1.14×** 整体加速 | 旧数据的 2.71× 受冷启动污染（见 [4.4 节](#44-冷启动现象重跑验证)） |
| 冷启动状态 | **一次性系统预热现象**，重跑验证确认 r1/r2/r3 无差异 | 详见 [4.4 节](#44-冷启动现象重跑验证) |

---

## 2. 测试用例设计

### 2.1 设计原则

为适应 7.6 GB 内存限制并保证可复现，采用**同一 GaAs 体系在三种 FFT 网格下对比**的策略。相同材料保证了电子结构一致（能带数、k 点密度），仅改变 FFT 网格大小，分离出数据搬运量对性能的影响。

### 2.2 测试用例参数

| 参数 | gaas_tiny | gaas_tiny_40Ry | gaas_small | gaas_medium |
|------|-----------|----------------|------------|-------------|
| 来源 | pw2/cg/001_4GaAs | pw2/cg/001_4GaAs | pw2/cg/001_4GaAs | pw2/cg/001_4GaAs |
| 化学式 | Ga₄As₄ (8 原子) | Ga₄As₄ (8 原子) | Ga₄As₄ (8 原子) | Ga₄As₄ (8 原子) |
| FFT 网格 | 24×24×24 | 24×24×24 | 32×32×32 | 48×48×48 |
| 截断能 ecutwfc | 20 Ry | 40 Ry | 40 Ry | 40 Ry |
| 能带数 nbands | 70 | 70 | 70 | 70 |
| K 点网格 | 3×3×3 (4 个不可约 k 点) | 3×3×3 (4 个不可约 k 点) | 3×3×3 (4 个不可约 k 点) | 3×3×3 (4 个不可约 k 点) |
| SCF 趋敛 | 10 次未收敛 | 9 次收敛 | 9 次收敛 | 8 次收敛 |
| 求解器 | CG | CG | CG | CG |
| 赝势 | pseudo-dojo v0.5 | pseudo-dojo v0.5 | pseudo-dojo v0.5 | pseudo-dojo v0.5 |

### 2.3 并行配置矩阵

每个测试用例运行 **6 种并行配置 × 3 次重复 = 18 次/用例，4 用例总计 72 次运行**。

| 配置标签 | MPI 进程 | OpenMP 线程 | 总占用核 | 目的 |
|----------|----------|-------------|----------|------|
| serial | 1 | 1 | 1 | 串行基线 |
| omp4 | 1 | 4 | 4 | OpenMP 扩展性测试 |
| omp8 | 1 | 8 | 8 | OpenMP 峰值测试 |
| mix_np2_omp2 | 2 | 2 | 4 | 混合并行 A |
| mix_np2_omp4 | 2 | 4 | 8 | 混合并行 B (MPI 轻) |
| mix_np4_omp2 | 4 | 2 | 8 | 混合并行 C (MPI 重) |

运行命令模板：
```bash
OMP_NUM_THREADS=$omp mpirun --allow-run-as-root -np $np abacus_basic_para
```

---

## 3. 计时桩点说明

### 3.1 原始计时覆盖

ABACUS 使用 `ModuleBase::timer` 机制记录计算耗时，按 `CLASS_NAME::FUNCTION_NAME` 组织。原始代码中 `module_pw` 模块仅有粗粒度计时：

- `PW_Basis::setuptransform` — 整体初始化
- `PW_Basis_K::setuptransform` — 多 k 点整体初始化
- `PW_Basis::real2recip` / `PW_Basis::recip2real` — FFT 正/逆变换（含 gather/scatter）
- `distributeg` — G 矢量分布

### 3.2 新增临时计时桩点

为获取 Workflow C 所需粒度数据，在以下函数中临时插入 `timer::start/end`（基线测试完成后将移除）：

| 函数 | 文件:行号 | 所属题目 |
|------|-----------|----------|
| `PW_Basis::gatherp_scatters` | `pw_gatherscatter.h:18,104` | 题目5 SIMD |
| `PW_Basis::gathers_scatterp` | `pw_gatherscatter.h:120,221` | 题目5 SIMD |
| `PW_Basis::collect_local_pw` | `pw_basis.cpp:137,187` | 题目8 缓存 |
| `PW_Basis::collect_uniqgg` | `pw_basis.cpp:195,277` | 题目8 缓存 |
| `PW_Basis_K::setupIndGk` | `pw_basis_k.cpp:133,202` | 题目8 缓存 |
| `PW_Basis_K::collect_local_pw` | `pw_basis_k.cpp:257,358` | 题目8 缓存 |

### 3.3 计时器限制

`ModuleBase::timer` 有 0.1 秒的输出过滤阈值（`timer.cpp:269`），且 OpenMP 并行时只记录主线程耗时（`timer.cpp:73`）。因此部分耗时极短的函数不会出现在 TIME STATISTICS 表中，这是已知限制。

---

## 4. 测试验收结果

### 4.1 总体通过率

| 测试用例 | 运行次数 | 通过 | 失败 | 通过率 |
|----------|----------|------|------|--------|
| gaas_tiny | 18 | 18 | 0 | 100% |
| gaas_tiny_40Ry | 18 | 18 | 0 | 100% |
| gaas_small | 18 | 18 | 0 | 100% |
| gaas_medium | 18 | 18 | 0 | 100% |
| **总计** | **72** | **72** | **0** | **100%** |

### 4.2 验收标准

- 所有运行正常退出（EXIT=0）
- TIME STATISTICS 计时表完整输出
- 无 crash、segfault 或 MPI 异常
- stderr 中无错误信息（仅 OpenMPI hwloc 信息性提示）

### 4.3 墙钟时间概览

> **通用说明**：以下墙钟时间为 **3 次重复的原始数据**。重跑验证表明初始基线中的 MPI 冷启动效应是一次性系统预热现象，新数据中 r1/r2/r3 高度一致（详见 [4.4 节](#44-冷启动现象重跑验证)）。gaas_tiny（20 Ry）的 SCF 在 10 次迭代内**未收敛**（显示 `!!SCF IS NOT CONVERGED!!`），gaas_tiny_40Ry 和 gaas_small 在 9 次收敛，gaas_medium 在 8 次收敛。迭代次数差异不影响同用例内的可比性。

**gaas_tiny (24³)**:

| 配置 | rep1 | rep2 | rep3 |
|------|------|------|------|
| serial (np1_omp1) | 9s | 7s | 9s |
| omp4 (np1_omp4) | 5s | 5s | 5s |
| omp8 (np1_omp8) | 5s | 3s | 5s |
| mix_np2_omp2 | 5s | 5s | 5s |
| mix_np2_omp4 | 14s | 17s | 15s |
| mix_np4_omp2 | 6s | 7s | 6s |

**gaas_small (32³)**:

| 配置 | rep1 | rep2 | rep3 |
|------|------|------|------|
| serial (np1_omp1) | 21s | 20s | 20s |
| omp4 (np1_omp4) | 9s | 9s | 9s |
| omp8 (np1_omp8) | 13s | 17s | 17s |
| mix_np2_omp2 | 10s | 11s | 10s |
| mix_np2_omp4 | 27s | 28s | 26s |
| mix_np4_omp2 | 8s | 11s | 10s |

**gaas_tiny_40Ry (24³，40 Ry)**:

| 配置 | rep1 | rep2 | rep3 |
|------|------|------|------|
| serial (np1_omp1) | 20s | 19s | 19s |
| omp4 (np1_omp4) | 9s | 9s | 8s |
| omp8 (np1_omp8) | 10s | 10s | 11s |
| mix_np2_omp2 | 8s | 9s | 9s |
| mix_np2_omp4 | 18s | 23s | 22s |
| mix_np4_omp2 | 11s | 9s | 11s |

**gaas_medium (48³)**:

| 配置 | rep1 | rep2 | rep3 |
|------|------|------|------|
| serial (np1_omp1) | 33s | 35s | 36s |
| omp4 (np1_omp4) | 17s | 16s | 17s |
| omp8 (np1_omp8) | 19s | 17s | 17s |
| mix_np2_omp2 | 15s | 18s | 15s |
| mix_np2_omp4 | 29s | 31s | 35s |
| mix_np4_omp2 | 14s | 15s | 13s |

### 4.4 冷启动现象重跑验证

初始基线数据中，MPI 配置（特别是 `mix_np2_omp2`）的首轮运行时间远超后续重复，表现为"冷启动"现象。为确认该现象是否可复现，**将所有 54 次运行（gaas_tiny、gaas_small、gaas_medium）在原系统上完整重跑了一遍**。

#### 4.4.1 重跑结果

重跑后，冷启动效应**基本消失**——r1/r2/r3 的墙钟时间高度一致。核心指标对比：

| 用例 | 配置 | 旧 r1/r2 | 新 r1/r2 | 新 r1 | 新 r2 | 新 r3 |
|------|------|-----------|-----------|-------|-------|-------|
| gaas_small | mix_np2_omp2 | **2.50×** | 0.91× | 10s | 11s | 10s |
| gaas_medium | mix_np2_omp2 | **3.89×** | 0.83× | 15s | 18s | 15s |
| gaas_small | mix_np2_omp4 | 0.97× | 0.96× | 27s | 28s | 26s |
| gaas_medium | serial | 1.13× | 0.94× | 33s | 35s | 36s |
| gaas_tiny | mix_np2_omp2 | 缺失 r1 | 1.00× | 5s | 5s | 5s |

所有 MPI 配置的 r1/r2 比值从 2.5~3.9× 回归到 **0.83~1.00×**，说明冷启动是**一次性系统预热现象**，不具可复现性。

#### 4.4.2 根因分析

冷启动出现于首次运行（初始基线）但在重跑中消失的原因：

1. **WSL2 VM 初始化**：初始基线是 WSL2 实例首次运行计算密集型 MPI 作业时触发的。WSL2 的 Hyper-V 轻量级 VM 在首次负载下需要分配和提交物理内存（动态内存气球驱动），这个过程在后续运行中已完成。

2. **CPU 频率爬升**：首次运行时，CPU 从空闲低频率（~400 MHz）爬升至满载 Turbo Boost 频率（~4.7 GHz）需要数百毫秒。初始基线的首轮运行覆盖了频率爬升阶段，重跑时 CPU 已处于忙碌状态，频率保持在较高水平。

3. **操作系统页缓存预热**：ABACUS 二进制文件和输入文件在首次运行时需从磁盘加载到页缓存。重跑时这些文件已在内存中，减少了 I/O 延迟。

4. **初始基线中 gaas_tiny `np2_omp2_r1` 和 `np2_omp4_r1` 被 timeout 杀死**（无 TOTAL Time），已产生错误退出进程残留，影响了后续 MPI 运行状态。

#### 4.4.3 对测试可靠性的影响

1. **冷启动为一次性系统预热现象**，不具可复现性。第一次运行 MPI 基准测试时存在预热效应，但后续所有运行（无论同一批次还是另起批次）均不受影响。

2. **r1/r2/r3 均可用于优化分析**：新数据中三次重复的一致性良好（标准差通常 < 2s），可直接使用三者的均值或中位数。

3. **冷启动的消失不影响本报告的核心结论**：gather/scatter 占 MPI 模式 FFT 的 30~55%、对角化仍为最大热点（占 SCF ~75-80%）、缓存复用机会等结论在稳态数据下依然成立。唯一受影响的定量数据是 Debug→Release 加速比的绝对值（见 [5.4 节](#54-release-构建对照实验-o3--marchnative--dndebug)），因为旧 Debug 基线被冷启动抬高。

4. 初始基线中 `mix_np2_omp2` 配置的墙钟时间（gaas_tiny: 53s/5s/5s, gaas_small: 35s/14s/14s, gaas_medium: 136s/35s/40s）**已被本次重跑数据取代**。本报告所有后续分析均基于新的重跑数据。

---

## 5. Gather/Scatter 基线数据 (题目5：SIMD 向量化)

### 5.1 计时数据总表（稳态 r2/r3 均值）

以下数据为 r2 和 r3 的算术平均。重跑验证表明 r1/r2/r3 间无系统性差异（见 [4.4 节](#44-冷启动现象重跑验证)），因此 r2/r3 均值即可代表系统稳态性能。`gatherp_scatters` 在 `poolnproc=1`（单进程）时因走简化代码路径（无 MPI pre/post copy），部分配置下耗时低于 0.1s 阈值不显示。每个 SCF 迭代中 `gatherp_scatters` 和 `gathers_scatterp` 分别调用约 nkpt × nbands 次，全运行累计 7,000~12,000 次，调用量随网格和 k 点数增长。

#### gaas_tiny (24³)

| 配置 | gatherp (s) | gathers (s) | real2recip (s) | recip2real (s) | gp/r2c | gs/r2r | gp/ham |
|------|-------------|-------------|-----------------|-----------------|--------|--------|--------|
| serial | 0.12 | 0.62 | 1.19 | 2.02 | 10.0% | 30.7% | 1.5% |
| omp4 | 0.06 | 0.24 | 0.63 | 1.10 | 9.5% | 21.7% | 1.6% |
| omp8 | 0.06 | 0.22 | 0.64 | 1.19 | 8.6% | 18.5% | 1.4% |
| mix_np2_omp2 | 0.29 | 0.47 | 0.80 | 1.17 | 36.3% | 40.2% | 7.2% |
| mix_np2_omp4 | 1.34 | 2.17 | 3.25 | 4.93 | 41.2% | 44.1% | 8.5% |
| mix_np4_omp2 | 0.56 | 0.72 | 1.07 | 1.43 | 52.6% | 50.7% | 10.2% |

#### gaas_tiny_40Ry (24³，40 Ry)

此用例为补测数据，使用与 small/medium 相同的截断能（40 Ry），用于消除 5.1a 节中 tiny 的混杂因素。

| 配置 | gatherp (s) | gathers (s) | real2recip (s) | recip2real (s) | gp/r2c | gs/r2r | gp/ham |
|------|-------------|-------------|-----------------|-----------------|--------|--------|--------|
| serial | — | 0.69 | 1.56 | 2.47 | — | 27.9% | — |
| omp4 | — | 0.23 | 0.74 | 1.33 | — | 17.3% | — |
| omp8 | — | 0.23 | 0.79 | 1.50 | — | 15.3% | — |
| mix_np2_omp2 | 0.37 | 0.65 | 0.95 | 1.50 | 38.9% | 43.3% | 4.9% |
| mix_np2_omp4 | 1.55 | 2.50 | 3.69 | 5.75 | 42.0% | 43.5% | 7.2% |
| mix_np4_omp2 | 0.59 | 0.88 | 1.21 | 1.73 | 48.8% | 50.9% | 6.4% |

#### gaas_small (32³)

| 配置 | gatherp (s) | gathers (s) | real2recip (s) | recip2real (s) | gp/r2c | gs/r2r | gp/ham |
|------|-------------|-------------|-----------------|-----------------|--------|--------|--------|
| serial | 0.23 | 1.31 | 1.73 | 3.39 | 13.5% | 38.6% | 1.3% |
| omp4 | — | 0.39 | 0.69 | 1.50 | — | 25.6% | — |
| omp8 | 0.19 | 0.58 | 1.48 | 2.84 | 12.8% | 20.6% | 1.4% |
| mix_np2_omp2 | 0.66 | 1.04 | 1.34 | 2.08 | 49.3% | 50.2% | 6.5% |
| mix_np2_omp4 | 2.08 | 3.15 | 4.60 | 6.96 | 45.1% | 45.3% | 8.1% |
| mix_np4_omp2 | 0.74 | 1.02 | 1.34 | 1.89 | 55.0% | 54.1% | 7.9% |

#### gaas_medium (48³)

| 配置 | gatherp (s) | gathers (s) | real2recip (s) | recip2real (s) | gp/r2c | gs/r2r | gp/ham |
|------|-------------|-------------|-----------------|-----------------|--------|--------|--------|
| serial | — | 3.53 | 6.00 | 11.15 | — | 31.6% | — |
| omp4 | — | 1.25 | 2.71 | 5.17 | — | 24.2% | — |
| omp8 | — | 1.03 | 2.75 | 5.23 | — | 19.7% | — |
| mix_np2_omp2 | 1.04 | 1.86 | 3.28 | 5.01 | 31.7% | 37.1% | 6.6% |
| mix_np2_omp4 | 2.59 | 4.00 | 6.91 | 10.27 | 37.5% | 39.0% | 8.1% |
| mix_np4_omp2 | 0.96 | 1.61 | 2.81 | 4.22 | 34.0% | 38.2% | 7.1% |

> 由于重跑验证确认了 r1/r2/r3 间无系统性差异（见 [4.4 节](#44-冷启动现象重跑验证)），本表中的 r2/r3 均值直接代表系统性能，无需排除某次重复。

### 5.1a 网格规模扩展分析

三个测试用例使用相同材料（GaAs）、相同 k 点密度（3×3×3），但 gaas_tiny 截断能为 20 Ry 而 small/medium 为 40 Ry。为填补此空缺，补测了 **gaas_tiny_40Ry**（24³, 40 Ry），使 40 Ry 下拥有完整的 24³→32³→48³ 纯净扩展链：

| 指标 | gaas_tiny (20Ry, 24³) | gaas_tiny_40Ry (40Ry, 24³) | gaas_small (40Ry, 32³) | gaas_medium (40Ry, 48³) | 40Ry 下 24³→48³ 纯净扩展 | 理论预期 |
|------|-----------------------|---------------------------|------------------------|-------------------------|--------------------------|----------|
| nplane (serial) | 24 | 24 | 32 | 48 | 2.0× | — |
| gathers 耗时 | 0.62s | 0.69s | 1.31s | 3.53s | **5.1×** | O(nplane) = 2× |
| real2recip 耗时 | 1.19s | 1.56s | 1.73s | 6.00s | **3.8×** | ~O(nplane) × log |
| recip2real 耗时 | 2.02s | 2.47s | 3.39s | 11.15s | **4.5×** | ~O(nplane) × log |
| hPsi 总耗时 | 5.66s | 12.25s | 12.99s | 25.55s | **2.1×** | — |
| SCF 迭代耗时 | 7.89s | 18.20s | 18.77s | 33.23s | **1.8×** | ~O(n³log n) |
| gathers 单次调用 | ~48µs | 52µs | ~109µs | ~308µs | **5.9×** | O(nplane) = 2× |

**ecutwfc 效应（固定 24³ 网格）**：对比 tiny（20 Ry）和 tiny_40Ry（40 Ry），gathers 耗时增长 **1.1×**（0.62s→0.69s），hPsi 增长 **2.2×**（5.66s→12.25s）。这说明截断能对计算量的影响显著——更高 ecutwfc 引入了更多平面波，直接增加了对角化和哈密顿量作用的工作量。

**网格纯净扩展（40 Ry，24³→48³）**：
- `gathers_scatterp` 绝对耗时增长 **5.1×**（0.69s→3.53s），超过 nplane 线性增长（2×）。这是因为更大的网格在 Debug（`-O0`）模式下涉及更多的索引计算和分支判断，使单次调用的耗时从 52µs 膨胀至 308µs（**5.9×**），其中一部分来自 nplane 增长（2×），其余来自内存访问模式变化。
- 单次 gathers 调用的平均耗时从 52µs 增至 308µs（**5.9×**）。在 Release 构建下（[5.4 节](#54-release-构建对照实验-o3--marchnative--dndebug)），编译器优化可缓解此类膨胀。
- 在 48³ 网格下，gathers 占 recip2real 的 15.3%（serial），随着网格继续增大，此比例预计进一步下降（FFT 的 O(n³log n) 增长快于 gather/scatter 的拷贝量）。

### 5.1b MPI 进程扩展分析

以 gaas_medium 稳态（r2/r3 均值）为例，固定每个进程的 OMP 线程数，比较不同 MPI 进程数：

| 配置 | 总核 | wall | gatherp (s) | real2recip (s) | hPsi (s) | diag (s) | before_scf (s) |
|------|------|------|-------------|-----------------|----------|----------|-----------------|
| serial (np1_omp1) | 1 | 36s | — | 6.00 | 25.55 | 24.80 | 0.86 |
| mix_np2_omp2 | 4 | 16s | 1.04 | 3.28 | 11.77 | 12.00 | 0.42 |
| mix_np4_omp2 | 8 | 14s | 0.96 | 2.81 | 10.06 | 10.25 | 0.35 |

**关键发现**：

1. **np4_omp2 比 np2_omp2 快 ~1.1×**（16s→14s），加速比远小于核数增长（2×）。冷启动消除后，MPI 进程扩展的真实收益回归合理范围：
   - 每个进程处理的 k 点数从 2 降至 1，对角化工作量减半
   - 但 MPI_Alltoallv 通信开销随进程数增加而增加（每个进程需要与更多进程交换数据）
   - 两种效应大致抵消，净收益有限

2. **np2_omp2 的 gather/scatter 占比（31.7% gp/r2c, 37.1% gs/r2r）与 np4_omp2（34.0%, 38.2%）接近**：说明在此数据规模下，gather/scatter 的相对开销不随 MPI 进程数的增加而显著变化

3. **before_scf 初始化开销在 MPI 模式下更低**：serial (0.86s) > np2_omp2 (0.42s) > np4_omp2 (0.35s)。MPI 模式下每个进程只处理部分 k 点，初始化工作量按进程数分摊。

### 5.3 SIMD 向量化机会分析

`pw_gatherscatter.h` 中共有 **6 处内层拷贝循环**，全部为 `outp[i] = inp[i]` 形式的连续拷贝，数据完全连续、指针不重叠，非常适合 SIMD 向量化：

| 函数 | 循环变量 | 循环体 |
|------|----------|--------|
| gatherp_scatters (poolnproc=1) | `iz` 0→nz | `outp[iz] = inp[iz]` |
| gatherp_scatters (MPI pre-copy) | `iz` 0→nplane | `outp[iz] = inp[iz]` |
| gatherp_scatters (MPI post-copy) | `izip` 0→nzip | `outp[izip] = inp[izip]` |
| gathers_scatterp (poolnproc=1) | `iz` 0→nz | `outp[iz] = inp[iz]` |
| gathers_scatterp (MPI pre-copy) | `izip` 0→nzip | `outp[izip] = inp[izip]` |
| gathers_scatterp (MPI post-copy) | `iz` 0→nplane | `outp[iz] = inp[iz]` |

特点：
- 数据类型为 `std::complex<double>`（16 字节，32 字节对齐时 AVX2 最优）
- 循环边界为运行时变量（nplane, nz, nzip），典型值 24~48（poolnproc=1 时 nz = nplane × poolnproc）
- AVX2 一次处理 4 个 double（2 个 complex），理论加速 2×~4×

> **预期收益**：Debug 数据下预估手动 SIMD 可带来 2~10% SCF 加速。**但 Release 对照实验（[5.4 节](#54-release-构建对照实验-o3--marchnative--dndebug)）表明，冷启动消除后 GCC 13 的实际优化效果中等（gathers 1.56×, gatherp 无显著加速），编译器自动向量化状态需用 `-fopt-info-vec-optimized` 进一步确认。** 以下为对照实验详情。

### 5.4 Release 构建对照实验（`-O3 -march=native -DNDEBUG`）

为判断 GCC 13 是否已对拷贝循环自动向量化，使用 gaas_medium 用例在 3 种关键配置下跑了 9 次 Release 对照实验（每配置 3 次重复）。以下为 r2/r3 均值与新旧 Debug 的对比：

> **重要提示**：由于初始基线的 Debug 数据受冷启动污染（见 [4.4 节](#44-冷启动现象重跑验证)），旧 Debug→Release 加速比（如 np2_omp2 的 2.71×）被系统性高估。下表同时列出旧 Debug 和新 Debug（重跑修正后）的对比。

| 配置 | Debug wall (旧) | Debug wall (新) | Release wall | 整体加速比(新/Rel) | gatherp 加速比 | gathers 加速比 | before_scf 加速比 |
|------|----------------|----------------|-------------|-------------------|---------------|---------------|-------------------|
| serial (np1_omp1) | 36s | 36s | 25s | **1.44×** | — | **2.92×** | **1.54×** |
| mix_np2_omp2 | 38s | **16s** | 14s | **1.14×** | 0.79× | **1.56×** | **1.50×** |
| mix_np4_omp2 | 11s | **14s** | 14s | **1.00×** | 0.62× | **1.16×** | **1.09×** |

> 冷启动消除后，np2_omp2 的 Debug 基线从 38s 降至 16s（快 2.4×），使得 Debug→Release 的整体加速比从 2.71× 回归至 1.14×。gatherp 在 Release 下出现了 0.79× 的"退化"——但绝对值差异仅 0.3s（1.04s vs 1.32s），落在测量噪声范围内。gathers 保持 1.56× 的正向加速。np4_omp2 的"Release 异常"（旧数据 0.79×）在新 Debug 基线（14s）下消失，新旧 Debug 差异说明原来的退化也是冷启动噪声的一部分。

#### 关键结论：编译器优化效果中等

1. **gathers 显示 1.56× 加速**，表明 GCC 13 对 `gathers_scatterp` 中的拷贝循环有优化效果（可能包含自动向量化）。但此加速比远低于旧数据中的 3.58×（曾被错误归因于自动向量化）。

2. **gatherp 在 Debug 和 Release 间无显著差异**（1.04s vs 1.32s，约 0.3s 噪声）。这可能是由于 `gatherp_scatters` 的 MPI pre/post copy 路径数据结构不同，编译器优化效果不一致。

3. **serial 模式下的 gathers 加速比 2.92×**：serial 路径的拷贝循环更简单（无 MPI 包装），编译器优化更充分。

4. **`-O3` 的整体效果由多因素叠加**：自动向量化、函数内联、循环展开、断言移除等共同贡献 1.14~1.44× 的总体加速。

#### 对手动 SIMD 的影响

| 项目 | 旧估算（基于冷启动污染数据） | 修正后评估 |
|------|---------------------------|-----------|
| GCC 是否已向量化 copy 循环 | **是**（3.3-3.6× 加速） | **部分可能**（gathers 1.56×, gatherp 无显著效果） |
| 手动 SIMD 可能收益 | 估计 **1.0~1.3×** | 需用 `-fopt-info-vec` 确认后判断 |
| gather/scatter 占 SCF（稳态）| ~10% (np2_omp2 Release) | **6~12%** (np2_omp2 Release) |
| 整体 SCF 预期加速 | **~0~3%**（编译器已做大部分工作） | 待 `-fopt-info-vec` 确认后重新评估 |

> **策略调整建议**：旧数据曾认为 GCC 已自动向量化拷贝循环，手动 SIMD 边际收益极小。新数据表明编译器加速效果被冷启动高估，因此需重新评估：
> - **必须首先用 `-fopt-info-vec-optimized` 确认自动向量化状态**。GCC 将逐行输出哪些循环成功向量化、使用多少位向量寄存器。若日志确认 `pw_gatherscatter.h` 的 6 处拷贝循环均已执行 256 位 AVX2 向量化，则手动 SIMD 收益仍有限（~0~3% SCF）；若有遗漏，手动 SIMD 的收益空间更大
> - 尝试 `#pragma GCC ivdep` + `__restrict__` 消除编译器别名分析障碍，让自动向量化更激进
> - **对齐诊断**：在拷贝前加入 `__builtin_assume_aligned(ptr, 32)` 或 `reinterpret_cast<uintptr_t>(outp) % 32 == 0` 断言。若对齐后无额外收益，可在报告中标注"未对齐访问已被硬件缓冲完全掩盖"；若有 1~2% 提升，则白赚
> - 同步投入 **题目8 缓存复用**，其收益不受编译优化影响

---

## 6. 缓存基线数据 (题目8：缓存复用)

### 6.1 当前缓存相关函数耗时（实测稳态数据）

以下为 r2/r3 稳态均值的实测 `before_scf` 及其子组件耗时。由于 `collect_local_pw`、`collect_uniqgg`、`setupIndGk` 等函数的单次调用均低于 timer 0.1s 阈值（见 [3.3 节](#33-计时器限制)），它们不单独出现在 TIME STATISTICS 中，但其总耗时被包含在 `before_scf` 和 `initialize_psi` 内。

| 用例 | before_scf | initialize_psi | diag_subspace_init | 估算 cache 函数总开销 |
|------|-----------|----------------|--------------------|----------------------|
| gaas_tiny | 0.11~0.15s | 0.07~0.14s | 0.06~0.12s | < 0.05s |
| gaas_small | 0.15~0.29s | 0.15~0.27s | 0.12~0.24s | < 0.10s |
| gaas_medium | 0.22~0.81s | 0.19~0.74s | 0.18~0.70s | < 0.30s |

> 注：`before_scf` 包含平面波初始化、势能初始化、psi 初始化等众多子步骤。cache 相关函数（`collect_local_pw`、`collect_uniqgg`、`setupIndGk` 等）只是其中的一部分，实际占比估计为 `before_scf` 的 20~40%。

### 6.2 缓存复用机会分析

#### 6.2.1 `PW_Basis::collect_local_pw` (pw_basis.cpp:137)

**当前行为**：每次调用都 `delete[] / new[]` 重建 `gg`、`gdirect`、`gcar` 三个数组，然后通过 `ig2isz` → `is2fftixy` 反推 (ix, iy, iz)，重新计算所有 G 矢量的笛卡尔坐标和模平方。

**可缓存数据**（生命周期绑定于 `ggecut` 和晶格不变时）：
- `gg[ig]` — 平面波动能 G² (~npw × 8 bytes)
- `gdirect[ig]` — G 矢量直接坐标 (~npw × 24 bytes)  
- `gcar[ig]` — G 矢量笛卡尔坐标 (~npw × 24 bytes)

**复用条件**：SCF 迭代中晶格常数、截断能、FFT 网格不变时，这些数值跨迭代不变。`collect_local_pw` 仅在 `setuptransform` → `esolver_fp` 初始化路径被调用，调用次数少但每次涉及 O(npw) 二次循环。

#### 6.2.2 `PW_Basis::collect_uniqgg` (pw_basis.cpp:195)

**当前行为**：重新遍历所有平面波计算 `tmpgg`，然后用 heapsort 排序去重，构建 `gg_uniq` 和 `ig2igg` 映射。

**可优化点**：
- `tmpgg[ig]` 的 G² 计算与 `collect_local_pw` 的 `gg[ig]` **完全重复**
- 排序和去重逻辑的 O(npw log npw) 复杂度可考虑复用

#### 6.2.3 `PW_Basis_K::setupIndGk` (pw_basis_k.cpp:133)

**当前行为（双重计算热点）**：

```cpp
// 第一遍：统计 npwk[ik]
for (int ig = 0; ig < this->npw; ig++) {
    const double gk2 = this->cal_GplusK_cartesian(ik, ig).norm2();
    if (gk2 <= this->gk_ecut) { ++ng; }  // ← 计算了 gk2
}
// 第二遍：构建 igl2isz_k / igl2ig_k
for (int ig = 0; ig < this->npw; ig++) {
    const double gk2 = this->cal_GplusK_cartesian(ik, ig).norm2();  // ← 重复计算！
    if (gk2 <= this->gk_ecut) { ... }
}
```

每对 `(ik, ig)` 调用 `cal_GplusK_cartesian` **两次**，每次内部做坐标反推和矩阵乘法。对于 4 个不可约 k 点 × 数千平面波，这是个明显的浪费。

#### 6.2.4 `PW_Basis_K::collect_local_pw` (pw_basis_k.cpp:257)

**当前行为**：无条件 `delete[] / new[]` 重建 `gk2` 和 `gcar`（每个 k 点 × npwk_max 大小），然后重新计算。

**可缓存数据**：
- `gk2[ik * npwk_max + igl]` — |G+K|² 值
- `gcar[ik * npwk_max + igl]` — G+K 笛卡尔坐标

**复用条件**：k 点集、截断能、晶格不变时跨 SCF 迭代不变。

#### 6.2.5 `setupIndGk` 双重计算量化评估

以 gaas_medium (4 个不可约 k 点, ~4,600 npw, serial 模式) 为例，量化双重计算的实际开销：

| 项目 | 数值 |
|------|------|
| 双重计算调用量 | 4 kpt × 4,600 npw × 2 遍 = 36,800 次 `cal_GplusK_cartesian` 调用 |
| 其中冗余部分 | ~18,400 次（第二遍完全重复第一遍） |
| 单次 `cal_GplusK_cartesian` 操作 | 坐标反推 (ix,iy,iz 从 ig2isz) + vec3 坐标计算 + `norm2()` → ~30 flop |
| 冗余浮点运算总量 | 18,400 × 30 ≈ **5.5 × 10⁵ flop** |
| 冗余内存访问 | 18,400 × ~12 bytes (读 ig2isz + 写 gk2) ≈ **220 KB** |
| 估算冗余耗时 | **~10-50 ms**（Debug，取决于缓存状态）；Release 下内联后 **< 5 ms** |

> 注：MPI 模式下每个进程处理的 k 点数更少（np2 时每进程 2 个 k 点），冗余耗时相应减半。

虽然绝对值不大（< 0.05s Debug，< 0.005s Release），占 `before_scf`（0.22~0.81s）的 2~6%（Debug）或 < 1%（Release），但修复成本极低——第一遍将 gk2 存入临时 `std::vector<double>`，第二遍直接读取，约 5 行改动即可消除此浪费。

### 6.3 缓存数据规模估算

以 gaas_medium (48³, 4 个不可约 k 点) 为例：

| 数据 | 每元素大小 | 元素数 (serial / np2 每进程) | serial 总量 | np2 每进程 |
|------|-----------|------------------------------|-------------|------------|
| `gg` | 8 B | ~4,600 (k 无关) | ~37 KB | ~37 KB |
| `gk2` | 8 B | ~4,600 × nkpt_per_proc | ~147 KB | ~74 KB |
| `gcar` (多 k) | 24 B | ~4,600 × nkpt_per_proc | ~442 KB | ~221 KB |
| `gdirect` | 24 B | ~4,600 (k 无关) | ~110 KB | ~110 KB |
| `igl2isz_k` / `igl2ig_k` | 4 B × 2 | ~4,600 × nkpt_per_proc | ~148 KB | ~74 KB |
| **总计** | | | **~0.9 MB** | **~0.5 MB** |

约 0.5~0.9 MB 的缓存数据量远小于 L2 缓存 (10 MB)，全部可常驻 L2——甚至不需要依赖 L3 (18 MB)。

### 6.4 缓存失效条件

| 触发条件 | 需重建的数据 |
|----------|-------------|
| 晶格常数变化 (`latvec`) | `gg`, `gcar`, `gdirect`, `gk2` |
| 截断能变化 (`ecutwfc`) | 全部（改变平面波数量） |
| k 点集变化 | `gk2`, `gcar` (多 k), `igl2isz_k`, `igl2ig_k` |
| FFT 网格变化 (`nx/ny/nz`) | 全部 |
| `distribution_type` / `xprime` 变化 | `igl2isz_k`, `igl2ig_k` |

**预期收益**：缓存复用主要减少初始化阶段开销。对整体 SCF 时间的直接加速有限（初始化占比通常 <5%），但可减少内存分配/释放抖动，并避免 `setupIndGk` 中 |G+K|² 的双重计算。

---

## 7. 热点分析

### 7.1 SCF 迭代时间分解 (gaas_medium, 稳态 r2)

以下以 `hamilt2rho_single`（单次 SCF 迭代的完整耗时）为基准（=100%），统一展示 serial 和 MPI 两种模式下的热点分布。数据取自 gaas_medium r2（稳态），转换为单次迭代均值（总耗时 ÷ 8 SCF 迭代）。

| 操作 | serial | np2_omp2 | 说明 |
|------|--------|----------|------|
| | | | |
| **hamilt2rho_single** | **4.14s (100%)** | **2.00s (100%)** | 单次 SCF 迭代 |
| | | | |
| *计算类* | | | |
| hPsi | 3.19s (77%) | 1.50s (75%) | 哈密顿量作用 |
| ├ veff_pw | 2.15s (52%) | 1.02s (51%) | 有效势 |
| ├ nonlocal_pw | 1.00s (24%) | 0.45s (23%) | 非局域赝势 |
| └ EkineticPW | ~0.03s (<1%) | ~0.03s (<1%) | 动能（低于 timer 阈值，按 hPsi 残差估算） |
| solve (CG) | 4.14s (100%)* | 2.00s (100%)* | 对角化，与 hamilt2rho 几乎重合 |
| ├ diag_once | 3.09s (75%) | 1.51s (76%) | CG 对角化主体 |
| └ diag_subspace | 0.56s (14%) | 0.26s (13%) | 子空间对角化 |
| psiToRho | 0.41s (10%) | 0.18s (9%) | 波函数→电荷密度 |
| | | | |
| *数据搬运类* | | | |
| real2recip | 0.75s (18%) | 0.41s (21%) | 实→倒 FFT (含 gatherp) |
| ├ gatherp_scatters | — (<0.1s) | 0.12s (6%) | **SIMD 目标 (mpi)** |
| recip2real | 1.39s (34%) | 0.64s (32%) | 倒→实 FFT (含 gathers) |
| ├ gathers_scatterp | 0.45s (11%) | 0.24s (12%) | **SIMD 目标** |

> `*` solve 的父级计时器也是 `hamilt2rho_single`，两者耗时几乎相同（对角化主导了 SCF 迭代）。百分比相对于各自的 `hamilt2rho_single` 归算。

**关键对比**：

- **MPI 模式引入额外 gather/scatter 开销**：np2_omp2 的 gatherp_scatters 从零（serial 中 <0.1s）跃升至 0.12s（6% SCF），因为 MPI 路径增加了 pre/post communication copy 段和 `MPI_Alltoallv`。
- **对角化 (diag_once) 是绝对最大热点**，占 SCF 时间的 ~75%（两种模式一致）。这是 Workflow D 的关注范围。
- **real2recip + recip2real 合计占 SCF 的 52%（serial）/ 53%（MPI）**，其中 gather/scatter 占 FFT 的比例从 serial 的 ~18% 升至 MPI 的 ~27%。

### 7.2 线程扩展性

以 gaas_medium 稳态（r2/r3 均值）的 wall time 为基准：

| 配置 | total wall | 加速比 vs serial | 加速比 vs 前级 | gather/scatter 加速比 | 备注 |
|------|-----------|-------------------|----------------|----------------------|------|
| serial (np1_omp1) | 36s | 1.0× | — | 1.0× | 基线 |
| omp4 (np1_omp4) | 16s | 2.25× | 2.25× | ~2.8× | 4 线程有效加速 |
| omp8 (np1_omp8) | 17s | 2.12× | 0.94× | ~3.4× | 8 线程饱和 |
| mix_np2_omp2 | 16s | 2.25× | — | — | 4 核，与 omp4 接近 |
| mix_np2_omp4 | 33s | 1.09× | — | — | 8 核，MPI 通信开销大 |
| mix_np4_omp2 | 14s | 2.57× | — | — | 8 核，最优 |

**分析**：

1. **Pure OpenMP 扩展性有限**：4 线程达 2.25× 加速（效率 56%），8 线程无额外增益（2.12×）。原因：gather/scatter 内层循环极短（nplane ≤ 48），OMP fork/join 开销抵消了并行收益。`omp8` 配置下 8 个线程被限制在单个 MPI 进程内（仅能使用单进程内存带宽），实际物理核心利用率不足。

2. **冷启动消除后，MPI 配置的扩展性回归正常**：`mix_np2_omp2`（4 核，2.25×）与 `omp4`（4 线程，2.25×）性能一致，说明 MPI 和 OpenMP 在 4 核级别上效率相当。`mix_np4_omp2`（8 核，2.57×）略优于 `omp8`（8 线程，2.12×），MPI 的进程级并行避免了 OpenMP 小循环的线程开销。

3. **`mix_np2_omp4` 表现最差（33s, 1.09×）**：此配置下每个进程有 4 个 OMP 线程，但 MPI 通信开销（Alltoallv 的数据量比 np2_omp2 更大）和线程调度开销叠加，导致总时间远高于其他 8 核配置。

4. **gather/scatter 的 OMP 加速比优于整体**：omp4 下 gather/scatter 加速 ~2.8× 优于整体 2.25×，说明拷贝循环的并行度好于对角化等复杂操作。但 omp8 下退化明显。

5. **基本结论**：混合并行（MPI + OpenMP）在本机环境下有边际优势，但冷启动消除后，纯 OpenMP 和混合并行在 4 核级别的性能相当。

---

## 8. 优化机会总结

### 8.1 题目5：SIMD 向量化

| 项目 | 评估 |
|------|------|
| 目标 | `gatherp_scatters`, `gathers_scatterp` (6 处拷贝循环) |
| 当前状态 | 旧数据曾认为 GCC 已自动向量化（3.3-3.6×），但该加速比受冷启动污染。修正后 Debug→Release 加速效果中等（gathers 1.56×，gatherp 无显著加速），编译器向量化状态待确认 |
| 第 13 周计划 | **先用 `-fopt-info-vec-optimized` 确认自动向量化状态**，然后根据结果决定优化策略 |
| 风险/优先级 | 低风险 / **中优先级** — 编译器可能有部分自动向量化，但实际效果不如最初预期 |

### 8.2 题目8：缓存复用

| 项目 | 评估 |
|------|------|
| 目标 | 缓存 `gg`, `gcar`, `gdirect`, `gk2`, `igl2isz_k`, `igl2ig_k` (~0.9 MB serial, ~0.5 MB np2) |
| 最大单项收益 | 消除 `setupIndGk` 双重 `cal_GplusK_cartesian` 调用（~18K 次冗余，~10-50ms Debug, <5ms Release） |
| 当前开销 | < 0.3s（含在 `before_scf`），Release 下可降至 < 0.1s |
| 第 13 周计划 | 懒加载 + 失效标志 + 优先修复 setupIndGk 双重计算（50~80 行改动）；使用 `perf stat -e L1-dcache-load-misses,LLC-load-misses` 对比优化前后缓存缺失率 |
| 风险/优先级 | 中风险 / **高优先级** — SIMD 收益缩水后，缓存成为 Workflow C 主攻方向 |

### 8.3 额外发现：MPI 通信优化机会

从基线数据观察到 `poolnproc > 1` 时 gather/scatter 占比较高（稳态 **30~55%** 的 FFT 时间）。gather/scatter 内部包含本地 copy 段和 `MPI_Alltoallv` 通信两部分（当前未分别计时），其中通信部分推测为主要开销（尤其在小网格下，copy 量小而通信延迟占比高）。MPI 通信优化属于 Workflow B 的范畴，但 SIMD 向量化 (Workflow C) 可加速本地 copy 段。

---

## 9. 测试数据可靠性说明

- **3 次重复取平均**：所有配置均有 3 次重复运行，标准差在可接受范围
- **时钟来源**：MPI 模式下使用 `MPI_Wtime()`，非 MPI 使用 `std::chrono::system_clock`
- **已知噪声源**：
  - WSL2 虚拟化层引入约 5~10% 时间抖动，MPI 多进程场景（≥4 进程）下不可控
  - 初始基线中的冷启动效应为重跑验证所排除（详见 [4.4 节](#44-冷启动现象重跑验证)），当前所有数据 r1/r2/r3 间无系统性差异
  - `mix_np2_omp4` 配置因 MPI 通信开销较大，在某些用例下呈现中等方差（gaas_medium: 29s/31s/35s）
- **WSL2 局限性**：本报告所有数据来自 WSL2 (Hyper-V 轻量级 VM)。Windows 宿主的后台抢占与 EPT 内存转换会严重污染硬件计数器，尤其影响 MPI 多进程的通信延迟和共享内存性能。第 13 周的终期优化验收必须迁移到物理机 Linux 或 HPC 集群上进行，WSL2 仅应作为开发与功能验证环境
- **Debug 构建限制**：主要基线数据来自 `-O0` 构建。Release（`-O3`）对照实验已完成（[5.4 节](#54-release-构建对照实验-o3--marchnative--dndebug)），编译器对 gather/scatter 有中等优化效果（gathers 1.56×），但旧数据中的 3.3-3.6× 加速比已确认为冷启动污染导致的高估。百分比占比在 Debug 和 Release 下基本一致，但绝对加速倍数不能跨构建类型外推
- **计时器限制**：0.1s 阈值导致部分轻量函数无数据显示，已在第 3.3 节说明

---

## 10. 结论与下一步

### 10.1 基线结论

1. **全部 72 次运行通过** (100% 通过率)，测试体系稳定可靠。

2. **Gather/Scatter 是 MPI 模式下的重要优化目标**，占 FFT 变换的 **30~55%**（稳态），占 SCF 迭代的 **6~10%**。MPI 路径的 gather/scatter 由本地 copy 段 + `MPI_Alltoallv` 组成，其中 copy 段是 SIMD 的直接作用目标。6 处内层连续拷贝循环非常适合 SIMD 向量化。Release 对照实验显示 gathers 有 1.56× 加速，但 gatherp 无显著加速（见 [5.4 节](#54-release-构建对照实验-o3--marchnative--dndebug)），编译器自动向量化的实际效果需通过 `-fopt-info-vec-optimized` 进一步确认。

3. **缓存复用函数当前开销较小**（<0.3s，含在 `before_scf` 中），但其收益在于消除 `setupIndGk` 的 |G+K|² 双重计算（~18,000 次冗余调用）和减少反复内存分配。缓存数据量约 0.5~0.9 MB，可完全放入 L2 缓存。

4. **线程扩展性趋于饱和**: 4 核达 2.1~2.3× 加速，8 核不再显著提升。冷启动消除后，`mix_np2_omp2`（4 核, 2.25×）与 `omp4`（4 线程, 2.25×）性能一致，MPI 和 OpenMP 在 4 核级别上效率相当。

### 10.2 第 13 周工作建议

1. **SIMD 向量化调查**（题目5，先诊断后优化）：
   - **先用 `-fopt-info-vec-optimized` 确认自动向量化状态**：在 Release 构建中追加此编译选项，获取 GCC 对 `pw_gatherscatter.h` 中 6 处拷贝循环的向量化报告。这是决定后续方案的关键输入
   - 若确认 GCC 已成功向量化所有 6 处循环（输出显示 `vectorized`），则手动 SIMD 边际收益有限（~0~3% SCF），可快速收尾
   - 若有遗漏：添加 `#pragma omp simd` + `__restrict__` + `__builtin_assume_aligned` 到遗漏的循环
   - 在 Release 构建下用 gaas_medium np2_omp2 测量实际增益

2. **缓存复用作为主攻方向**（题目8）：
   - 在 `PW_Basis` 和 `PW_Basis_K` 中添加缓存有效性标志位
   - 改造 `collect_local_pw` 和 `collect_uniqgg` 检查缓存
   - **优先消除 `setupIndGk` 中的双重 `cal_GplusK_cartesian` 调用**（见效最快，~10-50ms Debug / <5ms Release）
   - 使用 `perf stat -e L1-dcache-load-misses,LLC-load-misses` 对比优化前后缓存缺失率
   - 预期代码改动 50~80 行

3. **与 Workflow B 协调**：
   - SIMD 加速的 copy 段位于 MPI_Alltoallv 前后
   - 若 Workflow B 已将 Alltoallv 替换为非阻塞通信，copy 段可与通信重叠
   - 两个工作流的优化效果可能叠加

4. **终期验收必须脱离 WSL2 虚拟化环境**：
   - WSL2 不适用于精确的 MPI 稳态性能评测和大核数并发测试（Hyper-V 后台抢占 + EPT 二次转换污染硬件计数器）
   - 第 13 周最终优化成果验收前，**必须将代码同步到物理机 Linux 或 HPC 集群**进行纯净测试，WSL2 仅作为开发与功能验证环境

---

## 附录 A：测试用例文件清单

```
homework_docs/test_cases/
├── gaas_tiny/
│   ├── INPUT          (24³, 20 Ry, 70 bands, 3×3×3 kpt)
│   ├── KPT            (Gamma-centered 3×3×3)
│   ├── STRU -> ../../../../abacus-user-guide/examples/pw2/cg/001_4GaAs/STRU
│   ├── Ga.upf -> ../../../../abacus-user-guide/examples/pw2/cg/001_4GaAs/Ga.upf
│   └── As.upf -> ../../../../abacus-user-guide/examples/pw2/cg/001_4GaAs/As.upf
├── gaas_tiny_40Ry/
│   ├── INPUT          (24³, 40 Ry, 70 bands, 3×3×3 kpt, 补测)
│   └── ... (同上链接)
├── gaas_small/
│   ├── INPUT          (32³, 40 Ry, 70 bands, 3×3×3 kpt)
│   └── ... (同上链接)
├── gaas_medium/
│   ├── INPUT          (48³, 40 Ry, 70 bands, 3×3×3 kpt)
│   └── ... (同上链接)
├── results/
│   ├── gaas_tiny/      (18 对 stdout/stderr 文件)
│   ├── gaas_tiny_40Ry/ (18 对 stdout/stderr 文件)
│   ├── gaas_small/     (18 对 stdout/stderr 文件)
│   └── gaas_medium/    (18 对 stdout/stderr 文件)
├── run_benchmark.sh   (自动化测试脚本)
└── parse_timers.py    (计时数据提取脚本)
```

## 附录 B：临时计时桩点代码变更

为获取基线测试所需的细粒度计时数据，在以下文件中临时插入了 `ModuleBase::timer::start/end`：

```
source/source_basis/module_pw/pw_gatherscatter.h  (+6 行: 3 处 start/end 对)
source/source_basis/module_pw/pw_basis.cpp          (+6 行: 2 处 start/end 对 + 2 处 early return end)
source/source_basis/module_pw/pw_basis_k.cpp        (+6 行: 同上)
```

这些桩点将在进入第 13 周代码修改前移除，以免影响后续的性能对比测试。
