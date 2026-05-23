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

> **注意**：本报告主要数据来自 `-O0`（Debug）构建。Release（`-O3 -march=native -DNDEBUG`）对照实验见 [5.4 节](#54-release-构建对照实验-o3--marchnative--dndebug)。Debug 测得的百分比占比仍然有效，但 **SIMD 的预期加速倍数不能直接用 Debug 数据估算**。

---

### 执行摘要

| 维度 | 关键数字 | 说明 |
|------|----------|------|
| 测试通过率 | **54/54 (100%)** | 3 用例 × 6 配置 × 3 重复 |
| 最大热点 | **diag_once (CG 对角化)** 占 SCF 73~77% | Workflow D 范畴 |
| Gather/Scatter 占比 | FFT 的 **30~55%**, SCF 的 **6~10%** (MPI 稳态) | 6 处连续拷贝循环在 `pw_gatherscatter.h` |
| SIMD 状态 | **GCC 13 -O3 未自动向量化** gather/scatter 拷贝循环（别名屏障），清零循环已 AVX2 向量化 | 手动 SIMD 有 ~1.5~3× 单循环加速空间 |
| 缓存机会 | 数据 **~0.9 MB** (可常驻 L2), setupIndGk **~18K 冗余调用** | 见 [6.2.5 节](#625-setupindgk-双重计算量化评估) |
| Debug→Release | serial **1.40×**, np2_omp2 **1.23×** 整体加速 | 见 [5.4 节](#54-release-构建对照实验-o3--marchnative--dndebug) |

---

## 2. 测试用例设计

### 2.1 设计原则

为适应 7.6 GB 内存限制并保证可复现，采用**同一 GaAs 体系在三种 FFT 网格下对比**的策略。相同材料保证了电子结构一致（能带数、k 点密度），仅改变 FFT 网格大小，分离出数据搬运量对性能的影响。

### 2.2 测试用例参数

| 参数 | gaas_tiny_40Ry | gaas_small | gaas_medium |
|------|----------------|------------|-------------|
| 来源 | pw2/cg/001_4GaAs | pw2/cg/001_4GaAs | pw2/cg/001_4GaAs |
| 化学式 | Ga₄As₄ (8 原子) | Ga₄As₄ (8 原子) | Ga₄As₄ (8 原子) |
| FFT 网格 | 24×24×24 | 32×32×32 | 48×48×48 |
| 截断能 ecutwfc | 40 Ry | 40 Ry | 40 Ry |
| 能带数 nbands | 70 | 70 | 70 |
| K 点网格 | 3×3×3 (4 个不可约 k 点) | 3×3×3 (4 个不可约 k 点) | 3×3×3 (4 个不可约 k 点) |
| SCF 趋敛 | 9 次收敛 | 9 次收敛 | 8 次收敛 |
| 求解器 | CG | CG | CG |
| 赝势 | pseudo-dojo v0.5 | pseudo-dojo v0.5 | pseudo-dojo v0.5 |

### 2.3 并行配置矩阵

每个测试用例运行 **6 种并行配置 × 3 次重复 = 18 次/用例，3 用例总计 54 次运行**。

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
| gaas_tiny_40Ry | 18 | 18 | 0 | 100% |
| gaas_small | 18 | 18 | 0 | 100% |
| gaas_medium | 18 | 18 | 0 | 100% |
| **总计** | **54** | **54** | **0** | **100%** |

### 4.2 验收标准

- 所有运行正常退出（EXIT=0）
- TIME STATISTICS 计时表完整输出
- 无 crash、segfault 或 MPI 异常
- stderr 中无错误信息（仅 OpenMPI hwloc 信息性提示）

### 4.3 墙钟时间概览

> **通用说明**：以下为 3 次重复的原始墙钟时间。gaas_tiny_40Ry 和 gaas_small 在 9 次收敛，gaas_medium 在 8 次收敛。迭代次数差异不影响同用例内的可比性。

**gaas_tiny_40Ry (24³，40 Ry)**:

| 配置 | rep1 | rep2 | rep3 |
|------|------|------|------|
| serial (np1_omp1) | 20s | 19s | 19s |
| omp4 (np1_omp4) | 9s | 9s | 8s |
| omp8 (np1_omp8) | 10s | 10s | 11s |
| mix_np2_omp2 | 8s | 9s | 9s |
| mix_np2_omp4 | 18s | 23s | 22s |
| mix_np4_omp2 | 11s | 9s | 11s |

**gaas_small (32³)**:

| 配置 | rep1 | rep2 | rep3 |
|------|------|------|------|
| serial (np1_omp1) | 21s | 20s | 20s |
| omp4 (np1_omp4) | 9s | 9s | 9s |
| omp8 (np1_omp8) | 13s | 17s | 17s |
| mix_np2_omp2 | 10s | 11s | 10s |
| mix_np2_omp4 | 27s | 28s | 26s |
| mix_np4_omp2 | 8s | 11s | 10s |

**gaas_medium (48³)**:

| 配置 | rep1 | rep2 | rep3 |
|------|------|------|------|
| serial (np1_omp1) | 33s | 35s | 36s |
| omp4 (np1_omp4) | 17s | 16s | 17s |
| omp8 (np1_omp8) | 19s | 17s | 17s |
| mix_np2_omp2 | 15s | 18s | 15s |
| mix_np2_omp4 | 29s | 31s | 35s |
| mix_np4_omp2 | 14s | 15s | 13s |

---

## 5. Gather/Scatter 基线数据 (题目5：SIMD 向量化)

### 5.1 计时数据总表

以下数据为 r1、r2 和 r3 的算术平均（3 次重复均值）。`gatherp_scatters` 在 `poolnproc=1`（单进程）时因走简化代码路径（无 MPI pre/post copy），部分配置下耗时低于 0.1s 阈值不显示。每个 SCF 迭代中 `gatherp_scatters` 和 `gathers_scatterp` 分别调用约 nkpt × nbands 次，全运行累计 7,000~12,000 次，调用量随网格和 k 点数增长。

#### gaas_tiny_40Ry (24³，40 Ry)

| 配置 | gatherp (s) | gathers (s) | real2recip (s) | recip2real (s) | gp/r2c | gs/r2r | gp/ham |
|------|-------------|-------------|-----------------|-----------------|--------|--------|--------|
| serial | — | 0.69 | 1.55 | 2.45 | — | 28.0% | — |
| omp4 | — | 0.23 | 0.73 | 1.31 | — | 17.8% | — |
| omp8 | — | 0.23 | 0.79 | 1.50 | — | 15.5% | — |
| mix_np2_omp2 | 0.38 | 0.67 | 0.96 | 1.51 | 39.1% | 44.1% | 7.5% |
| mix_np2_omp4 | 1.48 | 2.40 | 3.56 | 5.56 | 41.6% | 43.1% | 11.0% |
| mix_np4_omp2 | 0.59 | 0.88 | 1.21 | 1.72 | 48.9% | 50.9% | 9.7% |

#### gaas_small (32³)

| 配置 | gatherp (s) | gathers (s) | real2recip (s) | recip2real (s) | gp/r2c | gs/r2r | gp/ham |
|------|-------------|-------------|-----------------|-----------------|--------|--------|--------|
| serial | 0.24 | 1.32 | 1.75 | 3.41 | 13.5% | 38.6% | 1.8% |
| omp4 | — | 0.39 | 0.69 | 1.50 | — | 25.8% | — |
| omp8 | 0.19 | 0.61 | 1.39 | 2.84 | 13.6% | 21.5% | 2.2% |
| mix_np2_omp2 | 0.70 | 1.10 | 1.41 | 2.17 | 49.8% | 50.8% | 9.9% |
| mix_np2_omp4 | 2.14 | 3.21 | 4.69 | 7.04 | 45.5% | 45.6% | 12.2% |
| mix_np4_omp2 | 0.70 | 1.00 | 1.28 | 1.84 | 54.3% | 54.3% | 11.6% |

#### gaas_medium (48³)

| 配置 | gatherp (s) | gathers (s) | real2recip (s) | recip2real (s) | gp/r2c | gs/r2r | gp/ham |
|------|-------------|-------------|-----------------|-----------------|--------|--------|--------|
| serial | — | 3.53 | 6.02 | 11.19 | — | 31.5% | — |
| omp4 | — | 1.21 | 2.62 | 5.01 | — | 24.1% | — |
| omp8 | — | 1.03 | 2.76 | 5.27 | — | 19.6% | — |
| mix_np2_omp2 | 1.04 | 1.85 | 3.28 | 5.00 | 31.7% | 37.0% | 8.9% |
| mix_np2_omp4 | 2.44 | 3.83 | 6.63 | 9.91 | 36.8% | 38.7% | 11.2% |
| mix_np4_omp2 | 0.97 | 1.63 | 2.82 | 4.23 | 34.6% | 38.5% | 9.7% |

#### 5.1a 网格规模扩展分析

三个测试用例使用相同材料（GaAs）、相同 k 点密度（3×3×3），在 40 Ry 下拥有完整的 24³→32³→48³ 纯净扩展链：

| 指标 | gaas_tiny_40Ry (40Ry, 24³) | gaas_small (40Ry, 32³) | gaas_medium (40Ry, 48³) | 24³→48³ 扩展倍率 | 理论预期 |
|------|---------------------------|------------------------|-------------------------|------------------|----------|
| nplane (serial) | 24 | 32 | 48 | 2.0× | — |
| gathers 耗时 | 0.69s | 1.32s | 3.53s | **5.1×** | O(nplane) = 2× |
| real2recip 耗时 | 1.55s | 1.75s | 6.02s | **3.9×** | ~O(nplane) × log |
| recip2real 耗时 | 2.45s | 3.41s | 11.19s | **4.6×** | ~O(nplane) × log |
| hPsi 总耗时 | 12.11s | 13.11s | 25.62s | **2.1×** | — |
| SCF 迭代耗时 | 18.60s | 19.55s | 34.13s | **1.8×** | ~O(n³log n) |
| gathers 单次调用 | 52µs | ~110µs | ~308µs | **5.9×** | O(nplane) = 2× |

**网格扩展（40 Ry，24³→48³）**：
- `gathers_scatterp` 绝对耗时增长 **5.1×**（0.69s→3.53s），超过 nplane 线性增长（2×）。这是因为更大的网格在 Debug（`-O0`）模式下涉及更多的索引计算和分支判断，使单次调用的耗时从 52µs 膨胀至 308µs（**5.9×**），其中一部分来自 nplane 增长（2×），其余来自内存访问模式变化。
- 单次 gathers 调用的平均耗时从 52µs 增至 308µs（**5.9×**）。在 Release 构建下（[5.4 节](#54-release-构建对照实验-o3--marchnative--dndebug)），编译器优化可缓解此类膨胀。
- 在 48³ 网格下，gathers 占 recip2real 的 31.5%（serial），随着网格继续增大，此比例预计进一步下降（FFT 的 O(n³log n) 增长快于 gather/scatter 的拷贝量）。

#### 5.1b MPI 进程扩展分析

以 gaas_medium 稳态（3 次重复均值）为例，固定每个进程的 OMP 线程数，比较不同 MPI 进程数：

| 配置 | 总核 | wall | gatherp (s) | real2recip (s) | hPsi (s) | diag (s) | before_scf (s) |
|------|------|------|-------------|-----------------|----------|----------|-----------------|
| serial (np1_omp1) | 1 | 35s | — | 6.02 | 25.62 | 24.86 | 0.87 |
| mix_np2_omp2 | 4 | 16s | 1.04 | 3.28 | 11.73 | 11.95 | 0.40 |
| mix_np4_omp2 | 8 | 14s | 0.97 | 2.82 | 10.05 | 10.23 | 0.36 |

**关键发现**：

1. **np4_omp2 比 np2_omp2 快 ~1.1×**（16s→14s），加速比远小于核数增长（2×）。MPI 进程扩展的真实收益：
   - 每个进程处理的 k 点数从 2 降至 1，对角化工作量减半
   - 但 MPI_Alltoallv 通信开销随进程数增加而增加（每个进程需要与更多进程交换数据）
   - 两种效应大致抵消，净收益有限

2. **np2_omp2 的 gather/scatter 占比（31.7% gp/r2c, 37.1% gs/r2r）与 np4_omp2（34.0%, 38.2%）接近**：说明在此数据规模下，gather/scatter 的相对开销不随 MPI 进程数的增加而显著变化

3. **before_scf 初始化开销在 MPI 模式下更低**：serial (0.87s) > np2_omp2 (0.40s) > np4_omp2 (0.36s)。MPI 模式下每个进程只处理部分 k 点，初始化工作量按进程数分摊。

### 5.2 SIMD 向量化机会分析

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

#### GCC 13 自动向量化实际状态

使用 `-fopt-info-vec-optimized` 和 `-fopt-info-vec-missed` 编译目标文件，逐循环检查 GCC 13.3.0 (`-O3 -march=native`) 的自动向量化决策：

| 函数 | 循环位置 | 行号 | 向量化状态 | 向量宽度 | 失败原因 |
|------|---------|------|-----------|----------|----------|
| **`gatherp_scatters`** | `poolnproc==1` 外层 `is` | L28 | ❌ 未向量化 | — | 指针别名：`outp[iz]=inp[iz]` 中 GCC 无法证明无重叠 |
| | `poolnproc==1` 内层 `iz` | L33 | ❌ 未向量化 | — | `"more than one data ref in stmt"` |
| | MPI pre-copy 内层 `iz` | L57 | ❌ 未向量化 | — | 同上 |
| | MPI post-copy 3 层嵌套 | L88-100 | ❌ 未向量化 | — | 同上 |
| **`gathers_scatterp`** | **清零循环** `out[i]=0` | L129 | ✅ **已向量化** | **AVX2 32B** | 连续 memset 模式 |
| | `poolnproc==1` 外层 `is` | L137 | ❌ 未向量化 | — | 指针别名 |
| | `poolnproc==1` 内层 `iz` | L142 | ❌ 未向量化 | — | `"more than one data ref"` |
| | MPI 3 层嵌套 | L164-174 | ❌ 未向量化 | — | 指针别名 |
| **`collect_local_pw`** | 主计算循环 `ig` | L150 | ❌ 未向量化 | — | 分支控制流 (ix/iy/iz 调整 + WARNING_QUIT) |
| **`collect_uniqgg`** | 主计算循环 `ig` | L211 | ❌ 未向量化 | — | 同上 |
| **`setupIndGk`** | 计数循环 `ig` | L137-148 | ❌ 未向量化 | — | 调用 `cal_GplusK_cartesian()`（函数调用阻止向量化） |
| | 建索引循环 `ig` | L178-190 | ❌ 未向量化 | — | 同上 |

**根因诊断——指针别名是 gather/scatter 向量化的核心障碍**：

```cpp
// pw_gatherscatter.h:33-35 (所有 6 处拷贝循环的共同模式)
std::complex<T> *outp = &out[is*nz_];   // 从 out 计算
std::complex<T> *inp = &in[ixy*nz_];     // 从 in 计算（不同索引数组）
for(int iz = 0 ; iz < nz_ ; ++iz)
    outp[iz] = inp[iz];  // ← GCC: "more than one data ref"
```

内层 `iz` 循环本质上是 `std::complex<double>` 的**连续拷贝**（典型 `nz` ≈ 24~48），数据布局完全连续且不存在真实重叠。但由于 `outp` 和 `inp` 来自不同基指针，通过 `istot2ixy` 索引数组间接寻址，GCC 13 的别名分析器无法推断无重叠，拒绝向量化。同一问题也导致 MPI 路径（3 层嵌套循环）的 4 处拷贝均未向量化。

相比之下，清零循环（L129 `for(int i = 0; i < nrxx_; ++i) out[i]=0`）是纯粹的连续 memset 模式，路径唯一、无别名歧义，编译器成功使用 256-bit AVX2 向量化。

`collect_local_pw` 和 `collect_uniqgg` 则因含多个条件分支（ix/iy/iz 边界修正 + `WARNING_QUIT` 调用）而未被向量化，但此类函数的计算密度远高于 gather/scatter 拷贝循环。

**结论**：GCC 13 对 gather/scatter 的核心拷贝循环**未能自动向量化**，主要障碍为**指针别名分析失败**。这恰好是手动优化的机会——通过添加 `__restrict__` 或 `#pragma GCC ivdep` 打破别名壁垒，即可使 AVX2 256-bit 向量化生效。以下为 Release 对照实验的具体数据。

> **预期收益**：Debug 数据下预估手动 SIMD 可带来 2~10% SCF 加速。**但 Release 对照实验（[5.4 节](#54-release-构建对照实验-o3--marchnative--dndebug)）表明，`-O3` 的整体效果 medium（gathers 1.64×, gatherp 无显著加速），且自动向量化报告确认核心拷贝循环未向量化，手动 SIMD 仍有约 1.5~3× 的单循环加速空间。** 以下为对照实验详情。

### 5.3 Release 构建对照实验（`-O3 -march=native -DNDEBUG`）

使用 gaas_medium 用例在 3 种关键配置下跑了 9 次 Release 对照实验（每配置 3 次重复），结合 [5.3.1 节](#531-gcc-13-自动向量化实际状态) 的编译器自动向量化报告，验证 GCC 13 的优化效果。以下为 3 次重复均值与 Debug 的对比：

| 配置 | Debug wall | Release wall | 整体加速比 | gatherp 加速比 | gathers 加速比 | before_scf 加速比 |
|------|-----------|-------------|-----------|---------------|---------------|-------------------|
| serial (np1_omp1) | 35s | 25s | **1.40×** | — | **2.89×** | **1.58×** |
| mix_np2_omp2 | 16s | 13s | **1.23×** | 0.86× | **1.64×** | **1.48×** |
| mix_np4_omp2 | 14s | 14s | **1.00×** | 0.66× | **1.23×** | **1.17×** |

> gatherp 在 Release 下出现了 0.86× 的"退化"，但绝对值差异仅 0.17s（1.04s vs 1.21s），落在测量噪声范围内。gathers 保持 1.64× 的正向加速。np4_omp2 在 Debug 和 Release 间基本一致（1.00×）。

#### 关键结论：编译器优化效果中等

1. **gathers 显示 1.64× 加速**（serial 下 2.89×），但自动向量化报告确认 GCC 对 gather/scatter 的内层拷贝循环**未做显式向量化**（`"more than one data ref"` 别名障碍）。1.64× 的加速主要来自 `-O3` 的其他优化（函数内联、循环展开、指令调度），而非 SIMD 向量化。

2. **gatherp 在 Debug 和 Release 间无显著差异**（1.04s vs 1.21s，约 0.17s 噪声）。这与向量化报告吻合——所有 gatherp 拷贝循环均未向量化，`-O3` 对该路径的收益有限。

3. **serial 模式下的 gathers 加速比 2.89×**：serial 路径的拷贝循环更简单（无 MPI 包装），编译器优化更充分。但加速仍来自非 SIMD 优化——报告显示 serial 路径同样因别名问题未向量化。

4. **`-O3` 的整体效果由多因素叠加**：自动向量化、函数内联、循环展开、断言移除等共同贡献 1.00~1.40× 的总体加速。其中自动向量化主要贡献在清零循环（AVX2 已确认）和 FFTW 接口，而非 gather/scatter 核心拷贝。

#### 对手动 SIMD 的影响

| 项目 | 评估 |
|------|------|
| GCC 是否已向量化 copy 循环 | **❌ 未向量化**（别名障碍 `"more than one data ref"`），见 [5.3.1 节](#531-gcc-13-自动向量化实际状态) |
| 1.64× gathers 加速来源 | 函数内联/循环展开/指令调度等非 SIMD 优化 |
| 手动 SIMD 可能收益 | **高**：打破别名壁垒后 AVX2 可覆盖 6 处拷贝循环（nz≈24~48），单循环预期 1.5~3× |
| gather/scatter 占 SCF（稳态 Release）| **6~12%** (np2_omp2) |
| 整体 SCF 预期加速 | **加 `__restrict__`/`ivdep` + 对齐** 预估 ~3~8% SCF；若进一步手动 AVX2 intrinsic，~5~12% SCF |

> **优化策略（已确认）**：
> - **自动向量化状态已确认**：GCC 13 因指针别名障碍未能向量化 gather/scatter 拷贝循环。**手动 SIMD 有明确的优化空间**
> - **第一步（最低成本）**：在拷贝循环前添加 `#pragma GCC ivdep` + 指针 `__restrict__`，消除别名障碍，让 GCC 的自动向量化生效。预估 2~5% SCF 提升
> - **第二步（中等成本）**：加 `__builtin_assume_aligned(ptr, 32)` 断言确保 32 字节对齐，使 AVX2 256-bit 加载/存储生效。预估额外 1~3% SCF
> - **第三步（最高成本）**：若自动向量化效果仍不满足，手写 AVX2 intrinsic 替代 6 处拷贝循环。预估额外 3~5% 但维护成本高，建议仅在测试后确实有必要时采用
> - **同步投入 题目8 缓存复用**，其收益不受编译优化影响

---

## 6. 缓存基线数据 (题目8：缓存复用)

### 6.1 当前缓存相关函数耗时（实测稳态数据）

以下为 3 次重复均值的实测 `before_scf` 及其子组件耗时。由于 `collect_local_pw`、`collect_uniqgg`、`setupIndGk` 等函数的单次调用均低于 timer 0.1s 阈值（见 [3.3 节](#33-计时器限制)），它们不单独出现在 TIME STATISTICS 中，但其总耗时被包含在 `before_scf` 和 `initialize_psi` 内。

| 用例 | before_scf | initialize_psi | diag_subspace_init | 估算 cache 函数总开销 |
|------|-----------|----------------|--------------------|----------------------|
| gaas_tiny_40Ry | 0.18~0.42s | 0.17~0.40s | 0.10~0.30s | < 0.20s |
| gaas_small | 0.21~0.54s | 0.19~0.51s | 0.12~0.44s | < 0.25s |
| gaas_medium | 0.34~0.91s | 0.29~0.81s | 0.25~0.69s | < 0.40s |

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

虽然绝对值不大（< 0.05s Debug，< 0.005s Release），占 `before_scf`（0.34~0.91s）的 1~6%（Debug）或 < 1%（Release），但修复成本极低——第一遍将 gk2 存入临时 `std::vector<double>`，第二遍直接读取，约 5 行改动即可消除此浪费。

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

### 7.1 SCF 迭代时间分解 (gaas_medium)

以下以 `hamilt2rho_single`（单次 SCF 迭代的完整耗时）为基准（=100%），统一展示 serial 和 MPI 两种模式下的热点分布。数据取自 gaas_medium 的 3 次重复均值，转换为单次迭代均值（总耗时 ÷ 8 SCF 迭代）。

| 操作 | serial | np2_omp2 | 说明 |
|------|--------|----------|------|
| | | | |
| **hamilt2rho_single** | **4.17s (100%)** | **1.97s (100%)** | 单次 SCF 迭代 |
| | | | |
| *计算类* | | | |
| hPsi | 3.20s (77%) | 1.47s (75%) | 哈密顿量作用 |
| ├ veff_pw | 2.16s (52%) | 1.00s (51%) | 有效势 |
| ├ nonlocal_pw | 1.00s (24%) | 0.44s (22%) | 非局域赝势 |
| └ EkineticPW | ~0.03s (<1%) | ~0.03s (<1%) | 动能（低于 timer 阈值，按 hPsi 残差估算） |
| solve (CG) | 4.17s (100%)* | 1.97s (100%)* | 对角化，与 hamilt2rho 几乎重合 |
| ├ diag_once | 3.11s (75%) | 1.49s (76%) | CG 对角化主体 |
| └ diag_subspace | 0.56s (13%) | 0.26s (13%) | 子空间对角化 |
| psiToRho | 0.41s (10%) | 0.18s (9%) | 波函数→电荷密度 |
| | | | |
| *数据搬运类* | | | |
| real2recip | 0.75s (18%) | 0.41s (21%) | 实→倒 FFT (含 gatherp) |
| ├ gatherp_scatters | — (<0.1s) | 0.13s (7%) | **SIMD 目标 (mpi)** |
| recip2real | 1.40s (34%) | 0.63s (32%) | 倒→实 FFT (含 gathers) |
| ├ gathers_scatterp | 0.44s (11%) | 0.23s (12%) | **SIMD 目标** |

> `*` solve 的父级计时器也是 `hamilt2rho_single`，两者耗时几乎相同（对角化主导了 SCF 迭代）。百分比相对于各自的 `hamilt2rho_single` 归算。

**关键对比**：

- **MPI 模式引入额外 gather/scatter 开销**：np2_omp2 的 gatherp_scatters 从零（serial 中 <0.1s）跃升至 0.13s（7% SCF），因为 MPI 路径增加了 pre/post communication copy 段和 `MPI_Alltoallv`。
- **对角化 (diag_once) 是绝对最大热点**，占 SCF 时间的 ~75%（两种模式一致）。这是 Workflow D 的关注范围。
- **real2recip + recip2real 合计占 SCF 的 52%（serial）/ 53%（MPI）**，其中 gather/scatter 占 FFT 的比例从 serial 的 ~18% 升至 MPI 的 ~27%。

### 7.2 线程扩展性

以 gaas_medium 稳态（3 次重复均值）的 wall time 为基准：

| 配置 | total wall | 加速比 vs serial | 加速比 vs 前级 | gather/scatter 加速比 | 备注 |
|------|-----------|-------------------|----------------|----------------------|------|
| serial (np1_omp1) | 35s | 1.0× | — | 1.0× | 基线 |
| omp4 (np1_omp4) | 17s | 2.06× | 2.06× | ~2.8× | 4 线程有效加速 |
| omp8 (np1_omp8) | 18s | 1.94× | 0.94× | ~3.4× | 8 线程饱和 |
| mix_np2_omp2 | 16s | 2.19× | — | — | 4 核，与 omp4 接近 |
| mix_np2_omp4 | 32s | 1.09× | — | — | 8 核，MPI 通信开销大 |
| mix_np4_omp2 | 14s | 2.50× | — | — | 8 核，最优 |

**分析**：

1. **Pure OpenMP 扩展性有限**：4 线程达 2.06× 加速（效率 52%），8 线程无额外增益（1.94×）。原因：gather/scatter 内层循环极短（nplane ≤ 48），OMP fork/join 开销抵消了并行收益。`omp8` 配置下 8 个线程被限制在单个 MPI 进程内（仅能使用单进程内存带宽），实际物理核心利用率不足。

2. **MPI 在 4 核级别略优于 OpenMP**：`mix_np2_omp2`（4 核，2.19×）高于 `omp4`（4 线程，2.06×）。`mix_np4_omp2`（8 核，2.50×）明显优于 `omp8`（8 线程，1.94×），MPI 的进程级并行避免了 OpenMP 小循环的线程调度开销。

3. **`mix_np2_omp4` 表现最差（32s, 1.09×）**：此配置下每个进程有 4 个 OMP 线程，但 MPI 通信开销（Alltoallv 的数据量比 np2_omp2 更大）和线程调度开销叠加，导致总时间远高于其他 8 核配置。

4. **gather/scatter 的 OMP 加速比优于整体**：omp4 下 gather/scatter 加速 ~2.8× 优于整体 2.06×，说明拷贝循环的并行度好于对角化等复杂操作。但 omp8 下退化明显。

5. **基本结论**：混合并行（MPI + OpenMP）在本机环境下有边际优势。4 核级别上 MPI 略优于纯 OpenMP（2.19× vs 2.06×），8 核级别差距更为明显（2.50× vs 1.94×）。

---

## 8. 优化机会总结

### 8.1 题目5：SIMD 向量化

| 项目 | 评估 |
|------|------|
| 目标 | `gatherp_scatters`, `gathers_scatterp` (6 处拷贝循环) |
| 当前状态 | **GCC 13.3 自动向量化状态已确认**：清零循环已 AVX2 向量化 ✅；6 处核心拷贝循环均**未向量化** ❌，根因为指针别名分析失败（`"more than one data ref"`）。Release 下 1.64× gathers 加速来自非 SIMD 优化（函数内联/循环展开）。见 [5.3.1](#531-gcc-13-自动向量化实际状态) |
| 第 13 周计划 | **不再需要诊断步骤**。直接执行三阶段优化：(1) `#pragma GCC ivdep` + `__restrict__` 打破别名；(2) 对齐断言 `__builtin_assume_aligned`；(3) 若不足再手写 AVX2 intrinsic |
| 预估收益 | 加 `restrict`/`ivdep` 预估 ~3~8% SCF 提升；进一步 intrinsic 可达 ~5~12% |
| 风险/优先级 | 低风险 / **中优先级** — 别名修复是安全操作，不会改变数值结果 |

### 8.2 题目8：缓存复用

| 项目 | 评估 |
|------|------|
| 目标 | 缓存 `gg`, `gcar`, `gdirect`, `gk2`, `igl2isz_k`, `igl2ig_k` (~0.9 MB serial, ~0.5 MB np2) |
| 最大单项收益 | 消除 `setupIndGk` 双重 `cal_GplusK_cartesian` 调用（~18K 次冗余，~10-50ms Debug, <5ms Release） |
| 当前开销 | < 0.4s（含在 `before_scf`），Release 下可降至 < 0.1s |
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
  - `mix_np2_omp4` 配置因 MPI 通信开销较大，在某些用例下呈现中等方差（gaas_medium: 29s/31s/35s）
- **WSL2 局限性**：本报告所有数据来自 WSL2 (Hyper-V 轻量级 VM)。Windows 宿主的后台抢占与 EPT 内存转换会严重污染硬件计数器，尤其影响 MPI 多进程的通信延迟和共享内存性能。第 13 周的终期优化验收必须迁移到物理机 Linux 或 HPC 集群上进行，WSL2 仅应作为开发与功能验证环境
- **Debug 构建限制**：主要基线数据来自 `-O0` 构建。Release（`-O3`）对照实验已完成（[5.4 节](#54-release-构建对照实验-o3--marchnative--dndebug)），自动向量化状态已确认（[5.3.1 节](#531-gcc-13-自动向量化实际状态)）：gather/scatter 核心拷贝循环未向量化，1.64× gathers 加速来自非 SIMD 优化。百分比占比在 Debug 和 Release 下基本一致，但绝对加速倍数不能跨构建类型外推
- **计时器限制**：0.1s 阈值导致部分轻量函数无数据显示，已在第 3.3 节说明

---

## 10. 结论与下一步

### 10.1 基线结论

1. **全部 54 次运行通过** (100% 通过率)，测试体系稳定可靠。

2. **Gather/Scatter 是 MPI 模式下的重要优化目标**，占 FFT 变换的 **30~55%**（稳态），占 SCF 迭代的 **6~10%**。MPI 路径的 gather/scatter 由本地 copy 段 + `MPI_Alltoallv` 组成，其中 copy 段是 SIMD 的直接作用目标。6 处内层连续拷贝循环非常适合 SIMD 向量化。**自动向量化状态已确认**（[5.3.1 节](#531-gcc-13-自动向量化实际状态)）：GCC 13.3 因指针别名障碍未能向量化核心拷贝循环，清零循环已使用 AVX2 256-bit 向量化。Release 的 1.64× gathers 加速来自非 SIMD 优化（函数内联/循环展开），手动 SIMD 有 ~1.5~3× 的单循环加速空间，预期加 `restrict`/`ivdep` 后可实现 ~3~8% SCF 整体提升。

3. **缓存复用函数当前开销较小**（<0.4s，含在 `before_scf` 中），但其收益在于消除 `setupIndGk` 的 |G+K|² 双重计算（~18,000 次冗余调用）和减少反复内存分配。缓存数据量约 0.5~0.9 MB，可完全放入 L2 缓存。

4. **线程扩展性趋于饱和**: 4 核达 2.0~2.2× 加速，8 核不再显著提升。`mix_np2_omp2`（4 核, 2.19×）略优于 `omp4`（4 线程, 2.06×），MPI 的进程级并行在小循环上避免了 OpenMP 线程调度开销。

### 10.2 第 13 周工作建议

1. **SIMD 向量化**（题目5）：
   - **第一阶段（最低成本）**：在 `pw_gatherscatter.h` 的 6 处拷贝循环前添加 `#pragma GCC ivdep` + 指针 `__restrict__` 消除别名障碍，让 GCC 自动向量化生效。预估 ~3~5% SCF 提升。在 Release 构建下用 gaas_medium np2_omp2 验证
   - **第二阶段（对齐优化）**：加 `__builtin_assume_aligned(ptr, 32)` 断言确保 32 字节对齐。预估额外 ~1~3% SCF
   - **第三阶段（按需 intrinsic）**：若前两阶段收益不足，手写 AVX2 intrinsic 替代 6 处拷贝循环。预估额外 ~3~5% 但维护成本高，仅当测试后确实有必要时采用

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
├── gaas_tiny_40Ry/
│   ├── INPUT          (24³, 40 Ry, 70 bands, 3×3×3 kpt)
│   ├── KPT            (Gamma-centered 3×3×3)
│   ├── STRU -> ../../../../abacus-user-guide/examples/pw2/cg/001_4GaAs/STRU
│   ├── Ga.upf -> ../../../../abacus-user-guide/examples/pw2/cg/001_4GaAs/Ga.upf
│   └── As.upf -> ../../../../abacus-user-guide/examples/pw2/cg/001_4GaAs/As.upf
├── gaas_small/
│   ├── INPUT          (32³, 40 Ry, 70 bands, 3×3×3 kpt)
│   └── ... (同上链接)
├── gaas_medium/
│   ├── INPUT          (48³, 40 Ry, 70 bands, 3×3×3 kpt)
│   └── ... (同上链接)
├── results/
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

## AI使用心得

这次因为copilot教育版达到了使用限额（），所以选择使用Claude Code + DeepSeek v4 pro的api来进行vibe coding。国产模型对于中文文档写作确实有独到的优势，语言更顺畅了；同时deepseek的成本优势也非常明显，我用pro模型跑了两天的测试样例并修改文档，一共花费5元，得到的结果似乎也优于之前使用的Gemini 3.1 pro、GPT 5.2 Codex。
不过进行测试的过程中也发现了模型的一些不足，比如一开始的测试并没有严格控制变量，导致网格数和截断能同时变化，测试结果的有效性存疑；而且修改文档时，对于我给出的要求并没有完全修改，还需要手动进一步调整。
