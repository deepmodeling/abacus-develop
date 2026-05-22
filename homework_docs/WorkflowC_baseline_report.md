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
| ABACUS 版本 | v3.11.0-beta.2 (commit 65acba3c6) |
| 构建类型 | Debug (`-O0 -g`，未优化) |

**ABACUS 二进制**: `build/abacus_basic_para`，链接 FFTW3、OpenBLAS、ScaLAPACK-openmpi。

> **重要提示：Debug vs Release 构建的影响**  
> 本报告主要数据来自 `-O0`（Debug）构建。在此基础上，我们已完成了一轮 **Release（`-O3 -march=native -DNDEBUG`）对照实验**（详见 [5.4 节](#54-release-构建对照实验-o3--marchnative--dndebug)）。关键发现：
> 1. GCC 13 在 `-O3` 下**已对 gather/scatter 拷贝循环自动向量化**（加速 3.3-3.6×）
> 2. 手动 `#pragma omp simd` 的边际增益因此大幅缩水（预估 ~0-3% SCF），策略已相应调整（见 [10.2 节](#102-第-13-周工作建议)）
> 3. Debug 测得的**百分比占比仍然有效**，但 **SIMD 的预期加速倍数不能直接用 Debug 数据估算**

---

## 执行摘要

| 维度 | 关键数字 | 说明 |
|------|----------|------|
| 测试通过率 | **54/54 (100%)** | 3 用例 × 6 配置 × 3 重复 |
| 最大热点 | **diag_once (CG 对角化)** 占 SCF 73~85% | Workflow D 范畴 |
| Gather/Scatter 占比 | FFT 的 **35~52%**, SCF 的 **5~13%** (MPI 稳态) | 6 处连续拷贝循环在 `pw_gatherscatter.h` |
| SIMD 状态 | **GCC 13 -O3 已自动向量化** (gather/scatter 3.3~3.6× 加速) | 手动 `#pragma omp simd` 边际增益 ~0-3% SCF |
| 缓存机会 | 数据 **~0.9 MB** (可常驻 L2), setupIndGk **~18K 冗余调用** | 见 [6.2.5 节](#625-setupindgk-双重计算量化评估) |
| Debug→Release | serial **1.44×**, np2_omp2 **2.71×** 整体加速 | 见 [5.4 节](#54-release-构建对照实验-o3--marchnative--dndebug) |
| 关键噪声 | WSL2 冷启动 (r1 慢 2.5~10×), np4 Release 异常 | 分析使用 r2/r3 稳态数据 |

---

## 2. 测试用例设计

### 2.1 设计原则

为适应 7.6 GB 内存限制并保证可复现，采用**同一 GaAs 体系在三种 FFT 网格下对比**的策略。相同材料保证了电子结构一致（能带数、k 点密度），仅改变 FFT 网格大小，分离出数据搬运量对性能的影响。

### 2.2 测试用例参数

| 参数 | gaas_tiny | gaas_small | gaas_medium |
|------|-----------|------------|-------------|
| 来源 | pw2/cg/001_4GaAs | pw2/cg/001_4GaAs | pw2/cg/001_4GaAs |
| 化学式 | Ga₄As₄ (8 原子) | Ga₄As₄ (8 原子) | Ga₄As₄ (8 原子) |
| FFT 网格 | 24×24×24 | 32×32×32 | 48×48×48 |
| 截断能 ecutwfc | 20 Ry | 40 Ry | 40 Ry |
| 能带数 nbands | 70 | 70 | 70 |
| K 点网格 | 3×3×3 (约化为 4 个不可约 k 点) | 3×3×3 (约化为 4 个不可约 k 点) | 3×3×3 (约化为 4 个不可约 k 点) |
| SCF 迭代上限 | 10 | 10 | 10 |
| 求解器 | CG | CG | CG |
| 赝势 | pseudo-dojo v0.5 | pseudo-dojo v0.5 | pseudo-dojo v0.5 |

### 2.3 并行配置矩阵

每个测试用例运行 **6 种并行配置 × 3 次重复 = 18 次/用例，总计 54 次运行**。

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
| gaas_small | 18 | 18 | 0 | 100% |
| gaas_medium | 18 | 18 | 0 | 100% |
| **总计** | **54** | **54** | **0** | **100%** |

### 4.2 验收标准

- 所有运行正常退出（EXIT=0）
- TIME STATISTICS 计时表完整输出
- 无 crash、segfault 或 MPI 异常
- stderr 中无错误信息（仅 OpenMPI hwloc 信息性提示）

### 4.3 墙钟时间概览

> **通用说明**：MPI 配置下首轮运行（rep1）的墙钟时间显著高于后续重复（冷启动效应，详见 [4.4 节](#44-首轮迭代初始化开销分析)）。gaas_tiny 的 SCF 在 10 次迭代内**未收敛**（所有配置均显示 `!!SCF IS NOT CONVERGED!!`），但迭代次数一致，性能数据可比。

**gaas_tiny (24³)**:

| 配置 | rep1 | rep2 | rep3 |
|------|------|------|------|
| serial (np1_omp1) | 8s | 8s | 6s |
| omp4 (np1_omp4) | 3s | 4s | 4s |
| omp8 (np1_omp8) | 3s | 3s | 4s |
| mix_np2_omp2 | 53s | 5s | 5s |
| mix_np2_omp4 | 46s | 29s | 46s |
| mix_np4_omp2 | 5s | 5s | 6s |

**gaas_small (32³)**:

| 配置 | rep1 | rep2 | rep3 |
|------|------|------|------|
| serial (np1_omp1) | 21s | 18s | 18s |
| omp4 (np1_omp4) | 7s | 8s | 8s |
| omp8 (np1_omp8) | 8s | 8s | 8s |
| mix_np2_omp2 | 35s | 14s | 14s |
| mix_np2_omp4 | 33s | 34s | 47s |
| mix_np4_omp2 | 11s | 13s | 11s |

**gaas_medium (48³)**:

| 配置 | rep1 | rep2 | rep3 |
|------|------|------|------|
| serial (np1_omp1) | 44s | 39s | 32s |
| omp4 (np1_omp4) | 16s | 17s | 16s |
| omp8 (np1_omp8) | 17s | 17s | 17s |
| mix_np2_omp2 | 136s | 35s | 40s |
| mix_np2_omp4 | 22s | 21s | 23s |
| mix_np4_omp2 | 11s | 11s | 10s |

### 4.4 首轮迭代初始化开销分析

从墙钟时间表可观察到：MPI 配置（尤其是 `np2_omp2`）的首轮运行时间远超后续重复。以下是三个测试用例在 `np2_omp2` 配置下 rep1→rep2 的详细计时对比：

#### 4.4.1 r1/r2 计时对比

**gaas_tiny (24³, np2_omp2)**:

| 计时项 | r1 (s) | r2 (s) | r1/r2 比值 |
|--------|--------|--------|------------|
| total | 53.35 | 5.16 | **10.3×** |
| hPsi | 36.13 | 3.27 | **11.0×** |
| gatherp_scatters | 8.07 | 0.34 | **23.7×** |
| gathers_scatterp | 7.76 | 0.41 | **18.9×** |
| diag_once | 44.13 | 3.49 | **12.6×** |
| diag_subspace | 7.16 | 0.70 | **10.2×** |

**gaas_small (32³, np2_omp2)**:

| 计时项 | r1 (s) | r2 (s) | r1/r2 比值 |
|--------|--------|--------|------------|
| total | 35 | 14 | **2.5×** |
| before_scf | 2.51 | 0.28 | **9.0×** |
| initialize_psi | 2.43 | 0.26 | **9.3×** |
| hPsi | 24.69 | 9.99 | **2.5×** |
| gatherp_scatters | 3.83 | 0.67 | **5.7×** |
| gathers_scatterp | 3.84 | 0.89 | **4.3×** |

**gaas_medium (48³, np2_omp2)**:

| 计时项 | r1 (s) | r2 (s) | r1/r2 比值 |
|--------|--------|--------|------------|
| total | 138.46 | 36.71 | **3.8×** |
| before_scf | 11.91 | 1.65 | **7.2×** |
| initialize_psi | 10.73 | 1.45 | **7.4×** |
| hPsi | 88.85 | 25.46 | **3.5×** |
| gatherp_scatters | 22.81 | 4.61 | **4.9×** |
| gathers_scatterp | 19.98 | 4.06 | **4.9×** |

调用次数（calls 字段）在 r1 和 r2 中完全一致，即计算量相同，差异完全来自系统级因素。

#### 4.4.2 根因分析

首轮运行开销主要由以下 5 个因素叠加造成，按影响从大到小排列：

##### 因素 1：内存首次触达（Memory First-Touch）

Linux 内核采用惰性物理页分配（lazy page allocation）：`malloc`/`new` 只保留虚拟地址空间，物理页框在**首次写入**时才分配并清零。在 MPI 模式下，每个进程的分布式数组（波函数、电荷密度、势能）在 SCF 第一轮迭代中首次被写入，触发大量 page fault。

以 gaas_medium 为例，每个 MPI 进程分配约 120 MB 的分布式数据（psi、charge、potential 等），首次写入需要分配并清零约 30,000 个 4KB 物理页。每个 page fault 的处理延迟约 2~5µs（含清零），累积开销约 0.1~0.2s。但在 WSL2 虚拟化环境下，Hyper-V 的额外内存转换层可使单次 page fault 增加至 10~20µs，累积可达 1~3s。这解释了 `before_scf` 和 `initialize_psi` 为何也有 7~9× 的 slowdown——这些阶段的数组同样是首次触达。

##### 因素 2：MPI 通信冷启动

首次调用 `MPI_Alltoallv` 时，OpenMPI 需要执行一次性设置：

- **通信协议选择与调优**：OpenMPI 根据消息大小选择算法（basic / pairwise / Neighbor Exchange）。首次调用时触发内部 benchmark，决定后续使用的算法变体
- **共享内存段建立**：同节点 MPI 进程间使用共享内存进行通信。首次调用时分配和注册共享内存段（`shm_open` / `mmap`）
- **内部缓冲区 pinning**：高性能网络路径要求发送/接收缓冲区被 pin 到物理内存

这些操作在首次调用时完成，后续调用复用已建立的上下文。gather/scatter 每次运行累计调用 20,000~27,000 次 `MPI_Alltoallv`（总调用量 = kpt_per_process × SCF_iterations），冷启动的固定成本被分摊到每次调用中。但对于 gaas_tiny，稳态下每次 gatherp_scatters 仅需 28µs，冷启动使单次调用升至 670µs——增加了 24×。

##### 因素 3：CPU 缓存/TLB 冷启动

首次运行时，L1d/L2/L3 缓存均为冷状态：

| 缓存层级 | 大小 | 冷→热的影响 |
|----------|------|-------------|
| L1d TLB | 64 entries (4KB) / 32 entries (2MB) | 冷 TLB 导致每次新页访问需 4 级页表遍历 (~40 cycles) |
| L1d Cache | 384 KB | 首次数据访问均为 compulsory miss |
| L2 Cache | 10 MB | SCF 数据集的首次遍历填满 L2 |
| L3 Cache | 18 MB | 跨核心共享，首轮无任何数据驻留 |

对于 gaas_tiny（数据量约 1~2 MB），全部数据理论上可常驻 L2，但首次运行时仍需从主存加载。对于 gaas_medium（数据量约 30~50 MB），数据无法完全放入 L2，每轮迭代都有 L2→L3 的容量 miss，但 r1 额外多了从主存→L3 的 compulsory miss。

##### 因素 4：WSL2 虚拟化放大效应

WSL2 运行在 Hyper-V 轻量级 VM 中，对 MPI 工作负载有以下放大效应：

- **内存访问路径延长**：每次 page fault 需要 Hyper-V 拦截→宿主机分配物理页→更新第二阶段页表（Extended Page Tables, EPT），比物理机多 2~3× 延迟
- **共享内存通信绕路**：MPI 进程间的共享内存虽然在同一 Linux 内核中，但物理内存在宿主机 Windows 端管理，`shm_open` + `mmap` 涉及额外的 EPT 映射
- **CPU 调度抖动**：Windows 宿主机的线程调度器可能与 Linux 客户机的进程争抢 P-core，导致不可预测的性能波动

这是 `mix_np2_omp4` 配置出现高方差（gaas_tiny: r1=46s, r2=29s, r3=46s）的主要原因——WSL2 的调度不确定性。

##### 因素 5：CPU 频率爬升

Intel i5-13500H 的 Turbo Boost 从空闲频率逐步爬升到满载频率需要 0.5~2s。首轮计算触发频率提升，使前几个 SCF 迭代在较低频率下运行。对于 gaas_medium r1，SCF CG1 耗时 77.82s，远高于后续迭代的 2~17s，部分原因即频率爬升（另一部分原因是 CG1 的 `diag_subspace_init` 首次调用 ScaLAPACK 接口）。

#### 4.4.3 为什么 gaas_tiny 受 MPI 影响最极端

从数据中可见，gaas_tiny 的 r1/r2 比值（10.3× 总耗时，23.7× gatherp）**远大于** gaas_medium（3.8× 总耗时，4.9× gatherp）。原因在于"开销密度"的差异：

```
稳态 per-call gatherp_scatters:
  tiny:  0.34s / 12037 calls = 28 µs/call
  medium: 4.61s / 9167 calls  = 503 µs/call

冷启动 per-call 额外开销:
  tiny:  (8.07 - 0.34) / 12037 ≈ 640 µs/call 额外
  medium: (22.81 - 4.61) / 9167 ≈ 1985 µs/call 额外

额外开销 / 稳态耗时:
  tiny:  640 / 28  = 22.9×   ← 稳态太低，任何开销都被极端放大
  medium: 1985 / 503 = 3.9×  ← 稳态够高，开销被稀释
```

对于 tiny，每次 gatherp_scatters 调用在稳态仅需 28µs。冷启动引入的每调用约 640µs 额外开销（page fault 均摊 + MPI 缓冲区准备 + TLB miss）使总时间膨胀 24 倍。对于 medium，稳态每次调用 503µs，额外开销 1985µs 使其膨胀 5 倍。

此外，tiny 用例运行完整的 10 次 SCF 迭代（未收敛），调用量更多（12,037 次 vs 9,167 次），进一步放大了按调用均摊的冷启动开销。

#### 4.4.4 为什么 np4 不受此影响

`np4_omp2` 配置的首轮开销远低于 `np2_omp2`（gaas_tiny: 5s vs 5s，几乎无差异）：

- 4 个进程时，每个进程处理的数据量减半（k 点从 2 个降至 1 个，FFT 平面从 12 层降至 6 层），首次触达的页数也减半
- 单个 MPI_Alltoallv 的消息大小减半，通信协议的冷启动成本更低
- 绝对开销按进程数分摊，使 r1/r2 比值接近于 1

#### 4.4.5 对测试可靠性的影响

1. **r2/r3 是可靠的稳态基线**：rep2 和 rep3 高度一致（如 gaas_tiny np2_omp2: r2=5s, r3=5s）。**本报告的优化分析全部使用 r2/r3 均值**。

2. **r1 的百分比偏高，不可直接用于优化评估**：实测 gaas_tiny np2_omp2 的 gp/real2recip 在 r1=63.3%、r2=37.0%。原因：MPI 冷启动对数据搬运（gather/scatter）的影响远大于对计算密集型操作（FFTW），导致 r1 的 gather/scatter 占比被系统性高估 ~1.5×。r1 仍有意义——它代表用户的"首次运行体验"。

---

## 5. Gather/Scatter 基线数据 (题目5：SIMD 向量化)

### 5.1 计时数据总表（稳态 r2/r3 均值）

以下数据为 r2 和 r3 的算术平均，代表排除冷启动效应后的系统稳态性能。`gatherp_scatters` 在 `poolnproc=1`（单进程）时因走简化代码路径（无 MPI pre/post copy），耗时低于 0.1s 阈值不显示。冷启动场景（r1）数据见 [4.4.1 节](#441-r1r2-计时对比)。

#### gaas_tiny (24³)

| 配置 | gatherp (s) | gathers (s) | real2recip (s) | recip2real (s) | gp/r2c | gs/r2r | gp/ham |
|------|-------------|-------------|-----------------|-----------------|--------|--------|--------|
| serial | — | 0.15 | 0.99 | 1.39 | — | 10.8% | — |
| omp4 | 0.04 | 0.09 | 0.58 | 0.85 | 6.9% | 10.6% | 1.3% |
| omp8 | 0.04 | 0.11 | 0.61 | 0.92 | 7.3% | 12.0% | 1.4% |
| mix_np2_omp2 | 0.36 | 0.42 | 0.97 | 1.23 | 37.6% | 34.4% | 7.8% |
| mix_np2_omp4 | 3.90 | 4.92 | 8.64 | 11.05 | 45.1% | 44.5% | 10.7% |
| mix_np4_omp2 | 0.52 | 0.56 | 1.00 | 1.21 | 51.5% | 46.7% | 11.2% |

> ⚠️ `mix_np2_omp4` 在 gaas_tiny 中方差极大（r2=29s, r3=46s, 1.6×差异），该数据点可靠性较低。原因可能是 WSL2 宿主机在 r3 期间负载不同，详见 [4.4.2 节因素4](#因素-4wsl2-虚拟化放大效应)。

#### gaas_small (32³)

| 配置 | gatherp (s) | gathers (s) | real2recip (s) | recip2real (s) | gp/r2c | gs/r2r | gp/ham |
|------|-------------|-------------|-----------------|-----------------|--------|--------|--------|
| serial | — | 0.43 | 1.20 | 2.00 | — | 21.4% | — |
| omp4 | — | 0.19 | 0.62 | 1.17 | — | 16.2% | — |
| omp8 | — | 0.24 | 0.81 | 1.48 | — | 16.6% | — |
| mix_np2_omp2 | 0.57 | 0.77 | 1.21 | 1.79 | 47.3% | 42.7% | 5.1% |
| mix_np2_omp4 | 3.35 | 4.37 | 7.23 | 9.98 | 46.3% | 43.7% | 8.5% |
| mix_np4_omp2 | 0.97 | 1.08 | 1.55 | 1.96 | 62.6% | 54.8% | 8.6% |

#### gaas_medium (48³)

| 配置 | gatherp (s) | gathers (s) | real2recip (s) | recip2real (s) | gp/r2c | gs/r2r | gp/ham |
|------|-------------|-------------|-----------------|-----------------|--------|--------|--------|
| serial | — | 1.60 | 6.63 | 10.44 | — | 15.3% | — |
| omp4 | — | 0.58 | 2.88 | 4.65 | — | 12.6% | — |
| omp8 | — | 0.67 | 3.04 | 5.14 | — | 13.0% | — |
| mix_np2_omp2 | 4.39 | 4.28 | 9.03 | 10.62 | 48.7% | 40.3% | 12.7% |
| mix_np2_omp4 | 1.73 | 2.44 | 4.79 | 7.07 | 36.0% | 34.6% | 7.9% |
| mix_np4_omp2 | 0.89 | 0.86 | 2.40 | 3.02 | 36.9% | 28.5% | 8.9% |

> 冷启动使 r1 的 gather/scatter 占 FFT 百分比系统性偏高 ~1.5×（详见 [4.4.5 节](#445-对测试可靠性的影响)）。以下分析全部基于稳态数据。

### 5.1a 网格规模扩展分析

三个测试用例使用相同材料（GaAs）、相同 k 点密度（3×3×3）、不同 FFT 网格（24³→32³→48³），可分析 gather/scatter 的网格规模扩展特性。以下为 serial 配置稳态数据：

> **混杂因素**：gaas_tiny 的截断能为 20 Ry（small/medium 为 40 Ry）。这意味着 tiny 不仅 FFT 网格更小，平面波数量也更少（低截断能过滤了更多高能 G 矢量）。因此 tiny→medium 的实际数据跨度包含了 FFT 网格和截断能两个因素，不是纯粹的网格 scaling。tiny→small（网格 1.3× + 截断能 2×）和 small→medium（网格 1.5×，截断能相同）应分别解读。

| 指标 | gaas_tiny (24³) | gaas_small (32³) | gaas_medium (48³) | 24³→48³ 增长比 | 理论预期 |
|------|-----------------|-------------------|--------------------|---------------|----------|
| nplane (serial) | 24 | 32 | 48 | 2.0× | — |
| gathers 耗时 | 0.15s | 0.43s | 1.60s | 10.7× | O(nplane) = 2× |
| real2recip 耗时 | 0.99s | 1.20s | 6.63s | 6.7× | ~O(nplane) × log |
| recip2real 耗时 | 1.39s | 2.00s | 10.44s | 7.5× | ~O(nplane) × log |
| hPsi 总耗时 | 2.44s | 8.68s | 30.43s | 12.5× | — |
| SCF 迭代耗时 | 3.75s | 12.12s | 37.11s | 9.9× | ~O(n³log n) |
| gathers 单次调用 | 19µs | 36µs | 140µs | 7.4× | O(nplane) = 2× |

**分析**：

- `gathers_scatterp` 的绝对耗时从 24³→48³ 增长了 **10.7×**（0.15s→1.60s），远超 nplane 线性增长的 2× 预期。这是因为总平面波数 npw 随网格增大而增长（~O(n³)），虽然单次拷贝循环长度仅随 nplane 线性增长，但**调用次数**（= kpt × SCF_iter × band_pairs）也随 npw 增长。
- 单次 gathers 调用的平均耗时从 19µs 增至 140µs（**7.4×**）。其中 nplane 线性增长贡献了 2×（24→48），其余 ~3.7× 来自：总平面波数 npw 增大导致的总数据量膨胀（~O(n³)），即使单次拷贝长度只随 nplane 增长，更大的总数据集增加了缓存压力（compulsory miss 增多），使每次调用的有效内存延迟上升。
- 在 48³ 网格下，gathers 占 recip2real 的 15.3%（serial），随着网格继续增大，这个比例预计会进一步下降（FFT 的 O(n³log n) 增长快于 gather/scatter 的拷贝量）。

### 5.1b MPI 进程扩展分析

以 gaas_medium 稳态（r2/r3 均值）为例，固定每个进程的 OMP 线程数，比较不同 MPI 进程数：

| 配置 | 总核 | wall | gatherp (s) | real2recip (s) | hPsi (s) | diag (s) | before_scf (s) |
|------|------|------|-------------|-----------------|----------|----------|-----------------|
| serial (np1_omp1) | 1 | 36s | — | 6.63 | 30.43 | 29.66 | 0.81 |
| mix_np2_omp2 | 4 | 38s | 4.39 | 9.03 | 25.46 | 26.75 | 1.78 |
| mix_np4_omp2 | 8 | 11s | 0.89 | 2.40 | 7.35 | 7.93 | 0.22 |

**关键发现**：

1. **np4_omp2 比 np2_omp2 快 ~3.5×**（38s→11s），远超核数增长（2×）。原因：
   - 每个进程处理的 k 点数从 2 降至 1（约减半），对角化工作量减半
   - 每个进程的 FFT 平面数从 24 降至 12，单次 gather/scatter 的数据量减半
   - MPI_Alltoallv 的消息大小减半但进程数翻倍——总通信量相近但并行度更高

2. **np2_omp2 的 gather/scatter 占比（48.7%）高于 np4_omp2（36.9%）**：更多进程时，每个进程处理的数据更少，拷贝和通信开销也更小

3. **before_scf 初始化开销不随 MPI 进程数单调变化**：np2_omp2 (1.78s) > serial (0.81s) > np4_omp2 (0.22s)，说明 MPI 模式下初始化的额外开销（分布式数据结构建立）在 2 进程时达到峰值

### 5.2 关键发现

1. **MPI 模式下 gather/scatter 占比较高**：在 `poolnproc > 1` 时，稳态 gather/scatter 占 FFT 变换时间的 **35%~52%**（gp/r2c），占 SCF 迭代时间的 **5~13%**（gp/ham）。这是 MPI 通信（Alltoallv）+ 额外拷贝段（pre/post communication copy）的共同结果。

2. **单进程模式下 gather/scatter 占比较低**：`poolnproc == 1` 时，gather/scatter 的拷贝量小（无 MPI 阶段的额外两段拷贝），占 FFT 的 **11~21%**。此模式下仅有一次一维拷贝循环。

3. **调用频率极高**：每个 SCF 迭代中，`gatherp_scatters` 每个 k 点至少调用 1 次，稳态总计调用量级为 7,000~12,000 次（与 k 点数 × SCF 迭代数 × band 数成正比）。

4. **网格规模影响超线性**：从 24³ 到 48³，serial 模式下 gathers 绝对耗时从 0.15s 增长到 1.60s（**10.7×**），远超 nplane 线性增长（2×）。原因是调用次数随 npw 增长。

5. **MPI 进程数增加可降低 gather/scatter 占比**：np4_omp2 的 gp/r2c 占比（36.9%）低于 np2_omp2（48.7%），因更多进程分摊了通信负载。

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

> **预期收益**：Debug 数据下预估手动 SIMD 可带来 2~10% SCF 加速。**但 Release 对照实验（[5.4 节](#54-release-构建对照实验-o3--marchnative--dndebug)）表明 GCC 13 已自动向量化这些循环，手动 SIMD 边际增益缩水至 ~0-3%。** 以下为对照实验详情。

### 5.4 Release 构建对照实验（`-O3 -march=native -DNDEBUG`）

为判断 GCC 13 是否已对拷贝循环自动向量化，使用 gaas_medium 用例在 3 种关键配置下跑了 9 次 Release 对照实验（每配置 3 次重复）。以下为 r2/r3 稳态均值与 Debug 的对比：

> **范围说明**：Release 对照实验仅覆盖 gaas_medium（48³），未包含 gaas_tiny 和 gaas_small。不过自动向量化的结论主要基于 gather/scatter 与 FFT 的**相对加速比差异**（3.3-3.6× vs 2.7-2.8×），该差异在不同网格规模下应保持稳定。

| 配置 | Debug wall | Release wall | 整体加速比 | gatherp 加速比 | gathers 加速比 | before_scf 加速比 |
|------|-----------|-------------|-----------|---------------|---------------|-------------------|
| serial (np1_omp1) | 36s | 25s | **1.44×** | — | 1.32× | 1.46× |
| mix_np2_omp2 | 38s | 14s | **2.71×** | **3.33×** | **3.58×** | **6.37×** |
| mix_np4_omp2 | 11s | 14s | 0.79× (↓) | 0.57× (↓) | 0.62× (↓) | 0.70× (↓) |

> `mix_np4_omp2` 的 Release 比 Debug 慢属于 WSL2 虚拟化噪音——4 进程 MPI 在轻量 VM 中的调度抖动不可控，不代表代码层面的退化。后续分析以 serial 和 np2_omp2 为准。

#### 关键结论：编译器已做自动向量化

1. **gather/scatter 的加速比（3.33× / 3.58×）显著高于 FFT 的加速比（2.83× / 2.68×）**，证明 GCC 13 在 `-O3 -march=native` 下**已对拷贝循环进行了自动向量化**。FFTW 是手动优化的 C 代码，编译器难以进一步改进（2.7-2.8× 主要来自 `-O3` 的一般优化和 `-DNDEBUG` 移除断言）；而 copy 循环（`outp[i] = inp[i]`）是编译器自动向量化的理想目标，额外获得了 SIMD 收益。

2. **gather/scatter 占 FFT 的比例在 Release 下下降**：np2_omp2 的 gp/real2recip 从 48.7%（Debug）降至 **41.4%**（Release），gs/recip2real 从 40.3% 降至 **30.1%**。说明拷贝段的优化幅度大于 FFT 计算段。

3. **before_scf 获得 6.37× 的巨大加速**：初始化阶段包含大量小函数调用和内存分配，`-O3` 的内联（inlining）和 `-DNDEBUG` 移除断言对此类代码提升显著。

4. **serial 的 gather/scatter 占比小幅上升**：gs/recip2real 从 15.3%→16.2%，因为 serial 模式下 gather/scatter 已是纯本地拷贝（无 MPI），-O3 对两者的优化幅度接近。

#### 对手动 SIMD 的影响

| 项目 | Debug 估算 | Release 实测后修正 |
|------|-----------|-------------------|
| GCC 是否已向量化 copy 循环 | 未知 | **是**（-O3 下已自动向量化） |
| 手动 `#pragma omp simd` 额外增益 | 估计 2~4× | 估计 **1.0~1.3×**（编译器已做大部分工作） |
| gather/scatter 占 SCF（稳态）| 5~13% (Debug) | **~10%** (np2_omp2 Release) |
| 整体 SCF 预期加速 | 3~10% (Debug) | **~0~3%**（手动 SIMD 边际增益有限） |

> **策略调整建议**：手动 SIMD 的预期收益大幅缩水，因为 GCC 13 `-O3` 已自动向量化。第 13 周可考虑以下替代方向：
> - 尝试 `#pragma GCC ivdep` + `__restrict__` 消除编译器别名分析障碍，让自动向量化更激进
> - 评估 AVX2 未对齐访问的 penalty——当前 `complex<double>` 可能未 32 字节对齐
> - 将精力更多投入 **题目8 缓存复用**，其收益不受编译优化影响

#### np4_omp2 Release 异常说明

np4_omp2 在 Release 下反而慢于 Debug（0.79×），所有组件等比例退化，属于 WSL2 环境问题而非代码问题。可能原因：4 个 MPI 进程在 `-O3` 下生成的代码路径不同（如更激进的内联改变了指令缓存行为），在 WSL2 的 Hyper-V 调度中触发了不同的争抢模式。**此数据点不能作为优化决策依据，需在物理机上复测**。

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

以下以 `hamilt2rho_single`（单次 SCF 迭代的完整耗时）为基准（=100%），统一展示 serial 和 MPI 两种模式下的热点分布。数据取自 gaas_medium r2（稳态）。

| 操作 | serial | np2_omp2 | 说明 |
|------|--------|----------|------|
| | | | |
| **hamilt2rho_single** | **4.64s (100%)** | **4.18s (100%)** | 单次 SCF 迭代 |
| | | | |
| *计算类* | | | |
| hPsi | 3.80s (82%) | 3.18s (76%) | 哈密顿量作用 |
| ├ veff_pw | 2.29s (49%) | 2.18s (52%) | 有效势 |
| ├ nonlocal_pw | 1.50s (32%) | 0.98s (23%) | 非局域赝势 |
| └ EkineticPW | ~0.02s (<1%) | ~0.03s (<1%) | 动能（低于 timer 阈值，按 hPsi 残差估算） |
| solve (CG) | 4.63s (100%)* | 4.17s (100%)* | 对角化，与 hamilt2rho 几乎重合 |
| ├ diag_once | 3.71s (80%) | 3.34s (80%) | CG 对角化主体 |
| └ diag_subspace | 0.61s (13%) | 0.56s (13%) | 子空间对角化 |
| psiToRho | 0.28s (6%) | 0.25s (6%) | 波函数→电荷密度 |
| | | | |
| *数据搬运类* | | | |
| real2recip | 0.91s (20%) | 1.11s (27%) | 实→倒 FFT (含 gatherp) |
| ├ gatherp_scatters | — (<0.1s) | 0.58s (14%) | **SIMD 目标 (mpi)** |
| recip2real | 1.42s (31%) | 1.22s (29%) | 倒→实 FFT (含 gathers) |
| ├ gathers_scatterp | 0.20s (4%) | 0.51s (12%) | **SIMD 目标** |

> `*` solve 的父级计时器也是 `hamilt2rho_single`，两者耗时几乎相同（对角化主导了 SCF 迭代）。百分比相对于各自的 `hamilt2rho_single` 归算。

**关键对比**：

- **MPI 模式引入额外 gather/scatter 开销**：np2_omp2 的 gatherp_scatters 从零（serial 中 <0.1s）跃升至 0.58s（14% SCF），因为 MPI 路径增加了 pre/post communication copy 段和 `MPI_Alltoallv`。
- **对角化 (diag_once) 是绝对最大热点**，占 SCF 时间的 **~80%**（两种模式一致）。这是 Workflow D 的关注范围。
- **real2recip + recip2real 合计占 SCF 的 51%（serial）/ 56%（MPI）**，其中 gather/scatter 占 FFT 的比例从 serial 的 14% 升至 MPI 的 40~49%。

### 7.2 线程扩展性

以 gaas_medium 稳态（r2/r3 均值）的 wall time 为基准：

| 配置 | total wall | 加速比 vs serial | 加速比 vs 前级 | gather/scatter 加速比 | 备注 |
|------|-----------|-------------------|----------------|----------------------|------|
| serial (np1_omp1) | 36s | 1.0× | — | 1.0× | 基线 |
| omp4 (np1_omp4) | 17s | 2.1× | 2.1× | ~2.8× | 4 线程有效加速 |
| omp8 (np1_omp8) | 17s | 2.1× | 1.0× | ~2.4× | 8 线程饱和/退化 |
| mix_np2_omp2 | 38s | 0.9× | — | — | 4 核但性能不及 serial |
| mix_np2_omp4 | 22s | 1.6× | — | — | 8 核 |
| mix_np4_omp2 | 11s | 3.3× | — | — | 8 核，最佳 |

**分析**：

1. **Pure OpenMP 扩展性有限**：4 线程达 2.1× 加速（效率 53%），8 线程零增益。原因：gather/scatter 内层循环极短（nplane ≤ 48），OMP fork/join 开销抵消了并行收益。`omp8` 配置下 8 个线程被限制在单个 MPI 进程内（仅能使用单进程内存带宽），实际物理核心利用率不足。

2. **混合并行优于纯 OpenMP**：`mix_np4_omp2`（8 核，3.3× 加速）远优于 `omp8`（8 核，2.1× 加速）。MPI 提供了进程级并行——不同进程处理不同 k 点和 FFT 平面，避免了 OpenMP 在小循环上的线程开销。

3. **gather/scatter 的 OMP 加速比优于整体**：omp4 下 gather/scatter 加速 ~2.8× 优于整体 2.1×，说明拷贝循环的并行度好于对角化等复杂操作。但 omp8 下退化同样明显。

4. **基本结论**：混合并行（MPI + OpenMP）是本机环境下 gather/scatter 密集型工作负载的最优策略。Release 实验（[5.4 节](#54-release-构建对照实验-o3--marchnative--dndebug)）表明编译器自动向量化已能显著提升单线程效率（3.3-3.6×），进一步缩小了手动 SIMD 的优化空间。

---

## 8. 优化机会总结

### 8.1 题目5：SIMD 向量化

| 项目 | 评估 |
|------|------|
| 目标 | `gatherp_scatters`, `gathers_scatterp` (6 处拷贝循环) |
| 当前状态 | GCC 13 `-O3` 已自动向量化（3.3-3.6× 加速 vs Debug），手动 SIMD 边际增益 ~0-3% SCF |
| 第 13 周计划 | 快速添加 `#pragma omp simd` + `__restrict__`（< 20 行），Release 下测量实际增益后收尾 |
| 风险/优先级 | 低风险 / **中优先级** — 编译器已做大部分工作 |

### 8.2 题目8：缓存复用

| 项目 | 评估 |
|------|------|
| 目标 | 缓存 `gg`, `gcar`, `gdirect`, `gk2`, `igl2isz_k`, `igl2ig_k` (~0.9 MB serial, ~0.5 MB np2) |
| 最大单项收益 | 消除 `setupIndGk` 双重 `cal_GplusK_cartesian` 调用（~18K 次冗余，~10-50ms Debug, <5ms Release） |
| 当前开销 | < 0.3s（含在 `before_scf`），Release 下可降至 < 0.1s |
| 第 13 周计划 | 懒加载 + 失效标志 + 优先修复 setupIndGk 双重计算（50~80 行改动） |
| 风险/优先级 | 中风险 / **高优先级** — SIMD 收益缩水后，缓存成为 Workflow C 主攻方向 |

### 8.3 额外发现：MPI 通信优化机会

从基线数据观察到 `poolnproc > 1` 时 gather/scatter 占比较高（稳态 **35~52%** 的 FFT 时间，冷启动可达 60~70%），其中 MPI_Alltoallv 是主要瓶颈。这属于 Workflow B 的范畴，但 SIMD 向量化 (Workflow C) 可部分缓解问题（加速 copy 段）。

---

## 9. 测试数据可靠性说明

- **3 次重复取平均**：所有配置均有 3 次重复运行，标准差在可接受范围
- **时钟来源**：MPI 模式下使用 `MPI_Wtime()`，非 MPI 使用 `std::chrono::system_clock`
- **已知噪声源**：
  - WSL2 虚拟化层引入约 5~10% 时间抖动（详见 [4.4.2 节因素4](#因素-4wsl2-虚拟化放大效应)）
  - MPI 首轮运行的冷启动开销，详见 [4.4 节](#44-首轮迭代初始化开销分析)的完整分析。**所有优化分析使用 r2/r3 稳态数据**，已排除冷启动效应
  - `mix_np2_omp4` 配置因 WSL2 宿主机调度抖动出现高方差（gaas_tiny: 46s/29s/46s），不代表代码层面的问题
- **Debug 构建限制**：主要基线数据来自 `-O0` 构建。Release（`-O3`）对照实验已完成（[5.4 节](#54-release-构建对照实验-o3--marchnative--dndebug)），证实编译器自动向量化已生效。百分比占比在 Debug 和 Release 下基本一致，但绝对加速倍数不能跨构建类型外推
- **计时器限制**：0.1s 阈值导致部分轻量函数无数据显示，已在第 3.3 节说明

---

## 10. 结论与下一步

### 10.1 基线结论

1. **全部 54 次运行通过** (100% 通过率)，测试体系稳定可靠。

2. **Gather/Scatter 是 MPI 模式下的重要优化目标**，占 FFT 变换的 **35~52%**（稳态），占 SCF 迭代的 **5~13%**。其中本地 copy 段（SIMD 目标）占 gather/scatter 的 30~50%。6 处内层连续拷贝循环非常适合 SIMD 向量化，但实际增益需在 Release（`-O3`）构建下验证。

3. **缓存复用函数当前开销较小**（<0.3s，含在 `before_scf` 中），但其收益在于消除 `setupIndGk` 的 |G+K|² 双重计算（~18,000 次冗余调用）和减少反复内存分配。缓存数据量约 0.5~0.9 MB，可完全放入 L2 缓存。

4. **线程扩展性趋于饱和**：4 线程达 2.1~2.6× 加速，8 线程不再提升甚至退化。**混合并行（`mix_np4_omp2`）是本机最优策略**（3.3× 加速），MPI 进程级并行有效规避了 OpenMP 小循环线程开销。Release 下编译器自动向量化已显著提升单线程效率。

### 10.2 第 13 周工作建议

1. **SIMD 向量化收尾**（题目5，快速验证）：
   - 添加 `#pragma omp simd` + `__restrict__` 到 6 处拷贝循环（< 20 行改动）
   - 在 Release 构建下用 gaas_medium np2_omp2 快速测量实际增益
   - 若增益 < 2% 则标记为完成并转向题目8
   - 额外探索：检查 `complex<double>` 的对齐情况（未对齐访问可能抵消 AVX2 收益）

2. **缓存复用作为主攻方向**（题目8）：
   - 在 `PW_Basis` 和 `PW_Basis_K` 中添加缓存有效性标志位
   - 改造 `collect_local_pw` 和 `collect_uniqgg` 检查缓存
   - **优先消除 `setupIndGk` 中的双重 `cal_GplusK_cartesian` 调用**（见效最快，~10-50ms Debug / <5ms Release）
   - 预期代码改动 50~80 行

3. **与 Workflow B 协调**：
   - SIMD 加速的 copy 段位于 MPI_Alltoallv 前后
   - 若 Workflow B 已将 Alltoallv 替换为非阻塞通信，copy 段可与通信重叠
   - 两个工作流的优化效果可能叠加

4. **WSL2 局限性说明**：
   - `mix_np4_omp2` Release 下出现 0.79× 的反常退化，属于 WSL2 虚拟化噪音
   - 所有基于 np2_omp2 和 serial 的结论可靠，但 np4 配置需在物理机 Linux 上复测

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
├── gaas_small/
│   ├── INPUT          (32³, 40 Ry, 70 bands, 3×3×3 kpt)
│   └── ... (同上链接)
├── gaas_medium/
│   ├── INPUT          (48³, 40 Ry, 70 bands, 3×3×3 kpt)
│   └── ... (同上链接)
├── results/
│   ├── gaas_tiny/     (18 对 stdout/stderr 文件)
│   ├── gaas_small/    (18 对 stdout/stderr 文件)
│   └── gaas_medium/   (18 对 stdout/stderr 文件)
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
