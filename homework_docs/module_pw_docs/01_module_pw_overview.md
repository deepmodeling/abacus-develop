# 第一部分：module_pw 模块的旧算法与实现方式详解

`module_pw` 是 ABACUS 核心模块之一，负责管理平面波（Plane Wave）基组与实/倒空间的快速傅里叶变换（FFT）。在理解后续优化工作流之前，必须深刻理解现有架构的设计逻辑与底层限制。

## 1. 核心架构与类体系

`module_pw` 的实现采用面向对象封装，并在关键路径上叠加 MPI + OpenMP 多级并行。与平面波相关的核心数据结构和通信映射均由 `PW_Basis` 系列类集中管理。

### 1.1 `PW_Basis`（基础单 k 点类）
**功能职责**：管理平面波基组、FFT 网格、分布式映射表，并提供 `real2recip / recip2real` 变换接口。

**关键初始化与参数含义**：
- `initgrids(lat0, latvec, gridecut | nx,ny,nz)`：确定 FFT 网格与晶格参数，计算 `G/GT/GGT` 与 `gridecut_lat` 等。
- `initparameters(gamma_only, pwecut, distribution_type, xprime)`：设置平面波截断 `ggecut`、并行划分方式与 Gamma-only 的半轴规则。
- `setuptransform()`：内部依次执行 `distribute_r()`、`distribute_g()`、`getstartgr()`，并调用 `fft_bundle.initfft()` 与 `setupFFT()` 完成 FFT 计划。

**G 矢量筛选与扫描范围**：
- 在 `count_pw_st()` 中遍历 $\mathbf{G}=(i_x,i_y,i_z)$ 的整数网格，并使用 $|\mathbf{G}|^2 = \mathbf{f}^T GGT \mathbf{f}$ 判定是否落入截断球。
- **Gamma-only + `xprime`** 时只扫描半轴（`ix >= 0` 或 `iy >= 0`），并相应缩小 `fftnx` 或 `fftny`。
- **full_pw** 模式会忽略 `ggecut`，扫描完整 FFT 盒并直接收集全部平面波。

**分布式映射表（核心数据结构）**：
- `istot2ixy`：全局 stick 索引到 `(ix,iy)`。
- `is2fftixy`：当前进程 stick 索引到 `(ix,iy)`。
- `ig2isz`：本地平面波 `ig` 到 `(is, iz)` 的线性索引。
- `fftixy2ip`：`(ix,iy)` 对应 stick 所属进程。
- `startz/numz`：每个进程负责的 Z-plane 区间。
- `numg/numr` 与 `startg/startr`：为 `MPI_Alltoallv` 构建收发 counts 与 offsets（stick->plane / plane->stick）。

**MPI 并行分配**：
- **倒空间（Reciprocal Space）**：以 stick 为单位分配，各 stick 对应固定 `(ix,iy)` 且 `iz` 连续。
- **实空间（Real Space）**：按 Z-planes 切片分配，`numz[ip]` 为每个进程持有的 plane 数。

### 1.2 `PW_Basis_K`（多 k 点扩展类）
**功能职责**：在 `PW_Basis` 基础上支持多 k 点，维护 k 点相关的平面波映射，并提供 k 点版本的 `real2recip / recip2real`。

**截断判据更新**：
- `gk_ecut` 使用 $|\mathbf{G}+\mathbf{k}|^2 \leq gk\_ecut$ 判定。
- `initparameters()` 会根据所有 k 点的最大 $|\mathbf{k}|$ 调整 `ggecut`，并在存在非零 k 时强制关闭 `gamma_only`（`kmaxmod > 0 -> gamma_only = false`）。

**索引与缓存结构**：
- `npwk[ik]`：第 `ik` 个 k 点的本地平面波数量；`npwk_max` 为 k 点中的最大值。
- `setupIndGk()` 会两次遍历 `ig`：先统计 `npwk`，再构建 `igl2isz_k/igl2ig_k`，并使用 `cal_GplusK_cartesian().norm2()` 判断是否落入截断球。
- `gk2` 是用于存储 $(G+K)^2$ 的数组，但在现有实现中仍存在重复计算路径（见后续 Workflow C）。

**局限性**：多 k 点场景下 `gamma_only` 被强制关闭，导致无法使用半轴 FFT 的内存与计算节省路径。

## 2. 变换机制（FFT Pipeline）

ABACUS 使用自研的 FFT 分解流程与手动通信，而非直接依赖全局 3D-FFT MPI 接口。`PW_Basis::real2recip / recip2real` 的关键数据流如下：

1. **XY 方向 FFT**：在本地 Z-planes 上先做 `fftxyfor/fftxybac` 或 Gamma-only 的 `fftxyr2c/fftxyc2r`。
2. **planes<->sticks 重排**：调用 `gatherp_scatters` 或 `gathers_scatterp` 将 `(nplane,fftxy)` 与 `(nz,nst)` 间转换。
3. **MPI AlltoAllv 转置**：进程间数据交换由 `MPI_Alltoallv` 负责，收发计数由 `numg/numr` 与 `startg/startr` 控制。
4. **Z 方向 FFT**：执行 `fftzfor/fftzbac` 完成完整 3D 变换。
5. **收敛与归一化**：`real2recip` 在输出前除以 `nxyz`；`recip2real` 在输出阶段应用系数 `factor`。

**Gamma-only 变体要点**：
- 仅在 `PW_Basis`（单 k 点）中启用，`PW_Basis_K` 遇到非 Gamma k 会强制关闭。
- 对实值输入使用 r2c/c2r 路径，并只保留半轴频域数据，`fftnx/fftny` 会缩减到 `nx/2+1` 或 `ny/2+1`。
