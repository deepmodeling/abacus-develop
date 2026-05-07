# module_pw 综合报告

本文将 module_pw 的概览与三个工作流文档整合为一份对齐代码实现的综合报告，包含依赖关系与术语表。内容反映当前实现，并标出热点路径与优化目标。

## 1. 范围与目标

- 归纳 module_pw 的架构与数据流。
- 总结热点与对应的代码位置。
- 提供从公开 API 到内部原语的依赖关系。
- 提供术语与数据结构词汇表。

## 2. 主要组件

### 2.1 核心类

- PW_Basis: 单 k 点平面波基、FFT 网格、MPI 映射与实/倒空间变换。
- PW_Basis_K: 多 k 点扩展，维护 k 点相关映射与 $(G+K)^2$ 逻辑。

### 2.2 核心文件（CPU 路径）

- source/source_basis/module_pw/pw_basis.h
- source/source_basis/module_pw/pw_basis.cpp
- source/source_basis/module_pw/pw_basis_k.h
- source/source_basis/module_pw/pw_basis_k.cpp
- source/source_basis/module_pw/pw_distributeg.cpp
- source/source_basis/module_pw/pw_gatherscatter.h
- source/source_basis/module_pw/pw_transform.cpp

### 2.3 职责摘要

- pw_distributeg.cpp: 构建 sticks 与平面波分布映射。
- pw_gatherscatter.h: 完成 sticks <-> planes 转置与 MPI AlltoAllv。
- pw_transform.cpp: FFT 流水线、打包、变换与归一化。

## 3. 依赖关系图（自顶向下）

本节连接公开 API 与内部原语。

### 3.1 PW_Basis

- initgrids(...)
  - 计算 FFT 网格与晶格相关量（G, GT, GGT, gridecut_lat）。
- initparameters(gamma_only, pwecut, distribution_type, xprime)
  - 设置 ggecut、Gamma-only 行为与分布方式。
- setuptransform()
  - distribute_r() -> distribute_g() -> getstartgr() -> fft_bundle.initfft() -> setupFFT().
- real2recip(...)
  - XY FFT -> gatherp_scatters -> Z FFT -> 1/nxyz 归一化 -> ig2isz 写出。
- recip2real(...)
  - 写入 g-space -> Z FFT -> gathers_scatterp -> XY FFT -> 输出。

### 3.2 PW_Basis_K

- initparameters(gamma_only, gk_ecut, nks, kvec_d, distribution_type, xprime)
  - 计算 kmaxmod，若存在非零 k 则关闭 gamma-only。
- setuptransform()
  - distribute_r() -> distribute_g() -> getstartgr() -> setupIndGk() -> fft_bundle init.
- setupIndGk()
  - 对每个 (ik, ig) 计算 $(G+K)^2$ 并构建 igl2isz_k/igl2ig_k。
- real2recip/recip2real（k 点版本）
  - 与 PW_Basis 相同流程，但索引由 igl2isz_k 与 npwk 控制。

## 4. 数据流（FFT Pipeline）

### 4.1 通用路径（非 gamma-only）

1. 本地 Z-planes 上执行 XY FFT。
2. 通过 gatherp_scatters 将 planes -> sticks。
3. 通过 MPI_Alltoallv 进行跨进程 sticks/planes 交换。
4. sticks 上执行 Z FFT。
5. 归一化并通过 ig2isz 写出。

### 4.2 Gamma-only（仅单 k 点）

- 输入为实数，FFT 路径使用 r2c/c2r。
- 只保留半轴频域数据。
- fftnx 或 fftny 缩减到 nx/2+1 或 ny/2+1（由 xprime 决定）。

## 5. 热点与瓶颈

### 5.1 count_pw_st（pw_distributeg.cpp）

- 三重循环扫描 ix/iy/iz。
- 内层计算 modulus = f * (GGT * f)。
- gamma_only 与 full_pw 改变扫描范围与接受规则。

影响：大体系或高 ggecut 时启动成本显著。

### 5.2 gatherp_scatters / gathers_scatterp（pw_gatherscatter.h）

- 多次 (nplane, fftxy) <-> (nz, nst) 的重排。
- MPI_Alltoallv 强制全局同步。
- 内层循环连续但未显式向量化。

影响：规模扩大时通信与内存拷贝占主导。

### 5.3 real2recip / recip2real（pw_transform.cpp）

- 多阶段拷贝/转置/FFT。
- gamma-only 路径存在额外打包与半轴处理。

影响：缓存压力与冗余搬运成本上升。

## 6. 优化主题（对应工作流）

### 6.1 Workflow A：初始化与 FFT 拷贝

- 使用 OpenMP reduction 或分块私有计数并行化 count_pw_st。
- 减少内层分支并复用 ix/iy 相关预计算。
- 调整拷贝循环顺序，改善缓存局部性与减少间接访问。

### 6.2 Workflow B：MPI 通信

- 用非阻塞等价实现替换阻塞式 MPI_Alltoallv。
- pack/unpack 与通信重叠。
- float/double 的 MPI datatype 保持与现有一致。

### 6.3 Workflow C：SIMD 与缓存复用

- 在 gather/scatter 的连续拷贝循环中加入 omp simd。
- 视调用约定使用 __restrict__ 等以便向量化。
- 缓存 gg/gdirect/gcar 与 gk2，避免重复计算。
- 在 setupIndGk 中合并重复的 $(G+K)^2$ 遍历。

## 7. 风险与正确性注意事项

- 多 k 场景下 gamma-only 当前设计为关闭，除非完整实现逻辑，否则避免开启。
- gather/scatter 重构必须严格保持索引与映射规则。
- MPI datatype 选择依赖模板类型 T，需保持一致。

## 8. 术语表

- Stick: 固定 (ix, iy) 且 iz 连续的一条平面波线。
- Plane: 实空间的一个 Z-plane，每个进程持有一段 plane 区间。
- ig2isz: 本地平面波索引到 (is, iz) 线性索引的映射。
- is2fftixy: 本地 stick 索引到 (ix, iy) 的映射。
- istot2ixy: 全局 stick 索引到 (ix, iy) 的映射。
- fftixy2ip: (ix, iy) 到所属 MPI 进程的映射。
- ggecut: PW_Basis 的平面波截断（倒空间）。
- gk_ecut: PW_Basis_K 的 $(G+K)^2$ 截断。
- gamma_only: Gamma 点实值输入的半谱 FFT。
- xprime: FFT 轴顺序，同时决定 gamma-only 的半轴规则。
- nplane: 当前进程拥有的 Z-plane 数量。
- nst/nstot: 当前进程/全局 stick 数量。
- npw/npwtot: 当前进程/全局平面波数量。
- nrxx: 当前进程的实空间网格点总数。

## 9. 参考说明

更完整的实现细节见 homework_docs/module_pw_docs 下的概览与工作流文档。