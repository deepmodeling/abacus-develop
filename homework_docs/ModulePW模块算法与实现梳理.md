# module_pw 模块算法与实现梳理

本文整理 module_pw 的核心算法、数据结构与实现路径，覆盖 CPU/MPI/GPU/DSP 版本的平面波基函数与 FFT 变换流程。重点放在平面波集合的构造、并行分配、FFT 流水线与不同设备后端的差异点。

## 1. 模块定位与核心目标

module_pw 负责把实空间函数与倒空间平面波系数互相变换，并在并行环境中分配平面波与实空间网格。其核心任务包括：

- 给定晶格与能量截断，构造 FFT 网格与可用平面波集合
- 将倒空间球面内的 G 向量按 sticks（固定 x,y，沿 z 方向）组织
- 在 MPI 池内分配 sticks 与 z-plane，建立收发索引
- 在 CPU/MPI/GPU/DSP 上执行 real2recip / recip2real 变换

数学定义（以 Gamma-only 或一般 k 点为例）：

- 平面波：$\langle r | g \rangle = \frac{1}{\sqrt{\Omega}} e^{i g\cdot r}$
- 展开：$f(r) = \frac{1}{\sqrt{\Omega}} \sum_g c(g) e^{i g\cdot r}$
- 变换：$c(g) = \int f(r) e^{-i g\cdot r} dr$

对于 k 点：

- $\langle r|g,k\rangle = \frac{1}{\sqrt{\Omega}} e^{i (g+k)\cdot r}$
- $f(r) = \frac{1}{\sqrt{\Omega}} \sum_g c(g,k) e^{i (g+k)\otimes r}$

实现中为简化 FFT，实际使用 $f'(r)=e^{-i k\cdot r} f(r)$ 的等价形式，变换只包含 $g$。这样所有 k 点共享相同的 FFT 网格与 G 向量集合，只是在每个 k 点上筛选满足 $(G+K)^2$ 截断的子集并建立映射。

参考：
- 平面波基类 [source/source_basis/module_pw/pw_basis.h](source/source_basis/module_pw/pw_basis.h)
- k 点扩展 [source/source_basis/module_pw/pw_basis_k.h](source/source_basis/module_pw/pw_basis_k.h)

## 2. 主要类与职责

### 2.1 PW_Basis

职责：

- 建立 FFT 网格与倒空间球面（根据 cutoff）
- 分配平面波 sticks 与实空间 z-plane（MPI）
- 提供 real2recip / recip2real 变换
- 维护 G 向量、$G^2$、映射表等基础数据结构

关键实现：

- 初始化网格与晶格矩阵：
  - 输入 lat0/latvec/gridecut，计算 $\Omega$、G、GT、GGT
  - 根据 gridecut 估算最小 FFT box，并扫描可用 G 点修正到最小覆盖球面
  - 选择 2/3/5/7 可分解的 FFT 维度，保证 FFT 性能
  - 见 [source/source_basis/module_pw/pw_init.cpp](source/source_basis/module_pw/pw_init.cpp)
- 平面波 sticks 分配：
  - 方法 1：按 plane waves 数量均衡
  - 方法 2：按 sticks 数量均衡
  - 见 [source/source_basis/module_pw/pw_distributeg.cpp](source/source_basis/module_pw/pw_distributeg.cpp)
  - 具体算法实现：
    - 方法 1 [source/source_basis/module_pw/pw_distributeg_method1.cpp](source/source_basis/module_pw/pw_distributeg_method1.cpp)
    - 方法 2 [source/source_basis/module_pw/pw_distributeg_method2.cpp](source/source_basis/module_pw/pw_distributeg_method2.cpp)
- 实空间 z-plane 分配：
  - 将 nz 按进程数均分，余数平均分配
  - 输出 startz/numz 与本进程的 nplane/startz_current
  - 见 [source/source_basis/module_pw/pw_distributer.cpp](source/source_basis/module_pw/pw_distributer.cpp)
- MPI gather/scatter：
  - sticks <-> planes 转换，内部使用 MPI_Alltoallv
  - numg/numr/startg/startr 决定收发布局
  - 见 [source/source_basis/module_pw/pw_gatherscatter.h](source/source/source_basis/module_pw/pw_gatherscatter.h)
- FFT 变换：
  - CPU 版本采用 XY FFT + MPI + Z FFT 组合
  - GPU 版本采用单机 3D FFT 并依赖 1D/3D 映射
  - 见 [source/source_basis/module_pw/pw_transform.cpp](source/source_basis/module_pw/pw_transform.cpp)
  - 见 [source/source_basis/module_pw/pw_transform_gpu.cpp](source/source_basis/module_pw/pw_transform_gpu.cpp)

### 2.2 PW_Basis_K

职责：

- 在 PW_Basis 基础上支持多 k 点
- 对每个 k 点构建 $G+K$ 截断与映射表
- 支持 k 点相关 real2recip / recip2real
- 维护每个 k 的 plane waves 子集与能量信息

关键实现：

- k 点输入与 $G+K$ cutoff：
  - kvec_d -> kvec_c 变换
  - 取 kmaxmod，计算 gk_ecut 与 ggecut
  - 非 Gamma k 点自动关闭 gamma_only
  - 见 [source/source_basis/module_pw/pw_basis_k.cpp](source/source_basis/module_pw/pw_basis_k.cpp)
- 每个 k 的 plane wave 选择与映射：
  - 对全局 npw 扫描，选择满足 $(G+K)^2 \le gk_ecut$ 的子集
  - igl2isz_k 保存 stick+z 映射，igl2ig_k 保存全局 ig
  - gk2 保存 $(G+K)^2$，并支持 erf 平滑调制
  - 见 [source/source_basis/module_pw/pw_basis_k.cpp](source/source_basis/module_pw/pw_basis_k.cpp)
- 变换逻辑：
  - 使用与 PW_Basis 相同的 FFT 流程，但在读写平面波时使用 igl2isz_k
  - DSP 路径使用 3D FFT 与 ig2ixyz_k_cpu 映射
  - 见 [source/source_basis/module_pw/pw_transform_k.cpp](source/source_basis/module_pw/pw_transform_k.cpp)
  - 见 [source/source_basis/module_pw/pw_transform_k_dsp.cpp](source/source_basis/module_pw/pw_transform_k_dsp.cpp)

### 2.3 PW_Basis_Sup

职责：

- dense grid（超网格）与 smooth grid sticks 保持一致
- 用于 USPP 密网格，避免不同网格之间 stick 对齐差导致的通信与重排成本

关键实现：

- 方法 3 分配 sticks，优先复用 smooth grid sticks
- 在 dense grid 上构造 ig2isz/is2fftixy 时，先沿用 smooth grid 的顺序，再补齐剩余 sticks
- 见 [source/source_basis/module_pw/pw_basis_sup.cpp](source/source_basis/module_pw/pw_basis_sup.cpp)

### 2.4 PW_Basis_Big / PW_Basis_K_Big

职责：

- 将多个 FFT 网格组合成更大的 box（历史兼容，主要用于 LCAO 集成）
- 在 distribute_r 中以 bx/by/bz 为块重新划分 z-plane

关键实现：

- [source/source_basis/module_pw/pw_basis_big.h](source/source_basis/module_pw/pw_basis_big.h)
- [source/source_basis/module_pw/pw_basis_k_big.h](source/source_basis/module_pw/pw_basis_k_big.h)

## 3. 核心数据结构

以下结构在 plane wave 分配与 FFT 变换中关键：

- `ig2isz`: 1D plane wave index -> (is, iz) 映射，存为 `isz = is * nz + iz`
- `is2fftixy`: stick index -> ixy（$iy + ix * fftny$）
- `istot2ixy`: 全局 stick 序号 -> ixy
- `fftixy2ip`: ixy -> MPI rank（该 stick 属于哪个进程）
- `startz/numz`: z-plane 分配区间
- `numg/numr/startg/startr`: MPI_Alltoallv 发送、接收布局
- `gdirect/gcar/gg`: 直接坐标、笛卡尔坐标与 $G^2$ 列表
- `ig2igg/gg_uniq`: $G^2$ 去重索引与唯一列表，用于按能量分组

具体定义在：

- [source/source_basis/module_pw/pw_basis.h](source/source_basis/module_pw/pw_basis.h)

## 4. 网格与 cutoff 的构造

### 4.1 FFT 网格尺寸

流程：

1. 根据 gridecut 估计每轴最小 $ibox$
2. 扫描所有可用 G 点，找出覆盖球面的最大半径
3. 选取最小可由 2/3/5/7 分解的 FFT 维度
4. 在 full_pw 模式下允许奇偶约束并扩展 FFT box

实现见：

- [source/source_basis/module_pw/pw_init.cpp](source/source_basis/module_pw/pw_init.cpp)

### 4.2 plane wave cutoff 与 gamma-only

- $ggecut = ecut / (2\pi/lat0)^2$
- gamma-only 时利用 $F(-k)=F(k)^*$，仅保存半个频率轴
- xprime 决定半轴在 x 还是 y 上，影响 fftnx/fftny 取值

实现见：

- [source/source_basis/module_pw/pw_init.cpp](source/source_basis/module_pw/pw_init.cpp)

## 5. plane wave sticks 分配算法

### 5.1 方法 1：按 plane waves 数量均衡

流程：

1. 统计每个 (x,y) stick 长度 `st_length2D`
2. 收集 sticks 并按长度降序排列（优先分配长 sticks）
3. 每次选择当前 plane waves 最少的核分配 stick
4. 形成 `fftixy2ip` 与 `istot2ixy`，再构造 `ig2isz`/`is2fftixy`

实现：

- 统计与初始化 [source/source_basis/module_pw/pw_distributeg.cpp](source/source_basis/module_pw/pw_distributeg.cpp)
- stick 收集与分配 [source/source_basis/module_pw/pw_distributeg_method1.cpp](source/source_basis/module_pw/pw_distributeg_method1.cpp)

### 5.2 方法 2：按 sticks 数量均衡

流程：

1. 统计 sticks 与 plane waves
2. 按 stick 数量平均分配（顺序 ixy）
3. 在分配过程中同步统计 npw_per
4. 构建 `fftixy2ip` 与 `istot2ixy`

实现：

- [source/source_basis/module_pw/pw_distributeg_method2.cpp](source/source_basis/module_pw/pw_distributeg_method2.cpp)

### 5.3 方法 3：dense grid 与 smooth grid 一致化

流程：

1. smooth grid 已分配 stick
2. dense grid 先继承 smooth sticks 的分配
3. 剩余 sticks 按 plane waves 数量最少优先分配
4. 生成 dense grid 的 ig2isz/is2fftixy 映射，并维护 stick 顺序一致

实现：

- [source/source_basis/module_pw/pw_basis_sup.cpp](source/source_basis/module_pw/pw_basis_sup.cpp)

## 6. 实空间网格分配（z-plane）

策略：

- 将 z 方向平面均匀分配给各进程
- `numz[ip]` 和 `startz[ip]` 记录每个进程管理的 z 范围
- 本地 nplane 与 startz_current 决定本进程实空间存储区

实现：

- [source/source_basis/module_pw/pw_distributer.cpp](source/source_basis/module_pw/pw_distributer.cpp)

## 7. FFT 变换流程

### 7.1 CPU 路径（PW_Basis）

real2recip（非 gamma-only 版本）：

1. 输入复制到 auxr
2. FFT x/y（在每个进程的局部 plane 上）
3. gather planes -> sticks（MPI_Alltoallv）
4. FFT z（在 stick 方向）
5. 按 ig2isz 提取结果并除以 $nxyz$

recip2real：

1. 初始化 auxg
2. 按 ig2isz 写入 plane waves
3. IFFT z
4. scatter sticks -> planes（MPI_Alltoallv）
5. IFFT x/y

实现：

- [source/source_basis/module_pw/pw_transform.cpp](source/source_basis/module_pw/pw_transform.cpp)

### 7.2 CPU 路径（PW_Basis_K）

流程同 PW_Basis，但提取与写入时使用 igl2isz_k 和每个 k 的 npwk。对于 Gamma-only 的实数据，使用 r2c/c2r 路径并保持半轴约束。

实现：

- [source/source_basis/module_pw/pw_transform_k.cpp](source/source_basis/module_pw/pw_transform_k.cpp)

### 7.3 DSP 路径（PW_Basis_K）

- 使用 3D FFT 实现
- 通过 `ig2ixyz_k_cpu` 将 1D plane wave 映射到 3D FFT box
- 适用于非 gamma-only，目标是降低 MPI gather/scatter 开销

实现：

- [source/source_basis/module_pw/pw_transform_k_dsp.cpp](source/source_basis/module_pw/pw_transform_k_dsp.cpp)

### 7.4 GPU 路径（PW_Basis / PW_Basis_K）

特点：

- 使用 3D FFT
- poolnproc 必须为 1（单 GPU）
- 依赖 `ig2ixyz_gpu` 或 `ig2ixyz_k`
- 通过 pw_op 的 set_3d_fft_box_op / set_real_to_recip_output_op 进行映射与归一化

实现：

- PW_Basis [source/source_basis/module_pw/pw_transform_gpu.cpp](source/source_basis/module_pw/pw_transform_gpu.cpp)
- PW_Basis_K [source/source_basis/module_pw/pw_transform_k.cpp](source/source_basis/module_pw/pw_transform_k.cpp)

### 7.5 gather / scatter 的 MPI 数据路径

- gatherp_scatters: planes -> sticks + MPI_Alltoallv
- gathers_scatterp: sticks -> planes + MPI_Alltoallv
- 通过 istot2ixy 将连续内存布局转换为 stick 顺序，减少非连续访问

实现：

- [source/source_basis/module_pw/pw_gatherscatter.h](source/source_basis/module_pw/pw_gatherscatter.h)

## 8. 辅助算子（pw_op）

这些算子用于 3D FFT 路径中的数据装箱与输出处理：

- `set_3d_fft_box_op`: 1D plane waves -> 3D FFT box
- `set_real_to_recip_output_op`: 3D FFT box -> 1D plane waves（归一化）
- `set_recip_to_real_output_op`: 3D FFT box -> 实空间数组
- CPU 版本为简单循环，GPU 版本在 CUDA/ROCm 内核中实现

实现：

- 接口定义 [source/source_basis/module_pw/kernels/pw_op.h](source/source_basis/module_pw/kernels/pw_op.h)
- CPU 实现 [source/source_basis/module_pw/kernels/pw_op.cpp](source/source_basis/module_pw/kernels/pw_op.cpp)
- GPU 实现见 CUDA/ROCm 子目录

## 9. 关键调用顺序（推荐理解路径）

1. initmpi
2. initgrids
3. initparameters
4. setuptransform
5. collect_local_pw / collect_uniqgg
6. real2recip / recip2real

对应实现：

- init / 参数 [source/source_basis/module_pw/pw_init.cpp](source/source_basis/module_pw/pw_init.cpp)
- setuptransform [source/source_basis/module_pw/pw_basis.cpp](source/source_basis/module_pw/pw_basis.cpp)
- collect_local_pw [source/source_basis/module_pw/pw_basis.cpp](source/source_basis/module_pw/pw_basis.cpp)

## 10. 常见注意点

- gamma-only：仅半轴存储，需要确保 xprime 一致
- ggecut 与 gridecut：实际 plane wave cutoff 不能超过 FFT box cutoff
- GPU 路径限制：目前仅支持 poolnproc == 1
- PW_Basis_Sup：仅用于 dense grid 并依赖 smooth grid sticks
- k 点路径：非 Gamma k 点自动关闭 gamma_only

如需补充测试用例、性能分析或与上层模块的接口关系，请告诉我具体范围。