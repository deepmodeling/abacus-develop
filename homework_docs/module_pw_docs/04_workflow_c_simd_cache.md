# Workflow C: 底层数据吞吐优化与缓存复用分析

当跨节点并行的效率在通信优化方面提升后，底层的单核运算效率、向量化执行吞吐率和寄存器换页频次就会成为新的核心瓶颈。Workflow C的主要方向是在不改变大逻辑和外部接口的前提下，对最底层的算力进行优化。

## 0. Workflow 涉及的文件与模块范围

本 workflow 主要覆盖两类低层热点：一类是平面波数据的 gather/scatter 拷贝重排，一类是平面波基底与$G+K$相关的只读/不变数据构建与复用。

- `source/source_basis/module_pw/pw_gatherscatter.h`
    - `gatherp_scatters`：进程内或进程间的平面波分块重排，内层包含连续拷贝循环。
    - `gathers_scatterp`：与上面互为逆操作的拷贝重排路径，同样包含大量标量拷贝循环。
- `source/source_basis/module_pw/pw_basis.cpp`
    - `PW_Basis::collect_local_pw()`：构建当前进程的局域平面波集合，生成 `gg/gdirect/gcar` 等数组。
    - `PW_Basis::collect_uniqgg()`：对全体平面波做筛选/去重并计算派生量（如$GGT \cdot f$）。
- `source/source_basis/module_pw/pw_basis_k.cpp`
    - `PW_Basis_K::setupIndGk()`：基于$G+K$构建索引映射，统计 `npwk` 并生成 `igl2isz_k/igl2ig_k` 等结构。

## 0.1 module_pw 内部调用链速览

以下调用链均可在 module_pw 的实现或测试/注释中直接对应到。

- **实空间/倒空间变换路径**（`pw_transform.cpp`）
    - `PW_Basis::real2recip()` -> `fftxyfor()` -> `gatherp_scatters()` -> `fftzfor()` -> 按 `ig2isz` 打包输出。
    - `PW_Basis::recip2real()` -> 按 `ig2isz` 解包 -> `fftzbac()` -> `gathers_scatterp()` -> `fftxybac()` -> 写回实空间。
- **平面波基底初始化路径**（`pw_basis_k.cpp` 与 `pw_init.cpp`）
    - `initgrids()` -> `initparameters()` -> `setuptransform()`
    - `setuptransform()` 内部顺序：`distribute_r()` -> `distribute_g()` -> `getstartgr()` -> `setupIndGk()` -> `fft_bundle.initfft()`。
- **平面波几何量构建路径**（`pw_basis.cpp` + `pw_basis_k.h` 使用示例）
    - `setuptransform()` 完成后，调用 `collect_local_pw()` 生成 `gg/gdirect/gcar`，再调用 `collect_uniqgg()` 构建 `gg_uniq/ig2igg`。

以上模块互相独立但共享同一类性能瓶颈：对连续内存的标量拷贝缺乏显式向量化提示，以及对不变数据的重复遍历与重建。

## 1. 痛点：缺乏 SIMD 优化的散列拷贝

### 1.1 模块实现逻辑
`gatherp_scatters` 与 `gathers_scatterp` 的核心逻辑是把一维连续的平面波数据块按进程/平面/索引重排。
典型实现是双层或多层循环，外层定位分块边界，内层按 `iz`/`izip` 逐元素拷贝。
分布在 `pw_gatherscatter.h` 中的拷贝重排操作大多是标量循环，典型模式如下：

```cpp
// 旧的逐个元素赋值循环
for (int iz = 0; iz < nplane; ++iz)
{
    outp[iz] = inp[iz];
}
```

这类循环在 `gatherp_scatters` 与 `gathers_scatterp` 中多处出现，贯穿单进程路径与 MPI 路径。其性能主要受制于内层连续拷贝的吞吐率，以及编译器对别名与对齐的保守判断。

从算法流程上看，MPI 路径包含三段明确的重排：

- 对于 gatherp_scatters（从 Planes 转换为 Sticks），重排包括：
    - 本地重排 1（打包提取）：提取有效网格点，把稀疏的平面数据 (nplane, fftnxy) 重排为紧凑的全局有效切片 (nplane, nstot)，便于一次性 MPI_Alltoallv。
    - 进程间交换（矩阵转置）：MPI_Alltoallv 完成不同进程之间的数据切片交换，将 (nplane, nstot) 交换为按来源进程分块的片段集合 (numz[ip], ns, poolnproc)。
    - 本地重排 2（解包拼接）：将按进程分块的碎片片段 (numz[ip], ns, poolnproc)，沿 Z 轴拼接重排为当前进程负责的完整柱状数据 (nz, ns)。

- 对于 gathers_scatterp（从 Sticks 还原为 Planes），重排包括：
    - 本地重排 1（切断打包）：将当前进程负责的完整柱状数据 (nz, ns) 按目标进程所需的厚度切断，打包重排为按目标进程分块的集合 (numz[ip], ns, poolnproc)，便于一次性 MPI_Alltoallv。
    - 进程间交换（矩阵转置）：MPI_Alltoallv 完成不同进程之间的数据片段交换，将发件缓冲区的 (numz[ip], ns, poolnproc) 交换为当前进程所辖楼层接收到的紧凑有效薄片 (nplane, nstot)。
    - 本地重排 2（补零散布）：在对输出网格进行全局清零后，把紧凑的有效薄片 (nplane, nstot) 重新散布并映射回包含无效真空区域的稀疏平面网格 (nplane, fftnxy)。

而当只有一个进程的时候，`gatherp_scatters` 里直接将当前进程分配到的z平面数设置为所有的xy切片数 `nz=nplane`；`gathers_scatterp` 里直接将当前进程分配到的stick数量设置为所有的stick数 `ns=nstot`。

此外在这段代码里，使用一维数组来存储二维的数据，保证了内存的连续性以及兼容MPI和OpenMP。首先第一步，对于一个xy平面，用`ixy`存储这个点在平面上的索引，使用公式$ixy = y*nx+x$。第二步，由于每个stick的内存地址是连续的，因此为了寻找空间坐标为`ixy`和`iz`的点，它的索引就是$ixy*nz+iz$。

### 1.2 逻辑局限性
- **自动向量化不稳定**：`outp` 与 `inp` 可能别名，编译器不敢做 aggressive vectorize。
- **连续访问未显式标注**：内层循环逻辑上连续，但缺少 `simd` 或 `restrict` 语义提示。
- **拷贝粒度过小**：单元素赋值使得 cache line 无法充分填充，拷贝带宽难以饱和。

### 1.3 初步解决方向
- 在 `iz/izip` 内层循环上显式添加 `#pragma omp simd`，减少编译器不确定性。
- 若调用约定允许，增加 `__restrict__` 或等价约束来声明不重叠。
- 对连续区域尝试 `memcpy` 或按 cache line 的块拷贝模式，以提高带宽利用率。

## 2. 痛点：循环内只读/不变数据的冗余重建

### 2.1 模块实现逻辑
`PW_Basis` 与 `PW_Basis_K` 负责构建平面波基底相关的几何量、筛选量以及$G+K$的索引映射。典型流程如下：

- `PW_Basis::collect_local_pw()`
    - 读取晶胞与截断参数，将可用平面波映射为 `gg/gdirect/gcar` 等数组。
    - 每次调用都会重新分配并计算这些数组。
- `PW_Basis::collect_uniqgg()`
    - 再次遍历全部平面波，构建去重后的 `gg` 与派生量（如$GGT \cdot f$）。
- `PW_Basis_K::setupIndGk()`
    - 对每个 `(ik, ig)` 计算$|G+K|^2$并统计 `npwk`，然后再次遍历生成 `igl2isz_k/igl2ig_k`。

其中 `collect_uniqgg()` 的现有算法要点是：
- 对每个 `ig` 计算$G^2 = f \cdot (GGT \cdot f)$并存入 `tmpgg`。
- 用 `heapsort` 对 `tmpgg` 排序并记录 `sortindex`。
- 线性扫描去重：对相近的$G^2$（阈值 `1.0e-8`）求平均，构建 `gg_uniq`，同时生成 `ig2igg` 映射。

`setupIndGk()` 的现有算法要点是：
- 第一次遍历 `(ik, ig)` 计算$|G+K|^2$并计数得到 `npwk[ik]`，同时用 `MPI_Allreduce` 检查全局是否有平面波。
- 第二次遍历填充 `igl2isz_k/igl2ig_k` 索引映射。
- 末尾调用 `get_ig2ixyz_k()` 完成后续索引派生。

这些数据在 SCF 或多次变换过程中往往保持不变，但却反复重建，导致不必要的算力消耗与内存抖动。

### 2.2 逻辑局限性
- **生命周期与不变性未绑定**：`gg/gdirect/gcar` 与 `gk2` 等应与 `latvec/ggecut/kvec` 绑定，但当前每次调用均重建。
- **双重遍历与重复计算**：`setupIndGk()` 对同一 `(ik, ig)` 重复计算$|G+K|^2$。
- **频繁分配导致抖动**：重复 `new/resize` 使得缓存失效，降低局域性。

### 2.3 初步解决方向
- **缓存复用**：将 `gg/gdirect/gcar` 与 `gk2` 的构建与生命周期绑定到 `ggecut` 或晶胞变化上，非必要不重建。
- **单次遍历复用结果**：在 `setupIndGk()` 中一次计算 `gk2` 后复用统计与索引构建，避免双重遍历。
- **失效规则明确化**：仅当 `latvec`、`ggecut` 或 `kvec` 变化时标记失效，其他路径直接复用。