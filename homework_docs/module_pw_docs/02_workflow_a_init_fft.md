# Workflow A: 平面波分布初始化与 FFT 数据流分析

本工作流旨在解决 `module_pw` 在初始化截断球和网格以及发生 FFT 数据拷贝时的严重串行挂起和缓存低效问题。

## 1. 痛点：`count_pw_st` 的串行循环惩罚

### 旧算法表现与实现形式
在 `pw_distributeg.cpp` 内的 `count_pw_st` 函数中，核心任务是扫描 FFT 盒并统计落入截断球的平面波数量，同时记录每个 stick 的长度与底部 `z` 坐标。

旧实现的特征如下：
- **三重深度嵌套的 `for` 循环**：对 `ix/iy/iz` 进行步进扫描，扫描边界由 `nx/ny/nz`、`gamma_only`、`xprime` 与 `full_pw` 共同决定。
- **沉重的内部计算**：在最内层循环中计算 `f * (GGT * f)`，并做 `modulus <= ggecut` 判定；若 `full_pw` 为真则直接接受所有点。
- **完全串行**：当前实现没有并发标识。对大体系 + 高 `ggecut` 的三重扫描会成为明显的启动瓶颈。

### 需要被重构的维度
- **引入 OpenMP 线程级并行**：`npwtot/nstot` 与每个 stick 的长度可通过 reduction 或分块计数合并，避免竞态。
- **边界与坐标计算优化**：当前每个 `(ix,iy,iz)` 都会重复构造 `f` 并做矩阵乘法，可在 `iz` 外层预计算 `ix/iy` 相关部分，或预筛选 `ix/iy` 范围减少不必要计算。
- **分支路径可控化**：`gamma_only` 与 `full_pw` 的边界判断可上提，减少内层循环的分支开销。

## 2. 痛点：FFT 变换前后的不连续内存拷贝

### 旧算法表现与实现形式
位于 `pw_transform.cpp` 中的 `real2recip` 与 `recip2real` 负责实/倒空间变换。核心 FFT 本身之外，还包含多个数据重排与拷贝阶段：

```cpp
// 典型的旧代码伪实现范式
for (int is = 0; is < nsticks; ++is) {
    for (int iz = 0; iz < nzt; ++iz) {
        // 利用诸如 igl2isz, ig2ixyz 等间接数组来获得实际索引
        out_buf[ mapped_index_a ] = in_buf[ mapped_index_b ]; 
    }
}
```

- **查表开销与内存离散**：`ig2isz` 与 `is2fftixy` 触发频繁索引映射，导致缓存命中率偏低。
- **多阶段拷贝**：`real2recip` 先将输入拷贝到 `auxr`，然后 `gatherp_scatters` 将 `(nplane,fftxy)` 重排为 `(nz,nst)`，`recip2real` 路径同理。
- **OpenMP 并行粒度偏粗**：拷贝与重排循环虽已使用 OpenMP，但缺乏缓存友好的块划分与循环变换。

### 优化的改进思路
- **明确数据布局**：区分 `nrxx`、`nplane`、`nst` 与 `nz` 维度，识别可连续拷贝的内层维度并优先使用 `memcpy`。
- **块化重排**：对 `istot` 或 `is` 维度进行分块，让 `outp` 与 `inp` 拷贝走连续 cache line。
- **与 FFT 计划对齐**：结合 `xprime` 与 `gamma_only` 的半轴规则，避免在无效区域做多余拷贝。