# WorkflowC 旧实现与改进工作清单

本文档对应 WorkflowC 要求的“SIMD 向量化 + G 矢量/缓存复用”部分，整理现有旧实现的位置与改进所需工作。

## 一、旧实现（当前代码路径）

### 1. gather/scatter 重排循环（暂无 SIMD 向量化）

现有代码在内层循环做连续拷贝，仅有 OpenMP 并行，没有 SIMD 提示。例如：

```cpp
for (int iz = 0; iz < nplane; ++iz)
{
    outp[iz] = inp[iz];
}
```

分布在 gatherp_scatters 与 gathers_scatterp 的多处拷贝循环中，见 [source/source_basis/module_pw/pw_gatherscatter.h](source/source_basis/module_pw/pw_gatherscatter.h#L14-L183)。

### 2. G 矢量及相关量的重复计算

`collect_local_pw()` 每次调用都会重新分配并计算 `gg`/`gdirect`/`gcar`，没有缓存复用机制，见 [source/source_basis/module_pw/pw_basis.cpp](source/source_basis/module_pw/pw_basis.cpp#L135-L185)。

`collect_uniqgg()` 同样会遍历所有平面波并重新计算 `GGT * f`，见 [source/source_basis/module_pw/pw_basis.cpp](source/source_basis/module_pw/pw_basis.cpp#L193-L259)。

### 3. 多 k 点 |G+K|^2 的重复计算

`setupIndGk()` 中对每个 `ik`、`ig` 先计算一遍 `cal_GplusK_cartesian(ik, ig).norm2()` 统计 `npwk`，再计算一遍用于构建 `igl2isz_k`/`igl2ig_k`，存在重复开销，见 [source/source_basis/module_pw/pw_basis_k.cpp](source/source_basis/module_pw/pw_basis_k.cpp#L131-L189)。

`gk2` 本身已在数据结构中定义，但并未作为统一缓存复用的来源，见 [source/source_basis/module_pw/pw_basis_k.h](source/source_basis/module_pw/pw_basis_k.h#L75-L106)。

## 二、实现 WorkflowC 所需的改进工作

### A. gather/scatter 重排循环 SIMD 向量化

1. 在 `gatherp_scatters` 和 `gathers_scatterp` 的内层拷贝循环上添加 `#pragma omp simd`（包括单进程分支与 MPI 分支中的 `iz`/`izip` 循环）。
2. 视内存对齐情况选择 `aligned` 子句，若无法保证对齐，保持基础 `omp simd` 提示即可。
3. 结合现有 OpenMP 外层并行，避免过度嵌套影响调度；必要时明确 `schedule(static)`。
4. 对于短循环，可通过阈值判断跳过 SIMD，避免负收益。

### B. G 矢量与 |G+K|^2 缓存复用

1. **PW_Basis 缓存机制**
   - 新增缓存状态标志，例如 `gvec_cache_valid`、`gg_cache_valid`。
   - `collect_local_pw()` 在缓存有效时直接返回，避免重复分配与计算。
   - 在 `initgrids`、`distribute_g`、`get_ig2isz_is2fftixy` 等会改变 `ig2isz`/`is2fftixy` 或晶格参数的路径上失效缓存。
   - 若 `collect_uniqgg()` 频繁使用，可让其复用 `gg`/`gdirect` 的计算结果，避免再次遍历。

2. **PW_Basis_K 的 gk2 复用**
   - 在 `setupIndGk()` 中一次性计算 `gk2` 并填入缓存数组，后续统计与构建索引均复用该数组，避免双重 `cal_GplusK_cartesian()`。
   - 为 `gk2` 增加有效标志，确保 `kvec`/`gk_ecut`/`npw` 变更时失效。
   - 结合现有 `getgk2()` 接口，统一对外提供缓存数据。

## 三、建议的交付物

- 在 gather/scatter 的 4 处内层循环加入 `omp simd` 与必要的对齐提示，并记录向量化开关条件。
- 为 PW_Basis 增加 G 矢量相关缓存复用与失效规则，并更新 `collect_local_pw()`/`collect_uniqgg()`。
- 为 PW_Basis_K 增加 `gk2` 缓存的构建与复用逻辑，避免 `setupIndGk()` 中的重复计算。
- 补充测试：对比 SIMD 前后与缓存前后的结果一致性，并记录性能变化。
