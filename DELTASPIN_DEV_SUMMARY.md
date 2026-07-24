# DeltaSpin LCAO 内存优化开发总结

## 一、背景

`lcao_PI_sub_save_` 存储每个 k-point 下每个约束原子的子空间投影矩阵 P_I_sub(k)，用于 DeltaSpin lambda 优化循环中的子空间加速。原始实现通过 `gather_sub_matrix_to_all` 将 2D 块循环分布的局部数据聚合成完整 nbands×nbands 矩阵，导致每个进程都持有完整副本——内存冗余为 NPROC_IN_POOL 倍。

对 52 原子、424 bands、72 k-points、16 进程的典型体系，每进程缓存达 ~10.7 GB，而分布式存储仅需 ~0.67 GB。

---

## 二、已实施优化

### Phase 1：消除 gather 全量冗余（方案 3.1 + 3.4）

#### 2.1 数据结构变更

| 成员 | 旧类型 | 新类型 |
|------|--------|--------|
| `lcao_PI_sub_save_` | `vector<vector<vector<complex<double>>>>` | `vector<map<int, vector<complex<double>>>>` |
| `lcao_PI_sub_diag_` | （不存在） | `vector<map<int, vector<double>>>` |

**变更要点**：
- 外层第二维从 `vector<complex<double>>[nat]`（所有原子）改为 `map<int, vector<complex<double>>>`（仅约束原子），消除了 nat - nat_constrained 个空 vector 的开销和 `[iat]` 随机访问语义
- 内层 vector 从 `nbands²` 元素（完整矩阵）改为 `nrow × ncol_bands` 元素（本地 2D 块），大小降至原来的 ~1/NPROC_IN_POOL
- 新增 `lcao_PI_sub_diag_` 缓存每个约束原子 P_I_sub 的对角线元素（double 类型），供一阶模式直接使用，无需从 complex 矩阵中提取

#### 2.2 `calculate_PI_sub_from_hr` 修改

**文件**：`lcao_subspace.cpp`

旧路径：`folding_HR → pzgemm → PI_sub_local[nloc_eij] → gather → PI_sub_full[nn]`

新路径：`folding_HR → pzgemm → PI_sub_local[nloc_eij] → memcpy 到输出指针`

- 删除了 `gather_sub_matrix_to_all` 调用
- fp32 路径同理：改为 `cast_up_to_double → memcpy`
- 串行路径不变（本身就是完整的）

#### 2.3 `calculate_delta_hcc_lcao` 修改

**文件**：`lcao_subspace.cpp`

新增两个重载：

**重载 1（本地块版本）**：
```cpp
void calculate_delta_hcc_lcao(
    complex<double>* h_sub_local,                        // 本地块 [nrow*ncol_bands]
    const map<int, vector<complex<double>>>& PI_sub_local, // 约束原子 → 本地块
    const Vector3<double>* lambda, int nbands, int ik,
    bool full_update, const Parallel_Orbitals* ParaV);
```
- 遍历 `map` 的 `for (const auto& [iat, pi_local] : PI_sub_local)` 替代原 `for (iat=0; iat<nat; iat++) { if (empty) continue; }`
- axpy 循环长度从 `nbands*nbands` 变为 `nloc_eij = nrow * ncol_bands`
- npol=2 目前仅使用 coeff0（z 分量），spin-flip 项尚未启用

**重载 2（全量矩阵版本，向后兼容）**：
```cpp
void calculate_delta_hcc_lcao(
    complex<double>* h_sub,                               // 完整矩阵 [nbands²]
    const vector<vector<complex<double>>>& PI_sub,        // 所有原子
    const Vector3<double>* lambda, int nbands, int ik,
    bool full_update);
```
- 保留给诊断函数（`run_lambda_local_diagnostic`、`run_lambda_scan_diagnostic`）和 `DiagonalizationEngine` 使用
- 逻辑与原来完全一致，不做改动

#### 2.4 Branch 1（子空间缓存构建）修改

**文件**：`cal_mw_from_lambda.cpp`，i_step==-2 分支

- 循环中对每个约束原子：
  1. 调用 `calculate_PI_sub_from_hr` 直接写入 `lcao_PI_sub_save_[ik][iat]`（本地块）
  2. 调用 `extract_diagonal_from_local_block` 提取对角线，写入 `lcao_PI_sub_diag_[ik][iat]`（含 MPI reduce）
- 不再调用 `resize(nat)` 和 `clear()`

#### 2.5 Branch 2（一阶模式）修改

**文件**：`cal_mw_from_lambda.cpp`

旧代码：
```cpp
for (int iat = 0; iat < nat; iat++) {
    if (lcao_PI_sub_save_[ik][iat].empty()) continue;
    double p_diag = lcao_PI_sub_save_[ik][iat][ib + ib*nbands].real();
    ...
}
```

新代码：
```cpp
for (const auto& [iat, diag] : lcao_PI_sub_diag_[ik]) {
    double p_diag = diag[ib];
    ...
}
```

- 直接从 double 缓存读取对角线，避免 complex→real 转换和 `nbands*nbands` 索引计算
- map 遍历天然跳过非约束原子

#### 2.6 Branch 2（子空间模式）修改

**文件**：`cal_mw_from_lambda.cpp`

流程变更：
1. `scatter_sub_matrix_to_local`：将缓存的 H₀_sub 从完整矩阵转为本地块
2. `calculate_delta_hcc_lcao`（本地块版本）：在本地块上累加 DeltaSpin 修正
3. `gather_sub_matrix_to_all`：将修正后的本地块聚合回完整矩阵
4. `diag_hegvd`：在完整矩阵上对角化（接口不变）

多了一次 scatter/gather，但消除了 P_I_sub 的 NPROC 倍冗余存储。

#### 2.7 诊断函数更新

**文件**：`lambda_loop.cpp`

所有诊断函数（`run_lambda_local_diagnostic`、`run_lambda_linear_scan`、`run_lambda_scan_diagnostic`、`run_trace_vs_dmr_diagnostic`）中的 `cal_PI_sub` 调用改为写入局部 `vector<vector<vector<>>>` 临时变量，然后手动 scatter 到 `lcao_PI_sub_save_` 的 map 结构中并提取对角线。`calculate_delta_hcc_lcao` 调用使用向后兼容的全量版本。

#### 2.8 工具函数可见性

以下原本 `static` 的函数改为外部可见（去除 `static`）以供 `cal_mw_from_lambda.cpp` 调用：

- `gather_sub_matrix_to_all`
- `scatter_sub_matrix_to_local`
- `extract_diagonal_from_local_block`

在 `cal_mw_from_lambda.cpp` 中通过 `namespace spinconstrain` 内的 `extern` 声明引用。

---

### Phase 2：消除 BFGS 循环内分配抖动（方案 3.5）

#### 2.9 成员缓冲区

**文件**：`spin_constrain.h`

新增成员变量：
```cpp
vector<complex<double>> h_sub_local_buf_;   // 本地块缓存
vector<complex<double>> h_tmp_buf_;         // 完整 H 矩阵临时空间
vector<complex<double>> s_tmp_buf_;         // 完整 S 矩阵临时空间
vector<complex<double>> vcc_buf_;           // 特征向量临时空间
vector<complex<double>> s_copy_buf_;        // S 副本（diag_hegvd 会修改 S）
vector<double> eigenvalues_buf_;            // 特征值临时空间
```

首次使用时按需分配（通过 size 检查），后续 BFGS 迭代中复用，避免每次迭代 `vector` 的 `malloc`/`free` 造成 RSS 虚高。

`free_lcao_subspace_cache()` 中同步清空这些缓冲区。

---

## 三、内存效果估算

以 52 原子、424 bands、72 k-points、16 进程（KPAR=1）为例：

| 项目 | 优化前 | 优化后 | 降幅 |
|------|--------|--------|------|
| `lcao_PI_sub_save_` | nk × nat × nbands² × 16 B = **10.7 GB** | nk × nat × (nbands²/16) × 16 B ≈ **0.67 GB** | **94%** |
| `lcao_PI_sub_diag_` | — | nk × nat × nbands × 8 B ≈ 0.019 GB | 新增（可忽略） |
| BFGS 临时 vector | 每次迭代 ~1.7 GB 重新分配 | 首次分配后复用 | RSS 不再虚高 |

**若 KPAR=4**（k-point 也分布）：`lcao_PI_sub_save_` 降至 ~0.17 GB，总降幅 98%。

---

## 四、尚未实施的优化方案

### 4.1 [方案 3.2] 子空间对角化改为分布式 `pzhegv`

**当前瓶颈**：Branch 2 子空间模式中，`diag_hegvd` 要求完整 nbands×nbands 矩阵（所有进程持有副本），因此需要 `gather_sub_matrix_to_all`。这引入了一次多余的 gather，且 S_sub 也需要每个进程持有完整副本。

**方案**：替换为 ScaLAPACK 的分布式广义特征值求解器 `pzhegv`，直接在 2D 块循环分布的本地块上求解。

**预期收益**：
- 消除 H_sub 的 gather 和 S_sub 的全量存储
- S_sub 本地块可按需存储（或与 H_sub 共享同一分配）
- 内存再降约 2 × nk × nbands² / NPROC = ~0.08 GB（16 进程场景下）

**难度**：中等。需要：
1. 将 `diag_hegvd`（LAPACK 的 `zhegv` 封装）替换为 `pzhegv`
2. 修改特征向量 `vcc` 的输出格式：从完整矩阵改为本地块
3. 下游的 `rotate_psi_subspace_lcao` 需要适配本地块格式的 `vcc_all`

**风险**：`pzhegv` 的数值稳定性与 `zhegv` 可能有微小差异，需要回归测试。

### 4.2 [方案 3.3] k-point 维度进程分布

**当前瓶颈**：pool 内所有进程持有相同的 k-point 子集数据，k-point 间无分布。

**方案**：将 pool 内进程按 k-point 维度分组，每个进程子组只负责部分 k-point 的缓存构建和使用。

**预期收益**：内存降至原来的 1/KPAR_EFF 倍。

**难度**：大。需要：
1. 设计 k-point 组内通信机制（含非本地 k-point 数据的跨进程获取）
2. 修改 `cal_mi_lcao_subspace`、`rotate_psi_subspace_lcao` 等函数的并行语义
3. 可能需要调整 MPI communicator 分层结构

**建议**：仅在单节点内存确实不足以容纳单个 k-point 缓存时考虑。对大多数实际体系，Phase 1 已足够。

### 4.3 [补充] Branch 2 subspace 模式中 S_sub 也改为本地块存储

**当前状态**：Phase 1 中 `lcao_sub_s_save` 仍为完整 nbands² 矩阵（所有进程持有副本），因为 `diag_hegvd` 需要完整的 S。

**方案**：
- 若实施 4.1（分布式 `pzhegv`），则 S_sub 自然也只需存本地块
- 若不实施 4.1，可单独将 S_sub 改为本地块 + 按 gather，省掉常驻内存

**预期收益**：nk × nbands² × 16 B / NPROC —— 对 424 bands / 16 进程约 0.08 GB，收益有限。

### 4.4 [补充] `DiagonalizationEngine` 内部缓存同步

**当前状态**：`SubspaceDiagonalizer` 和 `FirstOrderResponseEngine` 内部有自己的 `SubspaceCache`，存储格式为 `vector<vector<vector<>>>`（完整矩阵）。这套缓存与 `SpinConstrain` 的 `lcao_PI_sub_save_` 独立，存在同样的冗余问题。

**方案**：让 `DiagonalizationEngine` 直接使用 `SpinConstrain` 的本地块缓存，或将其内部 `SubspaceCache` 也改为分布式存储。

**难度**：中等。需要修改 `SubspaceCache::build` 和 `SubspaceCache::P_I_sub` 接口。

---

## 五、实施优先级建议

| 优先级 | 方案 | 收益 | 难度 |
|--------|------|------|------|
| ★★★ | 4.1 分布式 `pzhegv` | H_sub + S_sub 不再全量存储 | 中 |
| ★★ | 4.4 Engine 缓存同步 | 消除另一套冗余 | 中 |
| ★ | 4.3 S_sub 本地块化 | 收益有限 | 低 |
| ☆ | 4.2 k-point 分布 | 最极端场景才需要 | 大 |
