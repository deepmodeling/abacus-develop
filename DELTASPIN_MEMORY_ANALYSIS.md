# `lcao_PI_sub_save_` 并行冗余存储分析与优化方案

## 一、当前存储模式

### 1.1 数据结构

```cpp
// spin_constrain.h line 1010
std::vector<std::vector<std::complex<double>>> lcao_PI_sub_save_;
//     [ik]           [iat]              [ib * nbands + jb]
//   k-point索引    原子索引              nbands × nbands 矩阵元素
```

### 1.2 存储位置与时机

**Branch 1（i_step==-2）**: `cal_mw_from_lambda.cpp:624-658`
```cpp
this->lcao_PI_sub_save_.resize(nk);                              // nk 个 k-point
for (int ik = 0; ik < nk; ik++) {
    this->lcao_PI_sub_save_[ik].resize(nat);                     // nat 个原子
    for (int iat = 0; iat < nat; iat++) {
        if (!constrained[iat]) { clear(); continue; }
        this->lcao_PI_sub_save_[ik][iat].resize(nn, {0,0});     // nbands² 个元素
        this->calculate_PI_sub_from_hr(...);                     // 两轮 pzgemm + gather
    }
}
```

### 1.3 使用位置与时机

**Branch 2（每次 BFGS 迭代）**: `cal_mw_from_lambda.cpp:781-801`
```cpp
for (int ik = 0; ik < nk; ik++) {
    std::memcpy(h_tmp, lcao_sub_h_save + ik*nn, ...);
    std::memcpy(s_tmp, lcao_sub_s_save + ik*nn, ...);
    this->calculate_delta_hcc_lcao(h_tmp, this->lcao_PI_sub_save_[ik], ...);
    // 修改 h_tmp: H_sub(λ) = H₀_sub + (λ - λ_ref) · P_sub
    diag_hegvd(h_tmp, s_tmp, eigenvalues, vcc);                  // 子空间对角化
}
```

---

## 二、冗余分析

### 2.1 k-point 并行冗余

ABACUS 使用 k-point 并行（`KPAR` 参数），将 k-points 分配到不同的 MPI pool 中。
每个 pool 内的进程拥有**相同的 k-point 子集**，每个 pool 处理 `nk_pool = nkstot / KPAR` 个 k-points。

但代码中：
```cpp
int nk = psi_t->get_nk();  // 返回的是 LOCAL nk（当前 pool 的 k-point 数）
```

**问题**：虽然 `nk` 已经是 pool 内的本地 k-point 数，但 `lcao_PI_sub_save_` 在**每个 pool 的每个进程**上都完整存储了该 pool 所有 k-points 的 P_I_sub 矩阵。在 pool 内部的进程间，这是完全冗余的。

### 2.2 子空间矩阵全量冗余

更严重的问题是 `gather_sub_matrix_to_all`（lcao_subspace.cpp:61-89）：

```cpp
static void gather_sub_matrix_to_all(...) {
    // 每个 MPI 进程都持有局部 P_I_sub_local（2D 块循环分布）
    // 然后通过 reduce_pool 聚合到完整的 nbands × nbands 矩阵
    // 结果：每个进程都有一个完整的 P_I_sub 副本
    Parallel_Reduce::reduce_pool(full_data, nbands * nbands);
}
```

**每个进程都存储了完整的 nbands × nbands 矩阵**，而实际上只有在子空间对角化时才需要完整矩阵。在 LCAO 的 2D 块循环分布中，每个进程只需要存储自己的本地块。

### 2.3 原子维度冗余

```cpp
this->lcao_PI_sub_save_[ik].resize(nat);  // ALL atoms, not just constrained ones
for (int iat = 0; iat < nat; iat++) {
    if (!constrained[iat]) {
        this->lcao_PI_sub_save_[ik][iat].clear();  // 清空但仍占 vector 开销
        continue;
    }
    ...
}
```

即使只需要约束原子的 P_I_sub，也分配了所有 nat 个原子的 vector 容器（虽然非约束原子的 inner vector 被 clear，但 outer vector 的 nat 个 vector 元素本身有 24 bytes/table 开销）。

### 2.4 量化计算

对于典型 52 原子、424 bands、72 k-points、16 MPI 进程（KPAR=1）：

| 数据 | 每进程存储 | 说明 |
|------|-----------|------|
| P_I_sub 完整矩阵 × ik × iat | 72 × 52 × 424² × 16 B = **10.7 GB** | 全量冗余 |
| 若只存约束原子 | 72 × 52 × 424² × 16 B = **10.7 GB** | 用户选了全原子约束 |
| 若保持 2D 分布不 gather | 72 × 52 × (424²/16) × 16 B ≈ **0.67 GB** | 16进程分布 |
| 若再加 k-point 分布(KPAR=4) | 18 × 52 × (424²/16) × 16 B ≈ **0.17 GB** | k+distribution |

**结论：当前实现比最优分布式存储多占用 16-60 倍内存。**

---

## 三、优化方案

### 3.1 [Critical] 不要 gather，保持 2D 块循环分布存储

**当前路径**（每次 calculate_PI_sub_from_hr）：
```
pre_hr → folding_HR → PI_k_local[nloc] → pzgemm → PI_sub_local[nloc_eij] → gather → PI_sub_full[nn]
                                                                               ↑ 冗余！
```

**优化后路径**：
```
pre_hr → folding_HR → PI_k_local[nloc] → pzgemm → PI_sub_local[nloc_eij]  ← 直接缓存这个
```

将 `lcao_PI_sub_save_` 从 `vector<vector<vector<complex<double>>>>` 改为**每个 iat 存储本地 2D 块**：

```cpp
// 新数据结构
struct PI_sub_cache {
    int nbands;
    std::vector<std::complex<double>> local_block; // size = ParaV->nrow * ParaV->ncol_bands
};
// 按 [ik][iat_constrained] 存储
std::vector<std::vector<PI_sub_cache>> lcao_PI_sub_save_;
```

对应修改 Branch 2 中的 `calculate_delta_hcc_lcao`：不再对完整的 nn 矩阵做 `daxpy`，而是对 local_block 做 `pdaxpy`（Scalapack 的分布式向量加法）。

**内存节省**: 从 `nk × nat × nbands² × 16 B` 降至 `nk × nat × nloc_eij_local × 16 B`，约为原来的 `1/NPROC_IN_POOL`。

### 3.2 [Critical] 子空间对角化也改为分布式

当前 `diag_hegvd` 要求每个进程持有完整的 nn 矩阵。但 Scalapack 有 `pdsyev`/`pzhegv` 可以直接在 2D 分布矩阵上做对角化。

**优化**: 将子空间对角化从 `diag_hegvd(full_matrix)` 改为 `pzhegv(local_block)`，避免 gather 和 scatter。

### 3.3 [High] k-point 并行化存储

当前 `nk = psi_t->get_nk()` 已经是 pool 本地 k-point 数。但 pool 内的多进程仍然重复存储同一组 k-points 的数据。

**方案**: 如果进一步将 k-point 维度分布到 pool 内不同进程组（类似 PW 基组的做法），则每个进程只需存储 `nk / NPROC_IN_POOL` 个 k-points 的缓存。

但这需要改动较大，建议作为第二阶段优化。

### 3.4 [Medium] 只存储约束原子的 P_I_sub

将 `lcao_PI_sub_save_[ik]` 从 `vector<complex<double>>[nat]` 改为 `map<int, PI_sub_cache>`，只存约束原子的：

```cpp
std::vector<abacus::map<int, PI_sub_cache>> lcao_PI_sub_save_;
//     [ik]     iat → local_block
```

查找时用 `find(iat)` 代替 `[iat]`，空间从 O(nat) 降到 O(nat_constrained)。

### 3.5 [Medium] Branch 2 临时 vector 改为成员变量复用

当前每次 BFGS 迭代分配：
```cpp
std::vector<std::vector<std::complex<double>>> vcc_all(nk);           // nk × nn
std::vector<std::complex<double>> h_tmp(nn), s_tmp(nn);              // 2 × nn
std::vector<std::complex<double>> vcc(nn), eigenvalues(nbands);       // nn + nbands
```

改为成员变量，首次分配后复用：

```cpp
// 成员变量
std::vector<std::vector<std::complex<double>>> vcc_all_buf_;
std::vector<std::complex<double>> h_tmp_buf_, s_tmp_buf_, vcc_buf_;
bool subspace_buf_initialized_ = false;

// Branch 2 中
if (!subspace_buf_initialized_) {
    vcc_all_buf_.resize(nk);
    for (auto& v : vcc_all_buf_) v.resize(nn);
    h_tmp_buf_.resize(nn); s_tmp_buf_.resize(nn); vcc_buf_.resize(nn);
    subspace_buf_initialized_ = true;
}
```

**效果**: 消除 BFGS 循环内的内存抖动，RSS 不再因 malloc/free 虚高。

---

## 四、优化后的内存估算

对 52 原子、424 bands、72 k-points、16 进程：

| 方案 | 每进程内存 | 降幅 |
|------|-----------|------|
| 当前（全量 gather 存储） | ~10.7 GB | - |
| 3.1 + 3.4（2D分布 + 只存约束原子） | ~0.67 GB | **94%** |
| 3.1 + 3.3 + 3.4（+k分布） | ~0.17 GB | **98%** |
| 3.1-3.5 全部实施 | ~0.17 GB + 0 抖动 | **98%** |

即使只实施 3.1（2D 分布存储不 gather），也能将内存降至当前水平的 ~6%，使 52 原子体系的每进程内存从 10.7 GB 降至 ~0.67 GB，完全避免 OOM。

---

## 五、实施优先级

1. **3.1 + 3.4**: 最关键，消除 gather 全量冗余 + 只存约束原子。预计改动 `calculate_PI_sub_from_hr`、`calculate_delta_hcc_lcao`、`lcao_PI_sub_save_` 数据结构。
2. **3.5**: 简单改动，消除分配抖动。
3. **3.2**: 中等难度，需将 `diag_hegvd` 替换为分布式 `pzhegv`。
4. **3.3**: 较大改动，需调整 k-point 数据分布，作为第二阶段。
