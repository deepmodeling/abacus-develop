# DeltaSpin 内存优化方案：用 D_I 路径替代 pre_hr / P_I_sub 密集缓存

## 1. 问题背景

DeltaSpin 模块在 LCAO 基组下存在两个大内存数据结构：

- **`pre_hr`**：实空间投影算符 HContainer，per 约束原子，存储 `<φ_μ|α_lm><α_lm|φ_ν,R>`
- **`lcao_PI_sub_save_`**：子空间投影矩阵密集缓存，per k 点 per 约束原子，存储 `P_I_sub(k) = C†(k)·P_I(k)·C(k)`

当前优化代码在 BFGS Branch 2 阶段同时保留两者，内存开销大。需要分析是否可以只保留一个，或用更轻量的替代方案。

## 2. pre_hr 的所有使用场景

| 使用场景 | 函数 | 代码位置 | Branch 2 是否需要 |
|---|---|---|---|
| (a) contributeHR() | 构建全空间 H(k) 修正 `dHR += λ·pre_hr` | `dspin_lcao.cpp:140-174` | **不需要** — Branch 2 用 calculate_delta_hcc_lcao 替代 |
| (b) cal_moment() | 计算 Mi = Σ dmR·pre_hr | `dspin_lcao.cpp:471-515` | **需要** — 被 cal_mi_lcao() 调用 |
| (c) calculate_PI_sub_from_hr() | 从 pre_hr 生成 P_I_sub(k) | `lcao_subspace.cpp:376-466` | **不需要** — 仅在 cache build (Branch 1) 时调用 |
| (d) cal_mi_lcao_subspace() | 被 cal_moment() 间接使用 | `cal_mw.cpp:90,113` | **需要** — subspace Mi 内部调用 full-space path |

## 3. lcao_PI_sub_save_ 的所有使用场景

| 使用场景 | 函数 | 代码位置 | 频率 |
|---|---|---|---|
| (a) calculate_delta_hcc_lcao() | H_sub += Σ_I Δλ_I·P_I_sub_local | `lcao_subspace.cpp:497-520` | 每个 BFGS 步，每个 k 点 |
| (b) Branch 2 first-order | lcao_PI_sub_diag_[ik][iat][ib] 一阶修正 | `cal_mw_from_lambda.cpp:759-761` | 每个 BFGS 步，每个 k 点 |
| (c) update_psi_charge (优化路径) | 构建 H_sub → 对角化 → 旋转 psi | `cal_mw_from_lambda.cpp:1112` | 每次收敛后一次 |

## 4. 数据依赖关系

```
B_I_data (原始重叠 <φ|α>, O(N) 稀疏)
    ↓ cal_pre_HR() 中 cal_HR_IJR() 做乘积
pre_hr[iat] (实空间投影 Σ_lm |α_lm><α_lm|, O(N) 稀疏)
    ↓ folding_HR + C†·P_I(k)·C (pzgemm)
lcao_PI_sub_save_[ik][iat] (子空间投影 P_I_sub, O(N²) 密集)
```

lcao_PI_sub_save_ 是 pre_hr 的**压缩表示**：从 NLOCAL×NLOCAL 压缩到 nbands×nbands。但压缩后失去了稀疏性，变成密集存储，扩展性从 O(N) 退化为 O(N²)。

## 5. 内存扩展性分析

| 数据结构 | 存储方式 | 扩展性 | 52原子 72k 16进程 | 1000原子 1k 16进程 | 10000原子 1k 16进程 |
|---|---|---|---|---|---|
| pre_hr[iat] | 实空间稀疏 HContainer | O(N) | ~数百MB | ~数GB | ~数十GB |
| lcao_PI_sub_save_[ik][iat] | 密集 nbands×nbands | O(N²) | ~126MB | **~数十GB** | **不可行** |
| B_I_data[iat] | 原始重叠系数 | O(N) | ~10-50MB | ~数百MB | ~数GB |
| **D_I(k)** (新方案) | r × nbands 矩阵 | O(N) | ~25MB/rank | ~108MB/rank | ~1GB/rank |

以 r = max_l² ≈ 16（投影轨道数）估算 D_I(k) 大小：
- 单原子单 k 点：16 × nbands × 16B
- 52原子 424bands 72k 16rank：~25 MB/rank
- 1000原子 1000bands 1k 16rank：~16 MB/rank
- 10000原子 5000bands 1k 16rank：~400 MB/rank

**结论：D_I(k) 的 O(N) 扩展性使大尺度体系可行，而 P_I_sub 的 O(N²) 在约 500 原子时已不可行。**

## 6. 关键发现：代码中已有两条计算路径

### 路径 A (folding_HR) — 当前主要路径

```
pre_hr[iat] → folding_HR → P_I(k) [NLOCAL×NLOCAL, 稀疏]
    → C† · P_I(k) · C [pzgemm] → P_I_sub(k) [nbands×nbands, 密集]
```

- 需要 pre_hr
- 中间步骤产生密集的 P_I_sub，O(N²) 存储

### 路径 B (overlap) — 已存在于 cal_PI_sub()

```
B_I_data[iat] → exp(ik·R) · B_I(k) [r × NLOCAL, 稀疏行]
    → B_I(k) · C(k) = D_I(k) [r × nbands]
    → D_I†(k) · D_I(k) = P_I_sub(k) [nbands×nbands]
```

- 只需 B_I_data，**不需要 pre_hr**
- D_I(k) 维度 r × nbands，远小于 P_I_sub(k)

代码位置：`dspin_lcao.cpp:564-646` (`cal_PI_sub()` 函数)

## 7. 用 D_I 替代 P_I_sub 的核心推导

### 7.1 H_sub 修正

当前公式：
```
H_sub(k) += Σ_I Δλ_I · P_I_sub(k)       // P_I_sub: nbands×nbands，需存储
```

用 D_I 重写：
```
D_I(k) = B_I(k) · C(k)                   // r × nbands，on-the-fly
H_sub(k) += Δλ_I · D_I†(k) · D_I(k)     // (nbands×r)·(r×nbands) → nbands×nbands
```

ScaLAPACK 实现：
```
// D_I_local: 每个进程从 B_I_data 构建 r × nrow_local 的局部块
// 先做 D_I · C_local (pzgemm, r × nrow_local · nrow_local × ncol_bands)
// 得到 D_I_local: r × ncol_bands
// 再做 D_I† · D_I (pzgemm, ncol_bands × r · r × ncol_bands, accumulate into Eij block)
```

复杂度：O(r × nbands²) per atom，r ≈ 16，**远优于** O(N × nbands²)

### 7.2 Mi 计算 — 更妙，甚至不需要 P_I_sub

当前 DMR 路径（O(N) 但慢）：
```
M_I = Σ_{μν,R} dmR[μ,ν,R] · pre_hr[I][μ,ν,R]    // 稀疏，需要 pre_hr
```

用 D_I 的迹公式（O(nbands) 且极快）：

**nspin=2 collinear：**
```
M_I,z = Σ_k w_k · spin_sign(k) · Σ_a f_{k,a} · Σ_{lm} |D_I(k)[lm,a]|²
```

**nspin=4 non-collinear：**
```
occ[0] = Σ_a f_a · (|D_I_up|^2 + |D_I_dn|^2)    // Tr(rho · I)
occ[1] = Σ_a f_a · (D_I_up† · D_I_dn)            // Tr(rho · σ_x) + i·Tr(rho · σ_y)
occ[2] = Σ_a f_a · (D_I_dn† · D_I_up)            // conjugate of occ[1]
occ[3] = Σ_a f_a · (|D_I_up|^2 - |D_I_dn|^2)     // Tr(rho · σ_z)

M_I,x = (occ[1] + occ[2]).real()
M_I,y = (occ[1] - occ[2]).imag()
M_I,z = (occ[0] - occ[3]).real()
```

展开证明（nspin=2）：
```
M_I = Σ_k w_k · s_k · Tr[P_I_sub(k) · ρ_sub(k)]
    = Σ_k w_k · s_k · Σ_{ij} [D_I†·D_I]_{ij} · [V†·f·V]_{ji}
    = Σ_k w_k · s_k · Σ_a f_{k,a} · Σ_{lm} D_I*[lm,a] · D_I[lm,a]
    = Σ_k w_k · s_k · Σ_a f_{k,a} · Σ_{lm} |D_I[lm,a]|²
```

计算复杂度：**O(r × nbands) per atom per k-point** — 极其廉价，无需构建 DMR。

## 8. 单 k 点的特殊优化

Γ-only 时：

1. **B_I(Γ)** 的构建：无需相位因子，直接从 B_I_data 累加。B_I_data 本身只涉及近邻原子（~20-50 个），提取对应的 C 行时也是稀疏操作。

2. **P_I(Γ) 仍然稀疏**：从 B_I(Γ) 构建出的 P_I(Γ) 继承了实空间的截断半径稀疏性。但**根本不需要走 P_I 路径**——直接用 D_I(Γ) on-the-fly 计算即可。

3. **DMR 路径同样 O(N)**：`cal_moment(pre_hr, dmR)` 对稀疏的 dmR 和 pre_hr 做乘积，复杂度 O(N)。但在 subspace 路径里我们改用 D_I 迹公式，**连 DMR 都不需要构建**。

## 9. 完整方案

### 9.1 新增函数

```cpp
// 计算 D_I(k) = B_I(k) · C(k)，返回 r × nbands 的矩阵
// 输入: B_I_data, C(k) (2D block distributed), kvec_d
// 输出: D_I_local (r × ncol_bands，每个进程持有列块)
std::vector<std::complex<double>>
calculate_D_I_k(
    int iat,
    const std::vector<BI_AdjacentData>& bi_data,
    int r,  // B_I_nproj[iat]
    const std::complex<double>* psi_k,
    const Parallel_Orbitals* ParaV,
    const ModuleBase::Vector3<double>& kvec_d,
    int nbands);

// 用 D_I 计算 Mi（迹公式），替代 cal_moment
std::vector<double>
cal_mi_from_D_I_trace(
    const std::vector<std::vector<std::vector<std::complex<double>>>>& D_I_all,
    // [ik][iat][r * nbands]
    const Parallel_Orbitals* ParaV,
    int nbands, int nk, int npol, int nspin,
    const K_Vectors& kv,
    const elecstate::ElecState* pelec);
```

### 9.2 用 D_I 重写 calculate_delta_hcc_lcao

```cpp
// 当前: H_sub_local += Σ_I Δλ_I · P_I_sub_local(k)
// 改为: 构建 D_I_local(k), 然后 H_sub_local += Δλ_I · D_I†·D_I
void calculate_delta_hcc_lcao_from_D_I(
    std::complex<double>* h_sub_local,
    const std::vector<std::vector<std::complex<double>>>& D_I_local_map,
    // map[iat] -> D_I_local (r × ncol_bands)
    const ModuleBase::Vector3<double>* lambda,
    int nbands, int ik, bool full_update,
    const Parallel_Orbitals* ParaV);
```

### 9.3 释放时序

```
SCF iteration N:
  ┌─ lambda_loop() 开始
  │   acceleration_active_ = false
  │   free_lcao_subspace_cache()
  │   
  │   BFGS 步 (Branch 3): contributeHR() 懒初始化 pre_hr ✓
  │   ... 若干 BFGS 步，用 pre_hr 计算 Mi ...
  │   
  │   RMS < threshold → 激活加速
  │   cal_mw_from_lambda(-2): Branch 1 → cache build
  │     ├── 全空间对角化 (使用 pre_hr)
  │     ├── 用 B_I_data + C(k) 构建 D_I(k)，缓存到 D_I_cache_[ik][iat]
  │     ├── ★ 释放 pre_hr: delete pre_hr[]; initialized = false
  │     └── ★ 释放 lcao_PI_sub_save_ (不再需要)
  │   
  │   BFGS 步 (Branch 2):
  │     ├── 用 D_I_cache 做相乘计算 H_sub 修正 ✓
  │     ├── 用 D_I 迹公式计算 Mi ✓
  │     └── 不调用 contributeHR，不调用 cal_moment ✓
  │   
  │   收敛 → update_psi_charge (优化路径) → 不需要 pre_hr ✓
  └─ lambda_loop() 结束
  
  cal_mi_lcao_wrapper() → cal_mi_lcao() → cal_moment()
    → 懒重建 pre_hr (调用 cal_pre_HR + snap) ← 短暂内存峰值
  
  iter_finish(): 电荷混合、势能更新等
  下一 SCF 迭代开始...
```

### 9.4 额外优化：release B_I_data

Cache build 完成后，B_I_data 也不再需要（D_I_cache 已包含所有信息）。可额外节省 ~10-50 MB。

但注意：若后续 SCF 迭代需要重建 pre_hr（cal_mi_lcao_wrapper 会触发），则需要 B_I_data 来重新调用 snap + cal_HR_IJR。因此：
- 若 cal_mi_lcao_wrapper 也改用 D_I 迹公式 → 可释放 B_I_data
- 否则保留 B_I_data 作为 pre_hr 重建的依赖

### 9.5 内存节省汇总

| 场景 | 释放前 | 释放后 | 节省 |
|---|---|---|---|
| 52 原子 72k 16rank | ~数百MB (pre_hr) + 126MB (P_I_sub) | ~25MB/rank (D_I cache) | 数百 MB |
| 1000 原子 1k 16rank | ~数GB (pre_hr) + 数十GB (P_I_sub) | ~16MB/rank (D_I cache) | **数十 GB** |
| 10000 原子 1k 16rank | ~数十GB (pre_hr) + 不可行 (P_I_sub) | ~400MB/rank (D_I cache) | **从不可行变为可行** |

## 10. 实现优先级

| 步骤 | 内容 | 难度 | 收益 |
|---|---|---|---|
| 1 | 实现 cal_mi_from_D_I_trace() | 低 | 消除 Branch 2 对 pre_hr 的运行时依赖 |
| 2 | 用 D_I 重写 calculate_delta_hcc_lcao | 中 | 消除 lcao_PI_sub_save_ 存储 |
| 3 | Cache build 后释放 pre_hr | 低 | 数百 MB ~ 数 GB 内存峰值下降 |
| 4 | cal_mi_lcao_wrapper 也改用 D_I trace | 低 | 消除每次 SCF 的 pre_hr 重建开销 |
| 5 | 释放 B_I_data | 低 | ~10-50 MB 额外节省 |
| 6 | nspin=4 Pauli 分解的 D_I trace | 中 | 支持非共线计算 |

## 11. 风险评估

| 风险 | 严重性 | 缓解方案 |
|---|---|---|
| D_I 迹公式与 DMR 路径数值差异 | 低 | 数学等价，可先并行运行对比验证 |
| pzgemm 构建 D_I 的 ScaLAPACK 分布正确性 | 中 | 参考 cal_PI_sub() 中已有的 MPI_Allreduce + zgemm 实现 |
| cal_mi_lcao_wrapper 重建 pre_hr 开销 | 低 | 每次 SCF 仅一次，~秒级；长期方案改用 D_I trace |
| nspin=4 需额外 Pauli 分解 | 中 | 当前仅支持 nspin=2；nspin=4 用 occ[4] 公式 |
| 代码路径切换的鲁棒性 | 中 | 保留原有 pre_hr + P_I_sub 路径作为 fallback |

## 12. 总结

**pre_hr 和 lcao_PI_sub_save_ 都可以释放**，用 D_I(k) 路径替代：

1. **D_I(k) = B_I(k) · C(k)**，维度 r × nbands（r ≈ 16），内存 O(N)
2. **H_sub 修正**：Δλ · D_I†·D_I 代替 Δλ · P_I_sub
3. **Mi 计算**：Σ |D_I|² 迹公式代替 DMR + cal_moment
4. **扩展性**：从 O(N²) 回到 O(N)，万原子规模可行
5. **代码基础**：cal_PI_sub() 已有完整的 B_I → D_I → P_I_sub 实现，只需将 D_I 作为主要路径而非中间步骤
