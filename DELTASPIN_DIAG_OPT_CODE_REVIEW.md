# DeltaSpin Fast Mode 对角化优化 — 代码变更详细Review

## 一、修改范围

涉及 12 个 commit（`9c495033d..f4548f424`），主要变更文件：

| 文件 | 变更类型 |
|------|----------|
| `lambda_loop.cpp` | 核心算法流程重构 |
| `cal_mw_from_lambda.cpp` | 三条执行路径 + 分布式P_I_sub缓存 + 混合精度 |
| `lcao_subspace.cpp` | P_I_sub计算改用folding_HR + 本地块存储 |
| `spin_constrain.h` | 新增成员变量（分布式缓存 + BFGS缓冲区 + 精度控制） |
| `spin_constrain.cpp` | 缓存释放逻辑 |
| `diagonalization_engine.cpp/h` | 新增诊断/引擎策略类 |
| `dspin_lcao.cpp` | lambda符号修正 + sc_strategy参数映射 |
| `esolver_ks_lcao.cpp` | 混合精度精度同步 + 诊断触发 |
| `cal_dm_psi.cpp/h` | 新增fp32 GEMM路径 |
| `density_matrix.cpp/h` | 新增float特化的func_exp_mul_dmk |
| `scalapack_connector.h` | 新增complex<float>和float的gemm声明 |

---

## 二、逐Commit代码变更分析

### 1. `9c495033d` feat(deltaspin): opt-in subspace acceleration for LCAO spin-constrained DFT

**核心变更**：将子空间/一阶加速改为 opt-in，全对角化恢复为默认。

**🔴 新增参数**：
- **`sc_acceleration_mode`**（`off` / `first_order` / `subspace`）：子空间加速模式开关，默认 `off`
- **`sc_acceleration_rms_thr`**（double）：RMS 误差低于此值时激活加速，默认 -1.0（不激活）

**`cal_mw_from_lambda()` 三条执行路径**：

```
i_step == -2 (Branch 1):
  全对角化 → C_ref, ε_ref
  calculate_lcao_sub_hs → H0_sub, S_sub
  cal_PI_sub → P_I_sub（完整nbands²矩阵，所有进程持有）
  迹公式计算 Mi(V=I)

acceleration_active && subspace_built (Branch 2):
  first_order: ε_new = ε_old + dl * diag(P_I_sub)，DMR路径算Mi
  subspace: H_sub = H0_sub + δH, 本地zhegvd, rotate psi, DMR路径算Mi

default (Branch 3):
  HSolverLCAO::solve 全对角化
```

**问题**：Branch 1 中 `cal_PI_sub` 使用 `B_I_data`（onsite 重叠积分缓存）路径，后续证实存在邻域截断偏差。

---

### 2. `b2c9d0498` fix(deltaspin): fix LCAO subspace acceleration bugs and add diagnostic

**三个关键修复**：

**① 特征值偏移符号错误**（`cal_mw_from_lambda.cpp` 一阶路径）：

```cpp
// 修复前：
ekb_new = ekb_old + spin_sign * delta_lambda * P_diag;
// 修复后：
ekb_new = ekb_old - spin_sign * delta_lambda * P_diag;
```

符号应与 `calculate_delta_hcc_lcao` 一致：H_DS = +λ·σ，所以 ε 对 λ 的响应是 -spin_sign·Δλ·P_diag。

**② 跨SCF迭代的子空间缓存过期**：

在 `run_lambda_loop()` 入口处新增状态重置：
```cpp
this->acceleration_active_ = false;
this->acceleration_subspace_built_ = false;
this->free_lcao_subspace_cache();
```

问题：一个SCF迭代构建的子空间被复用到后续迭代时，H/S/C已变化，Mi会冻结不收敛。

**③ 收敛后缺少最终全对角化**：

在 `update_psi_charge()` 中，如果加速曾经激活过，在收敛λ处执行一次全对角化：
```cpp
if (this->acceleration_active_ && this->acceleration_subspace_built_
    && this->nspin_ == 2)
{
    hamilt::DeltaSpin<...>* dspin = ...;
    dspin->update_lambda();
    hsolver_t.solve(hamilt_t, psi_t[0], ...);
    // + calculate_weights + calEBand
}
this->pelec->psiToRho(*psi_t);
```

**新增诊断工具**：`run_trace_vs_dmr_diagnostic()` 对比全对角化、DMR子空间、迹公式三种Mi随Δλ的变化。

---

### 3. `e61fc8f37` fix(deltaspin): replace cal_PI_sub with calculate_PI_sub_from_hr

**核心问题**：旧 `cal_PI_sub` 通过 `B_I_data`（onsite积分缓存）构建 D_I，再算 P_I_sub = D_I†·D_I。这条路径的邻域截断与 `pre_hr` 不完全一致，导致 BCC Fe AFM 中原子1的Mi差 4.66 μB（符号反转）。

**新方案** `calculate_PI_sub_from_hr` 与 H_sub 计算完全对称：

```
folding_HR(pre_hr[iat]) → P_I(k)（2D-block 分布）
pzgemm('N','N'): temp = P_I(k) × C(k)
pzgemm('C','N'): P_I_sub = C†(k) × temp
gather_sub_matrix_to_all: 汇总到所有进程
```

**关键保证**：P_I(k) 来自 `pre_hr[iat]`，与 `cal_moment` 使用同一数据源，迹公式与DMR路径精确一致。

**新增访问器**：
```cpp
dspin_op->get_pre_hr(iat)
dspin_op->get_constraint_atom_list()
```

---

### 4. `5735ea673` feat(deltaspin): simplify execution strategy with sc_strategy parameter

**🔴 新增参数**：

| sc_strategy | 映射的 sc_acceleration_mode | 映射的 sc_acceleration_rms_thr | 行为 |
|-------------|---------------------------|-------------------------------|------|
| `normal`（默认） | subspace | 1e-2 | 全对角化起步，RMS < 阈值后自动切换 |
| `fast` | subspace | 1e10 | 第一步即激活 subspace |
| `accuracy` | off | -1.0 | 始终全对角化 |

**🔴 默认值变更**：
- **`nsc`**: 100 → **5**
- **`sc_scf_thr`**: 1e-3 → **10**
- **`sc_scf_thr_mode`**: "threshold" → **"immediate"**

**映射逻辑**（`deltaspin_lcao.cpp`）：
```cpp
if (inp.sc_strategy == "fast") {
    accel_mode = "subspace";
    accel_rms_thr = 1e10;
} else if (inp.sc_strategy == "accuracy") {
    accel_mode = "off";
    accel_rms_thr = -1.0;
} else { // normal
    accel_mode = "subspace";
    if (!user_overrode_rms_thr) accel_rms_thr = 1e-2;
}
```

用户显式设置 `sc_acceleration_mode` 或 `sc_acceleration_rms_thr` 可覆盖策略默认值。

---

### 5. `c9f3ff876` feat(deltaspin): implement LCAO subspace mixed-precision GEMM optimization

复用 `GintPrecisionController` 机制，在Gint切换到fp32时同步将子空间内GEMM降为fp32。

**变更清单**：

**`scalapack_connector.h`**：新增 `pcgemm_` 和 `psgemm_` 声明，以及 `complex<float>` 和 `float` 的 `gemm` 重载（全描述符和简化描述符两个API）。

**`SpinConstrain` 类**：
```cpp
ModuleGint::GintPrecision subspace_exec_precision_;
void set_subspace_exec_precision(ModuleGint::GintPrecision p);
```

**`esolver_ks_lcao.cpp`**：
- `iter_init()` 中：Gint精度切换后将精度同步到SpinConstrain
- `iter_finish()` 中：同上（处理迭代结束后的精度变化）

**`cal_mw_from_lambda.cpp` Branch 2 一阶路径**：
```cpp
if (this->subspace_exec_precision_ == ModuleGint::GintPrecision::fp32)
    elecstate::cal_dm_psi_mixed(this->ParaV, this->pelec->wg, *psi_t, *this->dm_);
else
    elecstate::cal_dm_psi(this->ParaV, this->pelec->wg, *psi_t, *this->dm_);
```

**`cal_dm_psi.cpp`**：新增 `cal_dm_psi_mixed`、`psiMulPsiMpiMixed`、`psiMulPsiMixed`：
```cpp
// 流程：cast-down → fp32 GEMM → cast-up
std::vector<std::complex<float>> psi1_f(nloc_psi);
cast_down_to_float(psi1.get_pointer(), psi1_f.data(), nloc_psi);
// ... pcgemm_ ...
cast_up_to_double(dm_f.data(), dm_out, nloc_dm);
```

**`density_matrix.cpp`**：新增 `func_exp_mul_dmk<float>` 特化：
```cpp
template <>
void DensityMatrix_Tools::func_exp_mul_dmk<float>(
    const std::complex<double> kphase,
    const std::vector<std::complex<double>>& DMK_mat_trans,
    float* target_DMR_mat)
{
    const float kr = static_cast<float>(kphase.real());
    const float ki = static_cast<float>(kphase.imag());
    for (std::size_t i = 0; i < mat_size; i++)
        target_DMR_mat[i] += kr * static_cast<float>(DMK_mat_trans[i].real())
                            - ki * static_cast<float>(DMK_mat_trans[i].imag());
}
```

---

### 6. `f4548f424` Perf(deltaspin): eliminate gather redundancy and reduce BFGS allocation jitter

**6.1 数据结构变更**（`spin_constrain.h`）：

| 成员 | 旧类型 | 新类型 |
|------|--------|--------|
| `lcao_PI_sub_save_` | `vector<vector<vector<complex<double>>>>` | `vector<map<int, vector<complex<double>>>>` |
| `lcao_PI_sub_diag_` | （不存在） | `vector<map<int, vector<double>>>` |

变更要点：
- 第二维从 `vector[nat]` 改为 `map<int, vector>`，只存约束原子
- 内层vector从 `nbands²` 改为 `nrow × ncol_bands`（本地2D块），大小降至 ~1/NPROC
- 新增 `lcao_PI_sub_diag_` 缓存对角线元素，供一阶路径直接使用

**6.2 `calculate_PI_sub_from_hr` 修改**（`lcao_subspace.cpp`）：

旧路径：`pzgemm → PI_sub_local[nloc_eij] → gather → PI_sub_full[nn]`

新路径：`pzgemm → PI_sub_local[nloc_eij] → memcpy到输出指针`

- 删除 `gather_sub_matrix_to_all` 调用
- fp32路径同理：`cast_up_to_double → memcpy`
- 串行路径不变

**6.3 `calculate_delta_hcc_lcao` 重载**（`lcao_subspace.cpp`）：

**重载1（本地块版本）**：
```cpp
void calculate_delta_hcc_lcao(
    complex<double>* h_sub_local,                        // [nrow*ncol_bands]
    const map<int, vector<complex<double>>>& PI_sub_local,
    const Vector3<double>* lambda, int nbands, int ik,
    bool full_update, const Parallel_Orbitals* ParaV);
```
- 遍历 `map`：`for (const auto& [iat, pi_local] : PI_sub_local)`
- axpy循环长度从 `nbands*nbands` 变为 `nloc_eij = nrow * ncol_bands`
- npol=2目前仅使用coeff0（z分量）

**重载2（全量矩阵版本，向后兼容）**：
```cpp
void calculate_delta_hcc_lcao(
    complex<double>* h_sub,                               // [nbands²]
    const vector<vector<complex<double>>>& PI_sub,
    const Vector3<double>* lambda, int nbands, int ik,
    bool full_update);
```
- 保留给诊断函数和DiagonalizationEngine使用

**6.4 Branch 1（子空间缓存构建）修改**（`cal_mw_from_lambda.cpp`）：

```cpp
this->lcao_PI_sub_save_.resize(nk);
this->lcao_PI_sub_diag_.resize(nk);

for (int ik = 0; ik < nk; ik++) {
    for (int iat = 0; iat < nat; iat++) {
        if (!dspin_op->get_constraint_atom_list()[iat]) continue;

        this->lcao_PI_sub_save_[ik][iat].resize(nloc_eij, {0,0});
        this->calculate_PI_sub_from_hr(..., this->lcao_PI_sub_save_[ik][iat].data(), ...);

        this->lcao_PI_sub_diag_[ik][iat].resize(nbands, 0.0);
        extract_diagonal_from_local_block(
            this->lcao_PI_sub_save_[ik][iat].data(), this->ParaV,
            this->lcao_PI_sub_diag_[ik][iat].data(), nbands);
    }
}
```

不再调用 `resize(nat)` 和 `clear()`；map 天然只存约束原子。

**6.5 Branch 2（一阶路径）修改**：

```cpp
// 旧代码：
for (int iat = 0; iat < nat; iat++) {
    if (lcao_PI_sub_save_[ik][iat].empty()) continue;
    double p_diag = lcao_PI_sub_save_[ik][iat][ib + ib*nbands].real();
    ...
}
// 新代码：
for (const auto& [iat, diag] : this->lcao_PI_sub_diag_[ik]) {
    double p_diag = diag[ib];
    ...
}
```

直接从 double 缓存读取对角线，避免 complex→real 转换和 `nbands*nbands` 索引。

**6.6 Branch 2（子空间路径）修改**：

旧路径：直接对完整的 `h_tmp[nn]` 数组做 `daxpy`

新路径：
```
1. scatter_sub_matrix_to_local: H0_sub(k) → h_sub_local_buf_[nloc_eij]
2. calculate_delta_hcc_lcao(本地块版本): h_sub_local += Σ_I Δλ_I · P_I_sub_local
3. gather_sub_matrix_to_all: h_sub_local → h_tmp_buf_[nn]（恢复全量用于对角化）
4. 同理 S_sub 复制到 s_tmp_buf_[nn]
5. diag_hegvd(h_tmp_buf_, s_copy_buf_, ...)
```

多了一次 scatter/gather，但消除了 P_I_sub 的 NPROC 倍内存冗余。

**6.7 BFGS临时变量提升为成员**（`spin_constrain.h`）：

```cpp
vector<complex<double>> h_sub_local_buf_;
vector<complex<double>> h_tmp_buf_;
vector<complex<double>> s_tmp_buf_;
vector<complex<double>> vcc_buf_;
vector<complex<double>> s_copy_buf_;
vector<double> eigenvalues_buf_;
```

首次使用时按需分配（通过 size 检查），后续 BFGS 迭代中复用。`free_lcao_subspace_cache()` 中同步清空。

**6.8 工具函数可见性变更**：

以下原本 `static` 的函数改为外部可见并移入 `namespace spinconstrain`：
- `gather_sub_matrix_to_all`
- `scatter_sub_matrix_to_local`
- `extract_diagonal_from_local_block`

在 `cal_mw_from_lambda.cpp` 中通过 `extern` 声明引用。

---

### 7. 符号一致性修复（`de96b2753`, `eb1a8b62f`, `eb9a0e410`）

**`dspin_lcao.cpp` `cal_coeff_lambda`**：
- nspin=2：系数从 `(-λ, +λ)` 修正为 `(+λ, -λ)`
- nspin=4：系数从 `(-λz, -λx-iλy, -λx+iλy, +λz)` 修正为 `(+λz, +λx+iλy, +λx-iλy, -λz)`

**`lcao_subspace.cpp` `calculate_delta_hcc_lcao`**（nspin=2路径）：
```cpp
// 修正前：
const complex<double> coeff(-effective_lambda[iat][2] * spin_sign, 0.0);
// 修正后：
const complex<double> coeff(effective_lambda[iat][2] * spin_sign, 0.0);
```

**E_lambda公式修正**：
```cpp
// 修正前：escon -= λ·Mi
// 修正后：escon -= λ·(Mi - Mtarget)
```

**pzgemm维度修正**：
```cpp
// 修正前：使用 this->ParaV->nrow（本地行数）
// 修正后：使用 this->ParaV->get_global_row_size()（全局NLOCAL）
```

---

## 三、Phase 3a/3b 优化实现在代码中的位置

### Phase 3a: 合并初始化与缓存构建

**位置**：`lambda_loop.cpp:515-556`

```cpp
if (i_step == -1) {
    const bool fast_mode = (this->sc_acceleration_mode_ != "off")
                           && (this->sc_acceleration_rms_thr_ >= 1e8);
    if (fast_mode) {
        this->acceleration_active_ = true;
        this->acceleration_subspace_built_ = false;
        this->free_lcao_subspace_cache();
        this->cal_mw_from_lambda(-2);       // 一次全对角化 + 缓存构建 + 初始Mi
        this->acceleration_subspace_built_ = true;
    } else {
        this->cal_mw_from_lambda(i_step);   // 仅计算初始Mi
    }
    spin = this->Mi_;
    ...
}
```

消除的冗余：
- 原始：i_step=-1 全对角化 → i_step=0 全对角化（加速未激活）→ RMS<threshold 触发 i_step=-2 全对角化+缓存
- 优化后：i_step=-1 一次全部完成，后续 BFGS 步直接走 Branch 2

### Phase 3b: 子空间旋转替代收敛后全对角化

**位置**：`cal_mw_from_lambda.cpp:1072-1137`

```cpp
if (this->acceleration_active_ && this->acceleration_subspace_built_
    && this->nspin_ == 2 && this->lcao_subspace_initialized_)
{
    // 1. scatter H0_sub → 本地块
    // 2. 累加修正 h_sub_local += Σ_I (λ_conv - λ_ref)_I · P_I_sub_local
    // 3. gather 回全量矩阵
    // 4. diag_hegvd → V, ε
    // 5. rotate_psi: ψ_new(k) = ψ_ref(k) · V(k)
    // 6. calculate_weights(ε)
    // 7. psiToRho(ψ_new)
}
```

nspin=4 或缓存未初始化时回退到全对角化。

---

## 四、内存优化效果估算

以52原子、424 bands、72 k-points、16进程（KPAR=1）为例：

| 项目 | 优化前 | 优化后 | 降幅 |
|------|--------|--------|------|
| `lcao_PI_sub_save_` | nk × nat × nbands² × 16B = **10.7 GB** | nk × nat_constrained × (nbands²/16) × 16B ≈ **0.67 GB** | **94%** |
| `lcao_PI_sub_diag_` | — | nk × nat_constrained × nbands × 8B ≈ 0.019 GB | 新增（可忽略） |
| BFGS临时vector | 每次迭代 ~1.7 GB重新分配 | 首次分配后复用 | RSS不再虚高 |

若 KPAR=4：`lcao_PI_sub_save_` 降至 ~0.17 GB，总降幅 98%。

---

## 五、尚未实施的优化方向

| 优先级 | 方案 | 收益 | 难度 |
|--------|------|------|------|
| ★★★ | Branch 2 子空间对角化改为分布式 `pzhegv` | 消除 H_sub 的 gather 和 S_sub 的全量存储 | 中 |
| ★★ | DiagonalizationEngine 内部缓存使用本地块格式 | 消除另一套冗余 | 中 |
| ★ | S_sub 也改为本地块存储 | nk × nbands² × 16B / NPROC，收益有限 | 低 |
| ☆ | 子空间缓存跨SCF迭代复用 | 可能降至0次全对角化/SCF | 高（风险大） |
