# cal_mw.cpp / cal_mw_from_lambda.cpp / cal_h_lambda.cpp 代码对比分析

## 对比版本
- **当前版本**: `feat/dftu-pw-port-v2` 分支 (`/root/abacus-dftu-pw-port/`)
- **基准版本**: `tmp` (zdy) 分支 (`/root/abacus-zdy-tmp/`), ABACUS v3.7.1

---

## 一、cal_mw.cpp

### 1.1 各自逻辑与算法

#### 当前版本

包含两个函数：

| 函数 | 用途 | 算法 |
|------|------|------|
| `cal_mi_lcao` | LCAO 基组下计算原子磁矩 | 通过 `dm_->get_DMR_pointer(1)` 获取密度矩阵实空间部分，调用 `p_operator->cal_moment()` 计算每个原子的磁矩投影 |
| `cal_mi_pw` | PW 基组下计算原子磁矩 | 遍历所有 k 点，用 `OnsiteProjector::tabulate_atomic()` + `overlap_proj_psi()` 计算 `<Psi|α><α|Psi>` 投影系数 (becp)，然后调用 `accumulate_Mi_from_becp()` 从 becp 累加磁矩 |

**`cal_mi_pw` 核心流程**：
```
zero_Mi() → for each ik: tabulate_atomic → overlap_proj_psi → accumulate_Mi_from_becp → reduce_allpool
```

#### tmp 版本

包含三个函数：

| 函数 | 用途 | 算法 |
|------|------|------|
| `cal_MW_k` | 按 k 点计算 MW（使用 DMK 密度矩阵） | 通过 S 矩阵与密度矩阵的 ScaLAPACK `pzgemm_` 矩阵乘，再调用 `collect_MW` 收集 |
| `cal_MW` | LCAO 基组下计算原子磁矩（主函数） | 同当前版本 `cal_mi_lcao`，使用 DMR 密度矩阵实空间部分 + `cal_moment()` |
| `cal_Mi_pw` | PW 基组下计算原子磁矩 | 遍历 k 点，在函数内部直接展开 band×atom 循环计算磁矩 |

### 1.2 差异点

#### 差异 1：函数命名与组织

| 项目 | 当前版本 | tmp 版本 |
|------|----------|----------|
| LCAO 磁矩计算 | `cal_mi_lcao` | `cal_MW` |
| PW 磁矩计算 | `cal_mi_pw` | `cal_Mi_pw` |
| k 点 MW 矩阵计算 | **已删除** | `cal_MW_k`（被注释掉的备选路径） |

当前版本 `cal_MW_k` 路径已被移除（在 tmp 中也已注释），改为纯 DMR 路径。

#### 差异 2：PW 磁矩计算实现方式 — **最关键差异**

**tmp 版本** (`cal_Mi_pw`): 在函数内直接展开磁矩计算循环：

```cpp
// npol==2 (nspin=4) 情况
for (int ib = 0; ib < nbands; ib++) {
    const double weight = this->pelec->wg(ik, ib);
    for (int iat = 0; iat < Mi_.size(); iat++) {
        std::complex<double> occ[4] = {...};
        for (int ih = 0; ih < nh; ih++) {
            const int index = ib*2*nkb + begin_ih + ih;
            occ[0] += conj(becp[index]) * becp[index];
            occ[1] += conj(becp[index]) * becp[index + nkb];
            occ[2] += conj(becp[index + nkb]) * becp[index];
            occ[3] += conj(becp[index + nkb]) * becp[index + nkb];
        }
        Mi_[iat].z += weight * (occ[0] - occ[3]).real();
        Mi_[iat].x += weight * (occ[1] + occ[2]).real();
        Mi_[iat].y += weight * (occ[1] - occ[2]).imag();
    }
}
```

**当前版本** (`cal_mi_pw`): 提取为公共函数 `accumulate_Mi_from_becp()`:

```cpp
// cal_mi_pw 中只负责获取 becp，然后调用：
this->accumulate_Mi_from_becp(becp, nkb, nbands, npol, ik, &this->pelec->wg(ik, 0), &onsite_p->get_nh(0));
```

`accumulate_Mi_from_becp()` 内部实现与 tmp 版本数学等价，但增加了：
- 对 npol==2 使用 `pauli_to_moment(occ, weight)` 辅助函数转换
- 对 npol==1 统一通过 `get_spin_sign(ik)` 获取自旋符号

#### 差异 3：`cal_mi_pw` 中 npol==1 路径的 sign 处理

| 项目 | 当前版本 | tmp 版本 |
|------|----------|----------|
| npol==1 符号来源 | 统一调用 `get_spin_sign(ik)` | 在 `cal_Mi_pw` 内联计算：`const int sign = (is == 0) ? 1 : -1` |
| npol==2 符号 | 无 sign（始终为 1） | 无 sign |

当前版本 `get_spin_sign()`:
```cpp
int get_spin_sign(int ik) const {
    if (this->npol_ == 2) return 1;
    return (this->pelec->klist->isk[ik] == 0) ? 1 : -1;
}
```

tmp 版本在 `cal_Mi_pw` 中 npol==1 路径：
```cpp
const int is = this->pelec->klist->isk[ik];
const int sign = (is == 0) ? 1 : -1;
```

两者**数学等价**。

#### 差异 4：Parallel_Reduce 的 KPAR 参数

| 项目 | 当前版本 | tmp 版本 |
|------|----------|----------|
| KPAR 来源 | `PARAM.inp.kpar` | `GlobalV::KPAR` |

```cpp
// 当前版本
Parallel_Reduce::reduce_double_allpool(PARAM.inp.kpar, GlobalV::NPROC_IN_POOL, ...);

// tmp 版本
Parallel_Reduce::reduce_double_allpool(GlobalV::KPAR, GlobalV::NPROC_IN_POOL, ...);
```

这是参数访问方式的改变，不影响逻辑。

#### 差异 5：get_npol 的访问方式

| 项目 | 当前版本 | tmp 版本 |
|------|----------|----------|
| 获取 npol | `psi_t->get_npol()` | `psi_t->npol`（直接访问成员） |

---

## 二、cal_mw_from_lambda.cpp

### 2.1 各自逻辑与算法

该文件包含 DeltaSpin lambda 迭代的核心计算流程：

1. **`calculate_delta_hcc`**: 根据 delta_lambda 和 becp 计算哈密顿量修正 delta_H
2. **`cal_mw_from_lambda`**: lambda 迭代主函数，更新波函数后重新计算磁矩
3. **`update_psi_charge`**: 用修正后的哈密顿量更新波函数和电荷密度

#### `calculate_delta_hcc` 算法

对每个原子 iat，构建 Pauli 矩阵系数：
```
nspin=4 (npol=2):
  coeff0 = (λz, 0)      coeff1 = (λx, λy)
  coeff2 = (λx, -λy)    coeff3 = (-λz, 0)
  
  ps[index]    += coeff0 * becp1 + coeff2 * becp2
  ps[index+nkb] += coeff1 * becp1 + coeff3 * becp2

nspin=2 (npol=1):
  coeff0 = λz * spin_sign
  ps[index] += coeff0 * becp
```

然后用 GEMM 更新哈密顿量: `H += becp† × ps`

#### `cal_mw_from_lambda` 算法

```
if LCAO:
  update_lambda() → hsolver.solve() → calculate_weights() → cal_dm_psi() → cal_DMR() → cal_mi_lcao()
else (PW):
  if 首次调用:
    updateHk() → cal_hs_subspace() → 保存 H_k, S_k, becp
  对所有 ik:
    拷贝 H_k, S_k → calculate_delta_hcc() 添加 delta_H
    → diag_responce() 得到新 becp_tmp 和 ekb
  calculate_weights() → accumulate_Mi_from_becp() → reduce_allpool
```

#### `update_psi_charge` 算法

```
if LCAO:
  psiToRho()
else (PW):
  对所有 ik:
    拷贝 H_k, S_k → calculate_delta_hcc() → diag_subspace_psi()
  if pw_solve:
    HSolver.solve()  // 全 PW 空间对角化
  else:
    calculate_weights() → psiToRho()  // 仅更新权重和电荷
```

### 2.2 差异点

#### 差异 1：`calculate_delta_hcc` 的 sign 参数传递

| 项目 | 当前版本 | tmp 版本 |
|------|----------|----------|
| 函数签名 | `(..., const int ik)` | `(..., const int sign)` |
| npol==1 系数 | `delta_lambda[iat][2] * this->get_spin_sign(ik)` | `delta_lambda[iat][2] * sign` |
| 调用方传参 | 传 `ik` | 传 `sign = (isk==0 ? 1 : -1)` |

**这是重要差异**：当前版本将 sign 的获取延迟到 `calculate_delta_hcc` 内部通过 `get_spin_sign(ik)` 计算，而 tmp 版本在调用方计算好 sign 后传入。两者数学等价，但当前版本的封装更好。

#### 差异 2：`calculate_delta_hcc` 的 Timer 和上下文对象

| 项目 | 当前版本 | tmp 版本 |
|------|----------|----------|
| Timer | 有 `ModuleBase::TITLE` 和 `timer::start/end` | 无 Timer |
| GPU 内存操作 | 使用默认 context 的 `resize_memory_op` | 显式创建 `ctx` 和 `cpu_ctx` 对象 |
| gemm_op 调用 | `ModuleBase::gemm_op<...>()` | `hsolver::gemm_op<...>(ctx, ...)` |

#### 差异 3：`cal_mw_from_lambda` LCAO 路径的 HSolver 调用

| 项目 | 当前版本 | tmp 版本 |
|------|----------|----------|
| HSolver 类型 | `hsolver::HSolverLCAO<TK> hsolver_t(ParaV, ks_solver)` | `static_cast<hsolver::HSolver<TK, CPU>*>(this->phsol)` |
| solve 调用 | `hsolver_t.solve(hamilt_t, psi_t[0], pelec, *dm_, *charge, nspin_, true)` | `hsolver_t->solve(hamilt_t, psi_t[0], pelec, KS_SOLVER, true)` |
| DM 计算 | `elecstate::cal_dm_psi(ParaV, wg, *psi_t, *dm_)` | 仅对特定 ks_solver: `"genelpa"|"scalapack_gvx"|"lapack"|"cg_in_lcao"` |
| 磁矩计算 | `cal_mi_lcao(i_step)` | `cal_MW(i_step)` |

**关键差异**：tmp 版本对 LCAO 的 `cal_dm_psi` 调用有条件判断（仅特定对角化器），而当前版本无条件调用。

#### 差异 4：`cal_mw_from_lambda` PW 路径的磁矩计算 — **最关键差异**

**tmp 版本**：磁矩计算直接内联展开，分 npol==2 和 npol==1 两个分支：

```cpp
// npol==2 分支 (第 339-369 行)
for (int ik = 0; ik < nk; ik++) {
    const std::complex<double>* becp = &becp_tmp[ik * size_becp];
    for (int ib = 0; ib < nbands; ib++) {
        const double weight = this->pelec->wg(ik, ib);
        for (int iat = 0; iat < Mi_.size(); iat++) {
            std::complex<double> occ[4] = {...};
            for (int ih = 0; ih < nh; ih++) {
                const int index = ib * npol * nkb + begin_ih + ih;
                // 计算 occ[0..3]
            }
            Mi_[iat].x += weight * (occ[1] + occ[2]).real();
            Mi_[iat].y += weight * (occ[1] - occ[2]).imag();
            Mi_[iat].z += weight * (occ[0] - occ[3]).real();
        }
    }
}

// npol==1 分支 (第 370-397 行)
for (int ik = 0; ik < nk; ik++) {
    const int sign = this->pelec->klist->isk[ik] == 0 ? 1 : -1;
    // ... 类似结构，Mi_[iat].z += weight * occ * sign
}
```

**当前版本**：统一调用 `accumulate_Mi_from_becp()`:

```cpp
for (int ik = 0; ik < nk; ik++) {
    const std::complex<double>* becp = &becp_tmp[ik * size_becp];
    this->accumulate_Mi_from_becp(becp, nkb, nbands, this->npol_, ik,
        &this->pelec->wg(ik, 0), nh_iat);
}
```

**数学等价性分析**：

对于 npol==2，`accumulate_Mi_from_becp` 中的索引计算：
```cpp
const int index = ib * 2 * nkb + begin_ih + ih;  // 即 ib * npol * nkb
```

tmp 版本内联代码中的索引计算：
```cpp
const int index = ib * npol * nkb + begin_ih + ih;  // npol==2 时为 ib * 2 * nkb
```

两者**完全等价**。

对于 npol==1，两者都使用 `get_spin_sign(ik)` 或等效的 sign 计算，也**等价**。

#### 差异 5：`update_psi_charge` 的函数拆分

| 项目 | 当前版本 | tmp 版本 |
|------|----------|----------|
| 函数结构 | `update_psi_charge` → `update_psi_charge_pw_cpu` / `update_psi_charge_pw_gpu` | 单一 `update_psi_charge` 函数内联处理 CPU/GPU |
| 代码行数 | 每个函数 ~80 行，职责清晰 | 单一函数 ~130 行 |

#### 差异 6：`update_psi_charge` PW 路径的全 PW 对角化

**当前版本**：
```cpp
hsolver::HSolverPW<std::complex<double>, DEVICE_CPU> hsolver_pw_obj(
    this->pw_wfc_, PARAM.inp.calculation, ..., PARAM.inp.use_k_continuity);
hsolver_pw_obj.solve(hamilt_t, psi_t[0], this->pelec, this->pelec->ekb.c, ...);
```

**tmp 版本**：
```cpp
hsolver::HSolver<std::complex<double>, DEVICE_CPU>* hsolver_t = 
    static_cast<hsolver::HSolver<...>*>(this->phsol);
hsolver_t->solve(hamilt_t, psi_t[0], this->pelec, this->KS_SOLVER, false);
```

**关键差异**：当前版本直接构造 `HSolverPW` 对象，而 tmp 版本通过成员指针 `phsol` 获取 HSolver。

#### 差异 7：`update_psi_charge` 非 pw_solve 路径

| 项目 | 当前版本 | tmp 版本 |
|------|----------|----------|
| npol==1 权重更新 | `calculate_weights()` + `psiToRho()` | 直接 `psiToRho()` |
| CPU 路径 | `reinterpret_cast<ElecStatePW*>(pelec)->psiToRho(*psi_t)` | CPU: `pelec->psiToRho(*psi_t)`, GPU: `reinterpret_cast<...>` |

当前版本在 CPU 和 GPU 路径都使用 `ElecStatePW::psiToRho`，而 tmp 版本 CPU 路径用基类 `pelec->psiToRho`，GPU 路径用 `reinterpret_cast`。

#### 差异 8：内存管理方式

| 操作 | 当前版本 | tmp 版本 |
|------|----------|----------|
| 分配 sub_h_save 等 | `new std::complex<double>[...]` | `new std::complex<double>[...]` |
| GPU 分配 sub_h_save 等 | `base_device::memory::resize_memory_op<..., GPU>()(ptr, size)` | `base_device::memory::resize_memory_op<..., GPU>()(ctx, ptr, size)` |
| 释放 sub_h_save 等 | `delete[]` / `delete_memory_op` | `delete[]` / `delete_memory_op` |

主要差异是 GPU 内存操作是否需要显式 context 参数。

---

## 三、cal_h_lambda.cpp

### 3.1 各自逻辑与算法

该文件仅用于 LCAO 基组，计算 DeltaSpin 约束项对哈密顿量的贡献。

**算法**：对每对原子 (iat1, iat2) 的轨道 (iwt1, iwt2)，计算：
```
lambda = (lambda_[iat1] + lambda_[iat2]) / 2.0
h_lambda[mu, nu] = -Sloc2[...] * lambda 分量 (按 Pauli 矩阵展开)
```

对于 nspin==2：
```
h_lambda[icc] = (isk==0) ? -Sloc2[icc]*λz : -Sloc2[icc]*(-λz)
```

对于 nspin==4，按 iwt1%2 和 iwt2%2 分四种情况，对应 Pauli 矩阵 σx, σy, σz 的 2×2 块。

### 3.2 差异点

#### 差异 1：nspin==4 时 Sloc2 的索引偏移 — **最显著差异**

**column_major=true** 分支：

| 条件 | 当前版本 | tmp 版本 |
|------|----------|----------|
| iwt1%2==0, iwt2%2==1 | `Sloc2[icc + pv->nrow]` | `Sloc2[icc + 1]` |
| iwt1%2==1, iwt2%2==0 | `Sloc2[icc + 1]` | `Sloc2[icc - 1]` |
| iwt1%2==1, iwt2%2==1 | `Sloc2[icc + 1 + pv->nrow]` | `Sloc2[icc]` |

**column_major=false** 分支：

| 条件 | 当前版本 | tmp 版本 |
|------|----------|----------|
| iwt1%2==0, iwt2%2==1 | `Sloc2[icc + 1]` | `Sloc2[icc - 1]` |
| iwt1%2==1, iwt2%2==0 | `Sloc2[icc + pv->ncol]` | `Sloc2[icc + 1]` |
| iwt1%2==1, iwt2%2==1 | `Sloc2[icc + 1 + pv->ncol]` | `Sloc2[icc]` |

**这是一个重大逻辑差异！** 两种版本对 Pauli 矩阵非对角元在 Sloc2 中的索引方式完全不同。当前版本使用 `+nrow`/`+ncol` 偏移（跨行/列），而 tmp 版本使用 `+1`/`-1` 偏移（相邻元素）。

#### 差异 2：Timer 调用不一致

当前版本 `cal_h_lambda.cpp:104`:
```cpp
ModuleBase::timer::start("SpinConstrain", "cal_h_lambda");  // 应该是 end
```

tmp 版本 `cal_h_lambda.cpp:104`:
```cpp
ModuleBase::timer::tick("SpinConstrain", "cal_h_lambda");
```

当前版本在函数末尾错误地使用了 `timer::start` 而非 `timer::end`/`tick`，这是一个 bug。

---

## 四、四个关键差异的深度调研

### 4.1 cal_mw_from_lambda.cpp: LCAO 路径 HSolver 调用方式变化

#### tmp 版本做法
```cpp
// spin_constrain.h: void* phsol = nullptr;
// esolver_ks_lcao.cpp 初始化:
this->phsol = new hsolver::HSolverLCAO<TK>(&(this->ParaV));
this->phsol->method = GlobalV::KS_SOLVER;

// cal_mw_from_lambda.cpp 调用:
hsolver::HSolver<std::complex<double>, base_device::DEVICE_CPU>* hsolver_t
    = static_cast<hsolver::HSolver<std::complex<double>, base_device::DEVICE_CPU>*>(this->phsol);
hsolver_t->solve(hamilt_t, psi_t[0], this->pelec, this->KS_SOLVER, true);
```

`HSolverLCAO` 继承自 `HSolver<T, Device>` 基类，`solve` 签名为：
```cpp
virtual void solve(hamilt::Hamilt<T>* pHamilt,
                   psi::Psi<T>& psi,
                   elecstate::ElecState* pes,
                   const std::string method_in,    // 运行时传入对角化方法
                   const bool skip_charge) override;
```

#### 当前版本做法
```cpp
// spin_constrain.h 中 phsol 成员已被删除
// cal_mw_from_lambda.cpp 调用:
hsolver::HSolverLCAO<std::complex<double>> hsolver_t(this->ParaV, PARAM.inp.ks_solver);
hsolver_t.solve(hamilt_t, psi_t[0], this->pelec, *this->dm_, *this->pelec->charge, this->nspin_, true);
```

当前版本 `HSolverLCAO` **不再继承** `HSolver` 基类，`solve` 签名变为：
```cpp
void solve(hamilt::Hamilt<TK>* pHamilt,
           psi::Psi<TK>& psi,
           elecstate::ElecState* pes,
           elecstate::DensityMatrix<TK, double>& dm,  // ← 新增
           Charge &chr,                                // ← 新增
           const int nspin,                           // ← 新增
           const bool skip_charge);
```

#### 行为差异分析

| 维度 | tmp 版本 | 当前版本 |
|------|----------|----------|
| solve 内部处理 DM | `skip_charge=true` 时直接 return，不调用 `psiToRho` 或 `cal_dm_psi` | `skip_charge=true` 时仍然调用 `cal_dm_psi(ParaV, pes->wg, psi, dm)` + `dm.cal_DMR()` |
| 外部 cal_dm_psi | 有条件调用（仅特定 ks_solver）| `cal_mw_from_lambda` 中无条件再调一次 `cal_dm_psi` |
| 结果 | DM 计算正确 | **DM 被计算了两次**（HSolver 内部一次 + cal_mw_from_lambda 外部一次），但数学上等价（覆盖写入） |

**结论：行为等价，但当前版本有冗余的 cal_dm_psi 调用。** 两个版本在 lambda 迭代中得到的密度矩阵结果是相同的。

---

### 4.2 cal_mw_from_lambda.cpp: `update_psi_charge` 全 PW 求解方式变化

#### tmp 版本做法
```cpp
hsolver::HSolver<std::complex<double>, base_device::DEVICE_CPU>* hsolver_t
    = static_cast<hsolver::HSolver<std::complex<double>, base_device::DEVICE_CPU>*>(this->phsol);
hsolver_t->solve(hamilt_t, psi_t[0], this->pelec, this->KS_SOLVER, false);
```

`HSolverPW` 构造函数简单：`HSolverPW(pw_wfc, wavefunc* pwf)`，持有 `wavefunc*` 成员用于 PAW/PSI 初始化。
`solve` 签名：
```cpp
void solve(hamilt::Hamilt<T, Device>* pHamilt,
           psi::Psi<T, Device>& psi,
           elecstate::ElecState* pes,
           const std::string method_in,
           const bool skip_charge) override;
```
内部使用临时 `eigenvalues` vector 存储能量，最后 cast 到 `pes->ekb.c`。

#### 当前版本做法
```cpp
hsolver::HSolverPW<std::complex<double>, base_device::DEVICE_CPU> hsolver_pw_obj(
    this->pw_wfc_,
    PARAM.inp.calculation,
    PARAM.inp.basis_type,
    PARAM.inp.ks_solver,
    PARAM.globalv.use_uspp,
    PARAM.inp.nspin,
    hsolver::DiagoIterAssist<std::complex<double>>::SCF_ITER,
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_NMAX,
    hsolver::DiagoIterAssist<std::complex<double>>::PW_DIAG_THR,
    hsolver::DiagoIterAssist<std::complex<double>>::need_subspace,
    PARAM.inp.use_k_continuity);

hsolver_pw_obj.solve(hamilt_t, psi_t[0], this->pelec, this->pelec->ekb.c,
    GlobalV::RANK_IN_POOL, GlobalV::NPROC_IN_POOL, false, this->tpiba, this->get_nat());
```

当前版本 `HSolverPW` **不需要 `wavefunc*`**，但需要 11 个构造参数。`solve` 签名变为：
```cpp
void solve(hamilt::Hamilt<T, Device>* pHamilt,
           psi::Psi<T, Device>& psi,
           elecstate::ElecState* pes,
           double* out_eigenvalues,         // ← 直接输出到外部数组
           const int rank_in_pool_in,       // ← 新增
           const int nproc_in_pool_in,      // ← 新增
           const bool skip_charge,
           const double tpiba,              // ← 新增
           const int nat);                  // ← 新增
```

#### 行为差异分析

| 维度 | tmp 版本 | 当前版本 |
|------|----------|----------|
| solver 生命周期 | `phsol` 由 ESolver 堆分配，跨 SCF step 复用 | 栈上临时构造，每次调用创建新对象 |
| eigenvalue 存储 | 内部临时 vector，最后 cast 到 ekb | 直接传入 `pelec->ekb.c` 指针 |
| 对角化参数 | 从 solver 成员读取（method, diag_thr, nmax 等） | 通过构造函数传入，使用 DiagoIterAssist 常量 |
| wavefunc | 持有 `wavefunc*` 成员 | 不需要，PW 信息从 `pw_wfc_` 获取 |

**结论：行为等价。** 两种方式的 PW 对角化算法相同（都是调用相同的 `davcio`/`cg` 等子程序），只是参数传递和内存管理方式不同。当前版本的 `skip_charge=false` 参数传递正确，`psiToRho` 会被调用。

---

### 4.3 cal_mw_from_lambda.cpp: 非 pw_solve 路径 psiToRho 调用差异

#### 代码对比

| 路径 | tmp 版本 | 当前版本 |
|------|----------|----------|
| PW CPU | `this->pelec->psiToRho(*psi_t)`（基类指针虚函数调用） | `reinterpret_cast<ElecStatePW<..., CPU>*>(this->pelec)->psiToRho(*psi_t)` |
| PW GPU | `reinterpret_cast<ElecStatePW<..., GPU>*>(this->pelec)->psiToRho(*psi_t)` | 相同 |
| 外部 prepare | 无 | `calculate_weights()` 在外部已调用 |

#### psiToRho 实现差异

| 维度 | tmp 版本 `ElecStatePW::psiToRho` | 当前版本 `ElecStatePW::psiToRho` |
|------|------|------|
| `calculate_weights()` | 内部调用 (L88) | **不调用** |
| `calEBand()` | 内部调用 (L90) | **不调用** |
| `init_rho_data()` | 仅 `!init_rho` 时调用 | 每次无条件调用 |
| `add_usrho()` | 仅 `use_uspp` 时调用 | 无条件调用 |

#### 行为分析

1. **reinterpret_cast vs 基类指针**: `psiToRho` 是 virtual 函数，C++ 虚函数 dispatch 机制保证两种调用方式**完全等价**。`ElecState::psiToRho()` 基类实现为空（直接 return），实际执行的都是 `ElecStatePW::psiToRho`。

2. **calculate_weights**: tmp 版本在 `psiToRho` 内部调用；当前版本 `psiToRho` 不调用，但**调用方在之前已经调用了** `calculate_weights`：
   ```cpp
   // cal_mw_from_lambda.cpp, update_psi_charge 的 subspace 对角化后:
   // pw_solve=false 路径 (当前版本):
   elecstate::calculate_weights(pelec->ekb, pelec->wg, ...);  // L228
   reinterpret_cast<ElecStatePW*>(pelec)->psiToRho(*psi_t);    // L235
   ```
   tmp 版本 `HSolver::solve(skip_charge=false)` 内部不调用 `calculate_weights`，然后 `psiToRho` 内部调用。两者最终都计算了权重。

3. **calEBand**: 当前版本缺少 `calEBand()` 调用。tmp 版本 `psiToRho` 内部调用，当前版本需要在外部补上。

**结论：psiToRho 的调用方式等价，但当前版本在非 pw_solve 路径缺少 `calEBand()` 调用。** 这可能导致 band energy 输出不完整，但不影响电荷密度和 SCF 收敛。

---

### 4.4 cal_h_lambda.cpp: nspin==4 时 Sloc2 索引偏移

#### Sloc2 来源

`cal_h_lambda` 接收的 `Sloc2` 参数是 k 空间的 overlap 矩阵 Sk（通过 `hsk->get_sk()` 获取），是 2D 块分布的 ScaLAPACK 矩阵。

#### 调用链调研

**tmp 版本**:
```
OperatorScLambda::contributeHk(ik)
  → sc.cal_h_lambda(&h_lambda[0], this->hsk->get_sk(), IS_COLUMN_MAJOR_KS_SOLVER(), isk[ik])
  → hk[irc] += h_lambda[irc]
```
`sc_lambda_lcao.cpp:25` — 在 LCAO SCF 主循环中，通过 `contributeHk` 调用 `cal_h_lambda` 将 DeltaSpin 约束项加到 Hk 上。

**当前版本**:
- `cal_h_lambda` 函数存在但**未被任何生产代码调用**。搜索结果显示仅在 `spin_constrain.h` 声明和 `cal_h_lambda.cpp` 定义中出现，无调用点。
- 当前版本 LCAO 路径的 DeltaSpin 哈密顿量贡献通过 `DeltaSpin<OperatorLCAO>::contributeHR()` 在实空间计算，框架负责傅里叶变换到 k 空间。
- PW 路径通过 `calculate_delta_hcc()` 基于 becp 投影器计算 delta_H，**不使用 cal_h_lambda**。

#### 索引差异详情

tmp 版本有单元测试 `cal_h_lambda_test.cpp` 验证了正确性（nrow=2, ncol=2, npol=2, nspin=4, lambda=(1,1,1)）：
```
期望 column_major 结果: h_lambda = [-1, -1+i, -1-i, +1]
```
这对应于 `-λ·σ = [[-λz, -(λx+iλy)], [-(λx-iλy), +λz]]`。

| 模式 | 条件 | tmp 版本索引 | 当前版本索引 | 哪个正确 |
|------|------|-------------|-------------|----------|
| column_major | (0,1) S↑↓ | `icc+1` | `icc+nrow` | **tmp 正确** |
| column_major | (1,0) S↓↑ | `icc-1` | `icc+1` | **tmp 正确** |
| column_major | (1,1) S↓↓ | `icc` | `icc+1+nrow` | **tmp 正确** |
| row_major | (0,1) S↑↓ | `icc-1` | `icc+1` | **当前正确** |
| row_major | (1,0) S↓↑ | `icc+1` | `icc+ncol` | **当前正确** |
| row_major | (1,1) S↓↓ | `icc` | `icc+1+ncol` | **当前正确** |

#### 结论

1. **Column Major**: tmp 版本索引正确（经测试验证），当前版本错误。
2. **Row Major**: 当前版本索引正确，tmp 版本错误。
3. **实际影响**: 该函数在当前版本的 PW DeltaSpin 计算中**不被调用**，因此不影响 PW 路径的正确性。
4. **LCAO 路径**: 如果未来需要启用 LCAO DeltaSpin，需要修复 `cal_h_lambda` 的 column_major 索引（当前版本错误）。

**总结：cal_h_lambda 的 Sloc2 索引差异不影响当前 PW 路径，因为该函数未被调用。但代码本身存在 bug，应修复以防未来使用。**

---

## 五、差异汇总与结论

### 5.1 四个重点差异的结论

| # | 差异 | 行为是否等价 | 严重程度 | 是否需要修复 |
|---|------|-------------|----------|-------------|
| 1 | LCAO HSolver 调用方式 | **等价**（DM 被重复计算但不影响结果） | 低 | 可选：移除 cal_mw_from_lambda 中冗余的 cal_dm_psi 调用 |
| 2 | PW 全求解 HSolver 构造 | **等价**（算法相同，参数传递方式不同） | 无 | 否 |
| 3 | 非 pw_solve 路径 psiToRho | **基本等价**（虚函数 dispatch），但缺 calEBand | 低 | 建议：在外部补 calEBand 调用 |
| 4 | cal_h_lambda Sloc2 索引 | **不等价**（column_major 当前版本错误），但**函数未被调用** | 中 | 建议：修复 column_major 索引以防未来使用 |

### 5.2 对 PW DeltaSpin 调试的指导

当前项目的核心问题是 PW 基组下 DeltaSpin 导致 ~1000 eV 能量偏差。基于以上调研：

1. **四个重点差异均不是 PW 能量偏差的直接原因**：
   - 差异 1、2、3 行为等价
   - 差异 4 的函数在 PW 路径中不被调用

2. **PW 路径的实际差异在于** `cal_mw.cpp` 中 `cal_mi_pw` 和 `cal_mw_from_lambda.cpp` 中 PW 子路径的磁矩计算：
   - 当前版本提取为 `accumulate_Mi_from_becp()` 公共函数
   - tmp 版本内联展开
   - **经索引分析，两者数学等价**

3. **其他可能的偏差来源**（需继续排查）：
   - `update_psi_charge_pw_cpu/gpu` 中 HSolverPW 的构造参数（diag_thr, nmax 等）与 tmp 版本默认值是否不同
   - `psiToRho` 中 `init_rho_data` 的调用策略变化（每次无条件 vs 仅首次）
   - `add_usrho` 的调用条件变化
   - SCF 重启逻辑和电荷混合参数的变化
