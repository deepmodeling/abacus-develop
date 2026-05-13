# 24_LCAO_DS_S2_Z 测试失败根因分析

## 问题现象

1. **内循环 lambda 约束失效**：Lambda 内循环跑满 100 步仍未收敛（RMS ~9.7），delta_lambda 的 z 分量持续增大（从 3e-05 → 16），lambda 值暴增至 216/169 eV/uB，远超合理范围。

2. **GE9/GE10 DRHO=0 导致计算异常退出**：
   - GE8: DRHO=9.3457e-04, ETOT=-6.75040239e+03
   - GE9: DRHO=0.0000e+00, ETOT=-2.54404582e+03 ← 能量突变！
   - GE10: DRHO=0.0000e+00, ETOT=-2.17240744e+03 ← 继续突变！

## 根因分析

### 问题 1：mag_converged 设置导致 HSolver 被跳过

在 `deltaspin_lcao.cpp:run_deltaspin_lambda_loop_lcao()` 中：

```cpp
if (!sc.mag_converged() && drho > 0 && drho < inp.sc_scf_thr)
{
    sc.run_lambda_loop(iter);
    sc.set_mag_converged(true);  // ← 设置 mag_converged = true
    skip_solve = true;           // ← 跳过 HSolver
}
else if (sc.mag_converged())
{
    sc.run_lambda_loop(iter);
    skip_solve = true;           // ← 继续跳过 HSolver
}
```

**执行流程**：
1. GE8: drho=9.3457e-04 < sc_scf_thr(~1e-3)，进入 lambda 循环
2. lambda 循环跑 100 步后退出（未收敛）
3. `set_mag_converged(true)` 被调用
4. GE9: mag_converged=true → `skip_solve=true` → **HSolver 被完全跳过**
5. GE10: 继续 `skip_solve=true` → HSolver 继续被跳过

### 问题 2：Lambda 循环内部 skip_charge=true 导致 dm2rho 未调用

在 `cal_mw_from_lambda.cpp` 的 LCAO 路径中：

```cpp
hsolver_t.solve(hamilt_t, psi_t[0], this->pelec, *this->dm_, *this->pelec->charge, this->nspin_, true);
// ↑ skip_charge=true 意味着不调用 dm2rho
```

`HSolverLCAO::solve` 内部：
```cpp
elecstate::cal_dm_psi(...);  // 计算密度矩阵 ✓
dm.cal_DMR();                // 更新 DMR ✓

if (!skip_charge) {
    LCAO_domain::dm2rho(...); // ← skip_charge=true 时不调用！
}
```

**后果**：
- Lambda 循环内部每一步都更新了**密度矩阵 DM**
- 但**实空间电荷密度 chr.rho 从未被更新**（因为 skip_charge=true）
- Lambda 循环退出时调用 `update_psi_charge`，但 LCAO 路径的 `update_psi_charge` 只调用 `pelec->psiToRho(*psi_t)`，这是基类的**空实现**（直接 return）

### 问题 3：DRHO=0 的死循环

GE9 开始：
1. `skip_solve=true` → HSolver 被跳过
2. 没有新的波函数/密度矩阵/电荷密度
3. `chr.rho` 和 `chr.rho_save` 都是旧值 → `get_drho()` 返回 0
4. `mag_converged` 仍为 true → 继续 `skip_solve=true`
5. **死循环直到最大迭代次数**

### 问题 4：Lambda 收敛失败的可能原因

Lambda 循环未收敛（RMS 从 0.34 增长到 9.7）的原因可能与之前删除的 `cal_dm_psi` + `dm->cal_DMR()` 调用有关：

虽然 `HSolverLCAO::solve` 内部已经调用了这两个函数，但 `cal_mi_lcao` 在计算磁矩前可能需要额外的 DMR 同步。删除后的冗余调用移除可能导致了某些边界情况下的状态不一致。

## 修复方案

### 方案 1：在 skip_solve=true 时补充 dm2rho 调用

在 `esolver_ks_lcao.cpp:iter_run()` 中，当 `skip_solve=true` 时手动调用 `dm2rho`：

```cpp
if (!skip_solve)
{
    hsolver::HSolverLCAO<TK> hsolver_lcao_obj(&(this->pv), PARAM.inp.ks_solver);
    hsolver_lcao_obj.solve(..., skip_charge);
}
else
{
    // Lambda loop has updated DM, but not real-space charge density
    // Need to sync charge density from DM
    LCAO_domain::dm2rho(this->dmat.dm->get_DMR_vector(), PARAM.inp.nspin, &this->chr);
}
```

### 方案 2：修复 update_psi_charge 的 LCAO 路径

在 `cal_mw_from_lambda.cpp` 中，LCAO 路径的 `update_psi_charge` 应该调用 `dm2rho`：

```cpp
#ifdef __LCAO
if (PARAM.inp.basis_type == "lcao")
{
    // Update real-space charge density from density matrix
    LCAO_domain::dm2rho(this->dm_->get_DMR_vector(), this->nspin_, this->pelec->charge);
    Symmetry_rho::symmetrize_rho(this->nspin_, *this->pelec->charge, ...);
}
```

### 方案 3：检查 cal_dm_psi 删除的影响

考虑是否应该恢复 `cal_mw_from_lambda` 中被删除的 `cal_dm_psi` + `dm->cal_DMR()` 调用，或者确认 `HSolverLCAO::solve` 内部的调用是否足够。

## 时间统计佐证

从日志的 TIME STATISTICS 可以看到：
- `dm2rho` 只被调用了 **8 次**（对应 GE1-GE8）
- `cal_dm_psi` 被调用了 406 次（包括 lambda 循环内部）
- `cal_mi_lcao` 被调用了 409 次

这证实了 GE9 和 GE10 **没有调用 dm2rho**，电荷密度没有被更新。

## 基准版本对比

在 tmp 版本 (`log2-tmp-1`) 中：
- GE10: DRHO=5.0241e-04, ETOT=-6.77782768e+03 ← 正常收敛
- `psiToRho` 被调用了 16 次

在新版本 (`log`) 中：
- GE10: DRHO=0.0000e+00, ETOT=-2.17240744e+03 ← 异常
- `dm2rho` 只被调用了 8 次
