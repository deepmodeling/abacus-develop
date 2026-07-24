# DeltaSpin Fast模式对角化次数优化方案

## 1. 问题分析

### 1.1 当前Fast模式下每次SCF迭代中的完整对角化次数

在fast模式（`sc_strategy="fast"`）下，每次外部SCF迭代的lambda循环中，
完整的HSolverLCAO对角化被调用 **4次**：

| 序号 | 触发位置 | i_step | 分支 | 目的 | 冗余原因 |
|------|----------|--------|------|------|----------|
| ① | 初始化 | -1 | Branch 3 | 计算初始Mi | 可与③合并 |
| ② | 首步优化 | 1 | Branch 3 | 加速未激活时的首次优化 | 加速此时已可激活，与③在同一lambda处做相同对角化 |
| ③ | 缓存构建 | -2 | Branch 1 | 完整对角化+构建子空间缓存 | 与②输入完全相同，仅额外缓存子空间数据 |
| ④ | 收敛后更新 | update_psi_charge | — | 将psi从参考λ旋转到收敛λ | 可用子空间旋转替代完整对角化 |

非DeltaSpin计算每次SCF仅需1次对角化 → fast模式开销约为 **4倍**。

### 1.2 冗余详解

**🔴 ①②③可合并**：初始化(i_step=-1)调用Branch 3仅计算Mi，
不做缓存。首步(i_step=1)又调用Branch 3，然后加速激活检查调用Branch 1再次对角化。
但②和③在**完全相同的lambda值**处进行完整对角化，区别仅在于③额外缓存了子空间数据。
如果我们预先激活加速，则：
- ①直接调用Branch 1（完整对角化 + 缓存 + Mi计算），一步完成
- ②直接走Branch 2（子空间对角化，O(nbands³)），无需完整对角化
- ③不再需要（缓存已在①中构建）

**🔴 ④可消除**：收敛后`update_psi_charge`做完整对角化仅为了将psi
从参考lambda旋转到收敛lambda。但子空间方法已经给出了旋转变换V，
可以直接应用 `ψ_new = ψ_ref · V`，代价仅 O(N × nbands²)。

## 2. 优化方案

### 2.1 Phase 3a: 合并初始化与缓存构建（节省2次完整对角化）

**核心思想**：在fast模式下，预先设置 `acceleration_active_ = true`，
让i_step=-1直接调用 `cal_mw_from_lambda(-2)` (Branch 1)，
即一次完整对角化 + 构建缓存 + 计算初始Mi。

**修改位置**：`lambda_loop.cpp` 的 `run_lambda_loop()` 初始化块

```cpp
// 修改前（当前代码）
if (i_step == -1)
{
    this->cal_mw_from_lambda(i_step);  // Branch 3: 仅计算Mi
    spin = this->Mi_;
    where_fill_scalar_else_2d(this->constrain_, 0, zero, this->lambda_, initial_lambda);
    i_step++;
}
```

```cpp
// 修改后
if (i_step == -1)
{
    if (this->sc_acceleration_mode_ != "off" && this->sc_acceleration_rms_thr_ >= 1e8)
    {
        // Fast模式：预激活加速，合并初始化与缓存构建
        this->acceleration_active_ = true;
        this->acceleration_subspace_built_ = false;
        this->cal_mw_from_lambda(-2);  // Branch 1: 完整对角化 + 缓存 + Mi
        this->acceleration_subspace_built_ = true;
    }
    else
    {
        this->cal_mw_from_lambda(i_step);  // Branch 3: 正常初始化
    }
    spin = this->Mi_;
    where_fill_scalar_else_2d(this->constrain_, 0, zero, this->lambda_, initial_lambda);
    i_step++;
}
```

**效果**：
- 合并①②③为1次完整对角化
- 首步优化(i_step=1)直接进入Branch 2（子空间方法）
- 从4次降至2次完整对角化/SCF

**fast_mode检测**：用 `sc_acceleration_rms_thr_ >= 1e8` 检测。
（fast模式设置 rms_thr = 1e10，远大于1e8；normal模式设置1e-2，远小于1e8）

### 2.2 Phase 3b: 用子空间旋转替代收敛后的完整对角化（节省1次完整对角化）

**核心思想**：收敛后在`update_psi_charge`中，
不再做完整对角化，而是：
1. 在收敛lambda处重新构建H_sub（从缓存H₀_sub + Δλ·P_I_sub）
2. 子空间对角化 → 得到V（O(nbands³)，远比完整对角化O(N²nbands)便宜）
3. 应用旋转 ψ_new = ψ_ref · V
4. 从旋转后的psi构造电荷密度

**修改位置**：`cal_mw_from_lambda.cpp` 的 `update_psi_charge()` LCAO路径

```cpp
// 修改后
if (this->acceleration_active_ && this->acceleration_subspace_built_
    && this->nspin_ == 2)
{
    psi::Psi<std::complex<double>>* psi_t = ...;
    
    // 1. 在收敛lambda处重新构建H_sub并子空间对角化
    //    （复用Branch 2的逻辑，但这次保留旋转后的psi）
    const int nk = psi_t->get_nk();
    const int nbands = PARAM.inp.nbands;
    const int nn = nbands * nbands;
    std::vector<std::vector<std::complex<double>>> vcc_conv(nk);
    
    for (int ik = 0; ik < nk; ik++)
    {
        // 构建H_sub(k) = H₀_sub(k) + Σ_I Δλ_I · P_I_sub(k)
        // （使用缓存的H₀_sub和P_I_sub，Δλ = lambda_conv - lambda_ref）
        // 子空间对角化 → V_conv, ε_conv
        // 保存V_conv[ik]
    }
    
    // 2. 应用旋转 ψ_new = ψ_ref · V_conv
    this->rotate_psi_subspace_lcao(*psi_t, this->ParaV, vcc_conv, nbands, nlocal, nk);
    
    // 3. 更新Fermi权重（使用子空间对角化得到的本征值）
    elecstate::calculate_weights(...);
    elecstate::calEBand(...);
    
    // 4. 从旋转后的psi构造电荷密度
    this->pelec->psiToRho(*psi_t);
}
else
{
    // 原有完整对角化路径（fallback：acceleration未激活或nspin=4）
    ...
}
```

**效果**：从2次降至1次完整对角化/SCF，与非DeltaSpin计算成本相同。

**风险与缓解**：
- **风险**：子空间旋转得到的psi可能与完整对角化结果有微小差异，
  电荷密度略有偏差 → 可能增加外层SCF迭代次数
- **缓解**：
  1. 外层SCF循环本身会修正电荷密度误差，不会导致发散
  2. 可设置 `sc_strategy = "normal"` 作为更稳健的备选方案
  3. 在收敛后期lambda变化很小，子空间近似精度极高，误差可忽略

### 2.3 远期展望：子空间缓存跨SCF迭代复用

**当前行为**：每次`run_lambda_loop`开始时重置加速状态并释放缓存。

**优化方向**：如果相邻SCF迭代的哈密顿量变化不大（电荷密度RMS 小于 阈值），
可以复用上一次的子空间缓存，省去缓存构建的完整对角化。

**效果**：在SCF收敛后期降至0次完整对角化/SCF（仅需子空间更新）。

**风险评估**：较高。使用过期的子空间数据可能导致收敛失败。
建议仅在fast模式+后期SCF迭代中启用，并设置严格的重置条件。
当前不建议实施。

## 3. 成本对比

| 模式 | 当前 | Phase 3a后 | Phase 3a+3b后 | 非DeltaSpin |
|------|------|-----------|--------------|-------------|
| fast | 4次完整diag | 2次 | 1次 | 1次 |
| normal | ~3次 | 不变 | 不变 | 1次 |
| accuracy | ~(2nsc+1)次 | 不变 | 不变 | 1次 |

## 4. 实现优先级

1. **Phase 3a**（高优先级）：简单、低风险、效果显著（4→2次）
2. **Phase 3b**（中优先级）：中等复杂度、需测试验证（2→1次），达到与非约束计算平价
3. **缓存跨迭代复用**（低优先级）：复杂、风险高，可作为未来增强

## 5. 测试计划

- 使用 `tests/fdf_fe_nspin2` 作为基准测试
- Phase 3a：验证收敛后磁矩Mi与原实现一致（误差 < 1e-4 uB）
- Phase 3b：比较SCF收敛速度（迭代次数）与非约束计算
- 并行测试（1核、4核、16核）
- 对比每SCF迭代的wall时间
