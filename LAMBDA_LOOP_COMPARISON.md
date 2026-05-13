# `run_lambda_loop` / `run_lambda_loop_lcao` 代码对比分析

## 一、代码位置

| 文件 | 当前版本 (port) | tmp 分支 (zdy/v3.7.1) |
|------|-----------------|------------------------|
| 主实现 | `source/source_lcao/module_deltaspin/lambda_loop.cpp` | `source/module_hamilt_lcao/module_deltaspin/lambda_loop.cpp` |
| 空模板实例化 | `source/source_lcao/module_deltaspin/template_helpers.cpp` | `source/module_hamilt_lcao/module_deltaspin/template_helpers.cpp` |
| 头文件 | `source/source_lcao/module_deltaspin/spin_constrain.h` | `source/module_hamilt_lcao/module_deltaspin/spin_constrain.h` |
| PW 层调用入口 | `source/source_pw/module_pwdft/deltaspin_pw.cpp` | 内联在 `source/module_esolver/esolver_ks_pw.cpp` |
| LCAO 层调用入口 | `source/source_esolver/esolver_ks_lcao.cpp` | `source/module_esolver/esolver_ks_lcao.cpp` |

---

## 二、`run_lambda_loop` (complex<double> 特化) — 逻辑与算法

### 2.1 共同的算法框架（两个版本相同）

两个版本共享同一个核心算法——**BFGS类共轭梯度法迭代优化 lambda Lagrange 乘子**：

```
目标: 找到 lambda_i 使得原子磁矩 M_i(lambda) = target_mag_i
方法: 迭代式优化，共轭梯度搜索 + 最优步长 alpha_opt
```

**主循环流程** (两个版本一致):

```
for i_step = -1 to nsc_:
  step -1 (初始化):
    1. cal_mw_from_lambda() → 从当前 lambda 计算原子磁矩 Mi
    2. 记录 initial_lambda (仅在 constrain!=0 的位置保留 lambda)
    3. 打印初始状态

  step >= 0 (迭代优化):
    4. lambda = initial_lambda + delta_lambda
    5. [direction_only] 投影: 去除 lambda 平行于 target_mag 的分量
    6. cal_mw_from_lambda(step, delta_lambda) → 新磁矩
    7. check_gradient_decay() → 检测梯度衰减, 若满足则提前退出
    8. delta_spin = spin - target_mag
    9. 计算 RMS error
    10. step==0 时设置 current_sc_thr = max(rms * sc_drop_thr, sc_thr)
    11. check_rms_stop() → 若收敛则:
        - update_psi_charge(dnu_last_step, rerun) → 更新波函数和电荷
        - [PW] cal_mi_pw() 二次检查, 若误差过大且 higher_mag_prec 则递归 rerun
        - break
    12. i_step>=2 时: beta = mean_error/mean_error_old, 共轭方向更新
    13. check_restriction() → 限制搜索方向幅度
    14. dnu = dnu + search * alpha_trial
    15. [direction_only] 投影: 去除 dnu 平行于 target_mag 的分量
    16. delta_lambda = dnu, 裁剪到 [-10, 10]
    17. cal_mw_from_lambda() → 试探磁矩 spin_plus
    18. cal_alpha_opt() → 最优步长
    19. check_restriction(search, alpha_opt)
    20. dnu += search * (alpha_opt - alpha_trial)  → 最终 dnu
    21. alpha_trial 自适应: g = 1.5 * |alpha_opt| / alpha_trial, alpha_trial *= g^0.7
```

**关键子函数:**
- `cal_mw_from_lambda()`: 在给定 lambda 下求解 Hamiltonian 得到原子磁矩
- `cal_alpha_opt()`: 线性插值求最优步长
  ```
  alpha_opt = alpha_trial * (target - spin) / (spin_plus - spin)
  ```
- `check_gradient_decay()`: 检测 |dM/dlambda| 是否衰减到阈值以下 (early exit)
- `check_rms_stop()`: 判断 RMS error 是否低于收敛阈值
- `check_restriction()`: 限制搜索方向的最大步幅
- `update_psi_charge()`: 用最终 lambda 更新波函数和电荷密度

### 2.2 差异点总结

| 差异项 | 当前版本 (port) | tmp 分支 (zdy) |
|--------|-----------------|----------------|
| **命名空间** | `spinconstrain::SpinConstrain` | `SpinConstrain` (全局) |
| **类模板参数** | `template <typename TK>` | `template <typename FPTYPE>` |
| **`ds_diag` 诊断代码** | 有！完整的 baseline 对比框架 (lines 105-328) | 无 |
| **`[DS-DIAG]` 打印** | `run_lambda_loop` 入口处打印 `lambda_` (lines 368-374) | 无 |
| **`ds_diag::check_step()` 调用** | 4 处: step=-1, RMS check, alpha_opt check, step=999 | 0 处 |
| **`cal_mi_pw()` vs `cal_Mi_pw()`** | `this->cal_mi_pw()` (下划线命名) | `this->cal_Mi_pw()` (驼峰命名) |
| **变量初始化** | `beta = 0.0, g = 0.0` 等内联初始化 | 分开声明后赋值 |
| **delta_lambda 裁剪** | 有！`abs > 10.0` 裁剪到 ±10 (lines 551-557) | 无 |
| **PW rerun 诊断打印** | 有 step=999 的 `ds_diag::check_step()` 调用 | 无 |
| **`run_lambda_loop_lcao`** | **存在！** 新增的 LCAO 专用优化路径 | **不存在** |
| **Lambda 策略系统** | 有 `LambdaStrategyType` 枚举和 `LambdaUpdateStrategy` 基类 | 无 |
| **`pw_wfc_` 成员** | 有 | 无 |
| **`phsol` 成员** | 无 (用 `pelec` 代替) | 有 `void* phsol` |
| **KS_SOLVER 成员** | 无 | 有 `std::string KS_SOLVER` |

---

## 三、`run_lambda_loop_lcao` — 当前版本新增函数

**仅在当前版本 (port) 中存在**, tmp 分支没有此函数。

### 3.1 设计目标

针对 LCAO nspin=2 体系的**高效** lambda 优化路径，避免重复求解全空间 Hamiltonian。

### 3.2 算法流程

```
Phase 1 — 全空间对角化:
  1. cal_mw_from_lambda(-1) → 获得初始 Mi, spin
  2. 计算 initial_lambda
  3. 检查是否已收敛 → 若是则直接更新电荷并返回

Phase 2 — 计算投影算符 P_I_sub:
  4. 对每个 k 点 ik: cal_PI_sub() → 获得 nbands×nbands 的 P_I_sub 矩阵
     P_I_sub 是 DeltaSpin 投影算符在现有波函数子空间中的表示

Phase 3 — 解析 Jacobian (磁化率):
  5. 对每个原子 iat: 计算 chi[iat] = dM_iat/dlambda_iat
     使用一阶微扰理论:
     chi = sum_{k,n,m} 2*(f_n - f_m) * |P_nm|² / (e_n - e_m)
     两个自旋通道贡献符号相同 (sign*sign=1)，叠加增强

Phase 4 — Newton 迭代 + 子空间验证:
  6. 对每次内迭代:
     a. Newton 步: delta_lambda = alpha_damp * (target - spin) / chi
        (alpha_damp = 0.8 为阻尼因子)
     b. 子空间对角化: 对每个 k 点构建 H_sub = diag(e_k) + sign * sum_I dlambda_I * P_I_sub
        用 LAPACK zheev 对角化得到新的本征值和本征向量
     c. 重新计算费米权重 calculate_weights()
     d. 从旋转后的子空间计算新磁矩 Mi_new:
        P_rotated = V^dag * P * V, M_i = sum_k sign * w_k * diag(P_rotated)
     e. 检查 RMS 收敛

Phase 5 — 波函数旋转和电荷更新:
  7. C_new = C * V (通过 pzgemm 矩阵乘法旋转波函数)
  8. 更新本征值 ekb, 权重 wg, 密度矩阵 dm, 电荷密度 rho
```

### 3.3 与主 `run_lambda_loop` 的对比

| 方面 | `run_lambda_loop` (通用) | `run_lambda_loop_lcao` (LCAO专用) |
|------|--------------------------|-----------------------------------|
| **适用体系** | 通用 (PW/LCAO 均可) | LCAO nspin=2 |
| **优化方法** | BFGS 类共轭梯度 + 试探步 | Newton 法 + 解析 Jacobian |
| **Hamiltonian 求解** | 每次迭代全空间对角化 | 仅 Phase 1 一次全空间对角化，后续子空间对角化 |
| **Jacobian** | 数值试探 (spin_plus) | 解析计算 (微扰理论) |
| **阻尼** | alpha_trial 自适应调整 | 固定 alpha_damp=0.8 |
| **共轭方向** | 有 (beta = mean_error/mean_error_old) | 无 |
| **收敛判据** | 动态 sc_thr + RMS | 动态 sc_thr + RMS |
| **成本** | 高 (每次迭代都需对角化) | 低 (子空间 zheev, O(nbands³)) |
| **方向约束** | 支持 direction_only | 仅优化 z 分量 |
| **波函数更新** | update_psi_charge() | 手动旋转 C*V + 更新 DM |

---

## 四、调用方差异

### 4.1 ESolver_KS_LCAO (LCAO 求解器)

**当前版本 (port):**
```cpp
// esolver_ks_lcao.cpp line 402-421
if (PARAM.inp.sc_mag_switch) {
    auto& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
    // [DS-DIAG] 额外打印调试信息
    if (!sc.mag_converged() && this->drho > 0 && this->drho < PARAM.inp.sc_scf_thr) {
        sc.run_lambda_loop(iter - 1);
        sc.set_mag_converged(true);
        skip_solve = true;
    } else if (sc.mag_converged()) {
        sc.run_lambda_loop(iter - 1);
        skip_solve = true;
    }
}
```

**tmp 分支 (zdy):**
```cpp
// esolver_ks_lcao.cpp line 728-744
if (PARAM.inp.sc_mag_switch) {
    SpinConstrain<TK>& sc = SpinConstrain<TK>::getScInstance();
    if(!sc.mag_converged() && this->drho>0 && this->drho < PARAM.inp.sc_scf_thr) {
        sc.run_lambda_loop(iter-1);
        sc.set_mag_converged(true);
        skip_solve = true;
    } else if(sc.mag_converged()) {
        sc.run_lambda_loop(iter-1);
        skip_solve = true;
    }
}
```

**差异:** 调用逻辑完全相同，仅命名空间和调试打印有差异。

### 4.2 PW 求解器

**当前版本 (port):** 抽取为独立函数 `run_deltaspin_lambda_loop()` 在 `deltaspin_pw.cpp`
```cpp
// deltaspin_pw.cpp — 封装成独立模块
bool run_deltaspin_lambda_loop(const int iter, const double drho, const Input_para& inp) {
    if (!inp.sc_mag_switch) return false;
    auto& sc = spinconstrain::SpinConstrain<complex<double>>::getScInstance();
    if (!sc.mag_converged() && drho > 0 && drho < inp.sc_scf_thr) {
        sc.run_lambda_loop(iter);  // 注意: 用 iter 而非 iter-1
        sc.set_mag_converged(true);
        return true;
    } else if (sc.mag_converged()) {
        sc.run_lambda_loop(iter);
        return true;
    }
    return false;
}
```

**tmp 分支 (zdy):** 直接内联在 `esolver_ks_pw.cpp`
```cpp
// esolver_ks_pw.cpp line 708-724 — 内联在 SCF 循环中
if (PARAM.inp.sc_mag_switch) {
    SpinConstrain<complex<double>>& sc = SpinConstrain<complex<double>>::getScInstance();
    if(!sc.mag_converged() && this->drho>0 && this->drho < PARAM.inp.sc_scf_thr) {
        sc.run_lambda_loop(iter-1);  // 注意: 用 iter-1
        sc.set_mag_converged(true);
        skip_solve = true;
    } else if(sc.mag_converged()) {
        sc.run_lambda_loop(iter-1);
        skip_solve = true;
    }
}
```

**差异:** 
1. port 版本用 `iter`，tmp 版本用 `iter-1` — 这可能导致 outer_step 索引差 1
2. port 版本封装为独立模块函数，tmp 版本内联

---

## 五、新增基础设施（仅当前版本）

### 5.1 `LambdaUpdateStrategy` 策略模式

当前版本引入了可扩展的 lambda 更新策略系统:

```cpp
enum class LambdaStrategyType { BFGS, LinearResponse, AugmentedLagrangian, HybridDelayed };

class LambdaUpdateStrategy {  // 抽象基类
    virtual LambdaUpdateResult update_lambda(...) = 0;
};

// 三个具体策略:
class LinearResponseUpdate      // 方案 B: 线性响应 + 混合
class AugmentedLagrangianUpdate  // 方案 C: 增广 Lagrangian
class HybridDelayedUpdate        // 方案 D: 混合延迟策略
```

### 5.2 头文件成员变量对比

| 成员 | port | tmp |
|------|------|-----|
| `direction_only_` | ✅ | ✅ |
| `strategy_type_` | ✅ (LambdaStrategyType) | ❌ |
| `strategy_` | ✅ (unique_ptr<LambdaUpdateStrategy>) | ❌ |
| `current_sc_thr_` | ✅ | ✅ |
| `pw_wfc_` | ✅ (PW_Basis_K*) | ❌ |
| `phsol` | ❌ | ✅ (void*) |
| `KS_SOLVER` | ❌ | ✅ (string) |
| `read_target_mag` | ❌ | ✅ (bool) |

---

## 六、潜在风险点

1. **outer_step 索引差异**: port 的 PW 路径使用 `iter`，tmp 使用 `iter-1`。如果外层 SCF 计数从 1 开始，这会导致诊断数值偏移 1。

2. **delta_lambda 裁剪**: port 新增了 `abs(delta_lambda) > 10.0` 的硬裁剪。在 tmp 分支能正常工作的情况下，此裁剪可能在某些场景下干扰收敛。

3. **ds_diag 诊断开销**: 大量 `check_step()` 调用在调试模式外可能产生不必要的 I/O 开销。

4. **`run_lambda_loop_lcao` 未集成到调用链**: 函数已实现但在 `esolver_ks_lcao.cpp` 中未被调用，仍使用通用的 `run_lambda_loop`。
