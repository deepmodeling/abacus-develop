# DeltaSpin Lambda Mixing 方案算法评估

## 1. 问题描述

在 DeltaSpin 的 SCF 迭代过程中，每个外层 SCF 步会运行一个 BFGS 内循环来优化 lambda（Lagrange 乘子），使得原子磁矩 Mi 逼近目标值 M_target。当使用 subspace 加速模式时，内循环收敛的 lambda 存在系统性偏差（因为 subspace 近似引入误差）。

当前的 lambda 传递方式是：
```
lambda_{iter+1} = lambda_BFGS_converged(iter)   // 直接使用上次内循环的收敛值
```

这导致：
- subspace 模式下 lambda 有偏差 → 电荷密度响应有偏差 → SCF 收敛困难或振荡
- 偏差方向可能一致（拖慢收敛）或振荡（正负交替叠加）

**提案：** 将 lambda 也纳入 charge_mixing 框架：
```
lambda_{iter+1} = (1 - β_λ) · lambda_old + β_λ · lambda_BFGS_converged(iter)
```
其中 β_λ 可以与电荷混合系数 β_rho 相同或有独立的比率。

## 2. 数学分析

### 2.1 自洽方程的结构

DeltaSpin 约束 DFT 的自洽问题是求解联立方程组：

```
ρ = F[λ, ρ]    (Kohn-Sham 自洽：给定 λ 和 ρ，求解得到新的 ρ)
Mi(λ, ρ) = M_target,i    (约束条件：给定 ρ，调整 λ 使磁矩满足约束)
```

这是一个**嵌套不动点问题**：
- 外循环：ρ 的不动点迭代（标准 SCF）
- 内循环：给定 ρ，求解 λ 使得 Mi ≈ M_target（BFGS 优化）

稳态条件是联立解 (ρ*, λ*) 满足：
```
ρ* = F[λ*, ρ*]
Mi(λ*, ρ*) = M_target,i
```

### 2.2 当前算法的收敛动力学

将问题写作不动点迭代：
```
ρ_{n+1} = F[λ_n, ρ_n]           ← 外层 SCF 步
λ_n = Λ[ρ_n]                    ← 内层 BFGS 优化的收敛值
```

其中 Λ[ρ] 是"给定电荷密度 ρ，BFGS 收敛得到的 λ"的映射。

展开到线性阶：
```
δρ_{n+1} = (∂F/∂ρ) δρ_n + (∂F/∂λ) δλ_n
δλ_n = (∂Λ/∂ρ) δρ_n
```

代入得：
```
δρ_{n+1} = [(∂F/∂ρ) + (∂F/∂λ)(∂Λ/∂ρ)] δρ_n
```

收敛条件是谱半径 ρ(J) < 1，其中 J = (∂F/∂ρ) + (∂F/∂λ)(∂Λ/∂ρ)。

**关键观察：** 交叉项 (∂F/∂λ)(∂Λ/∂ρ) 可能增大谱半径，导致收敛减慢或振荡。这是因为 λ 对 ρ 的响应 (∂Λ/∂ρ) 会放大 ρ 的误差，而 F 对 λ 的响应 (∂F/∂λ) 将这个放大的误差反馈到 ρ 中。

### 2.3 Lambda mixing 的数学效果

引入混合后：
```
λ_n = (1 - β_λ) λ_{n-1} + β_λ · Λ[ρ_n]
```

线性化：
```
δλ_n = (1 - β_λ) δλ_{n-1} + β_λ · (∂Λ/∂ρ) δρ_n
```

代入 ρ 迭代：
```
δρ_{n+1} = (∂F/∂ρ) δρ_n + (∂F/∂λ)[(1 - β_λ) δλ_{n-1} + β_λ (∂Λ/∂ρ) δρ_n]
```

这是一个**增广系统**：(δρ, δλ) 两维迭代。混合降低了 λ 的更新幅度，相当于在 δλ 方向上施加了一个**松弛因子**，减小交叉耦合的增益。

**定性结论：** λ mixing 在数学上是合理的，它可以：
1. **降低交叉耦合增益**：减小 (∂F/∂λ)(∂Λ/∂ρ) 项的有效贡献
2. **抑制振荡**：对 λ 的变化做低通滤波，平滑迭代轨迹
3. **补偿 subspace 偏差**：偏差被混合系数衰减而非全量传递

### 2.4 与标准电荷混合的类比

标准电荷混合处理的核心问题是：
```
ρ_{n+1} = (1 - β) ρ_n + β · Ψ[ρ_n]
```
其中 Ψ[ρ] 是"给定 ρ，求解 Kohn-Sham 方程得到的输出密度"。

Bell-Ramanan 定理表明：如果 Ψ 的 Jacobian 的谱半径 > 1（不混合不收敛），则存在 β ∈ (0,1) 使得混合后谱半径 < 1。

**DeltaSpin lambda 的关键区别：** λ 不是直接通过不动点迭代得到的——它是通过 BFGS 优化内循环精确求解的（在给定 ρ 的前提下）。因此 Λ[ρ] 的"输出"是精确满足约束的 λ，而非某个近似固定的映射输出。

这意味着：
- **如果 BFGS 内循环精确收敛**：Λ[ρ] 是 ρ 的确定性映射，λ mixing 类似于解决一个耦合不动点问题
- **如果 BFGS 内循环有误差**（subspace 模式）：Λ[ρ] = Λ_exact[ρ] + ε(ρ)，ε 是 subspace 近似误差。此时 λ mixing 更直接地起到误差衰减作用

## 3. 方案评估

### 3.1 方案 A：简单线性混合（推荐起步方案）

```
λ^{new}_i = (1 - β_λ) · λ^{old}_i + β_λ · λ^{BFGS}_i
```

参数选择：
- **β_λ = β_rho（与电荷混合系数相同）**：最简单的选择，保持 ρ 和 λ 的更新步长一致
- **β_λ = β_mag（使用磁矩混合系数）**：`mixing_beta_mag` 是 ABACUS 为磁性分量单独设置的混合系数（默认 2×β_rho），更激进
- **β_λ < β_rho（更保守）**：如 β_λ = 0.5 × β_rho，避免 λ 变化过快冲击电荷密度

**优点：**
- 实现极其简单（10 行代码）
- 与现有 charge_mixing 框架完全解耦，不侵入核心混合逻辑
- 物理直觉清晰：λ 是"磁场"，应与 ρ 的磁分量同比例更新

**缺点：**
- 不利用历史信息，收敛速度可能不如 Broyden
- β_λ 需要手动调参（但有默认值可覆盖大多数场景）

### 3.2 方案 B：Broyden 混合（推荐终态方案）

利用 `Base_Mixing::Broyden_Mixing` 或 `Pulay_Mixing` 框架：

```
residual_λ = λ^{BFGS} - λ^{old}
// Broyden 计算 c_1, c_2, ..., c_n
λ^{new} = λ^{old} + Σ c_k · residual_λ(k)
```

**实现方式：** 在 `Charge_Mixing` 中新增 `Mixing_Data lambda_mdata`，调用已有的 `mixing->push_data()` 和 `mixing->mix_data()`。

**优点：**
- 利用 SCF 历史信息，加速收敛
- 与 DFT+U 的 `mix_uom()` 模式完全类比（已在代码中实现）

**缺点：**
- 实现稍复杂（需要分配 `Mixing_Data`，在 `init_mixing()` 中初始化）
- 需要 `SpinConstrain` 与 `Charge_Mixing` 之间的数据传递接口

### 3.3 方案 C：耦合混合（高级方案，暂不推荐）

将 (ρ, λ) 拼接为一个联合向量，在一个统一的 Broyden 框架中混合：

```
x = [ρ_mag, λ]    // 联合变量
r = [ρ_out - ρ_in, λ^{BFGS} - λ^{old}]    // 联合残差
// 统一 Broyden 更新
```

**优点：** 最优的收敛速率（理论上）
**缺点：** 需要定义联合变量的内积，实现复杂度显著增加，且 ρ 和 λ 的量纲和尺度差异大，需要仔细调权

## 4. 关键问题与风险评估

### 4.1 内循环是否需要重复？

引入 λ mixing 后，BFGS 内循环的起始点变为混合后的 λ，而非上次收敛的 λ。这意味着：
- BFGS 起始偏离最优值 → 内循环可能需要更多步
- 但 λ 的变化被平滑 → 实际上偏离可能更小（尤其在振荡场景下）

**结论：** 内循环步数可能略微增加或减少，取决于是否振荡。总体影响应在实践中评估。

### 4.2 Lambda mixing 与 subspace 加速的交互

Subspace 模式下 lambda 有系统性偏差。设偏差为 δλ_bias，则：
- 无混合：λ_{n+1} = Λ[ρ_n] ≈ λ_exact + δλ_bias，偏差全量传递
- 有混合：λ_{n+1} = (1-β)λ_n + β(λ_exact + δλ_bias)，偏差被 β 衰减

如果 δλ_bias 在多次迭代中方向一致，混合后等效偏差 = β × δλ_bias（显著减小）。
如果 δλ_bias 方向振荡，混合后等效偏差 = β^{1/2} × δλ_bias（更显著减小）。

**结论：** λ mixing 对 subspace 偏差有直接的衰减效果，这正是期望的。

### 4.3 与 DFT+U 占据矩阵混合的对比

ABACUS 的 DFT+U 模块已经实现了 `mix_uom()`（`charge_mixing.cpp:276`），将 Hubbard U 的占据矩阵纳入 Broyden 混合。这是一个完全类比的方案：

| 特征 | DFT+U 占据矩阵 | DeltaSpin λ |
|---|---|---|
| 变量 | n^{I}_{mm'} (每个原子 2l+1 维) | λ_i (每个原子 3 维) |
| 混合方式 | Broyden (复用 mixing 对象) | 可以完全复用同样的框架 |
| 代码位置 | `Charge_Mixing::mix_uom()` | 新增 `Charge_Mixing::mix_lambda()` |
| 效果 | 显著改善 DFT+U SCF 收敛 | 预期改善 DeltaSpin SCF 收敛 |
| 初始化 | `allocate_mixing_uom()` | 新增 `allocate_mixing_lambda()` |

**结论：** DFT+U 的成功实践为 λ mixing 提供了直接的经验支持。

### 4.4 Lambda 的物理意义与混合的合理性

λ 在物理上是施加在原子上的"约束磁力"（constrained magnetic field），单位 eV/μB。在 SCF 过程中：
- 电荷密度 ρ 通过自洽迭代趋近稳态
- λ 也应趋近其稳态值 λ*
- 两者是**耦合的**：ρ 影响最优 λ，λ 影响最优 ρ

混合的物理直觉：如果上一步 λ 偏小，BFGS 优化会给出偏大的 λ（补偿 ρ 的变化），但 ρ 还没更新到与这个 λ 匹配的状态，所以不应全量采用新的 λ，而应取一个中间值。

**这正是混合的标准动机——防止过冲。**

### 4.5 潜在风险

| 风险 | 严重性 | 缓解方案 |
|---|---|---|
| β_λ 过大导致不收敛 | 中 | 默认 β_λ = β_rho，保守起步 |
| β_λ 过小导致收敛极慢 | 中 | 可设 β_λ = β_mag（≈2×β_rho） |
| 与 charge_mixing 的 Kerker 屏蔽不一致 | 低 | λ 无需 Kerker（不是长程量） |
| BFGS 内循环步数增加 | 低→中 | 实测评估；可设 β_λ 以小值起步 |
| 初始迭代（λ=0 阶段）用了 mixing 反而慢 | 低 | 仅在 lambda_loop 激活后启用 mixing |

## 5. 推荐实现方案

### 5.1 第一阶段：简单线性混合（快速验证）

**改动范围：** `esolver_ks_lcao.cpp` 中 lambda_loop 调用前后

```cpp
// 在 run_lambda_loop 返回后
auto lambda_bfgs = sc.get_sc_lambda();  // BFGS 收敛的 lambda

// 线性混合
double beta_lambda = p_chgmix->get_mixing_beta();  // 复用电荷混合系数
for (int ia = 0; ia < nat; ia++) {
    for (int ic = 0; ic < 3; ic++) {
        if (sc.get_constrain()[ia][ic] != 0) {
            double lambda_old = lambda_prev[ia][ic];
            double lambda_new = (1.0 - beta_lambda) * lambda_old 
                              + beta_lambda * lambda_bfgs[ia][ic];
            lambda_bfgs[ia][ic] = lambda_new;
        }
    }
}
sc.set_sc_lambda(lambda_bfgs);
```

**参数设计：**
- 新增参数 `sc_mixing_lambda_beta`（🔴 新增参数），默认 -1.0
  - -1.0 → 自动使用 `mixing_beta`（与电荷相同）
  - 0.0 → 禁用 lambda mixing（等价于当前行为）
  - >0 → 使用指定值

### 5.2 第二阶段：Broyden 混合（生产级别）

在 `Charge_Mixing` 中新增：

```cpp
// charge_mixing.h
class Charge_Mixing {
    ...
    Base_Mixing::Mixing_Data lambda_mdata;  // lambda mixing data
    
    void allocate_mixing_lambda(int nat_constrained);
    void mix_lambda(std::vector<ModuleBase::Vector3<double>>& lambda_in,
                    std::vector<ModuleBase::Vector3<double>>& lambda_save);
};
```

参照 `mix_uom()` 的实现模式，完全复用 `Base_Mixing::Broyden_Mixing`。

### 5.3 接口需求

`SpinConstrain` 需要暴露以下接口给 `Charge_Mixing`：

| 接口 | 方向 | 说明 |
|---|---|---|
| `get_sc_lambda()` | SC → CM | 获取 BFGS 收敛的 lambda |
| `set_sc_lambda()` | CM → SC | 设置混合后的 lambda |
| `get_constrain()` | SC → CM | 判断哪些分量有约束 |
| `get_nat_constrained()` | SC → CM | 约束原子总数（确定 mixing data 大小） |

## 6. 与 charge_mixing 的混合系数关系

### 6.1 推荐默认值

| 参数 | 推荐值 | 理由 |
|---|---|---|
| β_λ (nspin=2) | = β_rho | λ_z 与 ρ_mag 同量级耦合，同步更新 |
| β_λ (nspin=4) | = β_mag (= 2×β_rho) | 非共线有 3 个分量，可能需要更激进更新 |
| Kerker 屏蔽 | 不需要 | λ 是局域量，无长程成份 |
| mixing_ndim | 同 β_rho | 复用 Broyden 历史长度 |

### 6.2 特殊场景

- **direction_only 模式**：λ 的平行分量被投影掉，只有垂直分量参与更新。mixing 只对垂直分量生效。
- **linear_scan 模式**：λ 是扫描变量，不需要 mixing。
- **sc_scf_thr_mode = "off"**：λ 从 STRU 读取且不优化，不需要 mixing。

## 7. 结论

**Lambda mixing 在算法上是合理的，推荐实现。**

理由：
1. **数学上**：λ mixing 等价于对耦合不动点问题的松弛，降低交叉耦合增益，与标准电荷混合的动机完全一致
2. **物理上**：λ 是约束磁力，应与电荷密度协调更新，而非突变
3. **经验上**：DFT+U 的占据矩阵混合（`mix_uom`）已在 ABACUS 中成功实践，完全类比
4. **对 subspace 偏差**：直接衰减系统性偏差，是最自然的解决方案
5. **实现成本**：方案 A（线性混合）仅需约 20 行代码；方案 B（Broyden）约 80 行

**推荐路径：** 先实现方案 A 验证效果，确认有帮助后升级到方案 B。
