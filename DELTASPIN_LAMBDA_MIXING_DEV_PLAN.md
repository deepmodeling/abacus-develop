# DeltaSpin Lambda Mixing 详细开发文档

## 1. 背景与动机

### 1.1 问题

在 DeltaSpin 的 SCF 迭代过程中，每个外层 SCF 步运行一个 BFGS 内循环来优化 Lagrange 乘子 lambda，使原子磁矩 Mi 逼近目标值 M_target。当使用 subspace 加速模式时，内循环收敛的 lambda 存在系统性偏差，原因是 subspace 近似（H_sub = H0_sub + Σ Δλ·P_I_sub）仅在参考 lambda 附近有效。

当前 lambda 的传递方式：
```
lambda_{iter+1} = lambda_BFGS_converged(iter)   // 直接使用上次内循环的收敛值
```

这导致：
1. **偏差放大**：subspace 偏差 δλ_bias 全量传递到下一步，影响电荷密度计算
2. **SCF 振荡**：偏差方向可能交替变化，导致 (ρ, λ) 整体振荡
3. **收敛困难**：即使电荷密度本身已接近收敛，lambda 的跳动也能破坏收敛

### 1.2 提案

将 lambda 纳入 charge_mixing 框架：
```
lambda_{iter+1} = (1 - β_λ) · lambda_old + β_λ · lambda_BFGS_converged
```

### 1.3 合理性论证摘要

| 维度 | 论据 |
|---|---|
| 数学 | (ρ, λ) 构成耦合不动点问题，λ mixing 对交叉耦合项施加松弛因子，降低迭代矩阵谱半径 |
| 物理 | λ 是约束磁力 (eV/μB)，应与 ρ 的磁分量协调更新，而非突变 |
| 经验 | DFT+U 的 `mix_uom()` 在 ABACUS 中用同样框架混合占据矩阵，已成功实践 |
| 偏差衰减 | subspace 偏差 δλ_bias 经混合后残余 = β × δλ_bias（一致偏差）或 β^{1/2} × δλ_bias（振荡偏差） |

## 2. 设计方案

### 2.1 阶段一：简单线性混合（本期实现）

```
λ^{new}_i[ic] = (1 - β_λ) · λ^{prev}_i[ic] + β_λ · λ^{BFGS}_i[ic]    when constrain[i][ic] != 0
λ^{new}_i[ic] = 0                                                       when constrain[i][ic] == 0
```

参数：
- `sc_mixing_lambda_beta`：默认 -1.0（自动使用 mixing_beta），0.0 = 禁用，>0 = 指定值

适用条件：
- `sc_mag_switch == true` 且 `sc_scf_thr_mode != "off"` 且 `sc_lambda_strategy != "linear_scan"`
- 仅在有约束的分量上做混合；无约束分量始终为 0
- `direction_only` 模式下仅对垂直分量做混合

### 2.2 阶段二：Broyden 混合（未来升级）

在 `Charge_Mixing` 中新增 `Mixing_Data lambda_mdata`，复用 `Broyden_Mixing`，完全参照 `mix_uom()` 模式。

## 3. 实现细节

### 3.1 新增参数

| 参数名 | 类型 | 默认值 | 说明 |
|---|---|---|---|
| `sc_mixing_lambda_beta` | double | -1.0 | Lambda mixing beta。-1.0=自动(mixing_beta), 0.0=禁用, >0=指定值 |

在 `Input::read_sc()` 中解析，写入 `PARAM.inp`。

### 3.2 SpinConstrain 修改

在 `spin_constrain.h` 中新增：

```cpp
std::vector<ModuleBase::Vector3<double>> lambda_prev_;   // 上一步混合后的 lambda
bool lambda_mixing_enabled_ = false;                      // 是否启用 lambda mixing
double lambda_mixing_beta_ = 0.0;                         // 实际使用的 beta 值

void init_lambda_mixing(double beta);                    // 初始化
void mix_lambda();                                        // 执行混合
```

逻辑：
- `init_lambda_mixing(beta)`：设置 `lambda_mixing_beta_`，将 `lambda_prev_` 初始化为当前 `lambda_`
- `mix_lambda()`：
  1. 保存 `lambda_` 到 `lambda_bfgs`（BFGS 收敛值）
  2. 对每个约束分量做线性混合
  3. 更新 `lambda_` 为混合值
  4. 更新 `lambda_prev_` 为新的 `lambda_`

### 3.3 esolver_ks_lcao.cpp 修改

在 `iter_init()` 中，每次进入且 `sc_mag_switch` 开启时：
1. 首次调用前不做混合（没有历史数据）
2. 非首次调用前，调用 `sc.mix_lambda()` 对上一步的 BFGS 结果做混合

具体位置：在 `run_lambda_loop` 之前保存 lambda，之后执行混合。

### 3.4 特殊场景处理

| 场景 | 处理 |
|---|---|
| `sc_lambda_strategy = "linear_scan"` | 禁用 lambda mixing |
| `sc_scf_thr_mode = "off"` | 禁用 lambda mixing |
| `sc_mixing_lambda_beta = 0.0` | 禁用 lambda mixing |
| `direction_only` 模式 | 仅约束分量参与混合（自然由 constrain 数组控制） |
| Phase 1→2 过渡（direction_only 两阶段） | Phase 1 用刚得到的 lambda，Phase 2 衰减时自然受混合影响 |

### 3.5 数据流

```
SCF iteration N:
  1. 进入 iter_init()
  2. 如果 sc_mag_switch && 非首次:
     a. 执行 sc.mix_lambda()  ← 用 lambda_prev_ 和 lambda_ 做混合
     b. 混合后的 lambda_ 作为本轮 BFGS 的起始点
  3. sc.run_lambda_loop(iter-1)  → BFGS 优化 → lambda_ = λ^{BFGS}
  4. BFGS 结束后，lambda_prev_ = 混合前的值（已在 mix_lambda 中保存）
  
  下一次 SCF 迭代:
  1. mix_lambda(): λ_new = (1-β)·λ_prev + β·λ^{BFGS}
  2. 用 λ_new 作为起始点跑 BFGS
  ...

初始化:
  SCF iter 1: 不做 mixing（无历史）
  SCF iter 2 开始: 正常 mixing
```

### 3.6 对 BFGS 内循环的影响

混合后 lambda 可能略微偏离 BFGS 上次的收敛值，因此：
- BFGS 起始点可能不太理想 → 前几步搜索方向可能偏差
- 但 λ 的变化幅度被 β 衰减 → 实际偏离小于 |λ^{BFGS} - λ^{prev}|
- 对于振荡场景，混合后的起始点反而更接近真值（因为抑制了过冲）
- **建议**：BFGS 的 `alpha_trial` 不受影响，保持原有自适应逻辑

## 4. 文件修改清单

| 文件 | 修改内容 |
|---|---|
| `source/module_parameter/input.h` | 新增 `sc_mixing_lambda_beta` 成员 |
| `source/module_parameter/input.cpp` | 解析 `sc_mixing_lambda_beta`，默认 -1.0 |
| `source/module_parameter/parameter.h` | 无需修改（复用 existing 结构） |
| `source/module_parameter/parameter.cpp` | 无需修改 |
| `source/source_lcao/module_deltaspin/spin_constrain.h` | 新增 `lambda_prev_`, `lambda_mixing_beta_`, `lambda_mixing_enabled_`, 方法声明 |
| `source/source_lcao/module_deltaspin/spin_constrain.cpp` | 实现 `init_lambda_mixing()` 和 `mix_lambda()` |
| `source/source_esolver/esolver_ks_lcao.cpp` | 在 `iter_init()` 中集成 lambda mixing 调用 |
| `source/source_io/module_input/` 相关 | 输入文档更新（如需要） |

## 5. 测试计划

| 测试 | 预期结果 |
|---|---|
| `sc_mixing_lambda_beta = 0.0` | 行为与修改前完全一致（回归测试） |
| `sc_mixing_lambda_beta = -1.0`（默认） | 自动使用 mixing_beta，SCF 收敛有所改善或不劣化 |
| 振荡场景（subsace 模式下 Mi 徘徊） | 振荡幅度减小，收敛步数减少 |
| 非振荡场景 | 收敛步数相近或略增（β < 1 降低更新速率的代价） |
| `direction_only` 模式 | 正常运行，混合仅在约束分量生效 |

## 6. 后续扩展

1. **Broyden mixing for lambda**：在阶段一验证有效后，升级到 Broyden 框架
2. **自适应 β_λ**：根据 lambda 残差的历史动态调整 β_λ
3. **与 subspace 重建联动**：当 lambda mixing 导致 λ 显著偏离 λ_ref 时，触发 subspace cache 重建
