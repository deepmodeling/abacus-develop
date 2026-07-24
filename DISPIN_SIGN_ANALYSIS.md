# DeltaSpin 符号一致性分析与修复方案

## 1. 统一范式

```
E'[ρ] = E_DFT[ρ] + Σ_I λ_I · (M_I - M_target_I)
```

变分得哈密顿修正：

```
H_DS = +Σ_I λ_I · σ_I
```

### nspin=2 (collinear)

```
H_DS|↑> = +λ_z |↑>
H_DS|↓> = -λ_z |↓>
```

系数：spin-up = **+λ_z**，spin-down = **-λ_z**

### nspin=4 (non-collinear)

```
        | +λ_z       +λ_x + iλ_y |
H_DS =  |                             |
        | +λ_x - iλ_y    -λ_z       |
```

## 2. 代码现状

### PW 基组 — 全部符合统一范式 ✓

| 文件 | nspin=2 | nspin=4 |
|------|---------|---------|
| `op_pw_proj.cpp` (Hψ 主路径) | +λ_z × spin_sign ✓ | +λ_z, +λ_x+iλ_y, +λ_x-iλ_y, -λ_z ✓ |
| `cal_mw_from_lambda.cpp` calculate_delta_hcc | +λ_z × spin_sign ✓ | 同上 ✓ |

### LCAO 基组 — 全部与统一范式相反 ✗

| 文件 | nspin=2 | nspin=4 |
|------|---------|---------|
| `dspin_lcao.cpp` cal_coeff_lambda | -λ_z, +λ_z ✗ | -λ_z, -λ_x-iλ_y, -λ_x+iλ_y, +λ_z ✗ |
| `lcao_subspace.cpp` calculate_delta_hcc_lcao | -λ_z × spin_sign ✗ | +λ_z ✓ (已正确) |

注意：`lcao_subspace.cpp` 的 nspin=4 分支已经是正号（正确），但 nspin=2 分支是负号（错误）。而 `dspin_lcao.cpp` 两个自旋通道的符号都反了。

另外需检查：`lambda_loop.cpp` 和 `diagonalization_engine.cpp` 中是否有独立的系数构建。

## 3. escon 公式

统一范式下 H_DS = +λ·σ，额外贡献 +sum(λ·Mi) 进入 eband。

```
E_KS ≈ E_DFT + Σ λ·Mi
E_DFT = E_KS + escon
escon = -Σ λ·Mi
```

代码：`escon -= λ·Mi`，`etot += escon`。

## 4. 修复方案

### 修改 1: `dspin_lcao.cpp` `cal_coeff_lambda` — 翻转全部符号

**nspin=2** (第 85-86 行):
```cpp
// 修改前：
coefficients[0] = -current_lambda[0];
coefficients[1] =  current_lambda[0];
// 修改后：
coefficients[0] =  current_lambda[0];
coefficients[1] = -current_lambda[0];
```

**nspin=4** (第 98-101 行):
```cpp
// 修改前：
coefficients[0] = complex(-λ_z, 0);
coefficients[1] = complex(-λ_x, -λ_y);
coefficients[2] = complex(-λ_x, +λ_y);
coefficients[3] = complex(+λ_z, 0);
// 修改后（与 PW 完全一致）：
coefficients[0] = complex(+λ_z, 0);
coefficients[1] = complex(+λ_x, +λ_y);
coefficients[2] = complex(+λ_x, -λ_y);
coefficients[3] = complex(-λ_z, 0);
```

### 修改 2: `lcao_subspace.cpp` `calculate_delta_hcc_lcao` — 修复 nspin=2 符号

第 475 行:
```cpp
// 修改前：
const std::complex<double> coeff(-effective_lambda[iat][2] * spin_sign, 0.0);
// 修改后：
const std::complex<double> coeff(effective_lambda[iat][2] * spin_sign, 0.0);
```

nspin=4 分支（第 459 行）已经是 `+λ_z`，无需修改。

### 修改 3: `cal_escon()` — 移除 `is_Mi_converged` 门控，保持 `escon -= λ·Mi`

### 修改 4: 为 PW 基组 `ElecStatePW` 添加 `get_spin_constrain_energy()` override

### 修改 5: 同步更新 `dspin_lcao.cpp` 中的注释，使其与统一范式一致

## 5. 需要一并检查的文件

- `lambda_loop.cpp`: 是否有独立的 lambda 系数构建 →如有，同步翻转
- `diagonalization_engine.cpp`: 是否有独立的 lambda 系数构建 →如有，同步翻转
- `cal_mw_from_lambda.cpp`: `calculate_delta_hcc` 是 PW 路径，已正确 ✓
- 力和应力代码 (`onsite_proj_tools.cpp`): 使用独立布局，不在本次范围
