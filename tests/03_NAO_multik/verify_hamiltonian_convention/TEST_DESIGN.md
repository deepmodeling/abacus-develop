# LCAO 非共线哈密顿量复共轭错误修复

## 0. 问题描述

在 LCAO 基组下，nspin=4（非共线自旋）计算中，从 Pauli 基的实数有效势
`(V_0, B_x, B_y, B_z)` 构建旋量基复数哈密顿量的过程中存在两处错误：

1. **非对角元相位符号错误**：`H_{↑,↓}` 和 `H_{↓,↑}` 的虚部符号相反
2. **下三角填充缺少复共轭**：`H(-R)` 使用转置而非共轭转置，破坏 Hermiticity

**影响**：非共线哈密顿量与正确结果互为复共轭（仅对非对角元），导致自洽后的
波函数自旋纹理与磁化方向不一致。PW 代码路径不受影响。

**与 DM Fourier 变换符号修复的区别**：此 bug 位于哈密顿量构建环节
（grid integration → H(R)），与 DM 的 k→R Fourier 变换无关。DM 的 Fourier
符号修复（`density_matrix.cpp` 和 `dmr_complex.cpp` 中的 `sinp → -sinp`）
是另一个独立的问题。

---

## 1. 正确公式

### 1.1 Pauli 基 → 旋量基转换

有效势在 Pauli 基中为 4 个实数矩阵：
```
v(0) = V_0    (标量势，正比于电荷密度 n)
v(1) = B_x    (磁场 x 分量)
v(2) = B_y    (磁场 y 分量)
v(3) = B_z    (磁场 z 分量)
```

旋量基哈密顿量为 `H = V_0 * I + B_x * σ_x + B_y * σ_y + B_z * σ_z`，其中：

```
σ_x = [[0, 1], [1, 0]]
σ_y = [[0, -i], [i, 0]]     ← 关键：σ_y 的非对角元含 ±i
σ_z = [[1, 0], [0, -1]]
```

展开得 4 个旋量矩阵元：

```
H_{↑,↑}   = V_0 + B_z                  （纯实数）
H_{↑,↓}   = B_x - i*B_y                （虚部为负）
H_{↓,↑}   = B_x + i*B_y                （虚部为正）
H_{↓,↓}   = V_0 - B_z                  （纯实数）
```

### 1.2 实空间 Hermiticity

哈密顿量矩阵元 `H_{μν}(R) = <φ_{μ0}|H|φ_{νR}>` 需满足：

```
H_{μν}(R) = H*_{νμ}(-R)
```

即 `H(-R) = H†(R)`。在代码中，对 atom pair `(iat1, iat2)`，若已知
`(iat1, iat2, R)` 处的矩阵 `H_upper`，则 `(iat2, iat1, -R)` 处的矩阵为：

```
H_lower(-R) = H_upper†(R) = conj(H_upper^T(R))
```

---

## 2. 代码中的错误

### 2.1 错误位置

文件：`source/source_lcao/module_gint/gint_common.cpp`
函数：`merge_hr_part_to_hR()`

### 2.2 错误 1：非对角元相位

代码用系数 `(clx_i + i*clx_j)` 组合 Pauli 分量：

```cpp
// 错误代码：
std::vector<int> clx_i = {1, 0, 0, -1};
std::vector<int> clx_j = {0, 1, -1, 0};
//                         ↑   ↑
//                    is=1: +i  is=2: -i
//                    应为: -i  应为: +i
```

这导致：
```
H_{↑,↓} = B_x + i*B_y    （错误：虚部为正）
H_{↓,↑} = B_x - i*B_y    （错误：虚部为负）
```

与正确公式恰好互为复共轭。

### 2.3 错误 2：下三角填充

```cpp
// 错误代码：只转置，不取共轭
lower_mat->get_value(icol, irow) = upper_mat->get_value(irow, icol);
```

对实数矩阵 (is=0,3)，转置 = 共轭转置，无影响。
对复数矩阵 (is=1,2)，转置 ≠ 共轭转置，导致 `H(-R) ≠ H†(R)`。

---

## 3. 修复

```cpp
// 修复后：
std::vector<int> clx_j = {0, -1, 1, 0};
//                         ↑    ↑
//                    is=1: -i  is=2: +i  ← 正确

// 下三角填充：使用共轭转置
lower_mat->get_value(icol, irow) = std::conj(upper_mat->get_value(irow, icol));
```

---

## 4. 测试方案

### 4.1 测试体系

使用 bcc Fe 二聚体 + DFT+U，磁化方向沿 (1,1,1)，开启 SOC：

```
STRU:
  Fe 0.0  0.00  0.00  0.00  mag 1.0 1.0 1.0
  Fe 0.51 0.51  0.51  0.51  mag 1.0 1.0 1.0
```

该体系具有稳定的磁矩 (~3.47 μB/atom)，自洽后可产生显著的
非对角哈密顿量矩阵元。

### 4.2 验证条件

**条件 1：H(R=0) 的 Hermiticity**

```
max|H(R=0) - H†(R=0)| < 1e-10
```

修复前：~1e-2（因下三角缺少共轭）
修复后：~1e-19（机器精度）

**条件 2：非对角元相位符号**

对磁化方向 m = (mx, my, mz)，XC 势的非对角元为：
```
H_{↑,↓} ∝ mx - i*my
```

当 my > 0 时，Im(H_{↑,↓}) < 0。

测试中 m 沿 (1,1,1)，故 Im(H_{↑,↓}) 应为负值。

修复前：Im(H_{↑,↓}) > 0（错误）
修复后：Im(H_{↑,↓}) < 0（正确）

**条件 3：自旋纹理一致性（定性）**

对角化 H(k=0)，计算占据态的自旋期望值 `<σ>`。
修复后，`<σ>` 应与输入磁化方向 (1,1,1) 一致。
修复前，`<σ>` 指向共轭后的方向。

### 4.3 测试文件

| 文件 | 说明 |
|------|------|
| `tests/03_NAO_multik/verify_hamiltonian_convention/` | 测试目录 |
| `check_hamiltonian_convention.py` | Python 验证脚本 |
| `tests/03_NAO_multik/verify_dm_symmetry_soc/` | 复用已有的 SOC 测试算例 |

### 4.4 运行方式

```bash
# 运行 SCF 计算（需要 out_mat_hs2=1 输出 H(R)）
cd tests/03_NAO_multik/verify_dm_symmetry_soc
mpirun -np 4 abacus > log.txt 2>&1

# 运行验证脚本
python3 ../verify_hamiltonian_convention/check_hamiltonian_convention.py OUT.verify_soc/
```

预期输出：
```
TEST 1: H(R=0) Hermiticity: H = H^dagger
  max|H - H^dagger| = 1.73e-19
  [PASS]

TEST 2: Off-diagonal phase for m along y
  Mean Im(H_{up,down}) = -1.59e-03
  [PASS] Im(H_{up,down}) < 0, consistent with correct convention
```

---

## 5. PW vs LCAO 对比

| 代码路径 | 文件 | 状态 |
|---------|------|------|
| PW | `source/source_pw/module_pwdft/kernels/veff_op.cpp` | **正确** |
| LCAO | `source/source_lcao/module_gint/gint_common.cpp` | **已修复** |

PW 代码直接实现了正确的公式：
```cpp
// veff_op.cpp: 正确
sup   = out*(in[0]+in[3]) + out1*(in[1] - i*in[2]);  // V_0+V_z, V_x-i*V_y
sdown = out1*(in[0]-in[3]) + out*(in[1] + i*in[2]);  // V_0-V_z, V_x+i*V_y
```

LCAO 代码通过 `merge_hr_part_to_hR()` 间接转换，修复前符号错误。

---

## 6. 影响范围

- 所有 nspin=4 的 LCAO 计算（包括 SOC 和非共线无 SOC）
- CPU 和 GPU 路径均受影响（GPU 代码调用同一 `merge_hr_part_to_hR` 函数）
- regular 和 meta-GGA 均受影响

nspin=1,2（共线自旋）不受影响，因为此时非对角元为零。
