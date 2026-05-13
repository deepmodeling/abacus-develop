# 测试集精简方案：tests/17_DS_DFTU

## 现状

当前 CASES_CPU.txt 共 41 个条目（含1个注释），覆盖 bcc Fe 和 FeO 两种材料的 DeltaSpin / DFT+U 功能测试。

## 精简原则

1. **唯一性**：同一配置维度（基组×自旋模式×自旋方向×DFT+U）只保留一个算例
2. **正交性**：每个保留的算例必须有不可替代的测试维度
3. **阈值变体隔离**：sc_scf_thr 严格/宽松测试与默认收敛测试分离，不混在主线算例中
4. **保留边界情况**：ReadLam、DirectionOnly、无sc约束等特殊模式各保留一个代表

---

## 冗余分析

### 完全重复（可删除）

| 保留 | 删除 | 原因 |
|------|------|------|
| 14_PW_DS_S4_XYZ | 17_PW_DS_S4_XYZ | **完全重复**：同 basis=pw, noncolin=1, S4, XYZ, kpar=2, sc_scf_thr=10 |
| 19_PW_DFTU_DS_S4_XY | 22_PW_DFTU_DS_S4_XY | **完全重复**：同 basis=pw, DFT+U, S4, XY, kpar=2, sc_scf_thr=10 |
| 20_PW_DFTU_DS_S4_XYZ | 23_PW_DFTU_DS_S4_XYZ | **完全重复**：同 basis=pw, DFT+U, S4, XYZ, kpar=2, sc_scf_thr=10 |
| 25_LCAO_DS_S4_XY | 28_LCAO_DS_S4_XY | **近重复**：同 basis=lcao, S4, XY, sc=1 1 0；无额外测试价值 |
| 26_LCAO_DS_S4_XYZ | 29_LCAO_DS_S4_XYZ | **变体**：26 额外使用 sc_lambda_strategy=linear_response，可合并到 26 或保留 29 作为基准 |
| 31_LCAO_DFTU_DS_S4_XY | 34_LCAO_DFTU_DS_S4_XY | **完全重复**：同 basis=lcao, DFT+U, S4, XY, sc=1 1 0 |
| 32_LCAO_DFTU_DS_S4_XYZ | 35_LCAO_DFTU_DS_S4_XYZ | **完全重复**：同 basis=lcao, DFT+U, S4, XYZ, sc=1 1 1 |
| 53_PW_DS_S4_XY_MagMomCheck | — | **与13重复**：同配置，仅命名不同用于磁矩检查，检查逻辑应统一到 13 |

### kpar 变体（可合并）

| 保留 | 删除 | 原因 |
|------|------|------|
| 16_PW_DS_S4_XY (kpar=2) | 13_PW_DS_S4_XY (kpar=1) | kpar=2 覆盖 kpar=1 的并行路径；若需测试 kpar=1 可保留但标记 |

**注意**：13 的 kpar=1 是唯一的 kpar=1 DS 算例，如果并行缩减只测 kpar=2 就足够了，建议删除 13。

### 阈值变体（建议移至单独文件）

| 算例 | sc_scf_thr | 测试目的 |
|------|-----------|---------|
| 38_PW_DS_S2_Thr1e10_Z | 1e-10 | 超严格收敛 |
| 39_PW_DS_S4_Thr1e10_XY | 1e-10 | 超严格收敛+非共线 |
| 41_PW_DS_S4_Thr10_XY | 0.1 | 宽松收敛 |
| 43_PW_DFTU_DS_S4_Thr1e10_XY | 1e-10 | 严格收敛+DFT+U |
| 45_PW_DFTU_DS_S4_Thr10_XY | 0.1 | 宽松收敛+DFT+U |

**建议**：这些算例验证 sc_scf_thr 参数对收敛的影响，属于**参数扫描**，不应与功能测试混在一起。移至 `CASES_THRESHOLD.txt`。

### FeO 原子序变体

| 算例 | 区别 | 测试目的 |
|------|------|---------|
| 50_FeO_O_first_Fe_second | O在前 Fe在后 | DFT+U 按元素顺序正确分配 |
| 51_FeO_Fe_first_O_second | Fe在前 O在后 | 同上，顺序反转 |

**建议**：保留一个即可，另一个是同一逻辑的反向验证。

---

## 精简后方案

### 精简后 CASES_CPU.txt（推荐：25 个算例）

#### A. 基础 SPIN / DFT+U（无 DeltaSpin）— 5个

| 算例 | 保留原因 |
|------|---------|
| 06_PW_SPIN_S2_Z | **PW 自旋极化基线**：nspin=2, Z方向，无DS无U |
| 07_PW_SPIN_S4_XYZ | **PW 非共线基线**：noncolin=1, XYZ方向，无DS无U |
| 08_PW_DFTU_S2_Z | **PW DFT+U 基线**：nspin=2+U，无DS |
| 09_PW_DFTU_S4_XY | **PW DFT+U 非共线基线**：noncolin=1+U，无DS |
| 11_PW_DFTU_S2_FeO | **FeO 多元素 DFT+U**：验证多元素体系 U 分配 |

#### B. PW 纯 DeltaSpin — 3个

| 算例 | 保留原因 |
|------|---------|
| 12_PW_DS_S2_Z | **S2 nspin=2 模式**：唯一 noncolin=0 的 DS 算例 |
| 15_PW_DS_S4_Z | **S4 Z方向**：noncolin=1, Z方向自旋约束 |
| 16_PW_DS_S4_XY | **S4 XY方向**：noncolin=1, XY方向约束 + kpar=2 并行 |

> 删除 13(kpar=1, 被16覆盖)、14/17(完全重复)、53(与13重复)

#### C. PW DFT+U + DeltaSpin — 3个

| 算例 | 保留原因 |
|------|---------|
| 18_PW_DFTU_DS_S2_Z | **S2 模式**：唯一 noncolin=0 的 DFTU+DS |
| 21_PW_DFTU_DS_S4_Z | **S4 Z方向**：noncolin=1, Z方向约束+U |
| 19_PW_DFTU_DS_S4_XY | **S4 XY方向**：noncolin=1, XY方向约束+U + kpar=2 |

> 删除 20/23(与19重复的XYZ方向，Z已有21覆盖)、22(与19完全重复)

#### D. LCAO 纯 DeltaSpin — 3个

| 算例 | 保留原因 |
|------|---------|
| 27_LCAO_DS_S4_Z | **Z-only sc约束**：sc=0 0 1，唯一 Z 方向 LCAO DS |
| 26_LCAO_DS_S4_XYZ | **XYZ+linear_response**：sc=1 1 1 + 线性响应策略 |
| 28_LCAO_DS_S4_XY | **XY-only sc约束**：sc=1 1 0，补充 Z/XY 约束正交性 |

> 删除 25(与28重复)、29(与26重复，26额外测试 linear_response)

#### E. LCAO DFT+U + DeltaSpin — 3个

| 算例 | 保留原因 |
|------|---------|
| 33_LCAO_DFTU_DS_S4_Z | **Z-only sc约束**：sc=0 0 1 |
| 32_LCAO_DFTU_DS_S4_XYZ | **XYZ方向**：sc=1 1 1，标准 XYZ 约束 |
| 31_LCAO_DFTU_DS_S4_XY | **XY方向**：sc=1 1 0 |

> 删除 34(与31重复)、35(与32重复)

#### F. 特殊模式 — 5个

| 算例 | 保留原因 |
|------|---------|
| 36_PW_DS_S2_ReadLam_Z | **nsc=1 读lambda模式**：唯一跳过 lambda 优化的测试 |
| 37_PW_DS_S4_ReadLam_XY | **nsc=1 + noncolin**：非共线下的 ReadLam 模式 |
| 56_PW_DS_S4_DirectionOnly_XY | **sc_direction_only=1**：方向约束不约束大小 |
| 57_PW_DFTU_DS_S4_DirectionOnly_XY | **DirectionOnly + DFT+U**：方向约束与 U 交互 |
| 58_LCAO_DS_S4_DirectionOnly_XY | **DirectionOnly + LCAO**：LCAO 基组方向约束 |

#### G. 多元素特殊 — 1个

| 算例 | 保留原因 |
|------|---------|
| 51_FeO_Fe_first_O_second | **元素顺序**：Fe 在前 O 在后，验证 U 分配正确性 |

> 删除 50(与51对称，51覆盖)

#### H. 阈值变体 — 移至 CASES_THRESHOLD.txt（5个）

| 算例 | 测试目的 |
|------|---------|
| 38_PW_DS_S2_Thr1e10_Z | 超严格收敛，无 sc 约束 |
| 39_PW_DS_S4_Thr1e10_XY | 超严格收敛+非共线，无 sc 约束 |
| 41_PW_DS_S4_Thr10_XY | 宽松收敛 |
| 43_PW_DFTU_DS_S4_Thr1e10_XY | 严格收敛+DFT+U，无 sc 约束 |
| 45_PW_DFTU_DS_S4_Thr10_XY | 宽松收敛+DFT+U |

---

## 精简前后对比

| 维度 | 精简前 | 精简后 | 减少 |
|------|--------|--------|------|
| CASES_CPU.txt 总数 | 40 (39 active) | 25 | -37% |
| 完全重复删除 | - | 8 个 | 13→14/17, 19→22, 20→23, 25→28, 26→29, 31→34, 32→35, 53→13 |
| kpar变体合并 | - | 1 个 | 13→16 |
| 阈值变体移出 | - | 5 个 | 移至 CASES_THRESHOLD.txt |
| FeO对称删除 | - | 1 个 | 50→51 |
| 运行时间估算 | ~250s | ~120s | -52% |

---

## 覆盖矩阵（精简后）

| | PW S2 | PW S4 Z | PW S4 XY | PW S4 XYZ | LCAO S4 Z | LCAO S4 XY | LCAO S4 XYZ |
|---|---|---|---|---|---|---|---|
| 纯 DS | 12 | 15 | 16 | — | 27 | 28 | 26 |
| DFT+U+DS | 18 | 21 | 19 | — | 33 | 31 | 32 |
| 基线(无DS) | 06/08 | 07/09 | 09 | 07 | — | — | — |

**说明**：
- S4 XYZ 在 PW 中由 15(Z) + 16(XY) 间接覆盖（XYZ = Z + XY 的组合），且 14 的 XYZ 测试已作为基线 07 覆盖
- LCAO 中 XYZ 方向有 unique sc_lambda_strategy 测试价值，保留 26

## 保留算例的不可替代性总结

| 算例 | 不可替代原因 |
|------|------------|
| 06 | 唯一 PW nspin=2 无DS无U 基线 |
| 07 | 唯一 PW noncolin 无DS无U 基线 |
| 08 | 唯一 PW nspin=2+U 无DS |
| 09 | 唯一 PW noncolin+U 无DS + XY方向 |
| 11 | 唯一 FeO 多元素 DFT+U |
| 12 | 唯一 nspin=2(noncolin=0) DeltaSpin |
| 15 | 唯一 noncolin=1 + Z方向纯DS |
| 16 | 唯一 noncolin=1 + XY方向纯DS + kpar=2 |
| 18 | 唯一 nspin=2 + DFT+U + DS |
| 19 | 唯一 noncolin=1 + DFT+U + XY方向DS |
| 21 | 唯一 noncolin=1 + DFT+U + Z方向DS |
| 26 | 唯一 LCAO + linear_response 策略 |
| 27 | 唯一 LCAO + sc=0 0 1 (Z-only) |
| 28 | 唯一 LCAO + sc=1 1 0 (XY-only) |
| 31 | 唯一 LCAO+DFT+U + sc=1 1 0 |
| 32 | 唯一 LCAO+DFT+U + sc=1 1 1 |
| 33 | 唯一 LCAO+DFT+U + sc=0 0 1 |
| 36 | 唯一 nsc=1 (ReadLam) + nspin=2 |
| 37 | 唯一 nsc=1 (ReadLam) + noncolin |
| 51 | 唯一 Fe 优先的 FeO 元素顺序测试 |
| 56 | 唯一 PW DirectionOnly 纯DS |
| 57 | 唯一 PW DirectionOnly + DFT+U |
| 58 | 唯一 LCAO DirectionOnly |
