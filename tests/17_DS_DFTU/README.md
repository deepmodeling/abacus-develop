# 17_DS_DFTU — DeltaSpin & DFT+U 集成测试集

本目录包含 ABACUS 中 **DeltaSpin（自旋约束 DFT）** 和 **DFT+U** 功能的全部集成测试用例，
涵盖 LCAO 和 PW 基组、共线/非共线自旋、DFT+U、DeltaSpin 及其组合。

## 测试清单 (52 例)

### 一、LCAO Spin (01-02)

| # | 算例 | 说明 |
|---|------|------|
| 01 | LCAO_SPIN_S2_Z | 验证 LCAO 基组下共线自旋的基础 SCF 收敛性，作为 LCAO 磁性计算的基准对照 |
| 02 | LCAO_SPIN_S4_XYZ | 验证 LCAO 基组下非共线自旋的基础 SCF 收敛性，覆盖 LCAO 非共线计算路径 |

### 二、LCAO DFT+U (03-05)

| # | 算例 | 说明 |
|---|------|------|
| 03 | LCAO_DFTU_S2_Z | 验证 LCAO 基组下 DFT+U (U=5.0eV, l=2) 与共线自旋的耦合，确保 LCAO 路径的 DFT+U 占据矩阵计算正确 |
| 04 | LCAO_DFTU_S4_XY | 验证 LCAO 基组下 DFT+U 与非共线自旋 (XY 磁矩) 的耦合，覆盖 LCAO 路径中 nspin=4 的占据矩阵计算 |
| 05 | LCAO_DFTU_S4_XYZ | 验证 LCAO 基组下 DFT+U 与非共线自旋 (XYZ 磁矩) 的耦合，覆盖 LCAO 路径的最完整占据矩阵场景 |

### 三、PW Spin (06-07)

| # | 算例 | 说明 |
|---|------|------|
| 06 | PW_SPIN_S2_Z | 验证 PW 基组下共线自旋的基础 SCF 收敛性，作为 PW 磁性计算的基准对照 |
| 07 | PW_SPIN_S4_XYZ | 验证 PW 基组下非共线自旋的基础 SCF 收敛性，覆盖 PW 非共线计算路径 |

### 四、PW DFT+U (08-11)

| # | 算例 | 说明 |
|---|------|------|
| 08 | PW_DFTU_S2_Z | 验证 PW 基组下 DFT+U (U=5.0eV, l=2) 与共线自旋的耦合，确保 PW 路径的 DFT+U 有效势计算正确 |
| 09 | PW_DFTU_S4_XY | 验证 PW 基组下 DFT+U 与非共线自旋 (XY 磁矩) 的耦合，覆盖 PW 路径中 nspin=4 的 onsite 投影矩阵 |
| 10 | PW_DFTU_S4_XY | 与 09 相同参数但不同晶体结构，验证 PW DFT+U 非共线在不同晶格下的泛化能力 |
| 11 | PW_DFTU_S2_FeO | 验证 PW 基组下 DFT+U 在 FeO 体系上的正确性，确保 Fe-3d 轨道的 DFT+U 修正有效 |

### 五、PW DeltaSpin (12-17)

| # | 算例 | 说明 |
|---|------|------|
| 12 | PW_DS_S2_Z | 验证 PW 基组下 DeltaSpin 与共线自旋的耦合，确保 DeltaSpin 迭代优化磁矩到目标值的正确性 |
| 13 | PW_DS_S4_XY | 验证非共线 DeltaSpin 在 XY 磁矩约束下的迭代优化，覆盖 nspin=4 路径的 lambda 更新 |
| 14 | PW_DS_S4_XYZ | 验证非共线 DeltaSpin 在 XYZ 三方向磁矩约束下的迭代优化，覆盖最完整的自旋约束场景 |
| 15 | PW_DS_S4_Z | 验证非共线 DeltaSpin 仅约束 Z 方向磁矩时的行为，确保 noncolin=1 框架下单轴约束不引入非物理 XY 分量 |
| 16 | PW_DS_S4_XY | 与 13 相同参数但不同晶体结构，验证非共线 DeltaSpin XY 约束在不同晶格下的泛化能力 |
| 17 | PW_DS_S4_XYZ | 与 14 相同参数但不同晶体结构，验证非共线 DeltaSpin XYZ 约束在不同晶格下的泛化能力 |

### 六、PW DFT+U + DeltaSpin (18-23)

| # | 算例 | 说明 |
|---|------|------|
| 18 | PW_DFTU_DS_S2_Z | 验证 PW 基组下 DFT+U 与 DeltaSpin 联合 (共线自旋) 的耦合，确保 U 修正与磁矩约束不冲突 |
| 19 | PW_DFTU_DS_S4_XY | 验证非共线 DFT+U+DeltaSpin 联合在 XY 磁矩约束下的耦合，覆盖两种方法在 nspin=4 路径的联合迭代 |
| 20 | PW_DFTU_DS_S4_XYZ | 验证非共线 DFT+U+DeltaSpin 联合在 XYZ 三方向磁矩约束下的耦合，覆盖最完整的联合约束场景 |
| 21 | PW_DFTU_DS_S4_Z | 验证非共线 DFT+U+DeltaSpin 联合仅约束 Z 方向磁矩时的行为，确保单轴约束与 DFT+U 有效势的正确叠加 |
| 22 | PW_DFTU_DS_S4_XY | 与 19 相同参数但不同晶体结构，验证非共线 DFT+U+DeltaSpin 联合在不同晶格下的泛化能力 |
| 23 | PW_DFTU_DS_S4_XYZ | 与 20 相同参数但不同晶体结构，验证非共线 DFT+U+DeltaSpin 联合 XYZ 约束在不同晶格下的泛化能力 |

### 七、LCAO DeltaSpin (24-29)

| # | 算例 | 说明 |
|---|------|------|
| 24 | LCAO_DS_S2_Z | 验证 LCAO 基组下 DeltaSpin 与共线自旋的耦合，确保 LCAO 密度矩阵路径的自旋约束优化正确 |
| 25 | LCAO_DS_S4_XY | 验证 LCAO 基组下非共线 DeltaSpin 在 XY 磁矩约束下的迭代优化，覆盖 LCAO 路径中 nspin=4 的磁矩投影 |
| 26 | LCAO_DS_S4_XYZ | 验证 LCAO 基组下非共线 DeltaSpin 在 XYZ 三方向磁矩约束下的迭代优化，覆盖 LCAO 路径的最完整约束场景 |
| 27 | LCAO_DS_S4_Z | 验证 LCAO 基组下非共线 DeltaSpin 仅约束 Z 方向磁矩时的行为，确保 noncolin=1 框架下单轴约束的正确性 |
| 28 | LCAO_DS_S4_XY | 与 25 相同参数但不同晶体结构，验证 LCAO 非共线 DeltaSpin XY 约束在不同晶格下的泛化能力 |
| 29 | LCAO_DS_S4_XYZ | 与 26 相同参数但不同晶体结构，验证 LCAO 非共线 DeltaSpin XYZ 约束在不同晶格下的泛化能力 |

### 八、LCAO DFT+U + DeltaSpin (30-35)

| # | 算例 | 说明 |
|---|------|------|
| 30 | LCAO_DFTU_DS_S2_Z | 验证 LCAO 基组下 DFT+U 与 DeltaSpin 联合 (共线自旋) 的耦合，确保密度矩阵路径的 U 修正与磁矩约束不冲突 |
| 31 | LCAO_DFTU_DS_S4_XY | 验证 LCAO 基组下非共线 DFT+U+DeltaSpin 联合在 XY 磁矩约束下的耦合，覆盖 LCAO 密度矩阵路径的联合约束 |
| 32 | LCAO_DFTU_DS_S4_XYZ | 验证 LCAO 基组下非共线 DFT+U+DeltaSpin 联合在 XYZ 三方向磁矩约束下的耦合，覆盖 LCAO 路径的最完整联合场景 |
| 33 | LCAO_DFTU_DS_S4_Z | 验证 LCAO 基组下非共线 DFT+U+DeltaSpin 联合仅约束 Z 方向磁矩时的行为，确保单轴约束与 DFT+U 密度矩阵的正确叠加 |
| 34 | LCAO_DFTU_DS_S4_XY | 与 31 相同参数但不同晶体结构，验证 LCAO DFT+U+DeltaSpin 联合在不同晶格下的泛化能力 |
| 35 | LCAO_DFTU_DS_S4_XYZ | 与 32 相同参数但不同晶体结构，验证 LCAO DFT+U+DeltaSpin 联合 XYZ 约束在不同晶格下的泛化能力 |

### 九、PW DeltaSpin 特殊参数 (36-41)

| # | 算例 | 说明 |
|---|------|------|
| 36 | PW_DS_S2_ReadLam_Z | 验证 `nsc=1` 模式 (直接读取 lambda 文件不迭代优化) 的正确性，确保 DeltaSpin 在非自洽 lambda 模式下仍能正确计算磁矩 |
| 37 | PW_DS_S4_ReadLam_XY | 验证非共线 DeltaSpin 的 `nsc=1` 模式，覆盖 XY 磁矩约束下的非自洽 lambda 路径 |
| 38 | PW_DS_S2_Thr1e10_Z | 验证 DeltaSpin 在极严收敛阈值 (sc_scf_thr=1e-10) 下的稳定性，确保迭代优化能收敛到高精度解 |
| 39 | PW_DS_S4_Thr1e10_XY | 验证非共线 DeltaSpin 在极严收敛阈值 (sc_scf_thr=1e-10) 下的稳定性，覆盖 XY 磁矩约束场景 |
| 40 | PW_DS_S2_Thr10_Z | 验证 DeltaSpin 在极松收敛阈值 (sc_scf_thr=10) 下的行为，测试算法在低精度要求下的鲁棒性和 out_alllog 日志输出 |
| 41 | PW_DS_S4_Thr10_XY | 验证非共线 DeltaSpin 在极松收敛阈值 (sc_scf_thr=10) 下的行为，覆盖 XY 磁矩约束的低精度场景 |

### 十、PW DFT+U + DeltaSpin 特殊参数 (42-45)

| # | 算例 | 说明 |
|---|------|------|
| 42 | PW_DFTU_DS_S2_Thr1e10_Z | 验证 DFT+U 与 DeltaSpin 联合在极严收敛阈值 (sc_scf_thr=1e-10) 下的迭代稳定性，确保两种方法耦合时的收敛性 |
| 43 | PW_DFTU_DS_S4_Thr1e10_XY | 验证非共线 DFT+U+DeltaSpin 在极严收敛阈值 (sc_scf_thr=1e-10) 下的耦合稳定性，覆盖 XY 磁矩约束 |
| 44 | PW_DFTU_DS_S2_Thr10_Z | 验证 DFT+U 与 DeltaSpin 联合在极松收敛阈值 (sc_scf_thr=10) 下的行为，测试耦合算法在低精度要求下的鲁棒性 |
| 45 | PW_DFTU_DS_S4_Thr10_XY | 验证非共线 DFT+U+DeltaSpin 在极松收敛阈值 (sc_scf_thr=10) 下的行为，覆盖 XY 磁矩约束的低精度场景 |

### 十一、Relax 结构优化 (46-49)

| # | 算例 | 说明 |
|---|------|------|
| 46 | PW_DS_S2_Thr1e10_Z_bfgs | 验证 DeltaSpin 使用 BFGS 策略 (sc_lambda_strategy=bfgs) 的收敛行为，测试 BFGS 优化器在自旋约束 SCF 中的正确性 |
| 47 | PW_DS_S4_Thr1e10_XY_bfgs | 验证非共线 DeltaSpin 使用 BFGS 策略的收敛行为，覆盖 XY 磁矩约束下 BFGS 优化器的正确性 |
| 48 | PW_DFTU_DS_S2_Thr10_Z_bfgs | 验证 DFT+U 与 DeltaSpin 联合使用 BFGS 策略的收敛行为，测试 BFGS 在 DFT+U+DS 耦合场景中的正确性 |
| 49 | PW_DFTU_DS_S4_Thr10_XY_bfgs | 验证非共线 DFT+U+DeltaSpin 联合使用 BFGS 策略的收敛行为，覆盖 XY 磁矩约束下 BFGS 优化器的正确性 |

### 十二、FeO 原子顺序 (50-51)

| # | 算例 | 说明 |
|---|------|------|
| 50 | FeO_O_first_Fe_second | 验证 FeO 体系中 O 原子类型在前、Fe 在后的排序下 DFT+U 的正确性，确保原子类型顺序不影响 DFT+U 的 onsite 投影 |
| 51 | FeO_Fe_first_O_second | 验证 FeO 体系中 Fe 原子类型在前、O 在后的排序下 DFT+U 的正确性，与 50 对比确保 eff_pot_pw_index 索引计算与原子类型顺序无关 |

### 十三、SOC + DFT+U (52)

| # | 算例 | 说明 |
|---|------|------|
| 52 | PW_DFTU_SO | 验证 DFT+U 与自旋轨道耦合 (SOC) 同时开启时的兼容性，确保 DFT+U 的 onsite 投影与 SOC 的自旋混合正确耦合 |

## 运行方式

```bash
# 运行全部测试
cd tests/17_DS_DFTU
bash ../integrate/Autotest.sh -a <abacus路径> -n 4

# 运行单个测试
cd 08_PW_DFTU_S2_Z
bash ../../integrate/run_debug.sh ""
```

## 已知问题

- 19-23: PW DFT+U + DeltaSpin + 非共线 → port 和 zdy-tmp 均崩溃（上游 bug）

## 测试条件说明

- 09/10 (PW DFT+U + 非共线): 仅支持 **2 进程 MPI** 运行，已提供 `result.ref` 参考文件
