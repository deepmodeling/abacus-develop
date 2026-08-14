# LCAO 力/应力中 DMR 配对约定推导（DM/DMK/DMR 链路）

- 日期：2026-08-09
- 相关代码：`source/source_estate/module_dm/`（DMK/DMR）、
  `source/source_hamilt/module_hcontainer/func_folding.cpp`（folding_HR）、
  `source/source_lcao/module_operator_lcao/`（overlap/ekinetic/nonlocal 力与应力）
- 配套文档：`docs/dm_dmk_dmr_investigation.md`（排查）、
  `docs/dm_dmk_dmr_action_plan.md`（整改与测试方案）

本文档锁定 DMK→DMR 的傅里叶符号约定，以及力/应力计算中实空间密度矩阵
与两中心积分（含 nonlocal 的 R2−R1 配对）的收缩对应关系，作为代码注释
（`density_matrix.cpp` 中 "See Sec. 3 ..."）所指的推导依据。

## 1. 指标与 R 格矢约定

- 波函数系数 `C_{μn}(k)`：μ = 轨道（行），n = 能带（列），存储为
  `psi::Psi` 的 `wfc(ib, iw)`。
- 密度矩阵元素 `D_{μν}`：μ = 行（bra）轨道，ν = 列（ket）轨道。
- 格矢 R：`HContainer::AtomPair(iat1, iat2, R)` 的矩阵块为
  `O_{μν}(R)`，μ∈iat1（0 胞）、ν∈iat2（R 胞）；格矢满足
  `dtau = τ(iat2) + R·L − τ(iat1)`（`UnitCell::cal_dtau`）。
- k、R 均为**直接坐标**，相位 `k·R` 前乘 `2π`。

## 2. k 空间与实空间傅里叶对

- 正变换（算符，`folding_HR`，`func_folding.cpp`）：

  ```
  O(k) = Σ_R e^{+ikR} O(R)
  ```

- 逆变换（密度矩阵，`DensityMatrix_Tools::cal_DMR*`，`density_matrix.cpp`）：

  ```
  D(R) = Σ_k e^{-ikR} DMK(k)
  ```

- 归一化：k 点权重 `w_k` 已内嵌在 `DMK`（`wg = w_k·occ`，
  `occupy.cpp`），因此 `cal_DMR` **不**再乘 `1/Nk`。
- `DMK(μ,ν;k) = Σ_n f_nk C*_{μn}(k) C_{νn}(k) = (C f C†)^T = D_std^T`
  （`cal_dm_psi.cpp`，`zgemm('N','T')` + 先取共轭）。

## 3. 正逆变换互为逆（Sec. 3 约定）

对同一套 R 定义，`e^{+ikR}`（folding_HR）与 `e^{-ikR}`（cal_DMR）构成
自洽的正逆傅里叶对。若把 `density_matrix.cpp` 中的 `-sinp` 改回 `+sinp`，
则 `D(R)` 实际存的是正变换 `Σ_k e^{+ikR} DMK(k)`，与 `folding_HR` 不再
互为逆变换，非对角/非对称消费方（如 nspin=4 nonlocal 力）即出错。
该符号由单元测试 `T1_fourier_round_trip` 锁定（见
`docs/dm_dmk_dmr_action_plan.md` §T1）。

## 4. 力/应力的一般形式

Feynman–Hellmann 型受力（overlap/ekinetic/nonlocal 等）：

```
F = Σ_k Tr[ D(k) · ∂O(k)/∂τ ]
  = Σ_k Σ_{μν} D_{μν}(k) ∂O_{νμ}(k)/∂τ
```

代入 `∂O(k) = Σ_R e^{+ikR} ∂O(R)`：

```
F = Σ_R Σ_{μν} ∂O_{νμ}(R) D_{μν}(−R)          (*)
   其中 D(−R) := Σ_k e^{+ikR} DMK(k)
```

代码存储的是 `D(R) = Σ_k e^{-ikR} DMK(k)`，两者通过 §5 的闭合迹论证等价。

## 5. 闭合迹共轭保护（Sec. 5 配对）

代码在两中心积分路径中，对每个 `(iat1,iat2,R)` 同指标收缩：

```
W_code = Σ_{μνR} D(μ,ν;R) ∂O(μ,ν;R)
```

把 `D(R) = Σ_k e^{-ikR} DMK(k)` 代入：

```
W_code = Σ_k Σ_{μν} DMK(μ,ν;k) [Σ_R e^{-ikR} ∂O(μ,ν;R)]
       = Σ_k Tr( DMK(k) · ∂O'(−k) )
```

其中 `∂O'(−k) := Σ_R e^{-ikR} ∂O(R)` 是 folding 在 −k 处的取值。
若 `∂O` 在 HContainer 中**全方向成对**（对每个 `(iat1,iat2,R)` 存在
`(iat2,iat1,−R)`，且 `∂O(iat2,iat1,−R) = ∂O(iat1,iat2,R)†`），则
`∂O'(k) = Σ_R e^{+ikR} ∂O(R)` 逐 k 厄米，`∂O'(−k) = ∂O'(k)†`。于是

```
W_code = Σ_k w_k Tr( DMK(k) · ∂O'(k)† )
       = Σ_k w_k conj( Tr( DMK(k) ∂O'(k) ) )
       = Σ_k w_k Tr( DMK(k) ∂O'(k) )        [厄米矩阵对的迹为实数]
       = Σ_k w_k Re Tr( DMK(k) ∂O'(k) )
```

即**总和受共轭保护**，而非"逐项相等"。保护成立的三条件：

1. 收缩是**全轨道闭合 Frobenius 迹**（不按元素取用）；
2. 配对算符**逐 k 厄米**（HContainer 全方向成对存储）；
3. 只取**实部**。

破坏任一条件的消费方即暴露于 O(1) 配对误差。

### 5.1 overlap / ekinetic（两中心积分路径）

`operator_fs_utils.hpp` 的 `cal_force_stress_2center` 对
`(iat1,iat2,R)` 直接取 `D(μ,ν;R) ∂O(μ,ν;R)` 同指标收缩，且对
`(iat2,iat1,−R)` 同样各计一次（全 R 集迭代），故
`finalize_force_stress(force_factor = 1.0)`（无需因子 2）。

### 5.2 nonlocal（β 投影，R2−R1 配对）

`nonlocal_fs.cpp::cal_force_stress` 以中心原子 `iat0` 为基准枚举两个
邻原子：`iat1`（胞 `R1`）、`iat2`（胞 `R2`），取用的 DMR 块是

```
dmR->find_matrix(iat1, iat2, R2 − R1)
```

即配对结构为 `R_vector = R2 − R1`（`nonlocal_fs.cpp` 中
`R_index2 − R_index1`），不是简单同 R 闭合迹。该路径的受力公式为

```
F(iat0) = Σ_{iat1,iat2} <∂β_{iat1,R1}/∂τ_iat0 | D_{iat1,iat2}(R2−R1) | β_{iat2,R2}>
```

逐项按 `(iat1,iat2,R2−R1)` 取用复数 DMR 块（nspin=4 时由
`cal_DMR_full` 提供，不丢虚部）。该路径**不满足** §5 的同 R 闭合迹
保护条件，其正确性依赖：
1. `cal_DMR_full` 采用 `e^{-ikR}`（与 `folding_HR` 的 `e^{+ikR}` 互逆）；
2. DMR 块在 `(iat1,iat2,R2−R1)` 上的取值与 β 投影的指标顺序一致；
3. 数值上由 SOC 有限差分力测试锁定（见 `docs/dm_dmk_dmr_action_plan.md` §T4）。

若未来引入半集存储优化（只存 `(iat1,iat2,R)` 不存 `(iat2,iat1,−R)`），
§5 的保护论证与 `force_factor=1.0` 约定同时失效，须恢复 factor=2 并
重跑 T2/T3/T8 守卫。
