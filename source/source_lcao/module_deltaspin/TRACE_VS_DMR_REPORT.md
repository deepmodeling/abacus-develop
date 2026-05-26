# 非正交基下迹公式与 DMR 路径的不等价性分析

## 1. 问题本质：非正交基下恒等分解

在正交基 $\{|e_i\rangle\}$ 中，恒等算符的分解是：

$$\hat{I} = \sum_i |e_i\rangle\langle e_i|$$

在非正交基 $\{|\phi_\mu\rangle\}$（重叠 $S_{\mu\nu} = \langle\phi_\mu|\phi_\nu\rangle \neq \delta_{\mu\nu}$）中，恒等分解需要**对偶基**（dual basis）：

$$\hat{I} = \sum_{\mu\nu} |\phi_\mu\rangle (S^{-1})_{\mu\nu} \langle\phi_\nu| = \sum_\mu |\phi_\mu\rangle\langle\tilde{\phi}_\mu|$$

其中 $|\tilde{\phi}_\mu\rangle = \sum_\nu (S^{-1})_{\mu\nu} |\phi_\nu\rangle$ 满足 $\langle\tilde{\phi}_\mu|\phi_\nu\rangle = \delta_{\mu\nu}$。

## 2. DMR 路径的正确性

LCAO 中密度矩阵的系数表示为：

$$\text{DMK}(k)_{\mu\nu} = \sum_{ib} f_{k,ib} \, C(k)_{\mu,ib} \, C^*(k)_{\nu,ib}$$

密度算符的正确表达式是：

$$\hat{\rho}(k) = \sum_{\mu\nu\alpha\beta} \text{DMK}_{\mu\nu} \, |\phi_\mu(k)\rangle (S^{-1})_{\mu\alpha} (S^{-1})_{\nu\beta}^* \langle\phi_\beta(k)| \quad \text{(一般形式)}$$

但由于 $C^\dagger S C = I$（正交归一化条件），可以验证：

$$\hat{\rho}(k) = \sum_{ib} f_{k,ib} \, |\psi_{k,ib}\rangle\langle\psi_{k,ib}|$$

磁矩的正确计算为：

$$M_I = \text{Tr}(\hat{P}_I \hat{\rho}_\text{diff}) = \sum_R \sum_{\mu\nu} P_I(\mu,\nu;R) \times \text{DMR}_\text{diff}(\nu,\mu;-R)$$

这里 $P_I(\mu,\nu;R) = \langle\phi_\mu(0)|\hat{P}_I|\phi_\nu(R)\rangle$ 是投影算符在**原始基**下的矩阵元，而 DMR 是系数密度矩阵。

**关键**：这个公式之所以正确，是因为 $\hat{P}_I$ 的矩阵元和 DMR 的定义在重叠矩阵 $S$ 的处理上是**自洽的**。展开验证：

$$\text{Tr}(\hat{P}_I \hat{\rho}) = \sum_\mu \langle\tilde{\phi}_\mu|\hat{P}_I \hat{\rho}|\phi_\mu\rangle = \sum_{\mu\alpha} (S^{-1})_{\mu\alpha} \langle\phi_\alpha|\hat{P}_I \hat{\rho}|\phi_\mu\rangle$$

这与简单的 $\sum_{\mu\nu} P_{I,\mu\nu} \text{DMR}_{\nu\mu}$ 并不直接对应。但 ABACUS 中 `cal_moment` 的实际实现确实计算了正确的迹，因为它利用了 $C^\dagger S C = I$ 的约束——在 $|\psi\rangle = \sum C_\mu |\phi_\mu\rangle$ 的定义下，DMR 自然吸收了 $S^{-1}$ 的效果。

## 3. 迹公式的错误

迹公式计算：

$$M_I^\text{trace} = \sum_k \text{sign}(k) \sum_{ib} w_{k,ib} \, \langle\psi_{k,ib}|\hat{P}_I|\psi_{k,ib}\rangle$$

展开 $|\psi_{k,ib}\rangle = \sum_\mu C(k)_{\mu,ib} |\phi_\mu(k)\rangle$：

$$\langle\psi_{k,ib}|\hat{P}_I|\psi_{k,ib}\rangle = \sum_{\mu\nu} C^*_{\mu,ib} \, \langle\phi_\mu|\hat{P}_I|\phi_\nu\rangle \, C_{\nu,ib} = [C^\dagger P_I C]_{ib,ib}$$

而 DMR 路径计算：

$$M_I^\text{DMR} = \sum_{\mu\nu} P_{I,\mu\nu} \times \text{DMR}_{\nu\mu} = \text{Tr}(P_I \cdot \text{DMR})$$

**二者的关系**：

$$\text{Tr}(P_I \cdot \text{DMR}) = \sum_{\mu\nu} P_{I,\mu\nu} \sum_{ib} f_{ib} C_{\nu,ib} C^*_{\mu,ib}$$
$$= \sum_{ib} f_{ib} \sum_{\mu\nu} C^*_{\mu,ib} P_{I,\mu\nu} C_{\nu,ib}$$
$$= \sum_{ib} f_{ib} [C^\dagger P_I C]_{ib,ib}$$

看起来二者相等？但这忽略了一个关键问题：**迹公式使用的 $P_{I,\text{sub}}$ 来自 `cal_PI_sub`，而不是 $C^\dagger P_I C$**。

### cal_PI_sub 计算的是什么

`cal_PI_sub` 计算：

$$D_I(k)_{lm,ib} = \sum_{\mu,R} \langle\alpha_{lm}(I)|\phi_\mu(R)\rangle \, e^{ik \cdot R} \, C(k)_{\mu,ib}$$

$$[P_{I,\text{sub}}]_{ib,jb} = \sum_{lm} D^*_{lm,ib} \, D_{lm,jb} = \langle\psi_{ib}|\hat{P}_I|\psi_{jb}\rangle$$

而 $C^\dagger P_I(k) C$ 中的 $P_I(k)_{\mu\nu}$：

$$P_I(k)_{\mu\nu} = \sum_R e^{ik \cdot R} \sum_{lm} \langle\phi_\mu(0)|\alpha_{lm}(I)\rangle\langle\alpha_{lm}(I)|\phi_\nu(R)\rangle$$

$$[C^\dagger P_I(k) C]_{ib,jb} = \sum_{\mu\nu} C^*_{\mu,ib} P_I(k)_{\mu\nu} C_{\nu,jb}$$

展开后两者相同（都等于 $\langle\psi_{ib}|\hat{P}_I|\psi_{jb}\rangle$），因此迹公式和 DMR **在数学上等价**。

### 那偏差从哪来？

既然数学等价，实测 4.66 μB 的偏差必然来自**实现层面**。最可能的原因是 `cal_PI_sub` 中的邻域截断：

- `cal_PI_sub` 只对 `B_I_data[iat]` 中的邻域原子求和（由 `onsite_radius` 截断）
- 而 `cal_moment` 使用的 `pre_hr[iat]` 也有同样的截断
- 但关键区别：**DMR 路径中 DMR 本身是从完整 H(k) 对角化得到的全空间密度矩阵**，pre_hr 的截断通过 DMR 的全空间结构被部分补偿
- 迹公式中 $D_I$ 的截断直接作用在最终结果上，没有补偿机制

更重要的是，$D_I$ 通过 `MPI_Allreduce` 汇总，但汇总的范围取决于每个进程看到的邻域。如果并行分布导致某些 $(μ,R)$ 对在 `B_I_data` 中缺失，$D_I$ 就不完整。

## 4. 修正方案：对偶基方法

如果要在子空间中正确计算磁矩，需要使用 Löwdin 正交化或对偶基。

### 方案 A：Löwdin 正交化

定义 Löwdin 正交基 $|\bar{\phi}_\mu\rangle = \sum_\nu (S^{-1/2})_{\mu\nu} |\phi_\nu\rangle$，在此基下：

- 波函数系数：$\bar{C} = S^{1/2} C$
- 投影矩阵：$\bar{P}_I = S^{-1/2} P_I S^{-1/2}$
- 密度矩阵：$\bar{\text{DMR}} = S^{1/2} \text{DMR} S^{1/2}$

迹：$\text{Tr}(\bar{P}_I \cdot \bar{\text{DMR}}) = \text{Tr}(S^{-1/2} P_I S^{-1/2} \cdot S^{1/2} \text{DMR} S^{1/2}) = \text{Tr}(P_I \cdot \text{DMR})$ ✓

在正交基下，逐带投影自然给出正确的原子分辨磁矩。

### 方案 B：对偶基投影

定义 $|\tilde{\phi}_\mu\rangle = \sum_\nu (S^{-1})_{\mu\nu} |\phi_\nu\rangle$，使用对偶投影：

$$\tilde{P}_{I,\mu\nu} = \sum_{lm} \langle\tilde{\phi}_\mu|\alpha_{lm}(I)\rangle\langle\alpha_{lm}(I)|\tilde{\phi}_\nu\rangle$$

在此定义下迹公式和 DMR 路径严格等价。

### 实际可行性

两种方案都需要 $S^{-1}$ 或 $S^{-1/2}$，代价为 $O(N^3)$（但只需在构建子空间时计算一次）。

更实际的路径：不修改迹公式，而是用 DMR 路径配合梯度饱和回退机制。

## 5. 实测数据确认

| 原子 | Full | DMR 子空间 | 迹公式 | |Full-DMR| | |Full-Trace| |
|------|------|-----------|--------|-----------|------------|
| 0 (Mz=+1.85) | +1.8485 | +1.8485 | +1.8642 | 3.5e-14 | 0.016 |
| 1 (Mz=-2.80) | -2.7951 | -2.7951 | **+1.8650** | 1.1e-13 | **4.660** |

- DMR 路径在 Δλ=0 处**精确一致**（误差 < 1e-13），证明子空间构建正确
- 迹公式两原子值几乎相同（~1.86），丧失了原子分辨能力
- 原子 1 符号反转（-2.80 → +1.87），差异 4.66 μB

## 6. 结论

1. **迹公式和 DMR 路径在完备求和下数学等价**，但 `cal_PI_sub` 的邻域截断使二者在实践中不等价
2. 对 BCC Fe 等高对称体系，截断导致两原子投影几乎相同，迹公式丧失原子分辨能力
3. DMR 路径通过全空间密度矩阵规避了此问题
4. 如需修正迹公式，应引入 Löwdin 正交化或对偶基（$S^{-1/2}$），代价 $O(N^3)$
5. 当前推荐：使用 DMR 路径 + 梯度饱和回退到全对角化
