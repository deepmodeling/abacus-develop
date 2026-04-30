# SCANL 有限差分 Laplacian 实现

## 1. 问题

SCANL (deorbitalized SCAN) 将 SCAN 对动能密度 τ 的依赖替换为对密度 Laplacian ∇²ρ 的依赖。其 Kohn-Sham 势包含项

$$V_{xc}^{\nabla^2\rho}(\mathbf{r}) = e^2\,\nabla^2\!\left(\frac{\partial(\rho\varepsilon_{xc})}{\partial(\nabla^2\rho)}\right) \equiv e^2\,\nabla^2(\text{vlapl})$$

在平面波基组下，倒易空间中 ∇² → −|G|²，因此

$$\tilde{V}_{xc}^{\nabla^2\rho}(\mathbf{G}) = -e^2\,|\mathbf{G}|^2\,\widetilde{\text{vlapl}}(\mathbf{G})$$

|G|² 无界增长导致 SCF 发散：vlapl(G) 的高频数值噪声被 |G|² 放大后反馈到密度，形成正反馈。此发散不依赖 mixing 策略，是 Laplace 逆变换的数学病态性。

---

## 2. 有限差分 Laplacian

### 2.1 核心思想

FFT 网格上的密度是离散函数。对离散函数求 Laplacian，有限差分 (FD) 是自然且唯一自洽的方法。FD Laplacian 的 G 空间 kernel 在低 G 时与谱 kernel 一致，在高 G 时有界衰减，不会放大高频噪声。

### 2.2 坐标系与度规

设晶格矢量 $\mathbf{a}_\alpha$（α = 1,2,3），以 lat0 为单位存储在 `rhopw->latvec` 中。分数坐标

$$u^\alpha = n_\alpha / N_\alpha, \quad n_\alpha = 0, 1, \ldots, N_\alpha - 1$$

其中 $N_\alpha$ = `rhopw->nx, ny, nz`。实空间位置

$$\mathbf{r} = \text{lat0} \sum_\alpha u^\alpha \mathbf{a}_\alpha$$

Laplacian 在分数坐标下为

$$\nabla^2 = \sum_{\alpha,\beta} g^{\alpha\beta} \frac{\partial^2}{\partial u^\alpha \partial u^\beta}$$

其中逆变度规张量

$$g^{\alpha\beta} = \frac{\text{GGT}[\alpha][\beta]}{\text{lat0}^2}$$

GGT = G · Gᵀ 是倒格子度规矩阵（ABACUS 中的 `rhopw->GGT`），G = (latvec)⁻ᵀ。

### 2.3 平面波的谱 Laplacian

对平面波 $f(\mathbf{r}) = \exp(2\pi i \sum_\gamma m_\gamma u^\gamma)$，其中 $m_\gamma$ 为 Miller 指数（`rhopw->gdirect[ig]`）：

$$\frac{\partial^2 f}{\partial u^\alpha \partial u^\beta}\bigg|_{\text{spectral}} = -(2\pi)^2 m_\alpha m_\beta \, f$$

组合得谱 Laplacian：

$$\nabla^2_{\text{spectral}} f = -\underbrace{\left[\sum_{\alpha,\beta} \text{GGT}[\alpha][\beta] \cdot m_\alpha m_\beta\right]}_{\text{gg}_{\text{spectral}}} \cdot \text{tpiba2} \cdot f$$

其中 tpiba2 = (2π/lat0)²。这就是 ABACUS 中 `gg[ig] * tpiba2` 的来源。

### 2.4 平面波的 FD Laplacian

在分数坐标的均匀网格上对二阶导数做中心差分：

**对角项** (α = β)：

$$\frac{\partial^2 f}{\partial (u^\alpha)^2}\bigg|_{\text{FD}} = \frac{f(u^\alpha + 1/N_\alpha) + f(u^\alpha - 1/N_\alpha) - 2f(u^\alpha)}{(1/N_\alpha)^2}$$

对平面波：

$$= -2N_\alpha^2 \left(1 - \cos\frac{2\pi m_\alpha}{N_\alpha}\right) f$$

**交叉项** (α ≠ β)：

$$\frac{\partial^2 f}{\partial u^\alpha \partial u^\beta}\bigg|_{\text{FD}} = \frac{1}{4 N_\alpha N_\beta}\Big[f\Big(\frac{+1}{N_\alpha}, \frac{+1}{N_\beta}\Big) - f\Big(\frac{+1}{N_\alpha}, \frac{-1}{N_\beta}\Big) - f\Big(\frac{-1}{N_\alpha}, \frac{+1}{N_\beta}\Big) + f\Big(\frac{-1}{N_\alpha}, \frac{-1}{N_\beta}\Big)\Big]$$

对平面波：

$$= -N_\alpha N_\beta \sin\frac{2\pi m_\alpha}{N_\alpha} \sin\frac{2\pi m_\beta}{N_\beta} \, f$$

### 2.5 FD kernel

定义 FD 因子矩阵：

$$\text{FD}_{\alpha\alpha} = 2N_\alpha^2\left(1 - \cos\frac{2\pi m_\alpha}{N_\alpha}\right)$$

$$\text{FD}_{\alpha\beta} = N_\alpha N_\beta \sin\frac{2\pi m_\alpha}{N_\alpha} \sin\frac{2\pi m_\beta}{N_\beta} \quad (\alpha \neq \beta)$$

FD Laplacian 的等效 gg 为：

$$\text{gg}_{\text{FD}} = \frac{1}{(2\pi)^2} \sum_{\alpha,\beta} \text{GGT}[\alpha][\beta] \cdot \text{FD}_{\alpha\beta}$$

最终：

$$\nabla^2_{\text{FD}} f = -\text{gg}_{\text{FD}} \cdot \text{tpiba2} \cdot f$$

**与谱 kernel 的关系**：当 $m_\alpha / N_\alpha \to 0$ 时，

$$2(1-\cos\theta) \to \theta^2, \quad \sin\theta \to \theta$$

因此 $\text{FD}_{\alpha\alpha} \to (2\pi m_\alpha)^2$，$\text{FD}_{\alpha\beta} \to (2\pi)^2 m_\alpha m_\beta$，$\text{gg}_{\text{FD}} \to \text{gg}_{\text{spectral}}$。

**有界性**：gg_FD 的最大值有限（由 FFT 网格大小决定），而 gg_spectral 随 |m| 无界增长。

### 2.6 正交晶胞简化

对正交晶胞，GGT 为对角矩阵（GGT[αβ] = 0 for α ≠ β），交叉项消失：

$$\text{gg}_{\text{FD}} = \frac{1}{(2\pi)^2} \sum_\alpha \text{GGT}[\alpha][\alpha] \cdot 2N_\alpha^2\left(1 - \cos\frac{2\pi m_\alpha}{N_\alpha}\right)$$

每个方向独立贡献。

---

## 3. 实现

### 3.1 vlapl 势

在 `v_xc_meta` 函数中，对每个自旋通道 is：

1. **构造 vlapl·sgn**（实空间）：

$$\text{vlapl\_sgn}[ir] = \text{vlapl}[ir \cdot \text{nspin} + \text{is}] \cdot \text{sgn}[ir \cdot \text{nspin} + \text{is}]$$

2. **FFT 到 G 空间**：

$$\widetilde{\text{vlapl\_sgn}} = \text{FFT}(\text{vlapl\_sgn})$$

3. **应用 FD kernel**：对每个 G 向量 ig，Miller 指数 m = gdirect[ig]，

$$\widetilde{\nabla^2_{\text{FD}}(\text{vlapl\_sgn})}[ig] = -\text{gg}_{\text{FD}}[ig] \cdot \text{tpiba2} \cdot \widetilde{\text{vlapl\_sgn}}[ig]$$

其中 gg_FD[ig] 由 2.5 节公式计算。

4. **IFFT 回实空间**：

$$\nabla^2_{\text{FD}}(\text{vlapl\_sgn}) = \text{IFFT}(\widetilde{\nabla^2_{\text{FD}}(\text{vlapl\_sgn})})$$

5. **加入 V_xc**：

$$V_{xc}(\text{is}, ir) \mathrel{+}= e^2 \cdot \nabla^2_{\text{FD}}(\text{vlapl\_sgn})[ir]$$

对于 HSE 等 hybrid 泛函 (func_type == 5)，乘以 (1 − hybrid_alpha)。

6. **加入 vtxc**：

$$\text{vtxc} \mathrel{+}= e^2 \sum_{ir} \nabla^2_{\text{FD}}(\text{vlapl\_sgn})[ir] \cdot \rho[\text{is}][ir] \cdot \text{sgn}[ir \cdot \text{nspin} + \text{is}]$$

### 3.2 vlapl 应力

vlapl 对应力的贡献在 `gradcorr` 函数中计算，与谱 Laplacian 方案相同：

$$\sigma_{\alpha\beta}^{\nabla^2\rho} = \frac{2}{\Omega} \sum_{\mathbf{r}} \text{vlapl}(\mathbf{r}) \cdot H_{\alpha\beta}(\mathbf{r}) \cdot e^2$$

其中 $H_{\alpha\beta} = \partial^2\rho / \partial r_\alpha \partial r_\beta$ 是密度的 Hessian 矩阵（谱方法计算）。LibXC 返回的 vlapl = ∂(ρε)/∂(∇²ρ) 已包含 ρ 因子。

### 3.3 功能名注册

`dft_functional scanl` 映射为 LibXC 泛函 `XC_MGGA_X_SCANL` (ID=700) + `XC_MGGA_C_SCANL` (ID=702)。

---

## 4. FD kernel 数值验证

### 4.1 验证设置

FCC 原胞 (Si₂)，a = 10.2 Bohr。

晶格矢量（lat0 单位）：

$$\mathbf{a}_1 = (0, \tfrac{1}{2}, \tfrac{1}{2}), \quad \mathbf{a}_2 = (\tfrac{1}{2}, 0, \tfrac{1}{2}), \quad \mathbf{a}_3 = (\tfrac{1}{2}, \tfrac{1}{2}, 0)$$

度规矩阵：

$$\text{GGT} = \begin{pmatrix} 3 & -1 & -1 \\ -1 & 3 & -1 \\ -1 & -1 & 3 \end{pmatrix}$$

### 4.2 gg_FD vs gg_spectral（网格 36×36×36, ecutwfc=60 Ry）

| Miller (m₁,m₂,m₃) | gg_spectral | gg_FD | gg_FD / gg_spectral |
|---------------------|-------------|-------|---------------------|
| (0,0,0) | 0 | 0 | — |
| (1,0,0) | 3.000 | 2.992 | 0.9975 |
| (1,1,0) | 4.000 | 4.005 | 1.0013 |
| (2,0,0) | 12.00 | 11.88 | 0.9900 |
| (3,0,0) | 27.00 | 26.46 | 0.9800 |
| (5,0,0) | 75.00 | 70.36 | 0.9381 |
| (8,0,0) | 192.0 | 163.4 | 0.8511 |
| (10,0,0) | 300.0 | 231.2 | 0.7706 |
| (15,0,0) | 675.0 | 387.4 | 0.5740 |
| (18,0,0) | 972.0 | 393.9 | 0.4053 |

**解读**：

- |m|/N ≤ 0.1 时：FD 与谱偏差 < 1%，物理上等价
- |m|/N ≈ 0.2 时：偏差 ~6%，开始有差异
- |m|/N ≈ 0.5（Nyquist）：FD = 0.405 × 谱，大幅衰减

### 4.3 gg_FD vs gg_spectral（网格 45×45×45, ecutwfc=80 Ry）

| Miller (m₁,m₂,m₃) | gg_spectral | gg_FD | gg_FD / gg_spectral |
|---------------------|-------------|-------|---------------------|
| (1,0,0) | 3.000 | 2.994 | 0.9980 |
| (3,0,0) | 27.00 | 26.66 | 0.9875 |
| (5,0,0) | 75.00 | 72.55 | 0.9673 |
| (8,0,0) | 192.0 | 173.7 | 0.9048 |
| (9,0,0) | 243.0 | 211.7 | 0.8712 |

ecutwfc=80 对应的最大 |m| ≈ 8.4，此时 FD/谱 ≈ 0.89。即 PW 截断处 FD kernel 衰减约 11%。

### 4.4 有界性验证

Nyquist 频率 (m_α = N_α/2) 下 FD 因子的最大值：

$$\text{FD}_{\alpha\alpha}^{\max} = 2N_\alpha^2 \cdot 2 = 4N_\alpha^2$$

$$\text{gg}_{\text{FD}}^{\max} = \frac{1}{(2\pi)^2} \sum_\alpha \text{GGT}[\alpha][\alpha] \cdot 4N_\alpha^2 = \frac{4}{(2\pi)^2} \cdot 3 \cdot 3 N^2 = \frac{36 N^2}{4\pi^2}$$

对 N=36：gg_FD_max = 1181.8，对应 |G|_FD_max ≈ 21.2 Bohr⁻¹（有界）。
对 N=45：gg_FD_max = 1846.5。

而谱 Laplacian 在同一点：gg_spectral = 9 · (N/2)² → 随 N² 无界增长。

---

## 5. 有限差分应力验证

### 5.1 验证方法

对 Si₂ FCC 原胞（Γ 点）在 a₀ ± δ 下计算总能量 E，由差分得到应力：

$$\sigma_{\text{FD}} = -\frac{E(a_0+\delta) - E(a_0-\delta)}{V(a_0+\delta) - V(a_0-\delta)}$$

FCC 原胞体积 $V = |\det(\text{latvec})| \cdot a^3 / 4 = a^3 / 4$。

与 ABACUS 解析应力 σ_AB 比较。单位转换：

$$\sigma(\text{kbar}) = \sigma(\text{eV/Bohr}^3) \times \frac{147105}{13.6057}$$

其中 147105 kbar/(Ry/Bohr³)，13.6057 eV/Ry。

### 5.2 测试参数

| 参数 | 值 |
|------|-----|
| 体系 | Si₂ FCC 原胞 |
| a₀ | 10.2 Bohr |
| δ | 0.001 Bohr |
| 赝势 | Si_ONCV_PBE-1.0 (NCPP) |
| smearing | gauss, σ = 0.002 |
| mixing | Pulay, beta = 0.1–0.2, gg0 = 1.0 |
| mixing_tau | 1 |
| scf_thr | 1e-8 |
| pw_seed | 1 |

### 5.3 SCAN 基线（无 ∇²ρ 依赖）

SCAN 的 vlapl = 0，FD Laplacian 代码路径不激活，用于验证 FD 方法和应力代码。

| ecutwfc (Ry) | FFT grid | σ_FD (kbar) | σ_AB (kbar) | FD/AB | 误差 |
|-------------|----------|-------------|-------------|-------|------|
| 60 | 36³ | 434.29 | 434.30 | 1.0000 | 0.00% |
| 80 | 45³ | 434.47 | 434.48 | 1.0000 | 0.00% |
| 100 | 45³ | 434.54 | 434.55 | 1.0000 | 0.00% |

**结论**：SCAN 在所有 ecutwfc 下 FD/AB = 1.0000，验证了代码和方法的正确性。

### 5.4 SCANL（FD Laplacian 方案）

| ecutwfc (Ry) | FFT grid | σ_FD (kbar) | σ_AB (kbar) | FD/AB | 误差 |
|-------------|----------|-------------|-------------|-------|------|
| 60 | 36³ | 431.89 | 435.60 | 0.9915 | -0.85% |
| 80 | 45³ | 432.37 | 434.58 | 0.9949 | -0.51% |

ecutwfc=100 时 SCANL 的 SCF 收敛困难（能量振荡），FD 结果不可靠，未列入。

### 5.5 完整能量数据

ecutwfc = 80 Ry：

| | E(a₀+δ) (eV) | E(a₀−δ) (eV) | ΔE (eV) |
|------|-------------|-------------|---------|
| SCAN | −197.22465 | −197.21838 | −0.006271 |
| SCANL | −198.70596 | −198.69972 | −0.006241 |

ecutwfc = 60 Ry：

| | E(a₀+δ) (eV) | E(a₀−δ) (eV) | ΔE (eV) |
|------|-------------|-------------|---------|
| SCAN | −197.22344 | −197.21717 | −0.006269 |
| SCANL | −198.68669 | −198.68046 | −0.006234 |

SCAN 的 ΔE 跨 ecutwfc 稳定（−0.00627 eV），SCANL 的 ΔE 在 ecut=60/80 时接近（−0.00623/−0.00624 eV），与 SCAN 的差异约 0.5%，与 FD 应力误差一致。

### 5.6 与其他方案对比

| 方案 | V 含 vlapl | vtxc 含 vlapl | stress 含 vlapl | SCF | FD/AB (ecut=80) | 误差 |
|------|-----------|--------------|----------------|-----|-------|------|
| 严格谱 Laplacian | ✓ (谱) | ✓ (谱) | ✓ | **发散** | — | — |
| sigma-as-lapl | ✗ (假输入) | ✗ (假输入) | ✗ | 收敛 | 1.046 | +4.6% |
| V/vtxc 均省略 vlapl | ✗ | ✗ | ✓ | 收敛 | 0.906 | −9.5% |
| vtxc 含 vlapl | ✗ | ✓ | ✓ | 收敛 | 0.922 | −7.8% |
| 高斯衰减 | ✓ (衰减) | ✓ (衰减) | ✓ | 收敛 | ~0.9 | ~10% |
| 绝热 ramp | ✓ (渐增) | ✓ | ✓ | **发散** | — | — |
| **FD Laplacian** | **✓ (FD)** | **✓ (FD)** | **✓** | **收敛** | **0.995** | **−0.5%** |

---

## 6. 误差分析

### 6.1 残余 0.5% 应力误差来源

1. **FD kernel 与谱 kernel 在高 G 的偏差**：在 PW 截断处（|m|/N ≈ 0.37 for ecut=80），FD kernel 衰减到谱值的 ~89%。vlapl 势在高 G 分量上偏弱，密度与严格 SCANL 密度有微小差异。

2. **FD 应力本身的高阶误差**：δ = 0.001 Bohr 的中心差分为二阶精度，误差 O(δ²/a⁴) ≈ 10⁻⁸，对 0.5% 误差贡献可忽略。

3. **vlapl stress 项中 Hessian 仍为谱 Hessian**：gradcorr 中的 vlapl 应力项使用谱方法计算 Hessian，而势用 FD Laplacian。此不一致性贡献未知但量级较小。

### 6.2 ecutwfc=100 的 SCF 困难

ecutwfc=100 时 SCANL 的 ΔE = −0.00568 eV，显著偏离 SCAN 的 −0.00628 eV（偏差 ~10%），而 ecut=60/80 时偏差仅 ~0.5%。原因可能是：

- 45³ FFT 网格对 ecut=100 的 35749 个平面波不够大，某些高 G 分量的 FD kernel 在网格边界处行为异常
- SCF 收敛到与 ecut=60/80 不同的密度

**建议**：SCANL 使用 ecutwfc ≤ 80 Ry。若需更高截断，应确保 FFT 网格足够大。

### 6.3 与严格实现的关系

严格 SCANL 实现要求谱 Laplacian ∇² = −|G|²，但在任何实际 ecutwfc 下导致 SCF 发散。FD Laplacian 是在 FFT 网格上对离散函数求 Laplacian 的自洽方法——它是网格函数的"正确" Laplacian，但与连续 Laplacian 有差异。此差异随网格加密而减小（FD → spectral 当 N → ∞），但在实际网格上不可消除。

---

## 7. 使用建议

1. **dft_functional scanl**：已注册，映射到 XC_MGGA_X_SCANL + XC_MGGA_C_SCANL
2. **ecutwfc ≤ 80 Ry**：ecutwfc 过高时 SCF 收敛困难
3. **mixing_tau = 1**：meta-GGA 必需
4. **Pulay mixing, beta = 0.1–0.2**：对 SCANL 收敛效果较好
5. **优先考虑 r²SCAN 或 SCAN**：无 ∇²ρ 依赖，数值更稳定
