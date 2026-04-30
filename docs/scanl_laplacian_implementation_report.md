# SCANL Laplacian 依赖项实现与验证报告

## 1. 背景

SCANL（SCAN-L）是 SCAN 的 deorbitalized 版本，将 SCAN 中对动能密度 τ 的依赖替换为对密度 Laplacian ∇²ρ 的依赖 [Mejia-Rodriguez & Trickey, Phys. Rev. B **98**, 115161 (2018)]。这导致其 Kohn-Sham 势和应力公式中出现了 ∇²(∂ε/∂(∇²ρ)) 项，在平面波基组下引发严重的数值稳定性问题。

本文档记录了 ABACUS 中 SCANL 泛函 Laplacian 相关项的实现方案、有限差分（FD）应力验证结果，以及最终建议。

---

## 2. 理论公式

### 2.1 meta-GGA 能量泛函

meta-GGA 交换关联能量密度 ε_xc 依赖于四个变量：

$$E_{xc} = \int \rho(\mathbf{r})\,\varepsilon_{xc}(\rho, |\nabla\rho|^2, \nabla^2\rho, \tau)\, d\mathbf{r}$$

对于 SCANL，∂ε/∂(∇²ρ) ≠ 0（与 SCAN 不同，SCAN 的此项为零）。

### 2.2 Kohn-Sham 势

由变分原理，XC 势为：

$$V_{xc} = e^2 \frac{\delta E_{xc}}{\delta \rho} = e^2 \left[ \frac{\partial(\rho\varepsilon)}{\partial\rho} - \nabla \cdot \frac{\partial(\rho\varepsilon)}{\partial(\nabla\rho)} + \nabla^2 \frac{\partial(\rho\varepsilon)}{\partial(\nabla^2\rho)} + \frac{\partial(\rho\varepsilon)}{\partial\tau} \hat{\tau} \right]$$

其中 **vlapl 势项**为：

$$V_{xc}^{\nabla^2\rho} = e^2\,\nabla^2\!\left(\frac{\partial(\rho\varepsilon)}{\partial(\nabla^2\rho)}\right)$$

在倒易空间中，∇² → -|G|²，因此：

$$\tilde{V}_{xc}^{\nabla^2\rho}(\mathbf{G}) = -e^2\,|\mathbf{G}|^2\,\widetilde{\text{vlapl}}(\mathbf{G})$$

**关键问题**：|G|² 是高通滤波器，放大高频模式。vlapl(G) 中的数值噪声被 |G|² 放大后反馈到密度中，形成正反馈循环导致 SCF 发散。ecutwfc 越大，|G|_max 越大，放大效应越强，SCF 越容易发散。这并非 mixing 策略的问题，而是 Laplace 逆变换的数学病态性。

### 2.3 vtxc 与应力公式

vtxc（势对密度的积分）在应力计算中至关重要：

$$\text{vtxc} = \int V_{xc}(\mathbf{r})\,\rho(\mathbf{r})\, d\mathbf{r}$$

应力中 XC 对角贡献为：

$$\sigma_{ii}^{xc} = -\frac{E_{txc} - V_{txc}}{\Omega}$$

其中 $E_{txc} = \int e^2\,\rho\,\varepsilon_{xc}\, d\mathbf{r}$。

**vtxc 必须与实际使用的 V_xc 一致**。若 V 中不含 vlapl 而 vtxc 中含 vlapl，则 $-(E_{txc}-V_{txc})/\Omega$ 不正确。

### 2.4 vlapl 应力项

由能量泛函对应变的变分推导，vlapl 对应力的贡献为：

$$\sigma_{\alpha\beta}^{\nabla^2\rho} = \frac{2}{\Omega} \sum_{\mathbf{r}} \frac{\partial(\rho\varepsilon)}{\partial(\nabla^2\rho)}\, H_{\alpha\beta}(\mathbf{r})\, e^2$$

其中 $H_{\alpha\beta} = \partial^2\rho / \partial r_\alpha \partial r_\beta$ 是密度的 Hessian 矩阵。

**注意**：LibXC 返回的 vlapl = ∂(ρε)/∂(∇²ρ) 已包含 ρ 因子，无需再乘 ρ。

---

## 3. 实现方案与测试结果

### 3.1 验证方法

采用有限差分（FD）应力验证：对 Si2 FCC 原胞（Γ 点）在 a₀±δ 下计算总能量，由中心差分公式得到应力：

$$\sigma_{FD} = -\frac{E(a_0+\delta) - E(a_0-\delta)}{V(a_0+\delta) - V(a_0-\delta)}$$

与解析应力 σ_AB 比较，FD/AB 越接近 1.0 越好。

**单位转换**：ABACUS 总能量单位为 eV，应力原子单位为 Ry/Bohr³。FD 应力转换公式：

$$\sigma_{FD}(\text{kbar}) = \sigma_{FD}(\text{eV/Bohr}^3) \times \frac{2\,\text{Ry}}{27.211386\,\text{eV}} \times 147105\,\frac{\text{kbar}}{\text{Ry/Bohr}^3}$$

### 3.2 测试参数

| 参数 | 值 |
|------|-----|
| 体系 | Si₂ FCC 原胞 |
| a₀ | 10.2 Bohr |
| δ | 0.001 Bohr |
| ecutwfc | 60 / 80 / 100 / 400 Ry |
| smearing_sigma | 0.02 |
| ks_solver | dav_subspace |
| mixing_beta | 0.1 |
| mixing_tau | 1 |
| scf_thr | 1e-7 |
| pw_seed | 1 |

### 3.3 SCAN 基线验证

SCAN 不依赖 ∇²ρ（LibXC 返回 vlapl=0），用于验证 FD 方法和代码的正确性：

| ecutwfc (Ry) | ecutwfc (eV) | FD (kbar) | AB (kbar) | FD/AB | 误差 |
|-------------|-------------|-----------|-----------|-------|------|
| 60 | 816 | 265.69 | 265.66 | 1.0001 | 0.01% |
| 80 | 1088 | 266.11 | 266.11 | 1.0000 | 0.00% |
| 100 | 1361 | 266.11 | 266.12 | 0.9999 | -0.01% |
| 400 | 5442 | 266.11 | 266.09 | 1.0001 | 0.01% |

**结论**：SCAN 在所有 ecutwfc 下 FD/AB ≈ 1.000，验证了 FD 方法和应力代码的正确性。

### 3.4 各方案详情与结果

---

#### 方案 1：完整 e2·∇²(vlapl) 加入势和 vtxc（严格实现）

**实现**：
- V_xc：含 e2·∇²(vlapl·sgn)
- vtxc：含 e2·∫∇²(vlapl·sgn)·ρ·sgn
- stress：含 2·vlapl·Hess·e2

**结果**：**SCF 发散**。无论 mixing_beta（0.01–0.7）、Kerker mixing、绝热演化（50 步 ramp）均无法收敛。

**原因**：V(G) = -e2·|G|²·vlapl(G)，|G|² 高通滤波放大高频噪声，正反馈增益超过任何混合策略的阻尼能力。ecutwfc 越大，|G|_max 越大，发散越严重。

---

#### 方案 2：sigma 伪装 lapl（原始代码方案）

**实现**：
- `v_xc_libxc`：传 sigma（= |∇ρ|²）给 LibXC 的 Laplacian 参数位置
- `tau_xc`/`tau_xc_spin`：传 grho（= |∇ρ|²）给 LibXC 的 Laplacian 参数位置
- gradcorr stress：无 vlapl 应力项

**问题**：LibXC 在错误的 Laplacian 输入下计算 exc 和 vlapl，导致 etxc 和应力均不正确。

**FD 测试结果**：

| ecutwfc (Ry) | FD (kbar) | AB (kbar) | FD/AB | 误差 |
|-------------|-----------|-----------|-------|------|
| 60 | 249.90 | 240.91 | 1.037 | 3.73% |
| 80 | 252.32 | 241.11 | 1.047 | 4.65% |
| 100 | 252.39 | 241.12 | 1.047 | 4.68% |
| 400 | 252.39 | 241.24 | 1.046 | 4.62% |

FD/AB ≈ 1.04–1.05，存在系统性 4–5% 误差。SCANL 误差随 ecutwfc 趋于稳定（~4.6%），说明这是传假 Laplacian 导致的系统偏差，而非 ecutwfc 收敛问题。

---

#### 方案 3：V 和 vtxc 均不含 vlapl，stress 保留 vlapl 项（当前实现）

**实现**：
- `v_xc_libxc`：传真实 ∇²ρ 给 LibXC，获取正确的 exc、vrho、vsigma、vtau
- V_xc：不含 e2·∇²(vlapl)
- vtxc：不含 vlapl 贡献（与 V 保持一致）
- gradcorr stress：含 2·vlapl·Hess·e2

**物理逻辑**：vlapl 势无法加入 V（会发散），因此 vtxc 也不含 vlapl 以保持一致。etxc 中包含完整的 exc（含 Laplacian 贡献）。gradcorr 中的 vlapl 应力项从完整能量泛函严格推导。

**FD 测试结果**：

PW 基组：

| ecutwfc (Ry) | FD (kbar) | AB (kbar) | FD/AB | 误差 |
|-------------|-----------|-----------|-------|------|
| 60 | 260.22 | 283.88 | 0.917 | -8.3% |
| 80 | 254.26 | 280.80 | 0.906 | -9.5% |
| 100 | 270.06 | 278.70 | 0.969 | -3.1% |

LCAO 基组（ecutwfc=100 Ry）：

| 基组 | FD (kbar) | AB (kbar) | FD/AB | 误差 |
|------|-----------|-----------|-------|------|
| LCAO | 248.30 | 276.09 | 0.899 | -10.1% |

**误差分析**：AB 应力系统性大于 FD 应力约 8–10%。根本原因是 V 不含 vlapl → 密度在近似势下自洽（接近 SCAN 密度）→ etxc 用 SCANL 的 exc 评估（含 Laplacian 贡献）→ etxc-vtxc 不匹配 → -(etxc-vtxc)/Ω 偏大 → 对角应力偏大。这是省略 vlapl 势的根本性限制，无法通过调整 vlapl stress 项修正。

---

#### 方案 4：势中无 vlapl，vtxc 中加 e2·vlapl·∇²ρ

**实现**：
- V_xc：不含 vlapl
- vtxc：含 e2·∫vlapl·∇²ρ·sgn
- stress：含 2·vlapl·Hess·e2

**问题**：V 无 vlapl 但 vtxc 有 → vtxc 与实际 V 不一致 → $-(E_{txc}-V_{txc})/\Omega$ 不正确 → 对角应力偏小。此前测试 FD/AB = 0.922（7.78% 误差）。

---

#### 方案 5：高斯衰减

**实现**：
- V_xc：e2·∇²(vlapl·sgn·damp(G))，damp = exp(-α·G²/G_max²)
- vtxc：从衰减势计算或从未衰减 vlapl·∇²ρ 计算
- stress：含 2·vlapl·Hess·e2（无衰减）

**此前测试结果**（ecutwfc=200 Ry，仅参考趋势）：

| α | vtxc 来源 | FD/AB |
|---|----------|-------|
| 1.0 | 衰减势 | 0.789 |
| 2.0 | 衰减势 | 0.960 |
| 4.0 | 衰减势 | 0.918 |

**问题**：衰减改变了势的物理形状，密度在错误势下自洽。不同 α 给出不同结果，无法确定哪个是"正确的"。本质上是调参拟合，没有物理依据。

---

#### 方案 6：绝热演化（ramp）

**实现**：
- V_xc：e2·min(1, iter/ramp_steps)·∇²(vlapl·sgn)
- vtxc：e2·∫∇²(vlapl·sgn)·ρ·sgn（全强度）

**结果**：**SCF 发散**。即使 ramp_steps=50 + mixing_beta=0.03，当 ramp_factor → 1 时仍开始振荡。

---

#### 方案 7：有限差分 Laplacian（FD kernel in G-space）

**核心思想**：谱 Laplacian ∇²f(G) = -|G|²f(G) 中的 |G|² 是无界的，导致高频放大和 SCF 发散。有限差分（FD）Laplacian 在低 G 时与谱 Laplacian 一致，但在高 G 时自然衰减，最大放大因子有界。

**推导**：

在分数坐标 u^α = n_α/N_α（α = 1,2,3）下，Laplacian 表示为：

$$\nabla^2 = \sum_{\alpha,\beta} g^{\alpha\beta} \frac{\partial^2}{\partial u^\alpha \partial u^\beta}$$

其中 g^{αβ} = GGT[αβ]/lat0² 是逆变度规张量，GGT 是倒格子度规矩阵（ABACUS 中的 `rhopw->GGT`）。

对平面波 f(r) = exp(2πi Σ m_γ u^γ)，谱 Laplacian 给出：

$$\frac{\partial^2 f}{\partial u^\alpha \partial u^\beta} = -(2\pi)^2 m_\alpha m_\beta \, f$$

FD 近似给出：

- **对角项** (α = β)：中心差分

$$\left(\frac{\partial^2 f}{\partial (u^\alpha)^2}\right)_{\text{FD}} = -2N_\alpha^2 \left(1 - \cos\frac{2\pi m_\alpha}{N_\alpha}\right) \, f$$

- **交叉项** (α ≠ β)：四点差分

$$\left(\frac{\partial^2 f}{\partial u^\alpha \partial u^\beta}\right)_{\text{FD}} = -N_\alpha N_\beta \sin\frac{2\pi m_\alpha}{N_\alpha} \sin\frac{2\pi m_\beta}{N_\beta} \, f$$

对低频模式 (m_α ≪ N_α)，FD 退化为谱形式（1-cos(x) ≈ x²/2, sin(x) ≈ x）。对 Nyquist 频率 (m_α = N_α/2)，放大因子为有界值 N_α² · 4/(2π)² ≈ 0.405 N_α²，而非谱 Laplacian 的 (N_α/2)²。

组合所有项，FD Laplacian 的等效 G² 核为：

$$\text{gg}_{\text{FD}} = \frac{1}{(2\pi)^2} \sum_{\alpha,\beta} \text{GGT}[\alpha][\beta] \cdot \text{FD}_{\alpha\beta}$$

其中：

$$\text{FD}_{\alpha\alpha} = 2N_\alpha^2(1-\cos\frac{2\pi m_\alpha}{N_\alpha}), \quad \text{FD}_{\alpha\beta}^{(\alpha\neq\beta)} = N_\alpha N_\beta \sin\frac{2\pi m_\alpha}{N_\alpha} \sin\frac{2\pi m_\beta}{N_\beta}$$

最终，∇²_FD(vlapl·sgn) 在 G 空间中为：

$$\widetilde{\nabla^2_{\text{FD}}(\text{vlapl}\cdot\text{sgn})}(\mathbf{G}) = -\text{gg}_{\text{FD}}(\mathbf{G}) \cdot \text{tpiba2} \cdot \widetilde{\text{vlapl}\cdot\text{sgn}}(\mathbf{G})$$

**实现**：
1. 计算 vlapl_sgn = vlapl·sgn（实空间）
2. FFT: vlapl_sgn → vlapl_g
3. 应用 FD kernel: vlapl_g_lapl[ig] = -vlapl_g[ig] · gg_FD[ig] · tpiba2
4. IFFT: vlapl_g_lapl → vlapl_lapl（实空间）
5. 加入 V_xc 和 vtxc

**交叉项的物理意义**：对于非正交晶胞（如 FCC 原胞），度规张量 GGT 的非对角元不为零。此时，沿一个晶格方向的位移会影响其他方向的梯度分量。FD kernel 的交叉项 sin·sin 正确反映了这一点。对于正交晶胞（GGT 对角），交叉项消失，FD kernel 退化为三个方向的独立贡献。

**数值验证**（FCC 原胞, a=10.2 Bohr, nx=ny=nz=36）：

| Miller 指数 | gg_spectral | gg_FD | gg_FD/gg_spectral |
|-------------|-------------|-------|-------------------|
| (0,0,0) | 0 | 0 | 1.000 |
| (1,0,0) | 3.0 | 2.992 | 0.997 |
| (1,1,0) | 4.0 | 4.005 | 1.001 |
| (5,0,0) | 75.0 | 70.36 | 0.938 |
| (10,0,0) | 300.0 | 231.2 | 0.771 |
| (18,0,0) | 972.0 | 393.9 | 0.405 |

FD kernel 在低 G 时与谱 kernel 完全一致（<1% 偏差），在高 G 时衰减到谱值的 ~40%，有界且不会导致发散。

**SCF 收敛性**：使用 Pulay mixing (beta=0.2)，SCANL 在 ecutwfc=80 Ry 下 SCF 正常收敛。

**FD 应力验证**（Si₂ FCC 原胞, a₀=10.2 Bohr, δ=0.001, Pulay mixing beta=0.1–0.2）：

PW 基组：

| 泛函 | ecutwfc (Ry) | σ_FD (kbar) | σ_AB (kbar) | FD/AB | 误差 |
|------|-------------|-------------|-------------|-------|------|
| SCAN | 60 | 434.29 | 434.30 | 1.0000 | 0.00% |
| SCAN | 80 | 434.47 | 434.48 | 1.0000 | 0.00% |
| SCAN | 100 | 434.54 | 434.55 | 1.0000 | 0.00% |
| SCANL | 60 | 431.89 | 435.60 | 0.9915 | -0.85% |
| SCANL | 80 | 432.37 | 434.58 | 0.9949 | -0.51% |

**注意**：ecutwfc=100 时 SCANL 的 SCF 收敛较困难（能量振荡较大），FD 应力结果不稳定，未列入上表。建议 SCANL 使用 ecutwfc ≤ 80 Ry。

**与方案 3 的对比**：

| 方案 | V 含 vlapl | vtxc 含 vlapl | SCF | FD/AB | 误差 |
|------|-----------|--------------|-----|-------|------|
| 方案 3 (无 vlapl 势) | ✗ | ✗ | 收敛 | 0.906 | -9.5% |
| 方案 7 (FD Laplacian) | ✓ (FD) | ✓ (FD) | 收敛 | 0.995 | -0.51% |

**FD Laplacian 方案将应力误差从 9.5% 降至 0.51%，改善了约 20 倍。**

残余 0.5% 误差的可能来源：
1. FD kernel 在高 G 区与谱 kernel 的偏差（已衰减但仍有差异）
2. 有限差分应力本身的高阶误差（δ=0.001 的中心差分为二阶精度）
3. vlapl stress 项中使用的 Hessian 仍为谱 Hessian

---

| 方案 | V 含 vlapl | vtxc 含 vlapl | stress 含 vlapl | SCF | FD/AB | 误差 |
|------|-----------|--------------|----------------|-----|-------|------|
| SCAN 基线 | N/A | N/A | N/A | 收敛 | 1.000 | 0.00% |
| 方案 1 (严格) | ✓ | ✓ | ✓ | **发散** | — | — |
| 方案 2 (sigma-as-lapl) | ✗ (假) | ✗ (假) | ✗ | 收敛 | 1.046 | 4.6% |
| 方案 3 (严格应力) | ✗ | ✗ | ✓ | 收敛 | 0.906 | -9.5% |
| 方案 4 (vtxc含vlapl) | ✗ | ✓ | ✓ | 收敛 | 0.922 | -7.8% |
| 方案 5 (衰减) | ✓(衰减) | ✓(衰减) | ✓ | 收敛 | ~0.9 | ~10% |
| 方案 6 (ramp) | ✓(渐增) | ✓ | ✓ | **发散** | — | — |
| **方案 7 (FD Laplacian)** | **✓ (FD)** | **✓ (FD)** | **✓** | **收敛** | **0.995** | **-0.5%** |

**结论**：方案 7（FD Laplacian）是唯一同时满足 SCF 收敛和应力精度的方案。它通过用有界的 FD kernel 替代无界的 |G|² 谱 Laplacian，既解决了 SCF 发散问题，又将 vlapl 势纳入 V_xc 使得密度-势自洽，从而大幅改善应力精度（从 9.5% 误差降至 0.5%）。

---

## 4. 模守恒赝势与 SCANL 的根本矛盾

### 4.1 ∇²(vlapl) 发散的 ecutwfc 依赖

vlapl 势在倒易空间为 V(G) = -e2·|G|²·vlapl(G)。|G|_max 由 ecutwfc 决定：

$$|G|_{max} = \sqrt{2 \cdot E_{cut}}$$

|G|_max² 随 ecutwfc 线性增长，意味着 vlapl 势的高频放大随 ecutwfc 增强而加剧：

- **低 ecutwfc**：|G|_max 小，放大范围有限，SCF 可能收敛，但基组不完备导致其他误差
- **高 ecutwfc**：|G|_max 大，高频被剧烈放大，SCF 几乎必然发散

**严格实现方案 1 在任何实际使用的 ecutwfc 下都无法收敛**。

### 4.2 模守恒赝势需要高 ecutwfc

模守恒赝势（NCPP）比 PAW 方法需要更高的 ecutwfc 才能收敛：

| 赝势类型 | 典型 ecutwfc (Ry) | 典型 ecutwfc (eV) |
|---------|------------------|------------------|
| PAW/USPP | 30–50 | 400–700 |
| NCPP | 60–100+ | 800–1400+ |

ABACUS 目前使用模守恒赝势。NCPP 要求的高 ecutwfc 与 vlapl 势的 |G|² 发散直接矛盾：

- **低 ecutwfc**（< 50 Ry）：SCF 可能收敛，但 NCPP 的赝势不完备，结果不收敛
- **高 ecutwfc**（> 60 Ry）：NCPP 收敛，但 vlapl 势发散

**这是一个不可调和的矛盾**。在 NCPP + 平面波框架下，不存在一个 ecutwfc 值能同时满足赝势收敛和 vlapl 势稳定。

### 4.3 VASP 的做法

VASP wiki [METAGGA 页面] 明确指出：

> The results obtained with the meta-GGA functionals that depend on the Laplacian of the density ∇²n (e.g., SCAN-L) may not be reliable for large values of the energy cutoff ENCUT due to numerical instability. According to some tests, it is not recommended to use values of ENCUT above 800 eV.

VASP 使用 PAW 方法，ecutwfc 要求较低（~400–700 eV），可以在 800 eV 以下的窗口中获得结果。但即便如此，VASP 也承认高截断能下结果不可靠。对于 NCPP，800 eV（约 60 Ry）可能刚达到赝势收敛的门槛，而更高截断又会触发不稳定——几乎不存在可用的 ecutwfc 窗口。

---

## 5. 结论与建议

### 5.1 推荐方案：FD Laplacian（方案 7）

**FD Laplacian 方案是当前最佳方案**，它：

1. **SCF 收敛**：FD kernel 有界，不会导致 |G|² 发散
2. **密度-势自洽**：V 包含 vlapl 贡献（通过 FD Laplacian），密度在近似自洽势下产生
3. **应力精度高**：FD/AB = 0.995，误差仅 0.5%（远优于方案 3 的 9.5%）
4. **适用于非正交晶胞**：通过度规张量 GGT 正确处理交叉项
5. **无需调参**：不依赖衰减因子或 ramp 步数等人为参数

**残余 0.5% 误差**是 FD kernel 与谱 kernel 在高 G 区差异的代价，对于实际应用完全可以接受。

### 5.2 使用注意事项

1. **FD Laplacian 是对谱 Laplacian 的近似**：在 FFT 网格上，它是"正确"的 Laplacian（对网格函数的有限差分），但与连续 Laplacian 有差异
2. **需要 mixing_tau = 1**：meta-GGA 泛函需要混合动能密度
3. **建议使用 Pulay mixing**：mixing_beta=0.2，对 SCANL 效果较好
4. **仍建议优先考虑 r²SCAN**：r²SCAN 无 ∇²ρ 依赖，数值更稳定

### 5.3 代码实现

SCANL 已注册为 `dft_functional scanl`，内部映射为 `XC_MGGA_X_SCANL + XC_MGGA_C_SCANL`。

vlapl 势项实现于 `xc_functional_libxc_vxc.cpp` 中的 `v_xc_meta` 函数，使用 G 空间 FD kernel 计算 ∇²(vlapl·sgn)。

vlapl 应力项实现于 `xc_functional_gradcorr.cpp` 中的 `gradcorr` 函数。

### 5.4 替代泛函

| 泛函 | 依赖变量 | 特点 |
|------|---------|------|
| **r²SCAN** | ρ, ∇ρ, τ | SCAN 的正则化版本，数值稳定，推荐首选 |
| **SCAN** | ρ, ∇ρ, τ | 无 ∇²ρ 依赖，数值稳定 |
| **SCANL** | ρ, ∇ρ, ∇²ρ | 需 FD Laplacian 近似，0.5% 应力误差 |
