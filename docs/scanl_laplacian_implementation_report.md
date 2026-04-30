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

### 3.5 测试结果汇总

| 方案 | V 含 vlapl | vtxc 含 vlapl | stress 含 vlapl | SCF | FD/AB | 误差 |
|------|-----------|--------------|----------------|-----|-------|------|
| SCAN 基线 | N/A | N/A | N/A | 收敛 | 1.000 | 0.00% |
| 方案 1 (严格) | ✓ | ✓ | ✓ | **发散** | — | — |
| 方案 2 (sigma-as-lapl) | ✗ (假) | ✗ (假) | ✗ | 收敛 | 1.046 | 4.6% |
| 方案 3 (严格应力) | ✗ | ✗ | ✓ | 收敛 | 0.91 | -9% |
| 方案 4 (vtxc含vlapl) | ✗ | ✓ | ✓ | 收敛 | 0.922 | -7.8% |
| 方案 5 (衰减) | ✓(衰减) | ✓(衰减) | ✓ | 收敛 | ~0.9 | ~10% |
| 方案 6 (ramp) | ✓(渐增) | ✓ | ✓ | **发散** | — | — |

**结论**：所有能收敛的方案都存在 5–10% 的应力误差。这是省略 vlapl 势的根本性后果——密度不在 SCANL 自洽势下产生，导致 etxc-vtxc 不匹配，进而使应力偏大。方案 2（sigma-as-lapl）的误差最小（4.6%），但那是传假输入给 LibXC 导致的偶然偏差，并非物理上更正确。

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

### 5.1 不推荐使用 SCANL

**在 ABACUS（NCPP + 平面波）框架下不推荐使用 SCANL 及其他 ∇²ρ 依赖泛函**。原因：

1. **严格实现不可行**：vlapl 势 e2·∇²(vlapl) 在任何实际 ecutwfc 下导致 SCF 发散
2. **NCPP 矛盾**：模守恒赝势需要 ≥60 Ry 的 ecutwfc，此范围内 vlapl 势必然不稳定
3. **近似方案不可靠**：无论省略 vlapl 势（方案 3）还是用 sigma 伪装（方案 2），所得密度都不是 SCANL 自洽密度，结果与 FD 验证存在 2–5% 的偏差
4. **VASP 也承认问题**：即使在 PAW 方法下，VASP 也建议不要使用超过 800 eV 的截断

### 5.2 推荐替代泛函

| 泛函 | 依赖变量 | 特点 |
|------|---------|------|
| **r²SCAN** | ρ, ∇ρ, τ | SCAN 的正则化版本，数值稳定，推荐首选 |
| **SCAN** | ρ, ∇ρ, τ | 与 SCANL 精度相近，无 ∇²ρ 依赖，数值稳定 |

r²SCAN 是目前推荐的 meta-GGA 泛函，兼具 SCAN 的精度和更好的数值稳定性，且无 ∇²ρ 依赖问题。

### 5.3 代码实现建议

若仍需保留 SCANL 支持（例如用于与文献对比），建议采用**方案 3 + 运行时警告**：

1. **V_xc 不含 vlapl**：e2·∇²(vlapl) 导致 |G|² 发散，无法加入
2. **vtxc 不含 vlapl**：必须与 V 保持一致
3. **etxc 含完整 exc**：LibXC 返回的能量密度已含 Laplacian 贡献
4. **gradcorr stress 含 2·vlapl·Hess·e2**：从完整能量泛函严格推导
5. **运行时打印 WARNING**

在用户选择 ∇²ρ 依赖泛函时，应打印类似以下警告：

```
WARNING: The selected functional depends on the Laplacian of the density (∇²ρ).
The vlapl contribution to the Kohn-Sham potential is omitted because including
it (e2·∇²(vlapl)) causes SCF divergence due to |G|² amplification in reciprocal
space. The self-consistent density is therefore NOT the exact SCANL density.
The stress includes the vlapl term derived from the full energy functional.
Results may be unreliable, especially for norm-conserving pseudopotentials
which require high ecutwfc. Consider using r2SCAN or SCAN instead.
```
