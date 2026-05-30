# `source_md` 分子动力学模块算法说明

本文档说明 ABACUS 中 `source/source_md/` 目录下离子分子动力学（MD）的数学思想、时间积分流程与代码对应关系。物理量默认采用原子单位（a.u.），部分输入量在构造时换算为 a.u.（见各文件中的 `ModuleBase::` 换算因子）。

本文档存储路径：https://github.com/chencx04/abacus-develop/blob/courseHW/source/source_md/ALGORITHMS.md

---

## 1. 模块结构与文件职责

| 文件 | 职责 |
|------|------|
| `run_md.{h,cpp}` | MD 主循环：按 `md_type` 构造积分器，每步调用 `first_half` → 力/应力 → `second_half`，并负责输出、dump、重启。 |
| `md_base.{h,cpp}` | 所有积分器的基类：分配 `pos`/`vel`/`force` 等，实现**半 kick–漂移–半 kick** 中的**位置更新**与**半步速度更新**（与力耦合），以及打印、简单重启。 |
| `md_func.{h,cpp}` | 公共工具：初速、高斯随机数、**`force_virial`（调用任意 `ESolver`，含 LJ/DP/NEP）**、动能/温度、应力张量、轨迹 dump、重启信息读取等。 |
| `verlet.{h,cpp}` | 基于 `MD_base` 的 velocity-Verlet；在 `second_half` 末尾根据 `md_thermostat` 施加简单恒温策略（NVE/NVT 等）。 |
| `nhchain.{h,cpp}` | **Nose–Hoover 链**（文献中常称 NHC）：`md_type == nvt` 且 `md_thermostat == nhc"` 或 `md_type == npt"` 时使用；支持 NPT（各向同性/各向异性/三角晶格模式）。 |
| `langevin.{h,cpp}` | **Langevin** 动力学：阻尼 + 随机力，在速度积分中与真实力合并为 `total_force`。 |
| `fire.{h,cpp}` | **FIRE** 结构弛豫：非真实热浴 MD，通过混合速度/力方向、自适应时间步使体系向能量极小演化。 |
| `msst.{h,cpp}` | **MSST** 冲击波多尺度技术：约束近似 Hugoniot/Rayleigh 线，耦合体积变化率与原子速度。 |
| `CMakeLists.txt` | 将 `md` 编译为 OBJECT 库；在 `BUILD_TESTING` 且 `ENABLE_MPI` 时加入 `test/`。 |

头文件中的 Doxygen 注释描述接口；本文侧重**算法与数据流**。

---

## 2. 主循环与时间步逻辑（`run_md.cpp`）

积分器统一实现三步语义（与具体算法子类有关，顺序可能插入 thermostat/barostat）：

1. **`setup`**：仅在 `step_ == 0` 时调用一次。读重启、初速、调用 `MD_func::force_virial` 得到初始势能、力、virial，并算应力。
2. **每个 MD 步**（从第二步起 `step_ >= 1`）：
   - `first_half`：通常包含至少一次**半步速度**（力已知）和**整步位置**更新（见 `MD_base`）。
   - `MD_func::force_virial`：在新结构下自洽求能量、力、（可选）virial。
   - `second_half`：再次**半步速度**，并在此阶段或前后完成恒温/恒压/MSST 等附加演化。
   - `MD_func::compute_stress` + `current_temp`：更新应力与瞬时温度用于输出。

循环条件为 `(step_ + step_rst_) < md_nstep` 且未 `stop`（FIRE 在力收敛时置 `stop`）。`step_rst_` 来自重启文件，使总步数计数连续。

---

## 3. 速度 Verlet 与 `MD_base`（`md_base.cpp`）

经典哈密顿动力学中，对坐标 \(\mathbf{r}\) 与速度 \(\mathbf{v}\) 有：

\[
\mathbf{v}(t+\tfrac{h}{2}) = \mathbf{v}(t) + \tfrac{h}{2}\,\mathbf{F}(t)/m,\quad
\mathbf{r}(t+h) = \mathbf{r}(t) + h\,\mathbf{v}(t+\tfrac{h}{2}),\quad
\mathbf{v}(t+h) = \mathbf{v}(t+\tfrac{h}{2}) + \tfrac{h}{2}\,\mathbf{F}(t+h)/m.
\]

本模块中：

- **`update_vel(force)`**：实现 \(\mathbf{v} \leftarrow \mathbf{v} + \tfrac{h}{2}\,\mathbf{F}/m\)（对可动自由度；固定方向速度增量为 0）。
- **`update_pos()`**：用当前 \(\mathbf{v}\) 计算**分数坐标位移增量** `pos`（先按 `ionmbl` 屏蔽固定分量，再经 `ucell.GT` 变换），调用 `unitcell::update_pos_taud` 写回晶胞内离子坐标。

`md_dt` 为时间步长（a.u.），由输入 fs 换算。MPI 下 `pos`、`vel` 在 rank 0 更新后 `MPI_Bcast`，保证所有进程晶胞一致。

**注意**：`MD_base::first_half` / `second_half` 对派生类可能被重写；`Verlet` 将 `first_half` 设为与基类相同的 kick–drift 顺序，恒温放在 `second_half` 末尾。

---

## 4. 公共子程序 `MD_func`（`md_func.cpp`）

### 4.1 初速度与温度

- **`gaussrand`**：Box–Muller 型变换，由均匀随机数生成标准正态随机数（成对产生，相位交替返回）。
- **`rand_vel`**：对每个可动分量 \(v_{i\alpha} \sim \mathcal{N}(0, \sigma_i^2)\)，\(\sigma_i=\sqrt{T/m_i}\)。若某晶格方向整体平动被“冻结”（`frozen` 计数），则对该方向减去质心速度，使总动量为零。最后 **`rescale_vel`** 使动能严格满足 \(\frac{1}{2}\sum m v^2 = \frac{1}{2}(3N - f_{\mathrm{frozen}})T\)（\(f_{\mathrm{frozen}}\) 为约束掉的自由度数目）。
- **`init_vel`**：从 STRU 读入初速或随机初速；处理重启温度、`temperature < 0` 时自动采用 STRU 动能温度等分支。

### 4.2 力、能量、应力

- **`force_virial`**：调用 `p_esolver->runner` 后取能量、力矩阵、virial；**将 Rydberg 单位转为 Hartree**（势能与力、virial 乘 `0.5`），与 MD 中其他 a.u. 一致。对接 **DeePMD / NEP** 时的数据准备、库调用均在各 ESolver 的 `runner` 内完成，见 **第 13 节**。
- **`compute_stress`**：\(\boldsymbol{\sigma} = \mathbf{W}/\Omega + \mathbf{K}/\Omega\)，其中 \(\mathbf{W}\) 为电子/离子 virial（由 ESolver），\(\mathbf{K}\) 由 **`temp_vector`** 给出的离子动能张量 \(\sum_i m_i \mathbf{v}_i\mathbf{v}_i^\top\)。

### 4.3 温度定义

- **`current_temp`**：\(T = 2K/(3N - f)\)，其中 \(K\) 为离子动能；若全部自由度冻结则温度为 0。该定义与经典能均分一致。

---

## 5. `Verlet`：NVE 与简单恒温（`verlet.cpp`）

在 `second_half` 末尾 **`apply_thermostat`**：

| `md_thermostat` | 行为概要 |
|-----------------|----------|
| （`md_type == nve`） | 不修改速度。 |
| `rescaling` | 每步若 \(|T_{\mathrm{tgt}}-T_{\mathrm{cur}}|\)（换算为 K）大于 `md_tolerance`，则整体 **`thermalize`**：`v *= sqrt(T_tgt/T_cur)`。 |
| `rescale_v` | 每隔 `md_nraise` 步做一次上述 rescaling。 |
| `anderson` | 以概率 \(1/\texttt{md_nraise}\) 对每个原子速度按高斯重抽样（方差与 `md_tlast`、质量相关），实现随机碰撞型恒温。 |
| `berendsen` | 弱耦合：`fac = sqrt(1 + (T_tgt/T_cur - 1)/md_nraise)`（`nraise>0` 时），否则与简单 rescaling 相同因子形式。 |

目标温度 **`target_temp`**：在 `md_nstep` 内从 `md_tfirst` 线性插值到 `md_tlast`。

---

## 6. `Nose_Hoover`：NHC 与 NPT（`nhchain.cpp`）

### 6.1 思想

引入一组虚构变量 \(\{\eta_m, v_{\eta m}\}\)（热浴链）与（NPT 时）晶格自由度 \(\omega\) 及其共轭动量，使长时间平均逼近 **NVT** 或 **NPT**。实现采用 **可逆参考系统传播算法（RESPA）风格** 的算子分裂：系数 **`w[0..6]`** 为 Yoshida/Suzuki 七阶系数，在 `particle_thermo` / `baro_thermo` 中对 \(\Delta t\) 做多次子步。

### 6.2 粒子热浴 `particle_thermo`

- 更新 `mass_eta[m]` 与 \(g_{\eta 0} \propto (2K - f T_{\mathrm{tgt}})\)。
- 在对称分裂循环中交替：指数因子中的 \(v_{\eta,m+1}\)、更新 \(v_\eta\)、更新 \(\eta\)、用 \(\exp(-v_{\eta0}\Delta/2)\) 累积 **标度因子 `scale`** 作用在原子速度上，并反馈更新 \(g_{\eta 0}\) 中使用的动能 `KE`。

### 6.3 NPT 相关

- **`md_pmode`**：`iso` / `aniso` / `tri` 控制哪些应力分量参与 barostat；`tri` 要求晶格为下三角形式（代码中检查 `latvec` 非对角元）。
- **`target_stress`**：各目标应力分量在 MD 进程中由 `pstart`→`pstop` 线性插值。
- **`couple_stress`**：根据 `md_pcouple` 将应力张量对角元做平均（如 `xyz` 为各向同性压强耦合）。
- **`update_baro`**：由当前应力与目标、晶胞体积等更新 `v_omega`（Voigt 六分量中由 `pflag` 选中的部分）。
- **`vel_baro`**：对原子速度施加与 \(\omega\) 相关的指数因子及剪切项（ MTK 类修正通过 `mtk_term` 进入）。
- **`update_volume`**：对晶格向量按 \(\exp(v_\omega \Delta t/2)\) 等分段更新体积与剪切分量，最后 **`unitcell::setup_cell_after_vc`** 同步晶胞与离子坐标。

NPT 下 `first_half` 中在力 kick 前后插入体积与 barostat 更新，具体顺序见源码（保证与力、位置更新可配合）。

---

## 7. `Langevin`（`langevin.cpp`）

有效力：

\[
\mathbf{F}_{\mathrm{tot}} = \mathbf{F}_{\mathrm{DFT}} - \frac{m}{\gamma}\mathbf{v} + \mathbf{R},
\]

其中 \(\gamma = \texttt{md_damp}\)（换算到 a.u.），\(\mathbf{R}\) 各分量采用均匀随机近似高斯白噪声（幅度 \(\propto \sqrt{24 T_{\mathrm{tgt}} m /(\gamma \Delta t)}\) 的离散化形式），使涨落–耗散与目标温度 \(T_{\mathrm{tgt}}(\mathrm{step})\) 一致。`first_half` / `second_half` 均对 `total_force` 做半 kick。

---

## 8. `FIRE`（`fire.cpp`）

参考文献：*Fast Inertial Relaxation Engine*（PRL 97, 170201 (2006)）。在标准 Verlet 步之间插入：

1. **`check_fire`**（在第一次半 kick 之后、位置更新之前）  
   - 计算 \(P = \sum_i \mathbf{v}_i\cdot\mathbf{F}_i\)、\(|\mathbf{F}|\)、\(|\mathbf{v}|\)。  
   - 混合速度：\(\mathbf{v} \leftarrow (1-\alpha)\mathbf{v} + \alpha \mathbf{F}/|\mathbf{F}|\cdot|\mathbf{v}|\)。  
   - 若 \(P>0\)（力与速度同向，沿下坡）：累积 `negative_count`，若 \(\ge n_{\min}\) 则增大 `md_dt`（上限 `dt_max`）、减小 \(\alpha\)。  
   - 若 \(P\le 0\)：减小时间步、速度清零、\(\alpha\) 重置为 `alpha_start`。

2. **`check_force`**：若最大力分量（Hartree/Bohr）满足 `2*max < force_thr`（注意阈值在构造函数中被设为 `1e-3` 覆盖输入的 `force_thr`），则 `stop = true`。

重启文件除步数与温度外还保存 \(\alpha\)、`negative_count`、`dt_max`、`md_dt`。

---

## 9. `MSST`（`msst.cpp`）

用于模拟冲击波条件下沿 **`msst_direction`** 方向的体积/速度演化：

- 维护 \(\dot{V}\)（代码中 `omega[sd]`）、初始 \(E_0,P_0,V_0\)、拉格朗日坐标 `lag_pos` 等。
- **`propagate_voldot`**：由当前应力、体积、人工粘度 `msst_vis`、质量参数 `msst_qmass` 与冲击速度 `msst_vel` 更新 \(\dot{V}\)（含防止体积非物理增长的钳制与 \(\dot{V}\approx 0\) 的泰勒展开）。
- **`propagate_vel`**：在已知 \(\dot{V}\)、\(v^2\) 和与 \(\omega^2\) 相关的阻尼下对速度做半解析指数积分或二阶泰勒。
- **`rescale`**：按目标体积比缩放晶格对角元（各向同性膨胀/压缩选定方向），并重标度该方向速度分量；调用 `setup_cell_after_vc`。

`setup` 中可选 `msst_tscale` 初始化 \(\dot{V}\) 并缩放初速。

---

## 10. 单位、MPI 与重启

- 输入时间、温度、压强等常在 `MD_base` 或各子类构造中通过 `ModuleBase` 常数转为 a.u.。
- 大量操作仅在 `my_rank == 0` 执行，再 `MPI_Bcast` 速度/位置/力等，保证并行一致性。
- 通用重启：`Restart_md.txt`（步数、温度等）；FIRE、NHC、MSST 扩展不同字段，**不可混用不同积分器的重启文件格式**。

---

## 11. 与测试目录的关系

`test/` 下 GTest 用例（如 `md_func_test.cpp`、`verlet_test.cpp` 等）覆盖 `MD_func` 与各类积分器在小型晶胞上的行为；构建需 `-DBUILD_TESTING=ON` 且 `ENABLE_MPI=ON`。若配置了 **DeePMD_DIR** / **NEP_DIR**，还会生成 **`md_dp_mdfunc_test.cpp`** / **`md_nep_mdfunc_test.cpp`**（见 `test/CMakeLists.txt`），用于验证 `MD_func::force_virial` 与 ML 势的联调。

---

## 12. 阅读代码的建议顺序

1. `run_md.cpp` → 把握全局步进。  
2. `md_base.cpp` → 理解 kick/drift 与 `pos` 含义。  
3. `md_func.cpp` → 力、应力、温度、初速。  
4. 按需深入 `verlet.cpp` / `nhchain.cpp` / `langevin.cpp` / `fire.cpp` / `msst.cpp`。

若需与 INPUT 关键词一一对应，请同时查阅 `source_io/module_parameter` 中与 `MD_para` / `mdp` 相关的定义。

---

## 13. DeePMD 与 NEP：势函数如何计算、代码在何处、如何接入 `source_md`

本节说明两个机器学习势在 **ABACUS 中的实现位置**（不在 `source_md` 内）、**单次“算力”调用链**，以及 **`source_md` 如何只通过统一接口使用它们**。

### 13.1 架构分工（重要）

| 层次 | 目录 / 文件 | 作用 |
|------|-------------|------|
| **势与力、virial 的实现** | `source/source_esolver/esolver_dp.{h,cpp}`、`esolver_nep.{h,cpp}` | 把 `UnitCell` 转成 DeePMD-kit / NEP 库需要的数组，调用 `dp.compute` / `nep.compute`，并把结果存成 ABACUS 内部使用的能量、力、virial 矩阵。 |
| **MD 与 ESolver 的胶水** | `source/source_md/md_func.cpp` 中的 **`MD_func::force_virial`** | 每步（或 setup）调用 `p_esolver->runner`，再 `cal_energy` / `cal_force` / `cal_stress`，并把 ESolver 常用的 **Ry 约定** 转为 MD 使用的 **Hartree（×0.5）**。 |
| **时间积分与主循环** | `source/source_md/run_md.cpp`、`md_base.cpp` 及各积分器 | 只持有 **`ModuleESolver::ESolver* p_esolver`**，不区分 LJ / DP / NEP；位置、速度更新与恒温恒压逻辑与势函数无关。 |
| **ESolver 的创建** | `source/source_main/driver_run.cpp`、`source/source_esolver/esolver.cpp` | 根据 `INPUT` 中 `esolver_type` 构造具体子类（`dp` → `ESolver_DP`，`nep` → `ESolver_NEP`），再把指针传入 `Run_MD::md_line`。 |

因此：**DeePMD/NEP 的“公式/网络推理”在 ESolver 与外部库里完成**；`source_md` **不包含** DeepPot 或 NEP 的推理代码，只通过虚接口驱动。

### 13.2 从 INPUT 到 `md_line`：谁创建 `ESolver_DP` / `ESolver_NEP`

1. 用户在 `INPUT` 中设置 **`esolver_type dp`** 或 **`nep`**，并设置 **`pot_file`**（势文件路径，见 `mdp.pot_file` 等参数说明）。  
2. `source/source_esolver/esolver.cpp` 在初始化阶段把字符串规范为内部类型名（如 `"dp"` → `"dp_pot"`，`"nep"` → `"nep_pot"`），再通过工厂函数 **`init_esolver`（同类逻辑）** 执行 `new ESolver_DP(PARAM.mdp.pot_file)` 或 `new ESolver_NEP(PARAM.mdp.pot_file)`。  
3. `source/source_main/driver_run.cpp` 在 `calculation == md` 等分支里调用 **`Run_MD::md_line(ucell, p_esolver, PARAM)`**，把已构造好的 **`p_esolver`** 传入 MD 模块。

此后 MD 模块**始终**把 `p_esolver` 当作抽象基类 `ESolver` 使用。

### 13.3 `source_md` 内唯一“势相关”的核心调用：`MD_func::force_virial`

实现见 `md_func.cpp` 中 **`MD_func::force_virial`**，逻辑可概括为：

1. **`p_esolver->runner(unit_in, istep)`**  
   要求各具体 ESolver 根据**当前** `UnitCell`（晶格、离子分数坐标/笛卡尔坐标在库里的约定见下）更新内部缓存。对 DP/NEP 而言，**真正调用 ML 库**发生在 `runner` 内。  
2. **`potential = p_esolver->cal_energy()`**  
3. **`p_esolver->cal_force(...)`** → 填 `force_temp`  
4. 若 `cal_stress`：**`p_esolver->cal_stress(...)`** → 填 `virial`（此处语义为晶格 virial，后续与离子动能项合成总应力见 `compute_stress`）  
5. **统一单位**：代码假定 ESolver 侧能量/力/virial 为 **Rydberg 相关约定**，MD 使用 **Hartree**，因此对 `potential`、`force_temp`、`virial` **乘以 `0.5`**，再拷入 `force[]`。

因此：**任意**与 MD 对接的 ESolver（LJ、DP、NEP、DFT…）都必须与上述 **Ry → Ha** 约定一致，或由 `force_virial` 统一折算。

### 13.4 DeePMD（`ESolver_DP`）：如何计算、数据如何组织

**源码位置**：`source/source_esolver/esolver_dp.cpp`（声明 `esolver_dp.h`）。仅在编译定义 **`__DPMD`** 且链接 DeePMD-kit 时可用（顶层 CMake 配置 `DeePMD_DIR`）。

**初始化 `before_all_runners`**：

- 分配 `dp_force`、`dp_virial`，扩展 `atype`；从 `Input_para::mdp` 读取 **`dp_rescaling`**、**`dp_fparam`**、**`dp_aparam`**。  
- **`type_map(ucell)`**：从 frozen graph 读取元素顺序，与 `UnitCell` 中 `atoms[it].label` 对应，为每个原子填整数 **`atype[iat]`**，供 `dp.compute` 使用。  
- 可选写出 `STRU.cif` 便于检查结构。

**每步 `runner`（核心）**：

1. **晶格 `cell[9]`（Å，行主序 3×3）**  
   \(h_{ij} = \texttt{latvec}(i,j) \times \texttt{lat0\_angstrom}\)，与 DeePMD-kit 文档中“以 Å 为长度的晶格矩阵”一致。  
2. **坐标 `coord[3*nat]`（Å，笛卡尔）**  
   对每个原子：\(\mathbf{r}_\text{Å} = \tau \times \texttt{lat0\_angstrom}\)，其中 \(\tau\) 为 ABACUS 中存储的笛卡尔坐标（与 `latvec` 约定一致，见 `tau` 的更新方式）。  
3. 调用 **`dp.compute(dp_potential, f, v, coord, atype, cell, fparam, aparam)`**（C++ 或 C API 由 `__DPMDC` 区分），得到模型原生输出。  
4. **换算到 ABACUS 内部 Ry 系**（再交给 `force_virial` 做 ×0.5）：  
   - `fact_e = rescaling / Ry_to_eV` 等，对能量、力分量、virial 分量分别缩放（力含 `ANGSTROM_AU`，virial 含 `1/omega` 与 `Ry_to_eV`）。  
5. **`cal_stress`** 中还可按 `INPUT` 的 **`press1/2/3`** 减去外压对角项（与 DFT 路径一致）。

**小结**：DeePMD 一步 = **结构 → Å 制 cell/coord + 整数类型 → `dp.compute` → 缩放到 Ry 系力/能/virial**；MD 再在 `force_virial` 中转 Hartree。

### 13.5 NEP（`ESolver_NEP`）：如何计算、与 DeePMD 的差异

**源码位置**：`source/source_esolver/esolver_nep.cpp`（声明 `esolver_nep.h`）。仅在 **`__NEP`** 且链接 NEP 库时可用（顶层 CMake 配置 `NEP_DIR`）。

**初始化 `before_all_runners`**：

- 分配力、virial 及 NEP 返回缓存 `_e`、`_f`、`_v`。  
- **`type_map(ucell)`**：根据 NEP 模型文件中的 **元素列表顺序**，为每个原子建立 `atype[iat]`。

**每步 `runner`（核心）**：

1. **晶格按列主序（column-major）填 9 元**（注释 *“NEP are column major, thus a transpose is needed”*）：  
   与 `ESolver_DP` 的行主序 `cell[0]=e11, cell[1]=e12, …` 不同，NEP 使用  
   `(e11,e21,e31, e12,e22,e32, e13,e23,e33) × lat0_angstrom`。  
2. **坐标同样为 Å**；但 **`coord` 的布局为 NEP 常用“按分量块存”**：  
   `coord[i] = x_i`，`coord[i+nat] = y_i`，`coord[i+2*nat] = z_i`（与 DP 的 `3*iat` 连续三元组不同）。  
3. 调用 **`nep.compute(atype, cell, coord, _e, _f, _v)`**。  
4. **能量**：对每个原子的 `_e[i]` 求和，再乘 `fact_e = 1/Ry_to_eV`，得到 **`nep_potential`（内部以 Ry 为能量单位存储，与 DP 路径在 `force_virial` 处对齐）**。  
5. **力**：从 `_f` 的块布局解出 `nep_force(i,0/1/2)`，乘 `fact_f`。  
6. **virial**：`_v` 含每原子对 virial 的贡献，先在 9 个分量上对所有原子求和得 `v_sum`，再写入 `nep_virial(i,j)` 并乘 `fact_v`。同样 **`cal_stress`** 会处理外压。

**小结**：NEP 一步 = **转列主序 cell + 块布局 coord → `nep.compute` → 缩放到 Ry 系**；之后与 DP 一样进入 `MD_func::force_virial` 的 ×0.5 流程。

### 13.6 `source_md` 里具体在哪些地方“跑到” ML 势

以下调用链对 **DP / NEP / LJ / …** 完全相同，仅 `p_esolver` 的动态类型不同：

| 阶段 | 文件与函数 | 说明 |
|------|------------|------|
| 进入 MD | `driver_run.cpp` → `Run_MD::md_line` | 传入已在别处构造的 `p_esolver`。 |
| 首次能量/力 | `md_base.cpp` **`MD_base::setup`** | 内部调用 **`MD_func::force_virial(p_esolver, …)`**，完成第一步自洽力（对 ML 势即一次前向推理 + 力/virial）。 |
| 每一 MD 步 | `run_md.cpp` 主循环 | `first_half` 更新坐标后再次 **`MD_func::force_virial`**，再 `second_half`。 |

积分器类（`Verlet`、`Nose_Hoover` 等）中的 **`update_vel(force)`** 所使用的 **`force`**，即上一步 `force_virial` 写入、**已处于 Hartree/Bohr 约定** 的数组。

### 13.7 与单元测试的对应关系（可选阅读）

- `source/source_md/test/md_dp_mdfunc_test.cpp`：在配置 **`DeePMD_DIR`** 时编译，直接构造 **`ESolver_DP`** + 小体系，调用 **`MD_func::force_virial`**，验证与 DeePMD 联调。  
- `source/source_md/test/md_nep_mdfunc_test.cpp`：在配置 **`NEP_DIR`** 时编译，构造 **`ESolver_NEP`** + 与集成测试一致的 STRU/模型路径，同样走 **`force_virial`**。

二者都**不**在 `source_md` 内实现势函数，只验证 **“MD 胶水层 + ESolver”** 整条路径。

### 13.8 编译与依赖（备忘）

- **DeePMD**：CMake `DeePMD_DIR` → 定义 **`__DPMD`**，并链接 `deepmd_c` 或 `deepmd_cc`（见 `cmake/FindDeePMD.cmake`）。  
- **NEP**：CMake `NEP_DIR` → 定义 **`__NEP`**，并链接 **`NEP::nep`**（见 `cmake/FindNEP.cmake`）。

未配置上述变量时，主程序仍可编译（对应 ESolver 在 `runner` 内 `WARNING_QUIT`），`source_md` 的 MD 逻辑本身不依赖 ML 库。
