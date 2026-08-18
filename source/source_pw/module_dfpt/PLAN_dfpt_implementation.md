# ABACUS DFPT 完整实施计划

> 本文件记录 DFPT（密度泛函微扰理论）落地计划。执行状态见末尾"进度"章节。

## 总原则

- 顺序：**U0（DFT+U 预留）→ C 物理主体（C0–C7）→ B 数据层收编 → A irrep 分解**。C4（Metal）仅留接口。
- 每个子任务交付 = **代码 + 单元测试 + 物理对照**；阶段结束先 git 提交再继续。
- 生产代码零新增 `GlobalV`/`PARAM` 依赖（module_dfpt 不读 PARAM，决策由 esolver 接线层做）；C++11、LF、新文件进 CMakeLists。
- 验证体系：金刚石 `stru_lib[0]`（O_h，a=1，1 个 C 原子，Γ 网格）为主；沙箱 OpenMPI 不作为回归基准。
- 构建：`cmake --build build --target <target> -j8`；ctest 回归过滤：`"MODULE_CELL_klist_test$|MODULE_CELL_reciprocal_grid_test|MODULE_CELL_qlist_test|MODULE_CELL_little_group_test|MODULE_DFPT"`。
- 治理：`python3 tools/03_code_analysis/agent_governance_check.py --base upstream/develop --head HEAD --format text`。

---

## U0 — DFT+U 接口预留（数据+管线+桩+测试）

前提（已核实）：PW DFT+U 存在且已接线（`esolver_ks_pw.cpp:203` iter_init_dftu_pw、`hamilt_pw.cpp:126` OnsiteProj 入链、`setup_pot.cpp:75` OnsiteProjector 初始化），**但依赖 LCAO 轨道文件**（纯 PW 无文件则 locale 未初始化、实际跑不了）。设计据此自洽 on/off。

1. **`dfpt_pw_data.h/.cpp`**：
   - 头内前向声明 `class Plus_U;`（保头文件依赖最小）。
   - `init(..., int nat, const Plus_U* dftu)` 末位加参数；新增 `with_u()`（= 指针非空）、`u_active()`（= 非空 **且 `dftu_->is_locale_initialized()`**，覆盖无轨道文件退化）、`get_dftu()`、`set_docc/get_docc`（每 q complex 向量，惰性分配）。
   - 更新调用点：`dfpt_pw.cpp:71`、`dfpt_irrep_data_test.cpp:172`。
2. **`dfpt_pw.h/.cpp`**：`init(ucell, psi, nelec, ecutwfc, const Plus_U* dftu)`；头内前向声明、Impl 存指针；新增 `get_with_u()/get_u_active()`；更新测试调用点（传 nullptr）+ README 示例 + esolver 注释行。
3. **`dfpt_rho.h/.cpp`**：新增 `cal_docc(psi, wg, q_idx, data)` 桩（`if(!data.with_u()) return;`，注明 C3 实装点）。
4. **`dfpt_pert.h/.cpp`**：`build_dv` 末尾 `if(data.with_u()) build_dv_u(...)`；私有桩 `build_dv_u` 空实现（注明 C1 实装 frozen 项）。
5. **`dfpt_phon.h/.cpp`**：新增 `dftu_onsite(q_idx, data)` 桩，`assemble` 中 when with_u 调用（C5 实装）。
6. **`dfpt_q0.h/.cpp`**：仅加非局域 `[r,V_U]` commutator 预留注释（C6 延后）。
7. **测试**：`dfpt_irrep_data_test` 更新 init + docc roundtrip + with_u=false 安全路径；`dfpt_pw_run_test` SOURCES 加 `../../../source_lcao/module_dftu/dftu.cpp`，新增 `Plus_U dftu;` 传非空 → `u_active()` 为 false（locale 未初始化）、`run()` 不崩、桩 safe（覆盖"无轨道文件"退化路径）。
8. 验证：构建 2 测试目标 + 6 目标回归 + 治理。提交。

**DFT+U 物理特殊处理清单（后续实装点）**：
- 一阶占据矩阵 docc：交叉项 `becp(k+q,dψ)·becp(k,ψ)`（C3，依赖 dψ）+ 冻结项 `becp(k,ψ)·dbecp_f(k,ψ)`（GS k，复用 `cal_dbecp_f`）。
- 一阶 U 势 dV_U（C1 frozen 实装）：占据响应 `|φ(k+q)⟩U(diag·δ−docc)⟨φ(k)|ψ⟩`（需 SCF 自洽）+ 冻结 `|∂φ(k+q)/∂τ⟩V_eff⟨φ(k)|ψ⟩` 等（`Onsite_Proj_tools` 用 DFPT 的 k+q 基初始化即复用，相因子自动正确）。
- Stern 的 H(k+q) 含零级 V_U：复用含 OnsiteProj 的 ops 链即自动覆盖。
- 动力学矩阵 U 项（C5 `dftu_onsite`）、Q0 非局域项（C6）。
- 治理：DFPT 内层 SCF 绝不调 `cal_occ_pw`（防覆盖 GS locale）；零级 V_U 经 `get_eff_pot_pw_spin` 只读借用。

---

## C 阶段物理主体

**C0 — k+q 平面波基枚举**
- k+q 基枚举 helper：复用 `pw_basis_k.h` 的 `npwk/ig2ixyz_k/getgpluskcar`；生成每 (ik,q) 的 k+q 波矢 G 列表。
- 单元测试（核对 G 集合与 Gamma 平移关系）。构建+测试后提交。

**C1 — DFPT_Pert 扰动构建**
- `dVloc_dtau`/`dVnl_dtau` 真实实现（q 相因子、USPP 投影子导数）。
- `build_dv(q_idx, atom_idx, dir, data)` 组装 dV；`apply_dv`；`build_efield`。
- **U0 实装点**：`build_dv_u` frozen 项（OnsiteProjector 已初始化时启用）。
- 测试：dV 数值核对。提交。

**C2 — DFPT_Stern Sternheimer 求解**
- 移位哈密顿作用经 `LinearOperator` 抽象注入（dimension/apply）：生产适配器复用 `ops->hPsi(hpsi_info)`（hsolver_pw.cpp:268-274 模式，C7 接线），单测注入解析算子。
- 无状态 `solve(aop, occ_kq, b, max_iter, conv_thr, dpsi, residual)`：投影 CG（初始 x=0、r=P_c b，每步搜索方向重投影，收敛判据 `||P_c r||/||P_c b||`）；`apply_pv` 双扫 MGS 投影（alias-safe）；方程 `(H(k+q)−ε_n)P_c|dψ_n⟩=−P_c dV|ψ_n⟩`。
- 测试：对角算子闭式解 + 稠密 Hermitian（U D U† 谱展开参考，eps 落占据带内验证投影）+ 正交性保持 + 退化 RHS。提交。

**C3 — DFPT_Rho 密度响应**
- `compute_drho`：交叉核 `A(r)=Σ_{kn occ} wg·u*·du`（K/rho 两条 FFT 路径同网格点乘，均无 k 相位 → Bloch 相位合并为单个 e^{iq·r}）；`drho_g` 为 rho 网格 q 移位系数 `A_Δ=Σ wg Σ_G c*_G d_{G+Δ}`；Δ=−q 落在倒格矢时投影为零（电荷守恒，q=Γ 必触发）；`drho_r=2Re[e^{iqr}A]` 由投影后系数重建。USPP 增广项与 nspin>1 留 WARNING_QUIT 守卫（设计期）。
- `mix_drho` 直接复用 `Base_Mixing::Plain_Mixing::plain_mix`（q 移位复空间混合 + 残差 `||out−in||/||out||`，首步 in=0），不套 `Charge::rho`/`Charge_Mixing`（头依赖由 charge_mixing.h 降为前向声明 + matrix3.h 值成员）。
- **U0 实装点**：`cal_docc` 交叉项需 k/k+q 双端 β 投影子（PW 侧 vkb 适配器），与 Plus_U 生产接线同落 C7/U1 窗口；纯 PW 路径 `u_active()` 恒 false，安全退化保持。
- 测试：G 空间暴力双和对照 + 实空间直接求和对照 + Γ 电荷守恒（ig0 置零 + 网格和≈0）+ 混合两步组合与残差公式。提交。

**C4 — DFPT_Metal（仅接口）**
- 本期不实现：`compute_drho`/occupation 响应留接口与设计说明（`is_metal_`/`dmu_` 数据已备）。

**C5 — DFPT_Phon 动力学矩阵**
- `assemble`：`ion_ion(q,dynmat)` + `electron(q_idx,data,dynmat)` 真实现（Ewald α 选取复用 `force_pw.cpp:479 cal_force_ew` 惯例：1.1 起步递降、upperbound<1e-6、跳 `ig_gge0`；仓库中 `H_Ewald_pw::rgen` 已不存在，以 `cal_force_ew` 为蓝本；相因子经 `symm.gtrans[48]`+`kgmatrix[48]`）。
- `electron` 拆两步：`accumulate_electron(q,atom,dir,psi,wg,data)`（2n+1 形式 `D_ab=(2/Nk)Σ_kn wg·⟨dψ^b|dV^a_ext|ψ⟩` 复数逐 k 累积，虚部经 k-star 配对共轭相消/assemble Hermitian 对称化吸收；run() 每方向 SCF 收敛后即时调用，dpsi 存储不加方向维度，避免 B 前翻搅数据层）+ `assemble` 合并 ion_ion + 累积电子项。
- **U0 实装点**：`dftu_onsite`/`dftu_lambda` 电子项。
- `diagonalize` 真实现（`LapackConnector::zheev`，`ω=sign(e)√|e|` 换算 cm⁻¹，质量因子 1/√(MM′)）；`add_loto` 非解析项 `(4π/Ω)(q·Z*_a)(q·Z*_b)/(q·ε∞·q)/√(MM′)`；`check_sum_rule` Γ 声学 3 零模 + 列和。
- 测试：ion_ion vs 小胞朴素双和；Γ ASR；accumulate_electron vs 注入 dpsi 闭式收缩；zheev vs 已知矩阵；loto 各向同性解析极限。提交。

**C6 — DFPT_Q0 介电/Born/LO-TO**
- 新增 `v_hartree_q`：`dV_H(G)=4π|G+q|²⁻¹·drho_g`，跳过 |G+q|=0（ig=−q），约定对齐 `source_estate/module_pot/h_hartree_pw.cpp:16`；实现为 `DFPT_Rho` 成员，**同时服务 C7 全 q 点 SCF 屏蔽势**。
- **XC 一阶核（回调注入复用 `PotXC_FDM`）**：库里已有 `elecstate::PotXC_FDM`（`source_estate/module_pot/pot_xc_fdm.cpp:39`，`δV_xc=V_xc[ρ₀+δρ]−V_xc[ρ₀]`，LCAO 侧 `veff_dh.cpp:417 cal_dH_hf_xc` 已验证同构物理）。module_dfpt 内只定义回调契约 `XC_First_Order`（抽象类，镜像 `DFPT_Stern::LinearOperator` 惯例），esolver 接线层写 `PotXC_FDM_Adapter` 注入；module_dfpt 不 include `pot_xc_fdm.h`（头依赖最小）。复 δρ 拆 Re/Im 两次调用再重组（线性叠加合法，误差 O(δρ²)）。`LR::KernelXC`（LIBXC 解析核，module_lr）列为后续优化，本轮不动。
- `pos_matrix`：不走病态位置算符，用 `[Ĥ_SCF,r]` 速度算符等价式 `⟨u_m|r|u_n⟩=(ε_m−ε_n)⁻¹⟨u_m|[V_nl,r]|u_n⟩`（m≠n）；`[V_nl,r]` 复用 C1 `build_vkb` 平移基列表求导。
- 非局域 `[r,V_U]` commutator 项记录为 U 预留。
- 测试：v_hartree_q 单 G 闭式；XC 回调 vs 解析 LDA 核 `(4/9)v_xc/ρ`；ε/Z* 金刚石对称性约束；loto 方向极限。提交。

**C7 — run() 接线 + ESolver/INPUT**
- `DFPT_PW::init` 扩签名：`(..., pw_rho, pw_wfc, sf, wg, eig, const XC_First_Order* xc)`（规则 5：不加默认参，全调用点更新）；`nrxx=pw_rho->nrxx`；`pert_/rho_/q0_/phon_` 真 init。
- Stern 生产适配器 `HamiltShiftAdapter : DFPT_Stern::LinearOperator` 包 `p_hamilt->ops->hPsi`（`hsolver_pw.cpp:271-273` hpsi_info 模式）；占据态 `occ_kq` 由 GS ψ 经 k+q 基投影（复用 `apply_pv` MGS）。
- `run()` 实装：mode basis 为空（A 前置占位）→ 回退遍历 3N 方向；SCF 内环 `dv_sc=dv_ext+v_hartree_q+xc_->apply(drho_in)` → `stern_.solve` → `compute_drho` → `mix_drho`，残差<conv_thr 收敛 → `accumulate_electron`；全方向后 `assemble+diagonalize`（q=0 且 loto → `add_loto`）。单 k 基 k+q 球覆盖 WARNING_QUIT 守卫（C1 遗留边界，`ecutrho>=4*ecutwfc` 写入文档）。
- `esolver_dfpt_pw.cpp`：解开 `init` 注释，实参 `*this->stp.psi_cpu` + `PARAM.inp.nelec` + `PARAM.inp.ecutwfc` + `this->pw_wfc`/`this->pw_rho`/`this->sf` + `dft_plus_u ? &this->dftu : nullptr`（`esolver_ks.h:63`）+ PotXC_FDM 适配器（持 GS `Charge`，复 δρ 拆 Re/Im）；内层 SCF 禁调 `cal_occ_pw`（U0 治理条目）；`esolver.cpp` 工厂加 `"dfpt"` 分支。
- INPUT 行为若变则同步 `docs/parameters.yaml` + `input-main.md`；验证 `./build/abacus -h esolver_type` 与 `--check-input`。
- 金刚石端到端对照：声学 3 零模（ASR）+ 光学支 LDA 文献区间 + 介电/Born；`./build/abacus --version` 记录身份。
- 全量构建 + 回归 + 治理。提交。

---

## B — 工程化收编（2026-08-18 修订：先全流程工程验证，对称性后移）

> 决策记录：INPUT 主文件参数（非 dfpt.in 子文件）；irrep 保留接口、优先全流程工程验证；
> 调试插桩暂缓清理（A 阶段验证后统一收尾）。执行序：B1 → B0 → B2 → B3 → B4，每节点
> 完成 = 代码 + 构建/回归/治理 + git 提交 + 本文档进度回写，再进入下一节点。

**B1 — INPUT 参数接线（解除硬编码/死代码）**
- 新增 dfpt 前缀 INPUT 参数：`dfpt_qmesh`(3 int) / `dfpt_qfile`(str，走 `QList::read_from_file`)
  / `dfpt_compute_q0`(bool) / `dfpt_loto`(bool) / `dfpt_conv_thr` / `dfpt_max_iter` / `dfpt_mix_beta`。
- `esolver_dfpt_pw.cpp` 删除硬编码（set_qmesh(1,1,1)/conv_thr/max_iter/空桩 set_parameters），
  改从 `inp` 显式传递（规则 1）；`set_compute_q0`/`set_loto` 死代码开关经此激活。
- `docs/parameters.yaml` + `input-main.md` 同步；`-h dfpt_*` 与 `--check-input` 验证。

**B0 — 全流程工程验证（数值验收基线）**
- 真实金刚石 Γ 点（a≈3.567 Å、NC PP、收敛 k 网格；玩具胞 742 cm⁻¹ 仅为 FD 锁定基线）：
  声子（ASR + 光学支 vs LDA 文献）、ε∞（各向同性 ≈5.3–5.7）、Z*、LO-TO 方向依赖。
- 非 Γ q：`dfpt_qmesh` + 公度 k 网格（`build_occ_kq` 公度守卫已备）验证 q≠0 全流程。
- MPI>1 rank 冒烟（当前打印注释自述 single-rank，未验证面）。
- `--version`/`-h`/`--check-input` 记录；结论回写本文档。

**B2 — ε∞/Z*/LO-TO 输出正式化**
- esolver design-phase std::cout 转正式输出：多 q 布局、LO-TO 修正后频率（每方向）。
- loto 方向经数据层传递，消除 run() 中 (1,1,1)/√3 硬编码。

**B3 — Kerker 型预条件混合**
- `DFPT_Rho` 内自实现 `|G+q|²/(|G+q|²+a²)` 预条件（不引 charge_mixing.h，module_base 无
  现成 Kerker 已核实）；mix_type 支持 plain/kerker。
- 验收：λ_A1≈−2.2 模型问题 β=0.7 收敛（β 上界 0.62 解除，JPROBE 复用为验收工具）；
  金刚石频率与 β 无关（固定点正确性）；`dfpt_mix_beta` 默认回调并文档记录。

**B4 — 数据层收编**
- 收敛台账（converged_/residuals_/current_iter_ 按 (q,irrep)）并入 `DFPT_PW_Data`；
  删除 `DFPT_IrrepData` 适配层与 `get_dpsi_obj` static dummy；测试迁移；
  **保留 (q,irrep) 接口形状**（为 irrep 实装留插槽）；run() 外层 while 记账语义梳理。

**暂缓项（接口保留，后续单独立项）**
- A irrep 分解：LittleGroup 占位（nirr≡1/空 basis）与 run() 3N 回退保持现状。
- 插桩清理（PTCHK/DYNCHK/MDBG/JPROBE/VKBCHK/DFPT_MIX_BETA）：B3 验收仍需 JPROBE。
- KVectorUtils 薄封装删除：随 A 阶段收尾一并处理。

## A — irrep 分解（最后，工程验证完成后立项）

- `module_symmetry/little_group.{h,cpp}`：完整不可约表示表 + 投影算子 → 真实 `get_nirr`/`get_mode_basis`（替换占位 =1/空）。
- 测试：金刚石/闪锌矿 Γ/X/L 点 irrep 分解与理论表核对。提交。

---

## 风险与注意事项

- `Plus_U dftu` 是对象成员，指针永远非空 → `with_u` 由 esolver 在 `dft_plus_u` 时传非空指针决定，语义干净。
- 无轨道文件 → `u_active()` 为 false 安全退化；测试显式覆盖该路径。
- `dftu.cpp` 加入测试链接依赖 LCAO=OFF 配置；若 CI 在 LCAO=ON 下编译该测试需补 hamilt 依赖（记录在案）。
- 沙箱 OpenMPI 警告（`opal_ifinit`）为环境产物，不作为失败判据。

## 进度

- [x] 计划定稿（U0/C/B/A 序列、DFT+U 依赖轨道文件的 on/off 自洽设计）
- [x] U0 DFT+U 接口预留 `3599973fa`
    - `Plus_U*` 经 `DFPT_PW::init`/`DFPT_PW_Data` 线程化（esolver 接线层决策，module_dfpt 不读全局输入）
    - `with_u()`/`u_active()`（locale 未初始化即无轨道文件 → 安全退化）+ 每 q `docc_` 槽位
    - 桩：`DFPT_Rho::cal_docc`、`DFPT_Pert::build_dv_u`、`DFPT_Phon::dftu_onsite`、Q0 `[r,V_U]` 注释
    - 测试 5+3 全过；`dftu_test_support.cpp` shim（免 LCAO 侧链接闭包）；6 目标回归 + `abacus_pw_para` 链接通过
- [x] C0 k+q 平面波基枚举
    - `DFPT_KQ_Basis`：复用 GS 复杂 k 基的共享 G 网格，仅做 k+q 平移中心再过滤（`|G+k+q|^2<=gk_ecut`），无需新建 FFT 网格/重分发；前置条件 gamma_only=false + 网格截断覆盖 k+q 球（`gridecut_lat >= (sqrt(gk_ecut)+max|k|+max|q|)^2`，ecutrho>=4*ecutwfc 满足）
    - `get_npwk/get_ig/get_ig2isz/get_gcar/get_gpluskq/get_gk2/get_kplusq` 访问器；gamma_only 守卫 WARNING_QUIT
    - 测试 5 项全过（Γ q=0 全等复现、偏心非对称球、k+q 平移不变性、非零 q 与全网格穷举对照、null/gamma_only 拒绝）；7 目标回归
- [x] C1 DFPT_Pert
    - `dVloc_dtau`：rho 网格系数 `i·tpiba·(Δ+q)_dir·Vloc(|Δ+q|)·e^{i2π(Δ+q)·τ}`（Δ+q=0 分量剔除）；`vloc_at_g` Coulomb 解析式（复刻 `vl_pw::vloc_coulomb`）+ numeric 径向 FT（复刻 `vloc_of_g` 含 erf 补偿）
    - `dVnl_dtau`：NC 分离算符两项恒等式 `i·tpiba·(k+q+G'')_dir·(Vnl|ψ⟩ − Vnl[i·tpiba·(k+G')_dir|ψ⟩`；`build_vkb`（(−i)^l·Y_lm·(4π/√Ω)∫β j_l r dr·e^{i2π·gk·τ}，GS 约定对齐）+ `radial_vq`（Simpson）+ `real_ylm`（l≤2）；USPP/SOC WARNING_QUIT 守卫
    - `build_dv`→`set_dv_recip_c`→recip2real→`set_dv_rc`；`apply_dv`（纯循环卷积，q 相位已并入系数）+ `build_efield`（−E·r）+ `build_dv_u`（u_active 守卫，C7 激活）
    - 串行测试目录 `test_serial/`（`__MPI` 整体关闭 + `dfpt_planewave_serial` OBJECT 库，ABI 一致）8 项全过：rho_gvec≡gcar、dVloc 有限差分（含 q≠0/双方向）、apply_dv 卷积 vs 解析矩阵元、efield 斜坡闭式 FT、build_vkb 独立 Simpson+τ 纯相位、dVnl 两项恒等式 vs 算符有限差分、USPP 拒绝、with_u/u_active 安全退化
    - 测试捕获并修复 3 处约定/实现错误：① 相位幅角 `tpiba·(w·τ)` → `TWO_PI·(w·τ)`（GS `stru_fac` 的 e^{i2π(g·τ)} 约定，tau 为 lat0 单位）；② 实空间布局 `ir=(ix·ny+iy)·nz+iz`（z 最快，冲击响应探针钉死；build_efield 原假设反向）；③ rho/wfc 棒表枚举不同 G 球 → isz 编码不可互换，`real_space_dv` 改经 FFT 胞 (ix,iy,iz) 三元组反查
    - 已知边界（C7 处理）：单 k 基时 k+q 球需 wfc G 列表含 `sqrt(gk_ecut)+|k+q|` 半径（k 网格覆盖或 inflate）；并行 pool 实空间布局
    - 8 目标回归全过（CELL 4 + DFPT 4）；`abacus_pw_para` 链接通过
- [x] C2 DFPT_Stern
    - 无状态投影 CG：`DFPT_Stern::solve`（x=0 起步，α=|r|²/(pᵀAp)、β=|r_new|²/|r_old|²，搜索方向每步 `P_c` 重投影；pAp≤0 时残差方向重启）；`apply_pv` 双扫 MGS（alias-safe）；收敛 `||P_c r||/||P_c b||`；末步解 hygiene 投影
    - `LinearOperator` 注入（dimension/apply）：生产 hPsi 适配器留 C7；金属/dmu 分支留 C4
    - 边界行为：b 全在占据子空间 / b=0 / 维数不匹配 → dpsi=0、residual=0、返回 0 次迭代
    - 测试 5 项全过（MPI 侧 `MODULE_DFPT_stern_test`）：对角算子 vs 闭式补空间解、稠密 Hermitian（Givens+相位酉 U，eps=1.7 落占据带内）vs 谱展开参考、解对随机占据集正交性 <1e-9、占据子空间退化 RHS、零 RHS
    - 9 目标回归全过（CELL 4 + DFPT 5）；`abacus_pw_para` 链接通过；治理仅既有两类豁免 WARNING（头文件值类型 include、设计期模块 docs-sync）
- [x] C3 DFPT_Rho
    - `compute_drho`：每 (q,k) 经 `DFPT_KQ_Basis` 重建 k+q 基，dpsi 系数经 (ix,iy,iz) 反查散布到 rho 网格（C1 模式）；`u`=K 基 recip2real、`du`=rho 网格 recip2real，同网格共轭积累加 `A(r)`；real2recip → `drho_g`（q 移位系数）；Δ=−q（Miller 逆解 + 舍入判定）投影零；`drho_r` 从投影后系数重建（双存储一致）；占据门 `wg<1e-8` 跳过
    - `mix_drho`：`Plain_Mixing::plain_mix` 复空间混合（首步 in=0 → mixed=β·out，残差=1），混合后重建 `drho_r`；`init` 增加 `recip_matrix`（G 矩阵，q_frac→cart），非 plain 混合 WARNING_QUIT；nspin≠1 WARNING_QUIT（自旋- k 排序未钉死，C7 定）
    - 数据层 `set/get_drho_r/g` 由桩转正式存储；irrep 包装测试同步翻转（round-trip 非空）
    - 测试捕获并修复测试侧 2 处参考错误（生产代码无 bug）：① `PW_Basis_K::gcar` 是逐 k 数组（`ik*npwk_max+igl`，pw_basis_k.cpp:261-286），按基球 ig 读是错的——参考列表改按 igl 直读；② 直接求和参考混用 cart G 与 frac r（相位差 lat0 倍）——改 `r_cart=frac·latvec` 后 `g·r_cart`
    - 串行测试 5 项全过（`MODULE_DFPT_rho_serial`）：G 空间 vs 暴力双和（<1e-10）、实空间 vs 直接求和（5 采样点 <1e-9）、Γ 电荷守恒（ig0 置零 + Σdrho_r/|max|/N <1e-12）、混合首步=β·out 且残差=1、第二步组合公式 + 残差
    - 10 目标回归全过（CELL 4 + DFPT 6）；`abacus_pw_para` 链接通过；治理仅既有豁免 WARNING（头文件净减 charge_mixing.h）
- [x] C4 DFPT_Metal（仅接口）
    - `dfdeps`/`compute_dmu`/`compute_drho_metal` 加 WARNING_QUIT 守卫（"not supported in the design phase"），设计期金属分支显式拒绝而非静默错值；`sigma_`/`smearing_type_` 与数据层 `is_metal_`/`dmu_` 槽位保留
- [x] C5 DFPT_Phon
    - `ion_ion`：G 空间（Poisson 对偶恒等式，`w=G+q` 核 `w_a w_b/w²·e^{-w²/4α}`）+ 实空间（erfc Hessian 双循环，`r_c=6/√α`）+ 自项相位差；对角元 phase-free 交叉原子累积（`-√(Mb/Ma)` 系数）+ 自镜像 `(e^{i2πq·L}−1)` 项；α 选取复用 `cal_force_ew` 惯例（1.1×0.9^n，upperbound<1e-6）；Γ ASR 由构造精确成立
    - `accumulate_electron`：2n+1 复数累积 `2Σwg⟨dψ^b|dV^a_ext|ψ⟩`（cross 项经 `apply_dv` 复用 C1 全部约定）+ 同原子非谐项 `Σwg⟨ψ|d²V_loc+d²V_nl|ψ⟩`（`d2vloc_r` rho 网格核 + `apply_d2vnl` 四项 β 恒等式，均已在 C1/C5 实现并测试）；dpsi 槽备份/恢复（apply_dv 复用槽位）
    - `diagonalize`：`LapackConnector::zheev`，`ω=sgn(e)√|e|` 换算 cm⁻¹（独立 CODATA 常数交叉验证）；`add_loto` `(4πe²/Ω)(q̂Z*_a)(q̂Z*_b)/(q̂ε∞q̂)/√(MM′)`；`check_sum_rule` Γ 行和
    - 测试捕获并修复 2 处错误：① 生产 cross 项 `dot.real()` 丢虚部——q≠0 时单 k 矩阵元复数（虚部 k-star 配对相消），assemble Hermitian 对称化依赖复数项，改复数累积；② 测试期望动量缺 `+q`（用 `gpluskq` 直接当动量）——dV 系数动量是 `Δ+q`（与 C1 pert 测试 `AnalyticDVloc(gpp+q_cart)` 一致），手算数值双向定位后修正为 `w=g+q_cart`，d2 期望同步补 `wg` 占据因子
    - 串行测试 `MODULE_DFPT_phon_serial` 7 项全过：Γ ASR（双原子破对称胞）、Γ 声学 3 零模、非公度 q vs 朴素偶极 Hessian 双和、accumulate_electron vs 注入 dpsi 闭式收缩、zheev vs 已知矩阵、loto 各向同性解析极限、Γ 求和规则
    - 11 目标回归全过（CELL 4 + DFPT 7）；`abacus_pw_para` 链接通过；治理仅既有两类豁免 WARNING
- [x] C6 DFPT_Q0
    - `v_hartree_q`（DFPT_Rho 成员）：`dV_H(G)=e²·4π/(tpiba²·|G+q|²)·drho_g`（w=gcar+q_cart，1/lat0 单位），跳过 |G+q|=0（对齐 `h_hartree_pw.cpp` 跳 ig_gge0 惯例）；同函数服务 C7 全 q 屏蔽势
    - `XC_First_Order` 抽象契约（`apply(drho_r, dvxc_r)`，module_dfpt 不 include pot_xc_fdm.h，镜像 Stern::LinearOperator 注入惯例）；PotXC_FDM 适配器（复 δρ 拆 Re/Im）落 C7 esolver 层，XC 核数值对照随 C7 适配器一并测试
    - `build_vkb_dk`（C1 build_vkb 的 k 导数，转 public 供 Q0 复用）：三链解析导数——原子相位 `i2πτ_dir`、径向 `vq'(g)·tpiba·ghat_dir`（radial_vq 中心差分 dg=1e-4）、实谐函数方向链 `(e_dir−ghat·ghat_dir)/|G|`（l≤2 `grad_real_ylm`）；G=0 处 l≥1 行方向链奇异（测度为零，仅相位项，与 QE 同处理）
    - `pos_matrix` 速度算符形式：`⟨u_m|r_d|u_n⟩=−i·⟨u_m|dH/dk_d|u_n⟩/(tpiba·(ε_m−ε_n))`（[H,r]=−i·dH/dk，k 取 2π/lat0 无量纲导数与 build_vkb_dk 一致，r 出 bohr）；dH/dk = 动能 `2tpiba²(k+G)_d` + 非局域 `|dvkb⟩D⟨vkb|+|vkb⟩D⟨dvkb|`（D=dion·m 选择规则，dVnl_dtau 布局）；V_loc 与 k 无关；严格简并对跳过（规范依赖）
    - `compute_eps`：`ε_ab=δ_ab+(8π/Ω)Σ_k wg Re[r_a r_b]/(ε_c−ε_v)/Nk`（长度规范分母，与振子强度和规则一致；绝对值标定 C7 金刚石端到端）；`compute_born`：`Z*_{k,ab}=Z_k δ_ab−(4/Nk)Σ wg Re[⟨v|dV_b|m⟩⟨m|r_a|v⟩]/(ε_m−ε_v)`（m 跑全部带含占据；`⟨v|dV|m⟩=conj(dv_mv)`；经 C1 apply_dv@q=0 复用全部约定，dpsi 槽备份/恢复仿 phon 模式；离子 Z 只加 (a==b) 对角，每原子单次 set_born）
    - 串行测试 `MODULE_DFPT_q0_serial` 5 项全过：build_vkb_dk vs build_vkb 中心差分（泛型 gk 列表，1e-5）、pos_matrix 动能项闭式（−i 因子/tpiba 标定/Hermitian/简对跳过）、非局域收缩 vs 算符有限差分（ψ(G=0) 列置零避开奇点）、compute_eps 二能级全系数链（复激发态敏感于 conj 位置）、compute_born vs 闭式 G 求和（含离子对角+dpsi 恢复）
    - `MODULE_DFPT_rho_serial` 增 v_hartree_q 3 检查（单 G 闭式、|G+q|=0 跳过、尺寸守卫清空），6 项全过
    - 12 目标回归全过（CELL 4 + DFPT 8）；`abacus_pw_para` 链接通过；治理仅既有豁免 WARNING（docs-sync）
- [x] C7 run() 接线 + ESolver/INPUT（模块层 + esolver 工厂接线完成；金刚石端到端数值对照随 B 前置验证补做）
    - 模块层（C7a）：
      - `DFPT_PW::init` 新签名 `(ucell, psi, pw_rho, pw_wfc, sf, veff_r, wg, eig, xc, nelec, ecutwfc, dftu)`（规则 5：不加默认参，全调用点更新，含 pw_run_test 骨架模式传空基）；Impl 持 GS 基/veff/wg/eig + `DFPT_HamiltShift* hamilt_` + `occ_kq_` 缓存
      - `DFPT_HamiltShift : DFPT_Stern::LinearOperator`（新文件 `dfpt_hamilt_shift.{h,cpp}`）：H(k+q) 不复用 GS HamiltPW 链（ik 索引绑定 gk2/vkb，不可平移）→ 自组装三部分——动能 `tpiba²·kq.get_gk2(igl)` 对角 + veff_r FFT 卷积（kq2rho_ 经 FFT-cell triple 映射，C1 惯例）+ 缓存 k+q vkb 的分离非局域（dion m 选择规则同 dVnl_dtau 布局）；`set_context(q_idx,k_idx)` 缓存投影 / `set_shift(eps)` 每 solve 更新对角
      - `DFPT_Pert::apply_vr`（public）：屏蔽响应势作用全带（v_sc_r 与 dv_rc 同约定：q 移位复周期振幅）；`real_space_dv` 重构为委托私有 `apply_vr_core`（FFT-cell triple 散射/收集核心共用）；`build_vkb/build_vkb_dk` 保持 public（Q0 复用）
      - `DFPT_Rho::reset_mixing(q_idx)`：清 drho_in_/residual_，每位移重开 SCF
      - `build_occ_kq(q_idx)`：k+q 折叠匹配 GS k 列表（`kq ≡ k' (mod G)` 容差 1e-8；不匹配 WARNING_QUIT，需 Monkhorst 网格）；占据态经共享 FFT-cell triple 从 ikq 的 G 球映射到 k+q 列表
      - `solve_displacement(q_idx,iat,idir)` 完整位移级 SCF 内环：`v_hartree_q(drho_g) + xc_->apply(drho_r)` 组 v_sc_r → `apply_dv + apply_vr` 组 RHS → `set_shift + stern_.solve` 每占据带 → `compute_drho + mix_drho` 残差收敛判据
      - `run()`：q=0 时 q0 响应（eps/Born/loto）；每 irrep 位移循环 + `accumulate_electron`；`assemble + diagonalize + add_loto`（loto 方向默认 (1,1,1)/√3，一般方向随 A 阶段 irrep 机制）；null 基保持骨架首迭代收敛退化（测试兼容）
    - esolver 层（C7b）：
      - `esolver_dfpt_pw.{h,cpp}` 重写：`before_all_runners` 只做静态配置 + 从 `inp` 捕获 nspin/nelec/ecutwfc/dft_plus_u（规则 1：显式传递，init_dfpt 不读全局记录）；`runner` 先 `run_gs`（复用 `ESolver_KS_PW::runner`）→ `init_dfpt` 真接线（GS 收敛后 veff/charge/psi 才存在）→ `dfpt_->run()`
      - `init_dfpt` 实参：`*this->stp.psi_cpu`、`this->pw_rho/pw_wfc`、`&this->sf`、`get_veff_smooth()` 行 0 展开（`update_from_charge` 每迭代调 `interpolate_vrs`，收敛后即当前值）、`pelec->wg/ekb`、`XC_First_Order_FDM` 适配器（Re/Im 拆分过 `PotXC_FDM` 有限差分核，线性重组精确到 O(|δρ|²)；持 GS Charge + scratch Charge）、`dft_plus_u ? &this->dftu : nullptr`；守卫：nspin≠1 / charge 不在 rho 网格（USPP）/ veff_smooth 网格不匹配 → WARNING_QUIT
      - `esolver.cpp` 工厂：`determine_type` pw 分支加 `"dfpt"→"dfpt_pw"` + `init_esolver` 分支（治理豁免：determine_type 既有 PARAM 读取惯例，1 行）
      - `read_inp_sys.cpp`：esolver_types 合法值加 `"dfpt"` + 注释/description 更新；`docs/parameters.yaml` + `docs/advanced/input_files/input-main.md` 同步
    - 验证：`cmake --build` esolver/abacus_pw_para/12 测试目标全绿；ctest 12/12（CELL 4 + DFPT 8）；`abacus_pw_para -h esolver_type` 显示 dfpt 条目；`--version` v3.11.0-beta8；治理仅 determine_type 工厂 1 处豁免 ERROR + 既有 header/docs WARNING
    - 待办（随 B/前置验证）：金刚石端到端声子/ε∞/Z* 对照、`--check-input` 从有效算例目录验证
- [x] 校准：屏蔽通道三处修复（金刚石 2 原子 smoke，24³ rho 网格，NC PP，Γ 点）
    - FD/Ewald 锁定基线：e11=+0.08056、e12=−0.08059 Ry/bohr²（预质量 0.0028685/−0.0028701），光学 ~742 cm⁻¹；GS 力 FD 交叉验证 dV 装配（F_x 偏差 0.02%）
    - 修复 1（dfpt_rho.cpp compute_drho）：q=0 Hermitian 完成的 in-place `drho_g[ig] += conj(drho_g[gm])` 逐点双重处理 ±G 对（第二次访问读到已更新的第一项）→ 结果破坏 Hermitian 性，实空间重建混入 Re a(r) 寄生分量（均匀 ~1.25 过冲、对称违反被放大 1.3-1.8%、A1 投影 3.7%）；最终实现 = 实空间 `2 Re a(r)` 预对称化后再 real2recip（实数组 FFT 本征 Hermitian，单边 stick（−G 不在球内）亦获正确完成值，替代 G 空间逐点镜像）
    - 修复 2（esolver_dfpt_pw.cpp XC_First_Order_FDM）：前向差分 `Vxc[ρ+δρ]−Vxc[ρ]` 的曲率项 ½Vxc″δρ²（T2⊗T2⊃A1）向 v_sc 泄漏寄生 A1（band0 ⟨dv_sc⟩=+0.0173 违反 A1⊗T2⊗A1 选择定则、占据三重态迹 +0.052）且二次非线性反馈使混合迭代 β=0.7 超指数暴走；改 η=1e-6 中心差分（Re/Im 各一对 cal_v_eff 探测）后泄漏 ~1e-11，默认 β=0.7 恢复收敛
    - 修复 3（dfpt_pw.cpp solve_displacement）：`reset_mixing` 只清混合器内部态，data 层 `drho_g` 残留上一位移响应（含发散残渣）泄漏进新位移首迭代 v_sc；进入位移时同步清零
    - 修复后（默认 β=0.7，~76 s）：光学 742.367×3（FD ~742）、声学 6.40×3（ASR：e11+e12=3.1e-6）、e11=0.00286804（目标 0.0028685）、e12=−0.00286494（目标 −0.0028701）、非 irrep 元 ~1e-11、收敛 drho 小群违反 0.000000/A1 投影 5e-6（对称性精确）；裸响应（β=0.001 dump）小群违反 ~0.1% 确认裸链（Sternheimer/dV/dψ）干净
    - 后期漂移根因（本轮确诊）：残差降至 5e-5 后指数增长（1.27×/iter）、|in| 恒定而 out 偏离 → 垃圾方向与物理分量正交、混合映射本征值 μ=1.2765 恒定（纯本征模）；本征模身份 = {200} 壳 6 矢等幅实系数 + {111} 壳 8 矢 ±π/4 相位的 Hermitian 实 A1 呼吸模（seed ~1e-6 舍入级）；均匀探针实验（DFPT_JPROBE：注入纯 A1 模 + rhs 去 dV_ext 单迭代直测线性映射）给出 λ_A1 = −2.229（Hartree-only −3.180，XC 削减到 −2.23）——非符号 bug，是最小 G 壳的 Coulomb 刚性（4π/G² 硬核）：plain mixing 收敛条件 −2/β+1<λ 要求 β<0.62，物理 T2 模 λ=−1.42（小 G 头部含量少）在 β=0.7 恰好可收敛，故固定点正确而 A1 通道发散；β=0.4 时 μ_A1=−0.29 稳定
    - 修复：默认 mix_beta 0.7→0.4（注释记录测得的 |λ|~2.2 与 β 上界 2/(1+|λ_min|)，留裕量至 |λ|~5）；DFPT_MIX_BETA env 旋钮保留；β=0.4 时 6 位移全部经收敛旗标退出（平均 ~38 iter，总 228），频率/ele 矩阵与 β 无关逐位一致（固定点正确性再验证），收敛 drho manifest 干净（|FD| 比率 0.99994、cos 0.9993、逐点相对差 3.8%）；后续正解是 Kerker 型预条件混合（随 B 阶段排期）
    - ε∞/Z* 打印为空（随 B 阶段）；调试插桩（PTCHK/DYNCHK/MDBG/JPROBE dump/VKBCHK/drho dump/DFPT_MIX_BETA env）收尾节点统一清理评审
- [ ] B 工程化收编（2026-08-18 修订：B1 INPUT 接线 → B0 全流程验证 → B2 输出正式化 → B3 Kerker → B4 数据层）
    - 修订依据的差距盘点（代码 vs 计划交叉核对，2026-08-18）：
      ① `set_compute_q0`/`set_loto` 全仓库无调用者（死代码，q0/loto 分支不可达）；
      ② `set_parameters("dfpt.in")` 空桩 + esolver 硬编码 `set_qmesh(1,1,1)`/conv_thr/max_iter（非 Γ q 无法从输入驱动）；
      ③ ε∞/Z* 打印为 design-phase 临时 std::cout（esolver_dfpt_pw.cpp:322，单 rank 假定）；
      ④ 混合仍 plain β=0.4（Kerker 预条件排期 B3）；
      ⑤ DFPT_IrrepData irrep 维度占位穿透、get_dpsi_obj 返回 static dummy；
      ⑥ run() 外层 while 形式化（无条件 set_converged(true)、LO-TO 方向硬编码 (1,1,1)/√3；
      ⑦ 真实晶格金刚石端到端/非 Γ q/MPI>1 rank 均未冒烟（C7 待办未做）。
    - [x] B1 INPUT 参数接线 `297a2b3a2`
        - `read_inp_dfpt.cpp` 7 项（qmesh/qfile/compute_q0/loto/conv_thr/max_iter/mix_beta，
          check_value：loto 需 compute_q0、qmesh≥1、mix_beta∈(0,1]）；CMake 两处接入
        - esolver 删硬编码与 `dfpt.in` 空桩，`set_parameters` 移除，全走显式 setter（规则 1）；
          `DFPT_PW` 新增 set_qfile/set_mix_beta/set_compute_q0/set_loto（q 文件优先于 MP 网格）
        - `QList::read_from_file` 改填占位 A1 irrep（nirr=1）而非清空——q 文件路径下 run()
          依赖 get_nirr≥1 才进 3N 回退求解
        - 修复（测试捕获）：空默认串参数回写 INPUT 再读回时 `str_values[0]` 越界 → dfpt_qfile
          采用 pseudo_dir 的 `get_size()==0` 守卫模式
        - docs/parameters.yaml 新类目 + input-main.md 经 docs/generate_input_main.py 再生；
          `-h dfpt_qmesh/dfpt_mix_beta` 验证；回归 14/14（CELL 4 + DFPT 8 + IO 2）；
          治理仅既有豁免 WARNING
    - [ ] B0 全流程工程验证（真实金刚石 Γ / 非 Γ q / MPI 冒烟）
        - [x] 多 k（nk>1）数值错误根因与修复 `（本轮，6 文件 +374/−57）`
            - 排除法完成：sym=0/1、v_sc 假设、LAPACK/BLAS、npw 不匹配、球大小不匹配（{L 392, X
              388} 完全正确）、δρ-cube FD 路线（k 采样噪声地板 >> 位移信号，判死）、d2/计算
              通道（d2 混合=孤立预测精确一致）
            - **根因 1（标签折叠）**：build_occ_kq 假设 k+q 球与 k(ikq) 球共享同 FFT 胞 G 标签；
              {L,−L} 时 −L 折叠到 L 标签（差 b1），标签错位 → 投影态垃圾 → 核检查失败+发散。
              修复 = 倒格矢整数三元组匹配 f+dn=f'（dn=k(ik)+q−k(ikq)），ikq 侧标签经
              `PW_Basis_K::getgcar` 读取——关键发现：`collect_local_pw(erf)` 把 gcar 重建为
              per-k 球布局 [ik*npwk_max+igl]，父类全局 ig 布局已毁（nk=1 曾靠堆残留"幸运"通过）
            - **根因 2（smearing 投影悬崖）**：wg<1e-8 绝对阈使投影器随 k 采样跳变——{Γ,L} 采样
              的 E_f 使 L 带占据带尾 w=5e-6 跨过阈值，进入 P_c 投影 → 其空态通道在 (H−ε)⁻¹
              中关闭 → 收敛响应差 ~10%（X03 +46%、ASR 行和违反 21%）。权重 1% Γ 实验证实损伤
              与 w_Γ 无关（结构性）；β=0 单迭代实验证实同 rhs/本征值下 |dψ| 差 8%（纯投影效应）。
              修复 = 共享 `dfpt_band_occupied()`：wg(ik,ib) > 0.5·wg(ik,0)（多数占据判据），
              一致应用于投影器/求解驱动/drho/2n+1 装配/q0 v-c 划分
            - FD 验证矩阵（sym=0 金刚石 2 原子，FD 模板 b0_si_k050_fd2/run_fd.py 派生）：
              单 Γ D00 0.0208553 vs FD 0.020854；单 L 0.0129282 vs FD 0.012927（**新增 FD
              基准** b0_si_kLL1_fd）；{L,−L} = 单 L 逐位一致（原发散）；{Γ,L} 0.0166416 vs
              FD 0.016642（原 0.0182462，+9.6%）；{L,X}（两不等价非 Γ 点）与权重偏斜
              {Γ,L} 变体全部自洽；ASR 行和全部 ~1e-6
            - 调试方法学沉淀：XB per-(ik,ib) 分解仅 iter-1（v_sc=0）可比但受 GS 采样差异混淆；
              ASR 行和 = 免 FD 的在跑检测器；跨采样对比仅 {L}vs{L,−L}（BZ 等价）合法
            - 回归 14/14（CELL 4 + DFPT 8 + IO 2）；治理仅既有 header/docs WARNING
            - 遗留：调试插桩（OCCCHK miss/集合计数、PTCHK/DYNCHK2/4/XB/MDBG/JPROBE）保留至
              B0 收尾统一清理（用户决策）；真正的金属分数占据 DFPT（de Gironcoli 成对方程）
              超出当前绝缘体范围，dfpt_metal 占位
        - [x] 验证梯队扩展（单 k 0.25 / k 网格 / 金属占据区界 / ε∞Z* / MPI 冒烟）`（本轮）`
            - **新增 FD 通过项**：单 k=(0.25,0,0)（Λ 点，392 球）D 行 0 实部 vs FD 全部 ~6e-7
              （0.0104952/0.00304483/−0.0104942 vs 0.010495/0.003044/−0.010495）；单边 k 采样的
              Hermitian 虚部（D(0,3) imag 0.0061）= X_ba≠X_ab 的预期产物，物理力常数取实部
              （FD 证实）；2×2×2 非移位网格（σ=0.005 绝缘区）D00 0.0127458 vs FD 0.012739
              （0.05%），非对角/ASR 精确
            - **金属占据区界定量**（2×2×2 网格，Γ VBM 带尾）：σ=0.015 → VBM 92% 占据，FD 比
              DFPT 软 2.83×（FD 含 dμ/dτ 响应、Sternheimer 流无此通道）；σ=0.007 → 99.92%，
              差 3.8%；σ=0.005 → 99.9996%，差 0.05%——残差随带尾权重缩小，dμ 通道缺失的干净
              指纹。**守卫**：DFPT_PW::init 扫描最终 wg，任一带相对占据落入 (1e-3, 1−1e-3) 即
              WARNING_QUIT（显式拒绝而非静默错值，与 C4 哲学一致；k222 σ=0.015 实测触发）
            - **symmetry=1 单 k 陷阱**：k050(0.5,0,0) sym=1 与 sym=0 GS 本征值差 2.6e-3 Ry
              （−0.234996 vs −0.237553）→ D 差 0.3-4.4%；sym=0 重跑 = kLL1 逐位一致。FD 基准
              全部 sym=0，跨 sym 比较非法
            - ε∞/Z*：打印链路通；单 Γ 值（105/55.9）= 长度规范简并分母伪影；8-k 值
              （4.75/6.68）= 采样受限；真验证需密网格（8×8×8 ~10h 串行，推迟为过夜项）
            - MPI 冒烟（-np 2，kG）：DFPT 相位 MPI_ERR_TRUNCATE 硬崩溃（分布式布局未支持，
              响亮失败无静默错值）；Stern CG 标量积无 Allreduce 等串行假设已知，MPI 支持另立
              工作项
            - 案例清单（/tmp/opencode/）：k025/k050（已 sym=0 修正）、k222s007/k222s005（+
              _fd 配对）、gamma_k222_nosym（守卫触发样本）、kG_mpi2（崩溃样本）
    - [ ] B0 残余：非 Γ q 文件路径冒烟（dfpt_qfile + QList 端到端从未运行）；8×8×8 过夜验证；
      插桩清理评审
    - [ ] B2 输出正式化
    - [ ] B3 Kerker 预条件混合
    - [ ] B4 数据层收编
- [ ] A irrep 分解（保留接口，工程验证完成后立项）
