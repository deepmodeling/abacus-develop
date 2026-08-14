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
- `assemble`：`ion_ion(q,dynmat)`（已真）+ `electron(q_idx,data,dynmat)`（已真，Ewald 复用 `force_pw.cpp:479 cal_force_ew` + `H_Ewald_pw::rgen`；相因子经 `symm.gtrans[48]`+`kgmatrix[48]`）。
- **U0 实装点**：`dftu_onsite`/`dftu_lambda` 电子项。
- `diagonalize` 换真实现（替换 `freq[i]=i` 伪值）；`add_loto`/`check_sum_rule` 实装或明示桩。
- 测试：金刚石 Γ 声子频率对照。提交。

**C6 — DFPT_Q0 介电/Born/LO-TO**
- 新增 `v_hartree_q`：`|G+q|²` Poisson 因子、跳过 ig=−q（参考 `h_hartree_pw.cpp:16-97`）。
- XC 一阶核：库中无 v_xc 一阶 API，用 LIBXC kernel 或有限差兜底。
- 非局域 `[r,V_U]` commutator 项记录为 U 预留。
- 测试：金刚石介电张量/Born 电荷对照。提交。

**C7 — run() 接线 + ESolver/INPUT**
- 修正签名不一致：注释中 `pert_.build_dv(q,irrep,...)` 与真实 `build_dv(int q_idx,int atom_idx,int dir,DFPT_PW_Data&)`；mode basis 为空时逐 irrep 先遍历 3N 方向，irrep 收敛为代表模。
- `esolver_dfpt_pw.cpp`：解开 `init` 注释，实参 `*this->stp.psi_cpu` + `PARAM.inp.nelec` + `PARAM.inp.ecutwfc`；`dft_plus_u` 为真时传 `&this->dftu`（否则 nullptr）。
- INPUT 行为若变则同步 `docs/parameters.yaml` + `input-main.md`。
- 金刚石端到端对照：声子频率 + 介电；`./build/abacus --version` 记录身份。
- 全量构建 + 回归 + 治理。提交。

---

## B — 数据层收编

- `DFPT_IrrepData` 下沉为 `DFPT_PW_Data` 正式数据层（收敛台账、irrep 元数据并入），迁移现有 dfpt 测试。
- 提交。

## A — irrep 分解（最后）

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
- [ ] C3 DFPT_Rho
- [x] C3 DFPT_Rho
    - `compute_drho`：每 (q,k) 经 `DFPT_KQ_Basis` 重建 k+q 基，dpsi 系数经 (ix,iy,iz) 反查散布到 rho 网格（C1 模式）；`u`=K 基 recip2real、`du`=rho 网格 recip2real，同网格共轭积累加 `A(r)`；real2recip → `drho_g`（q 移位系数）；Δ=−q（Miller 逆解 + 舍入判定）投影零；`drho_r` 从投影后系数重建（双存储一致）；占据门 `wg<1e-8` 跳过
    - `mix_drho`：`Plain_Mixing::plain_mix` 复空间混合（首步 in=0 → mixed=β·out，残差=1），混合后重建 `drho_r`；`init` 增加 `recip_matrix`（G 矩阵，q_frac→cart），非 plain 混合 WARNING_QUIT；nspin≠1 WARNING_QUIT（自旋- k 排序未钉死，C7 定）
    - 数据层 `set/get_drho_r/g` 由桩转正式存储；irrep 包装测试同步翻转（round-trip 非空）
    - 测试捕获并修复测试侧 2 处参考错误（生产代码无 bug）：① `PW_Basis_K::gcar` 是逐 k 数组（`ik*npwk_max+igl`，pw_basis_k.cpp:261-286），按基球 ig 读是错的——参考列表改按 igl 直读；② 直接求和参考混用 cart G 与 frac r（相位差 lat0 倍）——改 `r_cart=frac·latvec` 后 `g·r_cart`
    - 串行测试 5 项全过（`MODULE_DFPT_rho_serial`）：G 空间 vs 暴力双和（<1e-10）、实空间 vs 直接求和（5 采样点 <1e-9）、Γ 电荷守恒（ig0 置零 + Σdrho_r/|max|/N <1e-12）、混合首步=β·out 且残差=1、第二步组合公式 + 残差
    - 10 目标回归全过（CELL 4 + DFPT 6）；`abacus_pw_para` 链接通过；治理仅既有豁免 WARNING（头文件净减 charge_mixing.h）
- [ ] C4 DFPT_Metal（仅接口）
- [ ] C5 DFPT_Phon
- [ ] C6 DFPT_Q0
- [ ] C7 run() 接线 + ESolver/INPUT + 金刚石对照
- [ ] B 数据层收编
- [ ] A irrep 分解
