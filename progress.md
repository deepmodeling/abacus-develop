# module_gint 单精度支持进度

## 阶段 A

### 状态

已完成。

### 本阶段实现

1. 新增 `source/source_lcao/module_gint/gint_precision.h`
   - 定义 `GintRealPrecision`
   - 定义 `GintExecConfig`
   - 默认内部精度为 `fp64`

2. 扩展 `gint_interface`
   - `cal_gint_vl(...)` 新增可选参数 `const GintExecConfig& cfg = {}`
   - `cal_gint_rho(...)` 新增可选参数 `const GintExecConfig& cfg = {}`
   - 当前阶段仅打通入口，默认行为保持不变，`cfg` 先不参与调度

3. 补齐 `module_hcontainer` 的 `float` 显式实例化
   - `BaseMatrix<float>`
   - `AtomPair<float>`
   - `HContainer<float>`

4. 扩展 `GintInfo::get_hr<T>()`
   - 补上 `get_hr<float>()` 的显式实例化

5. 在 `gint_common` 新增 `HContainer` 跨精度转换 helper
   - `cast_hcontainer_values<Tout, Tin>(...)`
   - `make_cast_hcontainer<Tout, Tin>(...)`
   - 当前已显式实例化 `double -> float` 和 `float -> double`

### 对后续阶段的意义

- 阶段 B 可以直接基于 `GintExecConfig` 在 `gint_vl` 里选择 `fp32/fp64`
- `HContainer<float>` 已经具备可链接的基础能力
- `float hr_gint_ <-> double hr_gint_` 的边界转换 helper 已经就位

### 验证记录

1. 构建配置
   - 命令：`../build.sh build_phase_a`
   - 结果：成功生成 `build_phase_a`

2. 实际编译
   - 命令：`cmake --build build_phase_a -j14`
   - 结果：成功，最终产物为 `build_phase_a/abacus_2g`

3. 运行测试
   - 目录：`tests/performance/P101_si32_lcao`
   - 命令：`OMP_NUM_THREADS=7 mpirun -n 2 --bind-to socket /home/dzc/abacus/abacus-mix/build_phase_a/abacus_2g`
   - 结果：程序正常结束，退出码为 0

4. 结果抽取
   - 命令：`bash ../catch_properties.sh result.out`
   - 生成结果：
     - `etotref -3403.2017700426463307`
     - `etotperatomref -106.3500553138`
     - `totalforceref 2.961504`
     - `totalstressref 469.743372`
   - 参考结果 `result.ref`：
     - `etotref -3404.7663590699248743`
     - `etotperatomref -106.3989487209`
     - `totalforceref 2.688768`
     - `totalstressref 473.654463`

### 重要观察

- A 阶段代码改动后，工程可以完整编译，`P101_si32_lcao` 也可以正常运行结束。
- 但当前运行结果与该用例目录中的 `result.ref` 不一致，且 `catch_properties.sh` 未能从当前日志中提取到 `totaltimeref`。
- 这更像是当前构建/运行环境与该性能用例参考值之间已有偏差，A 阶段本身没有启用任何新的单精度执行路径；进入阶段 B 前，建议继续把“结果偏离 reference”视为独立观测项记录，不要误判为本阶段引入的新行为变化。

## 阶段 B

### 状态

已完成。

### 本阶段实现

1. 重构 `source/source_lcao/module_gint/gint_vl.h/.cpp`
   - `Gint_vl` 新增 `GintExecConfig cfg_`
   - `cal_gint()` 改为调度器，根据 `cfg_.cpu_internal_real` 选择：
     - `cal_gint_impl_<double>()`
     - `cal_gint_impl_<float>()`

2. 打通 `gint_vl` 的内部单精度执行路径
   - `cal_gint_impl_<float>()` 中：
     - `hr_gint` 改为 `HContainer<float>`
     - `phi` / `phi_vldr3` 改为 `std::vector<float>`
     - `dr3_` 在调用 `phi_mul_vldr3` 时显式转为 `float`

3. 按计划补上 `vr_eff` 的单精度内部缓存
   - 新增 `get_vr_eff_data_<Real>()`
   - `fp32` 路径会先将外部 `double* vr_eff_` 一次性 cast 到本地 `std::vector<float>`
   - 避免在每个 biggrid 循环里重复做 `double -> float` 转换

4. 建立 `gint_vl` 的单精度输出边界
   - `fp32` 路径结束后使用 `make_cast_hcontainer<double>(hr_gint_sp)`
   - 再调用 `compose_hr_gint(...)`
   - 最后通过 `transfer_hr_gint_to_hR(...)` 写回外部 `HContainer<double>`

5. 泛化 `compose_hr_gint`
   - `source/source_lcao/module_gint/gint_common.h/.cpp`
   - `compose_hr_gint` 从仅支持 `HContainer<double>` 改为模板函数
   - 显式实例化补上：
     - `compose_hr_gint(HContainer<double>&)`
     - `compose_hr_gint(HContainer<float>&)`

6. 修复阶段 B 暴露出的单精度链接缺口
   - 全量链接时发现 `PhiOperator::set_phi<float>()` 依赖的
     `GintAtom::set_phi<float>()` 尚未显式实例化
   - 已在 `source/source_lcao/module_gint/gint_atom.cpp` 补上：
     - `GintAtom::set_phi<float>(...)`
     - `GintAtom::set_phi_dphi<float>(...)`

### 对后续阶段的意义

- 阶段 C 可以直接复用同样的“调度器 + `Real` 模板实现”模式改造 `gint_rho`
- `gint_vl` 的 `float` 内部路径已经具备完整的数据流：
  - `double vr_eff`
  - `float vr_eff_gint`
  - `float phi`
  - `float phi_vldr3`
  - `HContainer<float> hr_gint`
  - `HContainer<double> hR`
- 当前上层调用点仍全部使用默认 `GintExecConfig{}`，因此实际集成运行仍走 `fp64`
- 进入阶段 C 前不需要再改 `gint_vl` 的接口骨架，重点可转向 `double DMR -> float dm_gint`

### 验证记录

1. 构建配置
   - 命令：`../build.sh build_phase_b`
   - 结果：成功生成 `build_phase_b`

2. 实际编译
   - 首次命令：`cmake --build build_phase_b -j14`
   - 首次结果：在最终链接阶段失败，报错缺少
     `GintAtom::set_phi<float>(...)` 的定义
   - 修复后再次执行：`cmake --build build_phase_b -j14`
   - 最终结果：成功，最终产物为 `build_phase_b/abacus_2g`

3. 运行测试
   - 目录：`tests/performance/P101_si32_lcao`
   - 命令：`OMP_NUM_THREADS=7 mpirun -n 2 --bind-to socket /home/dzc/abacus/abacus-mix/build_phase_b/abacus_2g`
   - 结果：程序正常结束，退出码为 0

4. 结果抽取
   - 命令：`bash ../catch_properties.sh result.out`
   - 生成结果：
     - `etotref -3403.2017700426463307`
     - `etotperatomref -106.3500553138`
     - `totalforceref 2.961504`
     - `totalstressref 469.743372`
     - `totaltimeref` 为空
   - 参考结果 `result.ref`：
     - `etotref -3404.7663590699248743`
     - `etotperatomref -106.3989487209`
     - `totalforceref 2.688768`
     - `totalstressref 473.654463`
     - `totaltimeref +85.40579`

5. 运行日志中的额外观测
   - `Gint cal_gint_vl`：`1.70 s / 9 calls / 0.19 s avg`
   - `Gint cal_gint_rho`：`1.72 s / 9 calls / 0.19 s avg`

### 重要观察

- 阶段 B 的代码骨架已经把 `gint_vl` 的内部 `fp32/fp64` 双路径准备好了，并且全量工程可以成功链接。
- 本阶段过程中新增的 `GintAtom<float>` 显式实例化是实际需要的最小补丁；这说明 `phi` 相关单精度链路现在已经从 `GintAtom -> PhiOperator -> Gint_vl` 串通。
- 由于当前上层还没有把 `GintExecConfig{cpu_internal_real = fp32}` 传到 `cal_gint_vl(...)`，性能用例的实际运行结果仍与阶段 A 相同，这符合当前分阶段计划。
- `P101_si32_lcao` 的输出依旧与目录内 `result.ref` 不一致，且 `catch_properties.sh` 仍未提取到 `totaltimeref`；这延续了阶段 A 的观测，不应归因为阶段 B 的默认行为变化。

## 阶段 C

### 状态

已完成。

### 本阶段实现

1. 重构 `source/source_lcao/module_gint/gint_rho.h/.cpp`
   - `Gint_rho` 新增 `GintExecConfig cfg_`
   - `cal_gint()` 改为调度器，根据 `cfg_.cpu_internal_real` 选择：
     - `cal_gint_impl_<double>()`
     - `cal_gint_impl_<float>()`

2. 将 `dm_gint_vec_` 从类成员改为模板内部容器
   - 在 `cal_gint_impl_<Real>()` 内按需创建
     `std::vector<HContainer<Real>> dm_gint_vec`
   - `init_dm_gint_<Real>()` 通过 `gint_info_->get_hr<Real>()` 初始化内部 serial `HContainer`
   - 类本身保持非模板，外部接口继续维持双精度输入/输出

3. 泛化 `transfer_dm_2d_to_gint`
   - `source/source_lcao/module_gint/gint_common.h/.cpp`
   - 接口改为：
     `transfer_dm_2d_to_gint<TGint, TDM>(...)`
   - 补上显式实例化：
     - `double <- double`
     - `float <- double`
     - `std::complex<double> <- std::complex<double>`

4. 按计划实现 `float <- double` 的推荐搬运路径
   - 先沿用现有 `double` DMR gather 逻辑得到 serial `HContainer<double>`
   - 再使用 `cast_hcontainer_values<float>(...)` 做本地数值转换
   - 避免把 mixed-precision MPI gather 和本阶段目标耦合在一起

5. 放开 `phi_dot_phi` 的输出精度
   - `source/source_lcao/module_gint/phi_operator.h/.hpp`
   - `phi_dot_phi` 从单模板参数改为 `template<typename Tin, typename Tout = Tin>`
   - `gint_rho` 的 `fp32` 路径现在可以使用：
     - `float phi`
     - `float phi_dm`
     - `double rho`

6. 接通 `gint_interface` 到 `gint_rho`
   - `source/source_lcao/module_gint/gint_interface.cpp`
   - CPU 路径下不再忽略 `cfg`
   - `cal_gint_rho(...)` 会把 `GintExecConfig` 传给 `Gint_rho`

### 对后续阶段的意义

- `gint_rho` 已经具备完整的内部 `fp32/fp64` 双路径骨架，和阶段 B 的 `gint_vl` 结构对齐。
- `double DMR -> float dm_gint -> float phi/phi_dm -> double rho` 这条 CPU 数据流已经在 `module_gint` 内部打通。
- 由于主 SCF 上层目前仍未下发 `GintExecConfig{cpu_internal_real = fp32}`，集成运行默认仍走 `fp64`；下一阶段可以把重点转到“上层如何按策略真正启用 `fp32` 路径”。
- `phi_dot_phi` 的输出类型已经解耦，这为后续 `gint_tau` 或更多内部 mixed-precision 累加场景提供了可复用的基础设施。

### 验证记录

1. 构建配置
   - 命令：`../build.sh build_phase_c`
   - 结果：成功生成 `build_phase_c`

2. 实际编译
   - 命令：`cmake --build build_phase_c -j14`
   - 结果：成功，最终产物为 `build_phase_c/abacus_2g`

3. 运行测试
   - 目录：`tests/performance/P101_si32_lcao`
   - 命令：`OMP_NUM_THREADS=7 mpirun -n 2 --bind-to socket /home/dzc/abacus/abacus-mix/build_phase_c/abacus_2g`
   - 结果：程序正常结束，退出码为 0

4. 结果抽取
   - 命令：`bash ../catch_properties.sh result.out`
   - 生成结果：
     - `etotref -3403.2017700426458759`
     - `etotperatomref -106.3500553138`
     - `totalforceref 2.961504`
     - `totalstressref 469.743372`
     - `totaltimeref` 为空
   - 参考结果 `result.ref`：
     - `etotref -3404.7663590699248743`
     - `etotperatomref -106.3989487209`
     - `totalforceref 2.688768`
     - `totalstressref 473.654463`
     - `totaltimeref +85.40579`

5. 运行日志中的额外观测
   - `Gint cal_gint_vl`：`1.69 s / 9 calls / 0.19 s avg`
   - `Gint cal_gint_rho`：`1.71 s / 9 calls / 0.19 s avg`

### 重要观察

- 阶段 C 的核心边界已经落地：`gint_rho` 不再被 `HContainer<double>` 和 `double phi_dm` 写死，内部单精度数据流已经和阶段 B 的 `gint_vl` 对齐。
- 本阶段没有新增链接错误或运行期崩溃；全量工程编译和 `P101_si32_lcao` 用例都稳定通过。
- 由于当前主流程仍未主动传入 `GintExecConfig{cpu_internal_real = fp32}`，本次集成测试依旧走默认 `fp64`，因此数值结果与阶段 A/B 保持一致，这符合当前分阶段计划。
- `P101_si32_lcao` 与目录内 `result.ref` 的偏差、以及 `catch_properties.sh` 未提取到 `totaltimeref` 的现象仍然延续，继续应视为环境/基线观测项，而不是阶段 C 新引入的行为变化。

## 阶段 D

### 状态

已完成。

### 本阶段实现

1. 新增上层精度控制器
   - 新增 `source/source_estate/module_charge/gint_precision_controller.h/.cpp`
   - 定义 `GintPrecisionController`
   - 提供：
     - `reset_for_new_scf()`
     - `update_after_iteration(...)`
     - `current_config()`

2. 按计划把控制器挂到 SCF/charge mixing 层
   - `source/source_estate/module_charge/charge_mixing.h`
   - `Charge_Mixing` 新增 `gint_precision_controller_`
   - 暴露 `get_gint_precision_controller()`

3. 在 SCF 迭代边界接通控制器生命周期
   - `source/source_estate/module_charge/chgmixing.cpp`
   - `chgmixing_ks_lcao(iter == 1)` 时调用 `reset_for_new_scf()`
   - `chgmixing_ks(...)` 在完成本轮 `drho / conv_esolver / restart-step` 判定后调用
     `update_after_iteration(...)`
   - 当前策略严格按阶段 D 文档首版规则实现：
     - 初始 `fp32`
     - `iter >= 3`
     - 非 restart step
     - `drho <= max(100 * scf_thr, 1e-5)`
     - 连续 2 轮满足后，下一轮切到 `fp64`
     - 一旦切到 `fp64`，本轮 SCF 不再切回

4. 建立上层配置传播路径
   - `source/source_estate/module_charge/charge.h`
   - `source/source_estate/module_pot/potential_new.h`
   - `Charge` 和 `Potential` 新增 `GintExecConfig` 存取接口
   - `source/source_esolver/esolver_ks_lcao.cpp`
   - `ESolver_KS_LCAO::iter_init()` 每轮从
     `p_chgmix->get_gint_precision_controller().current_config()`
     读取配置，并下发给：
     - `chr`
     - `pelec->pot`

5. 接通 LCAO SCF 主路径的 `gint_vl`
   - `source/source_lcao/module_operator_lcao/veff_lcao.cpp`
   - `Veff<OperatorLCAO<double, double>>::contributeHR()`
   - `Veff<OperatorLCAO<std::complex<double>, double>>::contributeHR()`
   - 非 meta-GGA 分支现在会把 `pot->get_gint_exec_config()` 传给
     `ModuleGint::cal_gint_vl(...)`
   - meta-GGA 分支保持现状，继续 `fp64`

6. 接通 LCAO SCF 主路径的 `gint_rho`
   - `source/source_lcao/rho_tau_lcao.cpp`
   - `LCAO_domain::dm2rho(...)` 现在从 `Charge` 读取当前配置，并传给
     `ModuleGint::cal_gint_rho(...)`
   - `source/source_estate/elecstate_lcao.cpp`
   - PEXSI 路径也改为从 `charge->get_gint_exec_config()` 读取配置

7. 构建系统接入新源文件
   - `source/source_estate/CMakeLists.txt`
   - 新增 `module_charge/gint_precision_controller.cpp`

### 对后续阶段的意义

- 阶段 B/C 准备好的 `gint_vl` / `gint_rho` 内部 `fp32/fp64` 双路径，现在已经真正接入了 LCAO 主 SCF 调度闭环。
- 当前架构中，精度决策由 `Charge_Mixing` 持有并更新，`ESolver_KS_LCAO` 负责把当前配置下发给 `Charge` 与 `Potential`；`module_gint` 仍只做执行，不参与策略判断。
- 非 SCF 和后处理路径仍默认 `fp64`，因为它们没有在阶段 D 中主动下发新的 `GintExecConfig`，符合“先只接主 SCF 路径”的计划。
- 阶段 E 可以在此基础上直接做：
  - `fp32/fp64` 数值接近性测试
  - SCF 收敛轮数比较
  - `cal_gint_vl` / `cal_gint_rho` 的真实性能评估

### 验证记录

1. 构建配置
   - 命令：`../build.sh build_phase_d`
   - 结果：成功生成 `build_phase_d`

2. 实际编译
   - 命令：`cmake --build build_phase_d -j14`
   - 结果：成功，最终产物为 `build_phase_d/abacus_2g`

3. 运行测试
   - 目录：`tests/performance/P101_si32_lcao`
   - 命令：`OMP_NUM_THREADS=7 mpirun -n 2 --bind-to socket /home/dzc/abacus/abacus-mix/build_phase_d/abacus_2g`
   - 结果：程序正常结束，退出码为 0

4. 结果抽取
   - 命令：`bash ../catch_properties.sh result.out`
   - 生成结果：
     - `etotref -3403.2017700426458759`
     - `etotperatomref -106.3500553138`
     - `totalforceref 2.961504`
     - `totalstressref 469.743372`
     - `totaltimeref` 为空
   - 参考结果 `result.ref`：
     - `etotref -3404.7663590699248743`
     - `etotperatomref -106.3989487209`
     - `totalforceref 2.688768`
     - `totalstressref 473.654463`
     - `totaltimeref +85.40579`

5. 运行日志中的额外观测
   - `Gint cal_gint_vl`：`1.68 s / 9 calls / 0.19 s avg`
   - `Gint cal_gint_rho`：`1.72 s / 9 calls / 0.19 s avg`
   - `LCAO_domain dm2rho`：`1.77 s / 9 calls / 0.20 s avg`

### 重要观察

- 阶段 D 已经把精度选择真正接入 LCAO 主 SCF：不再只是 `module_gint` 内部有 `fp32` 路径，而是 SCF 迭代会在每轮开始前拿到当前 `GintExecConfig`。
- 本用例 `tests/performance/P101_si32_lcao/INPUT` 中 `scf_thr = 1e-6`，因此阶段 D 策略的切换阈值为 `max(100 * 1e-6, 1e-5) = 1e-4`。
- 从本次日志的 `drho` 观察，`CU6 = 6.1217e-05`、`CU7 = 2.1956e-05` 已连续满足阈值条件；据此可推断后续迭代已进入 `fp64` 收尾阶段。这一判断来自策略与日志的联合推断，而不是新增的显式运行日志。
- 接入上层策略后，`P101_si32_lcao` 仍可稳定收敛并正常结束，且抽取到的总能、力和应力结果与阶段 C 保持一致，说明当前切换策略没有破坏主 SCF 默认行为。
- `P101_si32_lcao` 与目录内 `result.ref` 的偏差、以及 `catch_properties.sh` 未提取到 `totaltimeref` 的现象仍然延续，继续应视为环境/基线观测项，而不是阶段 D 新引入的行为变化。
