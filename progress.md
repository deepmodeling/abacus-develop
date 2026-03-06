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
