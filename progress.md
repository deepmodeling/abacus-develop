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
