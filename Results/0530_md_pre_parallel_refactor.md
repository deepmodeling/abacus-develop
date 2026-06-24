# ABACUS MD 代码修改与重构报告

## 目标边界

本轮重构只调整代码结构，为后续正式 OpenMP/CUDA 并行优化降低改动风险，不引入并行代码。

## 重构点清单

| 编号 | 重构点 | 涉及文件 | 原因 | 行为影响 |
| --- | --- | --- | --- | --- |
| R1 | 抽取 MD 积分器创建逻辑 | `source/source_md/run_md.cpp` | 将 `md_type` 到具体积分器类的选择从 `md_line()` 主循环中分离，减少后续修改主循环时的生命周期风险。 | 使用 `std::unique_ptr<MD_base>` 替代手动 `new/delete`，积分器选择条件和 MD 主循环顺序不变。 |
| R2 | 复用 DP 接口层临时缓冲区 | `source/source_esolver/esolver_dp.h`、`source/source_esolver/esolver_dp.cpp` | DP runner 每步重复创建 `cell`、`coord`、`f`、`v` 容器；将其改为对象成员，为后续数据布局和接口层并行优化做准备。 | 缓冲区生命周期从函数局部变为 `ESolver_DP` 成员；每步仍按原顺序填充坐标、调用 `dp.compute()`、回填能量/力/Virial。 |
| R3 | 复用 NEP 接口层临时缓冲区 | `source/source_esolver/esolver_nep.h`、`source/source_esolver/esolver_nep.cpp` | NEP runner 每步重复创建 `cell`、`coord`、`v_sum`；成员缓存便于后续对坐标转换、力回填和 Virial 求和做并行化。 | 缓冲区生命周期从函数局部变为 `ESolver_NEP` 成员；NEP column-major 数据布局和数值转换公式不变。 |
| R4 | 拆出 MD 统计纯计算状态 | `source/source_md/md_statistics.h`、`source/source_md/md_func.h`、`source/source_md/md_func.cpp` | 动能、温度、应力后续会涉及 reduction 和线程一致性验证；先拆出无副作用的结果结构，保留旧接口作为包装器。 | 新增 `MDKineticState`、`MDStressState`、`calc_kinetic_state()`、`calc_stress_state()`；`current_temp()` 和 `compute_stress()` 的外部接口和返回结果不变。 |
| R5 | 建立 MD 测试夹具 | `source/source_md/test/md_test_fixture.h` 及 MD 测试文件 | 多个测试重复创建 `UnitCell`、`Parameter`、`ESolver_LJ`、积分器和数组；抽成夹具后，后续增加线程一致性测试和 benchmark 辅助测试更直接。 | 测试断言主体保持不变；用 RAII 管理测试对象，消除手动释放和遗漏释放风险。 |
| R6 | 移除测试 CMake 重复依赖 | `source/source_md/test/CMakeLists.txt` | `global_variable.cpp` 在 `depend_files` 中重复出现，删除重复项减少测试目标链接输入噪声。 | 不改变测试目标语义。 |

## 保留旧接口的约束

- `Run_MD::md_line()` 的公开函数签名不变。
- `MD_func::current_temp()`、`MD_func::compute_stress()`、`MD_func::temp_vector()` 等旧接口保留。
- DP/NEP 的 `before_all_runners()`、`runner()`、`cal_energy()`、`cal_force()`、`cal_stress()` 接口不变。
- MD 测试中的数值断言未改写。

## 验证记录

生产代码构建：

```bash
cmake -S . -B build-prod -DBUILD_TESTING=OFF -DENABLE_MPI=OFF -DENABLE_LCAO=OFF -DUSE_ELPA=OFF
cmake --build build-prod -j2
```

结果：`abacus_pw_ser` 构建成功。

MPI 测试配置：

```bash
cmake -S . -B build-mpi-test -DBUILD_TESTING=ON
```

结果：配置成功。普通沙箱下 OpenMPI 探测会打印 `opal_ifinit: socket() failed with errno=1`，因此该配置在允许 MPI 探测的环境中运行。

MD 测试目标构建：

```bash
cmake --build build-mpi-test --target MODULE_MD_LJ_pot MODULE_MD_func MODULE_MD_fire MODULE_MD_verlet MODULE_MD_nhc MODULE_MD_msst MODULE_MD_lgv -j2
```

结果：7 个 MD 测试目标均构建成功。

MD 单元测试：

```bash
ctest --test-dir build-mpi-test -R 'MODULE_MD' --output-on-failure
```

结果：7/7 通过。

## 回溯方式

本轮修改按重构主题提交，若后续发现问题，可用 `git log --oneline` 找到对应提交，再用 `git revert <commit>` 回退单个重构点。
