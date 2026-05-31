# OpenMP 基本部分与 NEP 优化修改说明

日期：2026-05-31

本次根据 `Planners/0531/md_parallel_optimization_plan.md`，只实现了收益较明确、数据依赖简单的 OpenMP 优化点；暂未处理随机数相关循环、NHC/MSST 扩展项、DPMD 和 LJ。

## 修改范围

### 基本 MD 循环

- `source/source_md/md_base.cpp`
  - `MD_base::update_pos()`：对 rank 0 上的逐原子位置增量计算加入 `#pragma omp parallel for schedule(static)`。
  - `MD_base::update_vel()`：对 rank 0 上的逐原子速度半步更新加入 OpenMP 并行。
  - MPI 广播和 `unitcell::update_pos_taud()` 顺序保持不变。

- `source/source_md/md_func.cpp`
  - `kinetic_energy()`：使用 OpenMP `reduction(+:ke)` 并行动能求和。
  - `force_virial()`：并行 `force_temp` 到 `force` 的逐原子回填。
  - `temp_vector()`：用 9 个标量 reduction 并行温度张量累加，并显式写回 3x3 矩阵，避免依赖 `matrix::create()` 是否清零。

### NEP 接口层

- `source/source_esolver/esolver_nep.h`
  - 新增 `atom_type_index` 和 `atom_local_index` 缓存，用于把全局原子序号映射到 `UnitCell` 的元素类型和类型内局部序号。

- `source/source_esolver/esolver_nep.cpp`
  - `before_all_runners()`：初始化全局原子到 `UnitCell` 存储位置的索引缓存。
  - `runner()`：
    - 并行 NEP 坐标缓冲区填充；
    - 使用 OpenMP reduction 并行每原子能量求和；
    - 并行 NEP 力回填和单位转换；
    - 使用线程局部 9 分量数组并行 per-atom virial 求和。

## 并行策略

- 所有新增循环使用 `schedule(static)`。
- 对短循环使用 `if (nat >= 256)` 或同类阈值，降低小体系线程启动开销。
- 不改变 MD 时间步顺序、MPI 广播位置、NEP 外部库调用边界。
- 浮点归约可能带来末位差异，后续测试应使用容差比较。

## 未执行项

- 按要求未运行测试。
- 未并行 `rand_vel()`、Anderson thermostat、Langevin 随机力等随机数相关循环。
- 未改 DPMD、NHC/MSST 扩展优化、LJ benchmark 优化。

## 回溯方式

本次修改已按 git commit 分开记录：

- 代码优化提交：`optimize: add OpenMP to MD base loops and NEP interface`
- 文档提交：`docs: record OpenMP NEP and MD base changes`

如需回滚，可对对应提交执行 `git revert <commit>`。
