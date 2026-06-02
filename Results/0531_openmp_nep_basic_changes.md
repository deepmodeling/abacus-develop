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

- 本文件初版按当时要求未运行测试；2026-06-01 已补充独立 microbenchmark，见下方“性能测试记录”。
- 未并行 `rand_vel()`、Anderson thermostat、Langevin 随机力等随机数相关循环。
- 未改 DPMD、NHC/MSST 扩展优化、LJ benchmark 优化。

## 性能测试记录

### 测试文件与运行方式

测试文件已放在代码文档专用测试目录：

- `Test/openmp_nep_basic_benchmark.cpp`
- `Test/run_openmp_nep_basic_benchmark.sh`
- `Test/results/openmp_nep_basic_benchmark.csv`
- `Test/results/openmp_nep_basic_benchmark.log`

本次测试使用独立 C++ OpenMP microbenchmark，而不是完整 ABACUS 端到端 MD 输入。原因是 NEP 端到端测试依赖外部 NEP 模型文件和编译宏；本 benchmark 直接复现本次改动中真正发生变化的逐原子循环和 reduction 结构，可以隔离测量 MD 基础循环与 NEP 接口层后处理的 OpenMP 收益。

编译与运行命令：

```bash
cd /root/abacus-md-refactor
./Test/run_openmp_nep_basic_benchmark.sh
```

脚本内部编译命令：

```bash
g++ -O3 -std=c++17 -fopenmp Test/openmp_nep_basic_benchmark.cpp -o Test/build/openmp_nep_basic_benchmark
```

测试参数：

| 参数 | 值 |
| --- | --- |
| 编译器 | g++ 11.3.0 |
| 原子数 | 2,000,000 |
| 每项重复次数 | 5 |
| 线程数 | 1 / 2 / 4 / 8 |
| OpenMP 绑定 | `OMP_PROC_BIND=close`, `OMP_PLACES=cores` |

每个 kernel 同时运行串行版本和 OpenMP 版本，表中 `serial_ms` 与 `omp_ms` 均为单次调用平均耗时。`speedup = serial_ms / omp_ms`，`efficiency = speedup / threads`。

### 测试结果

| kernel | threads | serial_ms | omp_ms | speedup | efficiency | max_abs_diff |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| md_update_pos | 1 | 13.801657 | 13.688899 | 1.008237 | 1.008237 | 0 |
| md_update_pos | 2 | 13.631519 | 6.941512 | 1.963768 | 0.981884 | 0 |
| md_update_pos | 4 | 13.600572 | 3.442708 | 3.950545 | 0.987636 | 0 |
| md_update_pos | 8 | 13.564951 | 1.819345 | 7.455954 | 0.931994 | 0 |
| md_update_vel | 1 | 15.260198 | 15.150155 | 1.007263 | 1.007263 | 0 |
| md_update_vel | 2 | 15.062960 | 7.496135 | 2.009430 | 1.004715 | 0 |
| md_update_vel | 4 | 14.971725 | 3.858649 | 3.880043 | 0.970011 | 0 |
| md_update_vel | 8 | 14.953506 | 2.099304 | 7.123076 | 0.890385 | 0 |
| md_kinetic_energy | 1 | 7.883401 | 7.899684 | 0.997939 | 0.997939 | 0 |
| md_kinetic_energy | 2 | 7.545974 | 3.789572 | 1.991247 | 0.995623 | 2.036415e-10 |
| md_kinetic_energy | 4 | 7.450567 | 2.004642 | 3.716657 | 0.929164 | 2.242899e-10 |
| md_kinetic_energy | 8 | 7.387472 | 1.038517 | 7.113485 | 0.889186 | 1.775646e-10 |
| md_temp_vector | 1 | 9.095514 | 8.946827 | 1.016619 | 1.016619 | 0 |
| md_temp_vector | 2 | 9.008466 | 4.422689 | 2.036875 | 1.018438 | 1.600000e-10 |
| md_temp_vector | 4 | 8.810445 | 2.291496 | 3.844844 | 0.961211 | 1.489155e-10 |
| md_temp_vector | 8 | 8.711893 | 1.178607 | 7.391688 | 0.923961 | 1.680718e-10 |
| md_force_copy | 1 | 10.244239 | 10.152378 | 1.009048 | 1.009048 | 0 |
| md_force_copy | 2 | 10.280891 | 5.134739 | 2.002223 | 1.001111 | 0 |
| md_force_copy | 4 | 10.098843 | 2.623814 | 3.848918 | 0.962229 | 0 |
| md_force_copy | 8 | 10.094556 | 1.422335 | 7.097174 | 0.887147 | 0 |
| nep_coord_fill | 1 | 11.474562 | 11.355256 | 1.010507 | 1.010507 | 0 |
| nep_coord_fill | 2 | 11.422869 | 5.694671 | 2.005887 | 1.002944 | 0 |
| nep_coord_fill | 4 | 11.268058 | 2.930989 | 3.844456 | 0.961114 | 0 |
| nep_coord_fill | 8 | 11.123313 | 1.572314 | 7.074487 | 0.884311 | 0 |
| nep_energy_sum | 1 | 3.048584 | 3.015369 | 1.011015 | 1.011015 | 0 |
| nep_energy_sum | 2 | 3.080593 | 1.524430 | 2.020816 | 1.010408 | 8.811185e-09 |
| nep_energy_sum | 4 | 3.075540 | 0.761431 | 4.039160 | 1.009790 | 2.277375e-08 |
| nep_energy_sum | 8 | 3.103870 | 0.386627 | 8.028073 | 1.003509 | 1.861918e-08 |
| nep_force_fill | 1 | 9.623283 | 9.488546 | 1.014200 | 1.014200 | 0 |
| nep_force_fill | 2 | 9.581044 | 4.787944 | 2.001077 | 1.000538 | 0 |
| nep_force_fill | 4 | 9.346130 | 2.425673 | 3.853005 | 0.963251 | 0 |
| nep_force_fill | 8 | 9.329194 | 1.263027 | 7.386377 | 0.923297 | 0 |
| nep_virial_sum | 1 | 29.281022 | 13.105685 | 2.234223 | 2.234223 | 0 |
| nep_virial_sum | 2 | 29.324353 | 6.786994 | 4.320669 | 2.160334 | 9.997166e-09 |
| nep_virial_sum | 4 | 29.311881 | 3.565785 | 8.220315 | 2.055079 | 6.324626e-09 |
| nep_virial_sum | 8 | 29.241326 | 2.000102 | 14.619914 | 1.827489 | 8.529241e-09 |

### 8 线程优化程度汇总

| 优化点 | 对应源码位置 | 8 线程加速比 | 优化程度判断 |
| --- | --- | ---: | --- |
| 位置更新 `update_pos()` | `source/source_md/md_base.cpp` | 7.46x | 接近线性扩展 |
| 速度更新 `update_vel()` | `source/source_md/md_base.cpp` | 7.12x | 接近线性扩展 |
| 动能归约 `kinetic_energy()` | `source/source_md/md_func.cpp` | 7.11x | 接近线性扩展 |
| 温度张量 `temp_vector()` | `source/source_md/md_func.cpp` | 7.39x | 接近线性扩展 |
| 力回填 `force_virial()` | `source/source_md/md_func.cpp` | 7.10x | 接近线性扩展 |
| NEP 坐标填充 | `source/source_esolver/esolver_nep.cpp` | 7.07x | 接近线性扩展 |
| NEP 能量求和 | `source/source_esolver/esolver_nep.cpp` | 8.03x | 接近线性，轻微超线性来自计时波动和 reduction 分块缓存收益 |
| NEP 力回填 | `source/source_esolver/esolver_nep.cpp` | 7.39x | 接近线性扩展 |
| NEP virial 求和 | `source/source_esolver/esolver_nep.cpp` | 14.62x | 改进后表现很好，超线性主要来自算法写法减少循环开销和整数除法 |

### 分析

本次 OpenMP 优化收益最明确的是逐原子写入类循环，包括 `update_pos()`、`update_vel()`、`force_virial()` 力回填、NEP 坐标填充和 NEP 力回填。这些循环的迭代之间没有写冲突，每个线程处理独立原子编号，使用 `schedule(static)` 后负载均匀、同步开销低。因此在 8 线程下普遍达到约 7x 到 8x 的加速，适合作为优先合入的低风险优化点。

动能和温度张量归约也表现较好。`kinetic_energy()` 使用单个 `reduction(+:ke)`，8 线程加速约 7.11x；`temp_vector()` 使用 9 个标量 reduction，8 线程加速约 7.39x。二者均属于低风险、高收益的 reduction 优化。归约结果与串行存在 `1e-10` 量级差异，来源是并行归约改变了浮点求和顺序。

NEP 接口层的坐标填充、能量求和、力回填均有明显收益。坐标和力回填是典型内存带宽型循环，在 8 线程下分别达到约 7.07x 和 7.39x；能量求和是简单 reduction，8 线程约 8.03x。它们不改变 `nep.compute()` 外部库调用边界，只减少 ABACUS 侧数据准备和后处理的串行开销，适合大原子数 ML-MD。

NEP virial 求和已从最初的线程局部数组加 `critical` 合并，改为 9 个显式标量 reduction。新写法避免了 `critical` 同步，也避免了循环内 `index / nat` 整数除法；同时把 9 个 virial 分量放在同一个按原子编号的循环中处理，减少循环控制开销。复测后 8 线程加速约 14.62x，1 线程 OpenMP 版本也快于原串行写法。这里的“超线性”不能简单理解为线程扩展超过物理上限，而是因为 OpenMP 版本的循环组织本身比原串行双层循环更高效。

正确性方面，所有逐原子写入类 kernel 的 `max_abs_diff` 为 0；归约类 kernel 的误差在 `1e-8` 以内，属于浮点并行归约改变求和顺序导致的末位差异。该结果符合本次方案中“使用容差比较，不追求 bitwise 相同”的预期。

## 正式构建与单测验证

2026-06-02 针对 PR 准备补充了正式 ABACUS 构建和 MD 单测。microbenchmark 只用于隔离评估循环性能，不能替代以下正式验证。

### 生产构建

使用已有 `build-prod` 配置构建完整串行 PW 目标：

```bash
cmake --build build-prod -j2
```

结果：构建通过，最终目标 `abacus_pw_ser` 成功生成。构建过程中 `esolver` object library 也完成构建，并生成：

```text
build-prod/source/source_esolver/CMakeFiles/esolver.dir/esolver_nep.cpp.o
```

这说明当前 `source/source_esolver/esolver_nep.cpp` 在 ABACUS 正式构建配置中可编译通过。当前机器未安装 NEP 外部库，未找到 `nep.h` 或 `libnep`，因此本次无法验证定义 `__NEP` 后与真实 NEP 库的链接配置。

### MD 单元测试

先构建 MD 相关测试目标：

```bash
cmake --build build-mpi-test --target MODULE_MD_LJ_pot MODULE_MD_func MODULE_MD_fire MODULE_MD_verlet MODULE_MD_nhc MODULE_MD_msst MODULE_MD_lgv -j2
```

结果：7 个 MD 测试目标全部构建通过。

随后运行 MD 测试集合：

```bash
ctest --test-dir build-mpi-test -R 'MODULE_MD' --output-on-failure
```

结果：

```text
100% tests passed, 0 tests failed out of 7

MODULE_MD_LJ_pot  Passed
MODULE_MD_func    Passed
MODULE_MD_fire    Passed
MODULE_MD_verlet  Passed
MODULE_MD_nhc     Passed
MODULE_MD_msst    Passed
MODULE_MD_lgv     Passed
```

## 回溯方式

本次修改已按 git commit 分开记录：

- 代码优化提交：`optimize: add OpenMP to MD base loops and NEP interface`
- 文档提交：`docs: record OpenMP NEP and MD base changes`

如需回滚，可对对应提交执行 `git revert <commit>`。
