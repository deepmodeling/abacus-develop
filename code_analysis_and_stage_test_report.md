# DP/NEP 代码分析与阶段性测试报告

## 1. 代码分析结论

本阶段分析目标是确认 ABACUS 当前仓库中 DeePMD 和 NEP 机器学习势函数的真实接入位置、调用链路、外部依赖和阶段性验证结果，为后续 CUDA 加速设计提供依据。

### 1.1 代码组织结论

当前仓库中的机器学习势函数并不位于旧资料中提到的 `source/source_md/potential/ml/` 目录。实际入口在 `source/source_esolver` 下：

- `source/source_esolver/esolver_dp.h`
- `source/source_esolver/esolver_dp.cpp`
- `source/source_esolver/esolver_nep.h`
- `source/source_esolver/esolver_nep.cpp`

`source/source_md` 主要负责 MD 积分、时间步推进、温控/压控等流程；DP/NEP 的能量、力和应力计算通过 `ModuleESolver::ESolver` 多态接口接入。因此，本任务后续 CUDA 改造应围绕 `ESolver_DP` 和 `ESolver_NEP` 的调用链展开，而不是按不存在的 `source_md/potential/ml` 目录设计。

### 1.2 运行调用链结论

两个目标样例的共同主流程为：

```text
INPUT
  -> source_io/module_parameter 读取 calculation/esolver_type/pot_file/md 参数
  -> source_esolver/esolver.cpp 根据 esolver_type 创建 ESolver_DP 或 ESolver_NEP
  -> source_md/run_md.cpp 进入 Run_MD::md_line()
  -> source_md/md_func.cpp 中 MD_func::force_virial()
  -> p_esolver->runner()
  -> cal_energy() / cal_force() / cal_stress()
```

`MD_func::force_virial()` 是 MD 与势函数求解器之间的统一边界。DP/NEP 在自身 `runner()` 中完成外部模型调用和单位换算，上层 MD 流程只读取统一的能量、力和应力结果。

### 1.3 DPMD 接入结论

DPMD 由 `ModuleESolver::ESolver_DP` 实现。它的核心职责是：

- 读取 `pot_file` 指向的 DeePMD 模型；
- 将 ABACUS 的晶胞、坐标和原子类型转换为 DeePMD 所需格式；
- 调用 DeePMD-kit 的 `dp.compute()` 完成模型推理；
- 将 DeePMD 返回的 eV、eV/Angstrom 等量转换到 ABACUS 内部单位；
- 将势能、力和 virial/stress 写回 ESolver 接口。

该路径受编译宏 `__DPMD` 控制。若构建时未指定 `DeePMD_DIR` 并成功链接 `libdeepmd_c.so` 或 `libdeepmd_cc.so`，程序即使能识别 `esolver_type = dp`，也会在运行时提示重新编译并退出。

### 1.4 NEP 接入结论

NEP 由 `ModuleESolver::ESolver_NEP` 实现。它的核心职责与 DP 类似，但数据布局不同：

- NEP 的坐标按分量分块存储，即 `x0..xN, y0..yN, z0..zN`；
- 晶胞矩阵使用列主序；
- NEP 返回每原子能量、力和每原子 virial；
- ABACUS 侧负责求和、单位换算和写回应力矩阵。

该路径受编译宏 `__NEP` 控制。CMake 中只有在指定 `NEP_DIR` 且找到 `include/nep.h` 和 `lib/libnep.so` 时才会启用。当前仓库的 `FindNEP.cmake` 注明 NEP 接口目前只支持 CPU 版本，因此单纯在 ABACUS 外壳层加入 CUDA kernel 并不能加速 NEP 的核心模型推理。

### 1.5 CUDA 改造判断

当前 DP/NEP 都属于“ABACUS 外壳 + 外部库计算核心”的结构。真正耗时的模型推理主要发生在 DeePMD-kit 或 NEP 库内部；ABACUS 侧主要承担输入打包、类型映射、单位换算和结果写回。

因此，后续 CUDA 加速应分两层考虑：

- 对 DP：优先确认 DeePMD-kit 是否以 GPU 后端构建，并让 `dp.compute()` 本身走 GPU；ABACUS 侧可进一步减少每步 host vector 重建和不必要的数据拷贝。
- 对 NEP：当前接入的是 CPU NEP 接口。若要获得实质加速，需要扩展或替换 NEP 计算核心，使 `nep.compute()` 内部支持 GPU；仅加速外层求和和单位换算收益有限。

### 1.6 构建与验证结论

当前已完成一版面向 DP/NEP 的最小构建，并用该构建跑通了两个优先样例。由此可以确认：

- `ESolver_DP` 与 `ESolver_NEP` 的入口位置和调用链是正确的；
- DeePMD 和 NEP 的外部依赖接入方式与仓库中的 CMake 逻辑一致；
- 这两个样例可以在最小依赖配置下独立验证，不需要开启整套 LCAO、Libxc、测试框架等额外模块。

## 2. 阶段性测试记录

### 2.1 测试目标

先验证 `readmeplan.md` 中列出的两个优先样例：

- `tests/04_FF/50_DP_Al`
- `tests/04_FF/101_NEP_HfO2`

重点确认当前仓库里的可执行程序是否能走通 DP / NEP 的样例链路，并记录实际输出与阻塞点。

### 2.2 测试环境

- 可执行文件：`build_dp_nep_minimal/abacus_1s`
- 构建方式：最小依赖构建，启用 DeePMD 和 NEP
- 运行方式：单机本地，`OMP_NUM_THREADS=1`
- MPI 处理：`I_MPI_FABRICS=shm`

### 2.3 测试过程

#### 2.3.1 `tests/04_FF/50_DP_Al`

执行命令：

```bash
I_MPI_FABRICS=shm OMP_NUM_THREADS=1 /share/abacus-develop-3.9.0.27/build_dp_nep_minimal/abacus_1s > log.minimal_dp_nep.txt
```

结果：

- 程序正常完成 4 步 MD
- 退出码：`0`
- 成功生成 `OUT.autotest/`

关键输出：

```text
STEP OF MOLECULAR DYNAMICS: 4
...
TIME STATISTICS
```

运行时间：

- `time.json` 记录总耗时约 `1 s`

#### 2.3.2 `tests/04_FF/101_NEP_HfO2`

执行命令：

```bash
I_MPI_FABRICS=shm OMP_NUM_THREADS=1 /share/abacus-develop-3.9.0.27/build_dp_nep_minimal/abacus_1s > log.minimal_dp_nep.txt
```

结果：

- 程序正常完成 4 步 MD
- 退出码：`0`
- 成功生成 `OUT.autotest/`

关键输出：

```text
STEP OF MOLECULAR DYNAMICS: 4
...
TIME STATISTICS
```

运行时间：

- `time.json` 记录总耗时约 `1 s`

### 2.4 测试结论

1. `50_DP_Al` 在最小构建下可以正常跑通。
2. `101_NEP_HfO2` 在最小构建下可以正常跑通。
3. 这说明 `ESolver_DP`、`ESolver_NEP` 以及它们的外部依赖接入链路都可以在当前仓库里完成闭环验证。

### 2.5 结果对比

已将两个样例的 `OUT.autotest/running_md.log` 通过仓库自带的 `catch_properties.sh` 抽取为 `result.out`，并与各自的 `result.ref` 对比。

对比结论如下：

- `etotref` 一致
- `etotperatomref` 一致
- `totalforceref` 一致
- `totalstressref` 一致
- `totaltimeref` 随运行环境变化，不作为严格数值回归项

具体数值：

- `50_DP_Al`：`etotref`、`etotperatomref`、`totalforceref`、`totalstressref` 均一致；`totaltimeref` 为 `0.90`，参考值为 `1.57`
- `101_NEP_HfO2`：`etotref`、`etotperatomref`、`totalforceref`、`totalstressref` 均一致；`totaltimeref` 为 `0.22`，参考值为 `0.02`

### 2.6 运行形态验证

仓库中当前只有这两个 DP/NEP 专项样例，没有额外的 DP 或 NEP 同类目录可再扩展。因此，运行形态验证改为对同一批样例做线程数变化测试：

- `OMP_NUM_THREADS=1`
- `OMP_NUM_THREADS=2`

验证结果表明，两种线程设置下的物理量结果一致，差异仅体现在 `totaltimeref` 上。这说明 DP/NEP 接入链路在单线程和多线程下的数值行为保持稳定。

### 2.7 当前版本性能基线

为后续重构提速建立对照，已对两个优先样例各做 5 轮重复测试，并整理当前版本的平均耗时：

- `50_DP_Al`
  - `OMP=1`：`total_s` 平均 `1.129829s`，`ESolver_DP::runner` 平均 `0.655738s`
  - `OMP=2`：`total_s` 平均 `1.058435s`，`ESolver_DP::runner` 平均 `0.613125s`
- `101_NEP_HfO2`
  - `OMP=1`：`total_s` 平均 `0.160003s`，`ESolver_NEP::runner` 平均 `0.023816s`
  - `OMP=2`：`total_s` 平均 `0.160783s`，`ESolver_NEP::runner` 平均 `0.024155s`

基线结论如下：

- 两个样例的物理量结果都与 `result.ref` 保持一致；
- `OMP=1` 和 `OMP=2` 下的结果一致，说明当前版本数值行为稳定；
- 之后的重构提速评估应以这组平均耗时作为基准，重点比较 `total_s` 与 solver runner 段耗时是否下降。

### 2.8 重构后建议补测项

当前阶段性测试已经覆盖了“能跑通、跑对、能定基线”这三件事。等后续重构完成后，还建议补以下几类测试，作为提速和正确性回归的正式收口：

1. 单元级回归
   - `before_all_runners()`：确认缓存数组、力/virial 容器、`atype` 和势函数对象初始化正确。
   - `type_map()`：确认模型元素表与 `STRU` 标签映射一致，缺元素时仍然能明确报错。
   - `runner()` 错误分支：在未启用 `__DPMD` / `__NEP` 时保持明确退出，而不是静默失败。
   - `cal_energy()` / `cal_force()` / `cal_stress()`：确认结果回写、外压修正、矩阵写回没有被重构破坏。

2. 接口与布局回归
   - DP 的 AoS 行主序坐标打包不变。
   - NEP 的 SoA 列主序坐标打包不变。
   - 任何新加的 host/device 拷贝都要确认不改变原有数值顺序。

3. 集成回归
   - 继续跑 `tests/04_FF/50_DP_Al` 和 `tests/04_FF/101_NEP_HfO2`。
   - 仍以 `OUT.autotest/running_*.log`、`MD_dump`、`result.ref` 为对照。
   - 先保留 1 到 4 步短轨迹回归，不要求长轨迹逐步完全一致。

4. 数值一致性回归
   - 固定结构下比较总能量、每原子力、3x3 应力/virial。
   - 先和当前 CPU baseline 对齐，再对比 GPU 路径。
   - 若后续引入混合精度，再单独放宽阈值，但需保留误差统计。

5. 性能回归
   - 对 `ESolver_DP::runner()` 和 `ESolver_NEP::runner()` 继续拆分计时。
   - 至少保留 `pack_cell_coord`、`model_compute`、`postprocess`、`device_to_host` 等阶段。
   - 用当前版本的 `total_s` 和 solver runner 耗时作为重构前基线，重构后比较是否真实提速。

6. 规模回归
   - 小规模：几十到几百原子。
   - 中规模：几千原子。
   - 大规模：一万原子以上。
   - 重点看重构后性能收益是否随规模放大，同时检查结果漂移是否可控。
