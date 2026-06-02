# 机器学习分子动力学模块并行优化项目计划书

## 一、项目背景与拟选方向

本项目围绕 ABACUS 分子动力学模块的并行优化展开。分子动力学通过数值求解牛顿运动方程，根据原子位置计算势能和受力，再通过积分器更新原子的位置与速度。当前仓库中 MD 主循环位于 `source/source_md/run_md.cpp`，基础积分算子位于 `source/source_md/md_base.cpp`，并由 `Verlet`、`Nose_Hoover`、`Langevin`、`MSST` 和 `FIRE` 等类实现不同动力学流程；DPMD、NEP 等机器学习势函数通过 `source/source_esolver` 中的 ESolver 接口接入。课程说明中指出，机器学习势函数能够在接近第一性原理精度的同时显著提升计算效率，但势函数接口开销、并行效率、扩展性和测试覆盖仍是需要优化的问题。

本项目拟按以下顺序推进：

1. 完成分子动力学积分器的 OpenMP 并行优化；
2. 建立与积分器优化配套的正确性测试和性能基准测试；
3. 完成机器学习势函数 ABACUS 接口层的 OpenMP 并行优化；
4. 在 DPMD 或 NEP 中选择一个方向进行专项深入优化；
5. 有余力时尝试机器学习势函数局部计算的 CUDA 加速。

整体路线为：先完成积分器优化与测试闭环，再进入势函数 OpenMP 优化，随后选择 DPMD 或 NEP 深入分析，最后尝试 CUDA 局部加速。

------

## 二、项目目标

本项目目标是完成一条较完整的科学计算软件优化流程：阅读 ABACUS 分子动力学模块代码，理解 MD 主循环、积分器、势函数接口和机器学习势函数接入方式；在此基础上，对积分器原子循环、ABACUS 侧势函数数据转换、力和 Virial 后处理等部分进行并行优化；最后通过测试脚本验证计算正确性，并用 benchmark 数据评价优化效果。

| 目标                   | 主要内容                                                     |
| ---------------------- | ------------------------------------------------------------ |
| 积分器并行优化         | 对 `MD_base`/`Verlet` 中按原子更新位置、速度及相关统计量的循环进行 OpenMP 并行 |
| 测试与性能框架         | 建立正确性测试、线程一致性测试和 benchmark 脚本              |
| 势函数 OpenMP 优化     | 分析 ABACUS 侧 DPMD 或 NEP 接口中的坐标转换、单位转换和结果后处理循环，完成局部 OpenMP 加速 |
| 机器学习势函数专项优化 | 选择 DPMD 或 NEP 之一，深入分析 ABACUS 接口层与外部库调用边界上的并行瓶颈 |
| CUDA 扩展尝试          | 对可独立抽取的局部转换或后处理计算尝试 GPU 加速，并与 CPU/OpenMP 版本对比 |

------

## 三、技术路线

分子动力学的核心流程可以概括为：

当前位置
→ 计算势能与受力
→ 积分器更新位置和速度
→ 温度/压力控制
→ 输出结果
→ 进入下一时间步。

其中，原子受力满足：

$$
\mathbf{F}_i = -\nabla_{\mathbf{r}_i}U.
$$

以当前 `Verlet` 类采用的 Velocity-Verlet 半步形式为例，速度和位置更新可写为：

$$
\mathbf{v}_i(t+\Delta t/2)
=
\mathbf{v}_i(t)
+
\frac{\mathbf{F}_i(t)}{2m_i}\Delta t,
$$

$$
\mathbf{r}_i(t+\Delta t)
=
\mathbf{r}_i(t)
+
\mathbf{v}_i(t+\Delta t/2)\Delta t .
$$

这类计算在代码中通常表现为按原子编号循环。对于不同原子之间相互独立或数据依赖较弱的循环，可以采用 OpenMP 进行线程级并行。

本项目主要使用以下并行与测试方法：

| 技术方法              | 用途                                   |
| --------------------- | -------------------------------------- |
| OpenMP `parallel for` | 并行化积分器和势函数中的原子循环       |
| `private` 或局部变量  | 处理线程私有临时变量                   |
| `reduction`           | 处理总能量、动能、温度等全局累加量     |
| `OMP_NUM_THREADS`     | 控制不同线程数下的运行实验             |
| 运行时间统计          | 计算加速比和并行效率                   |
| CUDA kernel           | 有余力时尝试加速机器学习势函数局部计算 |

性能评价指标为：
$$
S_p = \frac{T_1}{T_p},
$$

$$
E_p = \frac{S_p}{p},
$$

其中 $T_1$ 表示单线程或原始版本运行时间，$T_p$ 表示 $p$ 线程运行时间，$S_p$ 为加速比，$E_p$ 为并行效率。

------

## 四、代码阅读与修改范围

本项目重点关注当前仓库中与 MD 主循环、积分器、ESolver 势函数接口、DPMD/NEP 接入、晶胞更新和输入参数相关的代码。代码阅读与修改范围以下列当前实际存在的文件为准。

| 类别               | 实际相关文件                                                 | 文件信息与阅读/修改重点                                      |
| ------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| MD 主流程          | `source/source_md/run_md.cpp`、`source/source_md/run_md.h`   | `Run_MD::md_line()` 根据 `md_type` 选择 `FIRE`、`Nose_Hoover`、`Verlet`、`Langevin`、`MSST`，并组织 `setup()`、`first_half()`、`MD_func::force_virial()`、`second_half()`、`compute_stress()`、输出与 restart。用于理解完整时间步依赖，不作为主要并行修改点。 |
| MD 公共基类        | `source/source_md/md_base.cpp`、`source/source_md/md_base.h` | 保存质量、速度、位置增量、力、应力、Virial 等公共状态；`update_pos()` 和 `update_vel()` 是积分器共用的原子循环，目前采用主进程计算后广播，是积分器 OpenMP 优化的首要修改对象。 |
| Verlet 积分器      | `source/source_md/verlet.cpp`、`source/source_md/verlet.h`   | 实现 NVE/NVT Velocity-Verlet 流程；`first_half()` 调用半步速度和位置更新，`second_half()` 调用半步速度和简单温控器。重点阅读温控器 `thermalize()`、`apply_thermostat()` 中的原子循环与随机数逻辑。 |
| Nose-Hoover 积分器 | `source/source_md/nhchain.cpp`、`source/source_md/nhchain.h` | 实现 NVT/NPT Nose-Hoover Chain、barostat、`update_baro()`、`vel_baro()`、`update_volume()` 等流程。作为可变胞和压力控制依赖阅读对象，避免积分器优化破坏速度、应力和晶胞更新顺序。 |
| Langevin 积分器    | `source/source_md/langevin.cpp`、`source/source_md/langevin.h` | 在真实力基础上构造阻尼力和随机力 `total_force`；`post_force()` 是按原子循环，但涉及随机数和 MPI 广播，作为并行化风险较高的候选对象。 |
| MSST 积分器        | `source/source_md/msst.cpp`、`source/source_md/msst.h`       | 实现冲击波约束动力学，包含 `propagate_vel()`、`propagate_voldot()`、`rescale()` 和晶胞重建。主要用于理解可变胞场景下 `update_pos()` 优化的边界条件。 |
| FIRE 弛豫算法      | `source/source_md/fire.cpp`、`source/source_md/fire.h`       | 基于 `v·F` 调整速度、时间步长和收敛状态；`check_force()`、`check_fire()` 含原子循环与全局归约，是可选并行优化对象。 |
| MD 工具函数        | `source/source_md/md_func.cpp`、`source/source_md/md_func.h` | 提供初始速度、动能、温度、应力、`force_virial()`、输出等工具函数；其中 `kinetic_energy()`、`current_temp()`、`temp_vector()`、`dump_info()`、`get_mass_mbl()` 等循环需要作为测试与后续归约优化参考。 |
| MD 单元测试        | `source/source_md/test/CMakeLists.txt`、`source/source_md/test/verlet_test.cpp`、`source/source_md/test/nhchain_test.cpp`、`source/source_md/test/langevin_test.cpp`、`source/source_md/test/msst_test.cpp`、`source/source_md/test/fire_test.cpp`、`source/source_md/test/md_func_test.cpp`、`source/source_md/test/lj_pot_test.cpp`、`source/source_md/test/setcell.h` | 当前仓库已有 MD 相关单元测试，可作为修改后的正确性验证入口，并补充线程一致性或性能测试脚本。 |
| ESolver 总入口     | `source/source_esolver/esolver.cpp`、`source/source_esolver/esolver.h` | `determine_type()` 与 `init_esolver()` 根据 `esolver_type` 创建 `ESolver_DP`、`ESolver_NEP`、`ESolver_LJ` 等对象，是 MD 与势函数实现之间的接口入口。 |
| NEP 接口           | `source/source_esolver/esolver_nep.cpp`、`source/source_esolver/esolver_nep.h` | ABACUS 侧负责 type map、坐标/晶格数组转换、调用 `nep.compute()`、能量/力/Virial 单位转换和后处理；核心 NEP 推理在外部 NEP 库中，ABACUS 侧主要优化坐标转换、力后处理和 Virial 求和。 |
| DPMD 接口          | `source/source_esolver/esolver_dp.cpp`、`source/source_esolver/esolver_dp.h` | ABACUS 侧负责读取 `dp_rescaling`、`dp_fparam`、`dp_aparam`、type map、坐标/晶格转换、调用 `dp.compute()`、力和 Virial 后处理；核心 DeePMD 推理在外部 DeePMD-kit 中。 |
| LJ 参考势          | `source/source_esolver/esolver_lj.cpp`、`source/source_esolver/esolver_lj.h` | 当前仓库包含 LJ 势函数和对应 MD 测试，可用于不依赖外部 DPMD/NEP 库的 MD 流程正确性与性能基准。 |
| 晶胞与位置更新依赖 | `source/source_cell/update_cell.cpp`、`source/source_cell/update_cell.h`、`source/source_cell/unitcell.cpp`、`source/source_cell/unitcell.h`、`source/source_cell/unitcell_data.h` | `MD_base::update_pos()`、NPT/MSST、restart 输出依赖这些函数更新直接坐标、速度和晶胞派生量。并行修改 MD 位置循环时必须保持 `UnitCell` 更新顺序。 |
| MD 输入参数        | `source/source_io/module_parameter/md_parameter.h`、`source/source_io/module_parameter/input_parameter.h`、`source/source_io/module_parameter/parameter.h`、`source/source_io/module_parameter/read_input_item_md.cpp` | 定义并读取 `md_type`、`md_thermostat`、`md_dt`、`md_nstep`、`pot_file`、`dp_rescaling` 等参数。用于确定测试输入和 DPMD/NEP 编译运行条件。 |
| 构建配置           | `source/source_md/CMakeLists.txt`、`source/source_esolver/CMakeLists.txt`、`source/source_cell/CMakeLists.txt`、`source/source_io/module_parameter/CMakeLists.txt`、`CMakeLists.txt` | 用于确认新增 OpenMP 代码、测试文件和 DPMD/NEP 编译宏 `__DPMD`、`__NEP` 的构建关系。 |
| 外部测试脚本       | 自建脚本                                                     | 在已有单元测试基础上，自动设置 `OMP_NUM_THREADS`、运行小体系 MD、比较关键输出并统计耗时。 |

代码阅读重点包括：

1. MD 主循环的执行顺序；
2. `MD_func::force_virial()` 调用 ESolver 的位置；
3. `MD_base::update_pos()`、`MD_base::update_vel()` 及各派生积分器的调用顺序；
4. 按原子编号循环的核心代码；
5. 每个循环中的读写关系；
6. 全局累加变量的处理方式；
7. 线程并行可能涉及的数据竞争；
8. MPI 广播、随机数、restart 和可变胞更新对并行修改的限制；
9. 可以用于正确性验证的输出量。

------

## 五、实现内容

### 5.1 积分器 OpenMP 并行优化

积分器部分作为项目的基础实现模块。计划阅读 `MD_base`、`Verlet` 及其它派生积分器代码，找到位置、速度、力后处理和统计量计算中的原子循环，并对适合并行的循环加入 OpenMP。

主要工作包括：

1. 梳理积分器计算流程；
2. 找到 `update_pos()`、`update_vel()` 和相关统计函数中的核心原子循环；
3. 分析循环迭代间的数据依赖；
4. 使用 OpenMP 并行化原子循环；
5. 对动能、温度、最大力等统计量使用归约处理；
6. 保持原有函数接口和调用方式；
7. 记录修改位置和优化逻辑。

对应测试包括：

| 测试内容       | 说明                               |
| -------------- | ---------------------------------- |
| 编译测试       | 确认修改后代码可以正常编译         |
| 运行测试       | 使用相同输入运行原始版本和并行版本 |
| 结果一致性测试 | 比较位置、速度、能量等关键输出     |
| 多线程测试     | 测试 1、2、4、8 线程               |
| 性能测试       | 记录不同线程数下运行时间           |
| 加速比分析     | 计算加速比和并行效率               |

------

### 5.2 测试与性能基准框架

测试框架贯穿整个项目，用于验证积分器优化、势函数优化和 CUDA 扩展结果。

主要工作包括：

1. 编写自动运行脚本；
2. 自动设置不同线程数；
3. 保存不同配置下的输出日志；
4. 提取运行时间；
5. 比较串行和并行结果；
6. 生成性能数据表；
7. 支持后续势函数和 CUDA 测试复用。

测试脚本重点支持以下配置：

| 配置               | 含义                       |
| ------------------ | -------------------------- |
| 单线程运行         | 作为基准结果               |
| 多线程 OpenMP 运行 | 测试线程级并行效果         |
| 不同体系规模运行   | 观察原子数变化对性能的影响 |
| CPU / GPU 对比     | 用于 CUDA 扩展阶段         |

------

### 5.3 机器学习势函数 OpenMP 并行优化

势函数部分是项目的核心优化内容。当前 ABACUS 仓库中的 DPMD 和 NEP 核心推理分别由外部 DeePMD-kit 和 NEP 库完成；ABACUS 侧主要负责 type map、坐标/晶格格式转换、单位转换、力和 Virial 后处理。计划在 DPMD 或 NEP 接口中选择一个主要对象，识别 ABACUS 侧以原子为单位的热点循环，并进行 OpenMP 并行优化。

主要工作包括：

1. 阅读 ESolver 势函数接口；
2. 选择 DPMD 或 NEP 作为主要分析对象；
3. 梳理 ABACUS 侧数据准备、外部库调用和结果回填流程；
4. 识别坐标转换、力后处理、Virial 求和等热点循环；
5. 分析共享变量和数组写入关系；
6. 使用 OpenMP 并行化合适的循环；
7. 对 NEP 能量或 Virial 求和等累加量使用 reduction；
8. 使用测试框架验证不同线程数下结果一致性；
9. 记录局部运行时间和完整 MD 运行时间。

测试内容包括：

| 测试内容   | 说明                             |
| ---------- | -------------------------------- |
| 能量一致性 | 比较串行与并行版本势能结果       |
| 力一致性   | 比较不同原子的受力结果           |
| 线程一致性 | 测试 1、2、4、8 线程输出是否稳定 |
| 局部耗时   | 统计势函数计算部分耗时           |
| 整体耗时   | 统计完整 MD 运行时间             |
| 性能分析   | 计算加速比和并行效率             |

------

### 5.4 DPMD 或 NEP 专项深入优化

在完成势函数 OpenMP 基础优化后，选择 DPMD 或 NEP 之一进行深入优化。若不修改外部库，专项优化范围主要限定在 ABACUS 接口层和调用边界。

若选择 DPMD，重点分析：

1. DeePMD type map 和输入参数传递；
2. `coord`、`cell`、`atype` 的组织方式；
3. `dp.compute()` 前后的单位转换和结果回填；
4. 临时数组与内存访问模式；
5. 不同体系规模下的性能表现。

若选择 NEP，重点分析：

1. NEP type map 和多元素系统处理；
2. `coord`、`cell`、`_e`、`_f`、`_v` 的组织方式；
3. 每原子能量与 Virial 求和；
4. 数组访问顺序和局部缓存；
5. 不同体系规模下的性能表现。

GPUMD 的 NEP 文档说明，NEP 使用描述符向量表示原子局部环境，并通过简单前馈神经网络表示单个原子的 site energy；其径向和角向描述符都涉及对邻居原子的求和，这为理解 NEP 中的原子循环和并行机会提供了参考。([GPUMD](https://gpumd.org/potentials/nep.html))

------

### 5.5 CUDA 加速机器学习势函数的扩展尝试

CUDA 部分作为项目扩展内容。计划在已经完成 OpenMP 分析的基础上，选择 ABACUS 接口层中可独立抽取的局部计算环节尝试 GPU 加速；若要加速 DPMD/NEP 核心神经网络或描述符计算，需要进入对应外部库，不作为本项目基本目标。

候选对象包括：

1. 坐标数组格式转换；
2. 力数组后处理；
3. NEP Virial 求和；
4. 单位缩放；
5. 大量原子上的相同结构循环。

主要工作包括：

1. 分析候选计算是否适合 GPU 并行；
2. 设计 CPU 与 GPU 数据传输方式；
3. 编写局部 CUDA kernel；
4. 对比 CPU/OpenMP 和 CUDA 版本结果；
5. 测试不同体系规模下的运行时间；
6. 分析 GPU 加速中的数据传输开销。

CUDA 部分的目标是完成一个局部、可测试、可对比的 GPU 加速尝试，为报告提供机器学习势函数异构并行优化的扩展内容。

------

## 六、测试体系与实验设计

本项目采用由小到大的测试体系，先验证正确性，再进行性能分析。课程说明中建议的测试体系包括水、硅、铝和蛋白质等规模由小到大的体系。

| 测试体系   | 原子数     | 用途                     |
| ---------- | ---------- | ------------------------ |
| 小型体系   | 几十个原子 | 快速编译和运行验证       |
| 水分子体系 | 64         | 初级正确性测试           |
| 硅晶体体系 | 128        | 基准性能测试             |
| 铝体系     | 256        | 多线程性能测试           |
| 更大体系   | 512 以上   | 深入优化与 CUDA 扩展测试 |

实验安排如下：

| 实验       | 内容                                          | 目的                        |
| ---------- | --------------------------------------------- | --------------------------- |
| 基准实验   | 运行原始版本                                  | 获得 baseline 时间          |
| 积分器实验 | 对积分器进行 OpenMP 优化                      | 测试基础并行效果            |
| 积分器测试 | 比较串行和并行结果                            | 验证位置、速度、能量一致性  |
| 势函数实验 | 对 DPMD 或 NEP 接口层热点循环进行 OpenMP 优化 | 测试接口层优化效果          |
| 势函数测试 | 比较能量和力                                  | 验证势函数并行正确性        |
| 专项实验   | 深入分析 DPMD 或 NEP                          | 展示机器学习势函数优化效果  |
| CUDA 实验  | 对局部计算尝试 GPU 加速                       | 比较 CPU/OpenMP 与 GPU 性能 |

性能结果记录表如下：

| 模块          | 体系     | 原子数 | 线程/设备 | 运行时间 / s | 加速比 | 并行效率 |
| ------------- | -------- | ------ | --------- | ------------ | ------ | -------- |
| 积分器        | 水       | 64     | 1 线程    | 待测         | 1.00   | 1.00     |
| 积分器        | 水       | 64     | 4 线程    | 待测         | 待测   | 待测     |
| 势函数        | 硅       | 128    | 1 线程    | 待测         | 1.00   | 1.00     |
| 势函数        | 硅       | 128    | 4 线程    | 待测         | 待测   | 待测     |
| DPMD / NEP    | 铝       | 256    | 8 线程    | 待测         | 待测   | 待测     |
| CUDA 局部计算 | 更大体系 | 512+   | GPU       | 待测         | 待测   | 待测     |

------

## 七、成员分工

本项目采用“共同推进基础流程 + 模块主责分工”的方式。

| 成员             | 主责方向             | 主要任务                                                     |
| ---------------- | -------------------- | ------------------------------------------------------------ |
| 成员 A（陈怀瑜） | 积分器并行优化       | 阅读积分器代码，完成 OpenMP 并行修改，记录优化逻辑和局部性能结果 |
| 成员 B（胡婧瑶） | 机器学习势函数优化   | 阅读 DPMD 或 NEP 接口代码，完成热点分析、OpenMP 优化和专项深入优化 |
| 成员 C（李嘉宁） | 测试框架与 CUDA 扩展 | 编写测试脚本和 benchmark 工具，整理性能数据；有余力时负责 CUDA 局部加速尝试 |

协作安排如下：

1. 三人共同完成 ABACUS 编译、baseline 运行和 MD 主流程阅读；
2. 成员 A 优先推进积分器并行，成员 C 同步完成积分器测试脚本；
3. 成员 B 同步阅读势函数代码，随后推进势函数 OpenMP 优化；
4. 成员 C 的测试框架服务于积分器、势函数和 CUDA 三个阶段；
5. 成员 B 与成员 C 共同确定 CUDA 候选计算环节；
6. 所有代码修改通过 Git 分支管理，合并前完成编译和基本运行测试；
7. 最终报告由三人共同整理，成员 C 汇总性能数据，成员 A 和 B 分别完成对应模块说明。

------

## 八、进度安排

结合课程大作业从计划书、算法文档、测试总结、代码修改、优化报告到最终报告的整体节奏，本项目按六个阶段推进。

| 阶段                                 | 时间范围               | 主要任务                                                     | 阶段成果                             |
| ------------------------------------ | ---------------------- | ------------------------------------------------------------ | ------------------------------------ |
| 阶段一：项目调研与环境准备           | 第 10 周前后           | 明确选题方向，阅读课程说明，配置代码环境，编译并运行 ABACUS 最小 MD 示例 | 项目计划书，baseline 运行记录        |
| 阶段二：代码结构与算法流程梳理       | 第 11 至 13 周前半     | 阅读 `run_md.cpp`、`md_base.cpp`、各积分器和 ESolver 势函数接口，画出 MD 主流程，确定积分器与势函数优化位置 | 算法流程文档，代码入口说明           |
| 阶段三：积分器 OpenMP 优化与测试     | 第 13 周后半至第 14 周 | 完成 `MD_base`/`Verlet` 相关原子循环的 OpenMP 并行，建立积分器正确性测试和性能测试 | 积分器并行代码，测试数据，初步加速比 |
| 阶段四：势函数 OpenMP 优化与专项分析 | 第 14 至 15 周         | 选择 DPMD 或 NEP，完成热点循环分析和局部 OpenMP 优化，整理能量/力一致性测试 | 势函数并行代码，热点分析，性能对比   |
| 阶段五：深入优化与 CUDA 扩展         | 第 15 至 16 周         | 深入 DPMD 或 NEP 的 ABACUS 接口流程，优化内存访问或循环结构；尝试局部 CUDA kernel | 专项优化报告，CUDA 尝试结果          |
| 阶段六：结果汇总与最终报告           | 期末前                 | 汇总代码修改、测试结果、性能图表和项目结论，完成最终报告和展示材料 | 最终报告，代码说明，完整测试结果     |

------

## 九、预期成果

本项目预期提交以下成果：

1. 项目计划书；
2. 算法流程与代码结构说明文档；
3. 积分器 OpenMP 并行优化代码；
4. 积分器正确性测试和性能测试结果；
5. 机器学习势函数 OpenMP 并行优化代码；
6. 势函数能量、力一致性测试结果；
7. DPMD 或 NEP 专项深入优化分析；
8. *CUDA 局部加速尝试代码和性能对比；
9. benchmark 脚本和运行说明；
10. 运行时间、加速比、并行效率表格；
11. 代码修改与重构说明；
12. 最终优化效果和总结报告。

------

## 十、风险控制与协作方式

| 风险                  | 控制方式                                                     |
| --------------------- | ------------------------------------------------------------ |
| 代码规模较大          | 以 `run_md.cpp`、`md_base.cpp`、积分器和 ESolver 势函数接口为阅读主线 |
| OpenMP 数据竞争       | 重点检查共享变量、全局累加量和数组写入                       |
| 浮点结果存在微小差异  | 使用误差阈值比较能量、力、位置和速度                         |
| 整体加速比受限        | 同时统计局部模块耗时和完整 MD 耗时                           |
| CUDA 数据传输开销较大 | 选择计算密集且循环结构清晰的局部计算作为加速对象             |
| 三人代码合并冲突      | 使用 GitHub 分支管理，少量多次提交，合并前编译测试           |
| 测试结果难复现        | 统一输入文件、线程数配置、日志格式和数据记录方式             |

------

## 十一、参考文献与参考资料

1. 课程项目说明文件：`06_md.md`。其中包含项目背景、ABACUS 分子动力学模块结构、建议优化方向、测试体系和提交要求。
2. ABACUS 官方代码仓库：
   https://github.com/deepmodeling/abacus-develop
   用于阅读分子动力学模块、积分器模块、势函数接口和机器学习势函数相关代码。课程说明中也将该仓库列为参考资源。
3. ABACUS Pull Request #6603：
   https://github.com/deepmodeling/abacus-develop/pull/6603
   该 PR 主题为 “Add NEP as esolver”，可作为了解 NEP 与 ABACUS 接口关系的参考。([GitHub](https://github.com/deepmodeling/abacus-develop/pull/6603))
4. GPUMD NEP 文档：
   https://gpumd.org/potentials/nep.html
   该文档介绍了 NEP 的神经网络模型、描述符、模型维度和优化过程，可用于理解 NEP 中原子局部环境描述符与并行计算结构。([GPUMD](https://gpumd.org/potentials/nep.html))
5. OpenMP 官方文档：
   https://www.openmp.org/
   用于查询 `parallel for`、`private`、`reduction`、线程数控制等 OpenMP 语法和并行编程规范。
6. CUDA 官方文档：
   https://docs.nvidia.com/cuda/
   用于后续 CUDA kernel、线程组织、内存传输和 GPU 性能分析的学习与实现参考。