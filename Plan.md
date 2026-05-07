# Plan.md: ABACUS 平面波并行优化实施计划

## 1. 项目定位

本项目面向 ABACUS 平面波基组的生成与应用流程做并行优化，重点关注 `source/source_basis/module_pw` 中的平面波分布、FFT 变换、Gather/Scatter 通信重排，以及 `PW_Basis_K` 中多 k 点场景反复使用的 G 向量和 `|G+k|^2` 数据。计划服务于小组内部协作，不把本分支直接作为官方内容提交 upstream；最终代码改动应从干净的 `develop` 新建功能分支完成。

主线采用稳妥可交付策略：优先完成题目 1、2、3、8，并把题目 4/5/6/7 与完整 Gamma Only 算法作为扩展路线。这样既覆盖作业的核心并行主题，也避免在第一阶段就陷入高风险的数据格式重构。小组同学的计划书按 8 个基础题和 1 个附加题展开，本版在保留完整题目视野的同时，根据实际运行结果重新排列实施优先级。

## 2. 当前基线与热点判断

当前已用本仓库源码重新编译并重跑 `examples/02_scf/01_pw_Si2`，结果可复现：

- ABACUS: `v3.11.0-beta.1`, commit `dce7387e8`
- 运行方式：4 MPI processes, `OMP_NUM_THREADS=1`
- 算例：Si2, PW, `ecutwfc = 60 Ry`, 4x4x4 Gamma-centered k mesh，经对称性约化为 8 个 k 点
- SCF：5 步收敛
- 最终能量：`FINAL_ETOT = -215.5056984233578 eV`
- 与旧日志 `-215.5056984233570 eV` 的差值约 `8.2e-13 eV`，可作为正确性基线

该小算例显示 `count_pw_st/distributeg` 不是端到端热点，而 `PW_Basis_K` 相关 FFT 变换占比明确：

- `HSolverPW::solve`: 0.95 s, 80.59%
- `Diago_DavSubspace::diag_once`: 0.83 s, 70.03%
- `Operator::veff_pw`: 0.67 s, 56.94%
- `PW_Basis_K::recip2real`: 0.43 s, 36.61%
- `PW_Basis_K::real2recip`: 0.34 s, 28.70%

因此后续性能优化不能只看平面波初始化阶段，应以 `real2recip/recip2real -> gather/scatter -> VeffPW/hPsi` 这条热路径为主要目标。`count_pw_st` 仍作为独立题目完成，但用更大 FFT 网格和更多 MPI rank 验证收益。

## 3. 题目映射与目标校准

同学计划书中的 8 个基础题和附加题可映射为以下实施队列：

- 题目 1 `count_pw_st` OpenMP 并行化：进入 P1，是低风险、边界清晰的必做项。
- 题目 2 `gatherp_scatters` 非阻塞通信：并入 P2，与 FFT 前后数据重排一起评估。
- 题目 3 FFT 数据拷贝/重排 OpenMP 优化：进入 P2，是当前小算例中最直接对应热点的必做项。
- 题目 4 多 k 点 Gamma Only、题目 6 Gamma Only 紧凑存储、附加题完整 Gamma Only：作为高风险扩展，先完成调研和接口设计，主线不承诺大范围数据布局重构。
- 题目 5 gather/scatter SIMD：作为 P2 稳定后的局部加速项，优先使用 `#pragma omp simd` 和连续内存循环形态。
- 题目 7 FFT 通信-计算重叠：作为题目 2 的第二阶段，在非阻塞通信正确后再尝试流水线。
- 题目 8 G 向量与 `|G+k|^2` 缓存：进入 P3，适合作为多 k 点场景下的稳定收益项。

目标也需要从“理想加速比”校准为“可验证收益”：同学计划书中提到的 4 线程 4x、通信 2x、Gamma Only 约 50% 存储节省可作为方向性上限，但不能直接写成验收承诺。实际验收以同一机器、同一构建、同一算例下的计时下降、能量误差和无数据竞争为准；若小算例没有可见加速，应扩大到 `P000_si16_pw/P001_si32_pw/P002_si64_pw` 后再下结论。

## 4. 实施范围与优先级

### P0: 正确性与基线

- 固定基线算例：`examples/02_scf/01_pw_Si2`，用于快速 smoke test。
- 固定性能算例：`tests/performance/P000_si16_pw`、`P001_si32_pw`、`P002_si64_pw`，必要时增加 `P003_si128_pw`、`P004_cu4_pw`、`P006_Ge4As8_pw`、`P009_32H2O_pw`，用于观察 FFT/GatherScatter 随体系、元素和 k 点设置变化后的占比。
- 每次优化后必须记录：ABACUS commit、编译选项、CPU 型号、MPI 进程数、OpenMP 线程数、最终能量、SCF 步数、总时间和热点计时。

### P1: OpenMP 并行化 `count_pw_st`

- 修改位置：`source/source_basis/module_pw/pw_distributeg.cpp`
- 目标：并行化 `count_pw_st` 的 `(ix, iy, iz)` 扫描，保持 `st_length2D`、`st_bottom2D`、`npwtot`、`nstot`、`lix/rix/liy/riy` 与串行结果一致。
- 实现策略：每个 `(ix, iy)` stick 只由一个线程写入，计数使用 OpenMP reduction，边界值使用 thread-local 结果归并。
- 验证重点：gamma/non-gamma、`xprime` true/false、distribution method 1/2、不同网格和截断能。

### P2: FFT 数据搬运与 Gather/Scatter 优化

- 修改位置：`pw_transform.cpp`、`pw_transform_k.cpp`、`pw_gatherscatter.h`
- 目标：优先降低 `PW_Basis_K::real2recip` 和 `PW_Basis_K::recip2real` 的数据搬运开销。
- 实现策略：缓存 `auxr/auxg/rspace` 指针，减少重复 getter；优化清零、拷贝、取数循环的 schedule 和内存访问；对 `gatherp_scatters`、`gathers_scatterp` 的重排循环做连续内存复制和 OpenMP 调度优化。
- 通信策略：第一步保持函数签名不变，用 correctness-first 的非阻塞封装替代阻塞 `MPI_Alltoallv`；若 `MPI_Ialltoallv` 可用，优先用它；否则实现 `MPI_Irecv`/`MPI_Isend`/`MPI_Waitall`。
- 验证重点：round-trip 误差、MPI 1/2/4/8 rank、线程数 1/2/4/8、不同 k 点数。

### P3: G 向量与 `|G+k|^2` 缓存

- 修改位置：`pw_basis.cpp`、`pw_basis_k.cpp`
- 目标：减少 `collect_uniqgg`、`setupIndGk`、`collect_local_pw` 中重复恢复 `(ix, iy, iz)` 和重复计算模长。
- 实现策略：提取本地 helper 统一 `ig2isz/is2fftixy -> G direct/cartesian` 转换；在 `PW_Basis_K` 中复用已生成的 `gcar/gdirect` 或建立只读缓存；结构、FFT 网格、k 点或 cutoff 改变时重建缓存。
- 接口边界：不新增 INPUT 参数，不改变公开输出和文件格式；缓存只作为内部实现细节。

## 5. 扩展路线

- SIMD 向量化：在 P2 的重排循环稳定后，再为连续拷贝段添加 `#pragma omp simd` 或编译器友好的循环形式，避免过早写平台相关 intrinsic。
- 通信计算重叠：在非阻塞通信正确后，再尝试按 stick 或 z-plane 分块，形成双缓冲流水线；若小算例无收益，报告中说明通信粒度不足。
- Gamma Only：作为高风险扩展，先分析 `PW_Basis` 和 `PW_Basis_K` 现有 gamma path，再设计紧凑半谱存储和 pack/unpack；主线阶段不重构波函数/电荷密度公共数据布局。多 k 点场景必须先判断是否满足 Gamma Only 的共轭对称条件，不能因为 k 点网格含 Gamma 点就默认启用半谱存储。

## 6. 协作分工

- 工作流 A：`count_pw_st` OpenMP 化与分布正确性测试。
- 工作流 B：FFT 数据搬运、Gather/Scatter 重排和非阻塞通信。
- 工作流 C：G/G+k 缓存、SIMD 尝试和 Gamma Only 扩展调研。
- 工作流 D：测试脚本、性能表格、阶段报告和最终答辩材料。

每个工作流以独立功能分支开发，合并前必须给出最小可复现命令、结果差异和性能数据。共享分支只用于整合，不直接在 `develop` 上提交作业代码。

## 7. 测试与验收

- 单元测试：扩展 `source/source_basis/module_pw/test`，覆盖 PW 分布、FFT round-trip、MPI gather/scatter 和 gamma 分支。
- 快速回归：`examples/02_scf/01_pw_Si2` 必须保持 `FINAL_ETOT` 与基线差异小于 `1e-8 eV`，SCF 步数不异常增加。
- 并行回归：`mpirun -np 1/2/4/8` 与 `OMP_NUM_THREADS=1/2/4/8` 组合抽测，避免数据竞争和 MPI 死锁。
- 性能验收：至少在一个中等 PW 算例中展示 `PW_Basis_K::real2recip/recip2real` 或 Gather/Scatter 计时下降；若某项优化无收益，必须给出瓶颈分析。
- 文档验收：最终报告包含基线表、优化后表、加速比、能量误差、硬件环境和失败/无收益尝试说明。

## 8. 阶段交付

- 4.24（第 10 周）：完成项目调研、基线重跑、热点确认和计划书。
- 5.15（第 13 周）：完成题目 1/2/3/8 的详细算法文档，明确数据结构、同步点和测试方案；题目 4/6/附加题给出 Gamma Only 可行性分析。
- 5.22（第 14 周）：完成当前算法测试与 baseline 性能报告，至少包含 Si2 smoke test 和一个 `tests/performance` PW 算例。
- 5.29（第 15 周）：完成主线代码修改与重构报告，列出改动文件、核心函数、正确性结果和已知风险。
- 6.5（第 16 周）：完成性能优化报告，包含 OpenMP、MPI、SIMD/缓存方向的对比；无收益项要解释瓶颈。
- 6.12（期末）：整合最终报告、答辩材料和可复现实验脚本。

## 9. 风险与约束

- `develop` 始终与官方 upstream 保持一致；内部资料和作业实现留在功能分支。
- 主线不改变 ABACUS 用户输入接口，不改变数值语义，不引入默认启用的新算法开关。
- 当前环境缺少系统 ScaLAPACK 开发包，已通过本地解包方式完成一次 MPI 构建；后续正式性能数据需固定同一构建环境或在报告中说明差异。
- 小 Si2 算例只用于正确性和 smoke test，不能单独代表性能收益；性能结论必须来自更大 PW 算例。
- 非阻塞通信如果只替换调用而没有形成真实重叠，可能只降低很少开销；报告中应区分“通信 API 改造”和“端到端性能收益”。
- Gamma Only 改动牵涉波函数、电荷密度、势和 FFT 接口，风险高于 OpenMP/MPI 局部优化；除非前置任务完成且测试充分，否则只提交设计和原型，不强行进入主线。
