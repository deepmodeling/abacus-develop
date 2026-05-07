# Plan.md：平面波并行优化实施计划

## Summary
- 在 `/home/yangxu/abacus-develop/Plan.md` 新增内部项目计划，服务于小组协作，不作为官方内容直接提交 upstream。
- 主线选择“稳妥可交付”：优先完成题目 1、2、3、8；题目 4/5/6/7 和附加 Gamma Only 作为扩展。
- 核心目标：保持数值结果不变，降低平面波分布、FFT 数据搬运、Gather/Scatter 和重复 G/G+k 计算的耗时，并形成可复现实验报告。

## Key Changes
- 分布优化：在 `module_pw/pw_distributeg.cpp` 中并行化 `count_pw_st`，使用 OpenMP reduction 处理 `npwtot/nstot/lix/rix/liy/riy`，确保 `st_length2D/st_bottom2D` 每个 `(ix,iy)` 唯一写入、无数据竞争。
- 通信优化：在 `module_pw/pw_gatherscatter.h` 中保留现有函数签名，新增内部非阻塞 all-to-all 包装；先实现 correctness-first 的 `MPI_Irecv` + `MPI_Isend` + `MPI_Waitall`，再评估是否做真正 overlap。
- FFT 数据搬运优化：在 `module_pw/pw_transform.cpp` 中缓存 aux 指针、减少重复 getter、优化清零/拷贝/取出循环；不改变 FFT 调用顺序和归一化语义。
- 预计算与缓存：复用/减少 `gg/gcar/gk2/igl2isz_k` 相关重复计算，优先优化 `PW_Basis::collect_uniqgg` 与 `PW_Basis_K::setupIndGk/collect_local_pw` 中重复坐标恢复和模长计算。
- 接口边界：不新增 INPUT 参数，不改变公开 API、文件格式、数值输出格式；只允许新增私有 helper、测试用派生类或本地 benchmark 脚本。

## Workstreams
- 工作流 A：`count_pw_st` OpenMP 化，负责正确性对齐、线程数 1/2/4/8 加速比和 race 检查。
- 工作流 B：Gather/Scatter 非阻塞通信，负责 MPI 1/2/4 进程正确性、通信耗时对比和回退方案。
- 工作流 C：FFT 数据搬运与缓存优化，负责 `real2recip/recip2real` 热点清理和 round-trip 误差验证。
- 工作流 D：测试、基准和报告，负责统一构建命令、性能表格、实验环境记录和 PR/作业说明整理。

## Test Plan
- 单元测试：扩展 `source/source_basis/module_pw/test`，覆盖 gamma/non-gamma、xprime true/false、distribution method 1/2、多 MPI rank。
- 回归测试：运行 `ctest -R MODULE_PW_pw_test -V`，并用 `mpirun -np 3/4 ./MODULE_PW_pw_test` 验证并行路径。
- 端到端测试：选 `tests/01_PW` 小算例和 `tests/performance/P000_si16_pw`、`P001_si32_pw` 做结果与性能验证。
- 性能记录：固定编译选项、CPU 型号、MPI 进程数、OMP 线程数，报告 baseline 与优化后 wall time、speedup、最大误差。
- 验收标准：所有现有 PW 测试通过；1 线程结果与 baseline 一致；主要热点在目标场景中有可测加速或给出瓶颈解释。

## Assumptions
- 小组按 4 个工作流认领任务，成员姓名不写死。
- `develop` 继续保持与官方 upstream 一致，所有实现从新功能分支开展。
- 主线不做高风险 Gamma Only 存储格式重构；若主线完成且时间允许，再单独开扩展分支处理。
