# PR #8317 第 B 批修复结果

## 1. 本批完成项

### B 类审阅意见

- `4628814647`：`make_local_metadata()` 改为先计算精确 `nmeta`，构造
  `std::vector<int>(nmeta)` 后按偏移赋值，不再 `reserve + push_back`。
- `4628880948`：保留 `int` 作为 BLACS/ScaLAPACK 边界类型；在
  `matrix2d_mpi.cpp::checked_mpi_count()` 的 `size_t -> int` 边界保留
  `FIXME`，没有机械改成 `size_t`。
- `4629051593`：公共内部 helper 使用 `make_byte_offsets()` 命名，保留显式
  scan，没有为缩短代码强改为 STL 累加。
- `4629082027`：`allgather_nbytes()` 保留两阶段流程；第一阶段通过临时
  `ParaWorld::make_mpi(tag, comm)` 调用现有
  `Parallel::allgather_int()` 收集各 rank 的 byte count，第二阶段继续使用
  `allgatherv_bytes()`。没有新增 `Parallel_Common` wrapper API。
- `4657484541`：为 collective/LAPACK/reduction/ScaLAPACK 后端头和新的 MPI
  helper 增加 `@internal` 文档，明确调用者应只 include
  `matrix2d.h`/`matrix2d_linalg.h`；未改 PIMPL。
- `4657491818`：`assert_compatible_layout()` 改为非模板私有成员，实现移入
  `matrix2d.cpp`；原有显式类实例化保留。
- `4657493928`：分布句柄初始化函数改为按值接收
  `std::shared_ptr<const Parallel_Blkcyclic>` 并 `std::move` 到成员；
  `Matrix2dView` 继续持有 `shared_ptr`。
- `4657511707` / `4657523602`：删除每个 rank 一份
  `vector<unsigned char>` 的复制；all-gather/gather 只在单个连续
  `recvbuf` 上工作，并通过轻量 `RankBufferView {data, size}` 反序列化。
- `4657536301` / `4657539879`：保留自定义 gather/scatter。新增文件级说明，
  明确真正的 block-cyclic owner map 由
  `Parallel_Blkcyclic`/`Cblacs_gridinit`/`numroc_`/`descinit_` 定义；scatter
  打包循环提取为 `pack_rank_values()`，这里只按已有 owner map 搬值。
- `4657554100`：LAPACK 的 `n`/`nrhs`/`lda` 仍为 `int`；增加 `n_` 及
  `nrhs_` 作为 `std::size_t`，统一用于 vector size 和下标。
- `4657566089`：linalg dispatcher 继续统一执行 shape/layout 前置校验；仅将
  `desc[...]` 比较提取为 ScaLAPACK 专用
  `matrix2d_impl::scalapack::assert_compatible_layout()`。
- `4657568182` / `4657568544`：`axpby_impl()` 与 `scale_impl()` 外层循环增加
  `_OPENMP` 守卫及
  `#pragma omp parallel for collapse(2) schedule(static)`。
- `4657580013`：新增 `matrix2d_mpi.{h,cpp}`，集中
  `checked_mpi_count()`、`make_byte_offsets()` 和 rank-count gather；
  collective 与 reduction 共用，MPI 聚合通过 `ParaWorld` 调用现有 wrapper。
- `4657597193`：`Matrix2dView` 保持 `T* data_` 与分布共享句柄，没有改为
  `weak_ptr`。
- 额外发现：修正 ScaLAPACK `p[dz]getri_()` 和 `ScalapackConnector::getri()`
  输出参数上的错误 `const`，包括 `A`、workspace、workspace size、`iwork`
  与 `info`。

### 用户拍板

- 拍板③：按值接收 `shared_ptr` 后 move，`Matrix2dView` 继续共享句柄。
- 拍板④：在 `size_t -> int` 检查处保留 `FIXME`，未改成普通约束说明。
- 拍板⑥：dispatcher 保留统一 shape/layout 检查；低层 descriptor 判断进入
  ScaLAPACK helper。
- 拍板⑦：公共 API 拆为 `has_distribution()`（含 serial 1x1 handle）和
  `is_active()`（非 serial 的 active MPI 分布）；内部继续使用
  `distribution_is_active()`/`distribution_is_serial()`。所有调用点已按语义
  迁移，例如 LAPACK 的 `trans` 依据 column-major handle 使用
  `has_distribution()`，算法分派使用 `is_active()`。

## 2. W_INT 核对

批 A 已使用 `W_INT = sizeof(int)`，没有使用 `NBYTES_INT`。collective 和
reduction 的定义处均已有注释，说明其是 native `int` 的 byte 数，因此本批
无需再改名。

## 3. 跳过项

- 未发现因 develop merge 而消失、需要单独跳过的 B 类符号。
- C 批的数值返回类型/API 项（`real_type`、`abs2`、`inv` 数值告警等）按用户
  要求未做。
- 工作树中未找到 `THOROUGH_review_plan.md`；本批依据用户消息中的 B 类清单、
  拍板项和只读 GitHub 审阅评论上下文执行，没有据此遗漏列表中的 B 类条目。

## 4. 验证

所有构建和测试均在 `module load abacus-dev/1.0` 后执行。

- `cmake --build build-matrix2d-cxx11 --target base -j2`
  - 退出码 0；编译命令为 `-std=gnu++11`。
- `cmake --build build-matrix2d-mpi-cxx11 --target base -j2`
  - 退出码 0；编译命令为 `-std=gnu++11`。
- `cmake --build build-matrix2d-serial --target MODULE_BASE_matrix2d MODULE_BASE_matrix2d_linalg -j2`
  - 退出码 0。
- `cmake --build build-matrix2d-mpi --target MODULE_BASE_matrix2d MODULE_BASE_matrix2d_linalg -j2`
  - 退出码 0。
- `OMP_NUM_THREADS=1 ctest --test-dir build-matrix2d-serial --output-on-failure -R 'MODULE_BASE_matrix2d'`
  - 退出码 0，2/2 通过。
- `OMP_NUM_THREADS=1 ctest --test-dir build-matrix2d-mpi --output-on-failure -R 'MODULE_BASE_matrix2d'`
  - 退出码 0，6/6 通过；包含 serial 2 项及 MPI 2-rank/4-rank 各 2 项。
- OMP 防竞争复跑使用 `OMP_NUM_THREADS=2`：
  - serial 2/2 通过。
  - MPI 6/6 通过。
- `git diff --check`
  - 退出码 0。
- `python3 tools/03_code_analysis/agent_governance_check.py --staged --format text`
  - 退出码 0；仅保留 `matrix2d_mpi.h` 的直接 include 警告和本地无 PR body 的
    documentation-sync 警告。
- `agent_governance_check.py --base fbaf69833 --head 0ac210ae3`
  - 退出码 0，无 blocker；warning 同上。
- `code_quality_score.py` 扫描 17 个改动 C++ 文件，15/17 达到 60 分。
  - 未达线历史结构项为既有 612 行 `scalapack_connector.h` 和 1612 行
    `test_matrix2d.cpp`；其问题是文件规模等既有债务，不是新增确定性错误。

## 5. Commit

- B 类实现 commit：`0ac210ae3`
  - `refactor(base): address review structure/MPI/OMP findings (part B)`
- 本结果文件使用紧随其后的独立 docs commit 记录。

## 6. 工作树

报告提交后工作树干净；未 push、未修改 auto-merge、未回复 PR thread。
