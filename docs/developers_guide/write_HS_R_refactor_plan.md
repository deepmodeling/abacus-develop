# write_HS_R.cpp 重构方案与执行计划

## 目标

本文档用于规划 `source/source_io/module_hs/write_HS_R.cpp` 及其相邻实空间矩阵输出代码的重构。重构目标是降低职责耦合、减少重复逻辑、提升可测试性，同时保证后向兼容。

本次重构要拆清楚以下职责：

- 矩阵计算；
- MPI 并行矩阵向串行矩阵收集；
- 输出文件名、目录、追加模式策略；
- CSR 文件头和矩阵主体序列化；
- 旧 sparse-map 输出格式兼容；
- 精度、离子步、MD 输出行为处理。

后向兼容是最高优先级。默认行为下，现有输出文件名、输出目录、文件头格式、CSR R 块顺序、步数编号、文本精度默认值、二进制布局都不能改变。任何可能改变用户可见输出的修改，都必须拆成独立兼容性变更，通过 opt-in 参数或明确文档说明后再做。

## 目录级 review 更新

2026-06-12 对 `source/source_io/module_hs/` 做了一轮静态 review，范围包括：

- `write_HS_R.cpp/.h`;
- `write_HS_sparse.cpp/.h`;
- `single_R_io.cpp/.h`;
- `write_HS.hpp/.h`;
- `write_vxc.hpp`;
- `write_vxc_lip.hpp`;
- `write_vxc_r.hpp`;
- `output_mat_sparse.cpp/.h`;
- `cal_r_overlap_R.cpp/.h`;
- `cal_pLpR.cpp/.h`。

这次 review 的结论是：`module_hs` 不只是 `write_HS_R.cpp` 需要拆职责，整个目录都存在输出格式逻辑分散、MPI rank0 写文件逻辑重复、全局状态依赖重、手写内存管理较多的问题。因此后续重构应把 `write_HS_R.cpp` 作为切入口，但不要只在单文件内做局部整理。

新增的高优先级风险如下：

| 优先级 | 文件/位置 | 问题 | 建议处理阶段 |
| --- | --- | --- | --- |
| P0 | `write_HS_sparse.cpp`, `save_sparse()` binary header | `ofs.write(reinterpret_cast<char*>(0), sizeof(int))` 从空指针写数据，是未定义行为 | Phase 1 |
| P1 | `write_HS.hpp`, 非 MPI 文本输出分支 | `std::ofstream out_matrix(...)` 在 `if/else` 内重新声明，外层 stream 未打开 | Phase 1 |
| P1 | `write_vxc_r.hpp`, `write_Vxc_R()` | `elecstate::Potential*` 和 `vxcs_op_ao` 使用 `new` 后没有释放 | Phase 1 |
| P1 | `cal_pLpR.h/.cpp`, `AngularMomentumCalculator` | `ptr_log` 默认允许 `nullptr`，构造函数直接解引用 | Phase 1 |
| P1 | `output_mat_sparse.cpp/.h` | `p_ham`、`p_dftu` 参数未使用，接口语义漂移 | Phase 2 或独立 API 清理 |
| P2 | `single_R_io.cpp` | 为 CSR indices 使用 rank 临时文件，且 `line` 使用 `new[]/delete[]` | Phase 4 |
| P2 | `cal_r_overlap_R.cpp` | orbital 构造和 `out_rR/out_rR_other` 大量重复，函数过长 | Phase 4 或独立拆分 |
| P2 | `write_vxc*.hpp` | 业务逻辑放在模板头文件中，内存管理和输出逻辑混杂 | Phase 6 之后单独治理 |
| P2 | 多个头文件 | include guard 使用 `__...` 保留标识符风格，注释中有过多日期/作者型信息 | 风格清理阶段 |

这些问题中，P0/P1 属于“修复已有缺陷”，不应等到大规模结构重构后再处理。但修复时仍必须遵守兼容原则：除修复未定义行为或无效写入外，不顺手改变文件名、目录、文本格式、精度或步号。

## 当前状态梳理

当前代码中存在两套主要输出路径。

### 新 CSR 路径：HContainer/DMR 风格

相关代码：

- `ModuleIO::write_hsr()`：`source/source_io/module_hs/write_HS_R.cpp`
- `ModuleIO::write_hcontainer_csr()`：`source/source_io/module_hs/write_HS_R.cpp`
- `ModuleIO::write_dmr()` 和 `ModuleIO::write_dmr_csr()`：`source/source_io/module_dm/write_dmr.cpp`
- `hamilt::Output_HContainer`：`source/source_lcao/module_hcontainer/output_hcontainer.cpp`

典型文件头：

```text
 --- Ionic Step 1 ---
 # print H matrix in real space H(R)
 1 # number of spin directions
 1 # spin index
 8 # number of localized basis
 13 # number of Bravais lattice vector R

 [UnitCell information]

 #----------------------------------------------------------------------#
 #                               CSR Format                             #
 ...
```

目前 `out_mat_hs2` 走这条路径，`out_dmr` 也基本是同类格式。

### 旧 CSR 路径：legacy sparse-map 风格

相关代码：

- `ModuleIO::output_SR()`
- `ModuleIO::output_TR()`
- `ModuleIO::output_dHR()`
- `ModuleIO::output_dSR()`
- `ModuleIO::save_sparse()`
- `ModuleIO::save_dH_sparse()`
- `ModuleIO::output_single_R()`
- `cal_r_overlap_R::out_rR()`

典型文件头：

```text
STEP: 0
Matrix Dimension of T(R): 13
Matrix number of T(R): 7
```

目前 `calculation=get_s`、`out_mat_t`、`out_mat_dh`、`out_mat_ds`、`out_mat_r` 仍然使用这条路径。

### 需要特别保护的兼容点

- `out_mat_hs2` 当前输出 `--- Ionic Step` 风格的新 CSR 格式。
- `calculation=get_s` 当前通过旧路径输出 `sr_nao.csr`，集成测试引用文件也是旧 `STEP:` 风格。
- `out_mat_t`、`out_mat_dh`、`out_mat_ds`、`out_mat_r` 当前都是旧 `STEP:` 风格。
- `S_async` 当前由 `Overlap::output_SR_async_csr()` 输出，文件头是独立的 `IONIC_STEP:` 风格。
- `write_hsr()` 当前始终写到 `PARAM.globalv.global_out_dir`。虽然旧路径中部分 MD 非追加输出会写到 `global_matrix_dir`，但本重构不能默认改变 H/S 的输出位置。
- `single_R_io.cpp` 当前文本输出精度硬编码为 16 位。即使输入参数文档中有精度项，也不能在结构性重构中顺手改变这个行为。
- `write_HS_R.cpp` 中匿名命名空间的 `sparse_map_to_hcontainer()` 未被调用，并且返回裸指针。应在测试护栏建立后作为清理项移除。

## 非目标

第一轮重构不做以下事情：

- 不改变默认输出格式；
- 不改变文件名；
- 不改变输出目录；
- 不改变 MD append 行为；
- 不改变旧路径默认精度；
- 不把下游可见的旧格式文件默认替换为新格式文件；
- 不删除当前已有调用方仍在使用的函数；
- 不更新数值参考文件，除非只是新增覆盖现有行为的测试引用。

以下问题可以记录，但不建议混入本轮结构性重构：

- 旧 sparse-map 路径是否应该真正使用 `out_mat_t[1]`、`out_mat_dh[1]`、`out_mat_ds[1]`、`out_mat_r[1]` 的精度；
- MD 且 `out_app_flag=false` 时，H/S 是否应和 T/dH/dS/r 一样写入 `OUT.*/matrix/`；
- `get_s` 是否应支持新 HContainer CSR 文件头；
- T/dH/dS/r 是否应支持带 UnitCell 的新 CSR 格式。

这些都属于用户可见行为变化，应单独设计兼容开关或迁移策略。

## 后向兼容契约

动公共 writer 之前，必须先建立兼容性矩阵并用测试锁住。

| 功能 | 当前文件 | 当前文件头 | 当前目录 | 默认行为是否必须保持 |
| --- | --- | --- | --- | --- |
| `out_mat_hs2`, `nspin=1` | `hrs1_nao.csr`, `srs1_nao.csr` | `--- Ionic Step` | `OUT.${suffix}/` | 是 |
| `out_mat_hs2`, `nspin=2` | `hrs1_nao.csr`, `hrs2_nao.csr`, `srs1_nao.csr` | `--- Ionic Step` | `OUT.${suffix}/` | 是 |
| `out_mat_hs2`, `nspin=4` | 当前生成的 H/S 文件 | `--- Ionic Step` | `OUT.${suffix}/` | 是 |
| `calculation=get_s` | `sr_nao.csr` | `STEP:` | `OUT.${suffix}/` | 是 |
| `out_mat_t` | `trs1_nao.csr` 或 MD 步相关文件 | `STEP:` | 当前行为 | 是 |
| `out_mat_dh` | `dhrxs*_nao.csr` 等 | `STEP:` | 当前行为 | 是 |
| `out_mat_ds` | `dsrxs*_nao.csr` 等 | `STEP:` | 当前行为 | 是 |
| `out_mat_r` | `rr.csr` 或 MD 步相关文件 | `STEP:` | 当前行为 | 是 |
| `cal_syns` / async S | `syns_nao.csr` | `IONIC_STEP:` | `OUT.${suffix}/` | 是 |
| `out_dmr` | `dmrs*_nao.csr` | `--- Ionic Step` | `OUT.${suffix}/` | 是 |

凡是会改变上表任意一格的改动，都不能作为“纯重构”合入，应拆成独立行为变更。

## 总体执行顺序

推荐顺序如下：

1. 先冻结现有行为，补 characterization tests。
2. 修复 review 中确认的 P0/P1 缺陷，并补针对性回归测试。
3. 做不改变 API 和输出的低风险清理。
4. 抽取 HContainer CSR 公共 writer。
5. 抽取文件名、目录、append 策略。
6. 在兼容层内重构旧 sparse-map writer。
7. 优化 `Output_HContainer` 的 R 遍历内部实现。
8. 最后再考虑 opt-in 行为改进。

不能倒过来做。尤其不能先改文件名策略、不能先改旧格式精度、不能先改 R 块遍历顺序。

## Phase 0：冻结行为并补测试

### 目的

先让现有行为变成可观察、可回归的契约。后续所有重构都以这些测试为护栏。

### 具体任务

1. 增加文件名 helper 的单元测试：
   - `hsr_gen_fname()`;
   - `dhr_gen_fname()`;
   - 如果后续会共享 DMR 命名逻辑，也覆盖 `dmr_gen_fname()`。
2. 使用小型 synthetic `HContainer` 增加 writer header 测试：
   - H 文件头；
   - S 文件头；
   - DM 文件头；
   - `istep == 0` 创建或覆盖；
   - `istep > 0` append。
3. 增加不同文件头的最小测试：
   - `--- Ionic Step`;
   - `STEP:`;
   - `IONIC_STEP:`。
4. 增加或扩展现有集成参考覆盖：
   - `tests/03_NAO_multik/scf_out_hsr`;
   - `tests/03_NAO_multik/scf_out_hsr_spin4`;
   - `tests/03_NAO_multik/nscf_out_hsr_tr_rr`;
   - `tests/03_NAO_multik/scf_out_dh_t`;
   - `tests/03_NAO_multik/get_s`;
   - `tests/03_NAO_multik/md_out_syns`;
   - `tests/03_NAO_multik/scf_out_dmr_dmk`。
5. 对 MD `out_app_flag=true/false` 的当前路径和文件名行为做一份最小 fixture 或测试说明。

### 验证命令

先跑小范围：

```bash
ctest --test-dir build -N
ctest --test-dir build -R hcontainer -V
ctest --test-dir build -R source_io -V
ctest --test-dir build -R scf_out_hsr -V --timeout 1700
ctest --test-dir build -R get_s -V --timeout 1700
```

### 合入门禁

这个阶段不能改变任何现有输出。如果测试暴露已有问题，只记录问题，不在本阶段修。

### 当前完成状态

已新增 Phase 0 的最小单元测试护栏：

- `source/source_io/test/write_hs_r_compat_test.cpp`;
- `source/source_io/test/CMakeLists.txt` 中的 `MODULE_IO_write_hs_r_compat_test`。

当前覆盖内容：

- `hsr_gen_fname()`、`dhr_gen_fname()`、`dmr_gen_fname()` 文件名契约；
- `write_hcontainer_csr()` 的 H/S 文件头和 append 行为；
- `write_dmr_csr()` 的 DMR 文件头；
- legacy `save_sparse()` 的 `STEP:` 文件头；
- `--- Ionic Step`、`STEP:`、`IONIC_STEP:` 三类文件头样例区分。

已执行并通过：

```bash
g++ -std=c++14 -Isource -Isource/source_base/module_container -c source/source_io/test/write_hs_r_compat_test.cpp -o /tmp/write_hs_r_compat_test.o
cmake -B build-phase0 -DBUILD_TESTING=ON -DCMAKE_CXX_COMPILER=mpicxx
cmake --build build-phase0 --target MODULE_IO_write_hs_r_compat_test -j8
ctest --test-dir build-phase0 -R MODULE_IO_write_hs_r_compat_test -V
```

后续仍建议补充集成级 characterization tests，尤其是 MD append/non-append、binary legacy sparse、`out_mat_r`、`write_HS.hpp` 非 MPI 分支。

## Phase 1：修复明确缺陷与低风险清理，不改输出契约

### 目的

先修复 review 中确认的 P0/P1 级缺陷，再清掉明显无用代码和局部重复，为后续抽公共 writer 降低噪声。

本阶段允许修复未定义行为、无效写入和资源泄漏，但不允许借机改变用户可见输出契约。修复后如果输出从“无法正确写出”变为“按既有格式正确写出”，必须用新增测试说明这是缺陷修复，而不是格式迁移。

### 具体任务

1. 修复 `ModuleIO::save_sparse()` 的 binary header 未定义行为：
   - 将空指针写入改为写入明确的 step 变量；
   - 补 binary sparse 最小单测，锁定 `{step, nlocal, output_R_number}` 布局；
   - 文本 `STEP:` 行保持不变。
2. 修复 `ModuleIO::save_mat()` 非 MPI 文本分支 stream 遮蔽问题：
   - 使用外层 `out_matrix.open(...)`；
   - 补 `ENABLE_MPI=OFF` 或可编译覆盖；
   - 不改变现有文本排列、精度和三角矩阵行为。
3. 修复 `write_Vxc_R()` 的资源泄漏：
   - 优先用栈对象或 `std::unique_ptr`；
   - `vxcs_op_ao` 使用 `std::vector<std::unique_ptr<...>>` 或局部对象；
   - 不改变 Vxc(R) 文件名和 `save_sparse()` 调用行为。
4. 修复 `AngularMomentumCalculator` 日志指针空值风险：
   - 推荐将 `std::ofstream* ptr_log = nullptr` 改为必须传入引用，或在构造函数内部提供 null-safe 路径；
   - 如果改签名影响调用方，可先保持签名但加显式空值检查和失败提示；
   - `kernel()` 中也应先判空再访问 `ofs->is_open()`。
5. 删除 `write_HS_R.cpp` 匿名命名空间中未使用的 `sparse_map_to_hcontainer()`，或移到后续专门分支。
6. 保持现有函数签名不变：
   - 暂不删除未使用的 `kv` 参数；
   - 暂不删除 `binary` 参数；
   - 暂不删除 `output_mat_sparse()` 中的 `p_ham`、`p_dftu`，先记录为 API 清理项；
   - 不改变默认参数。
7. 在安全前提下降低 include 耦合：
   - 能放到 `.cpp` 的重 include 不放在 `.h`；
   - 但不能破坏当前调用方编译。
8. 增加小型内部 helper，但必须保证输出字节级一致：
   - `open_text_or_append_file()`；
   - `should_append_matrix_step(istep, append)`；
   - `serial_output_rank()`。

### 验证命令

```bash
cmake --build build -j8
ctest --test-dir build -R hcontainer -V
ctest --test-dir build -R source_io -V
ctest --test-dir build -R MODULE_IO_write_hs_r_compat_test -V
```

如果修改 `save_mat()` 非 MPI 分支，还应补一次非 MPI 配置构建：

```bash
cmake -B build-nompi -DBUILD_TESTING=ON -DENABLE_MPI=OFF
cmake --build build-nompi --target <target-covering-save_mat> -j8
```

### 合入门禁

所有 characterization tests 通过，测试覆盖到的输出必须字节级一致。

### 当前执行状态

已完成第一批 Phase 1 修复：

- 修复 `ModuleIO::save_sparse()` binary header 从空指针写入的问题，改为写入明确的 step 值；
- 修复 `ModuleIO::save_mat()` 非 MPI 文本分支中 `std::ofstream` 局部变量遮蔽导致外层 stream 未打开的问题；
- 修复 `write_Vxc_R()` 中 `elecstate::Potential` 和 `Veff` 对象的资源泄漏，改为栈对象管理；
- 修复 `AngularMomentumCalculator` 日志指针允许为空但直接解引用的问题；
- 删除 `write_HS_R.cpp` 中未使用且返回裸指针的 `sparse_map_to_hcontainer()` helper；
- 为 legacy binary sparse 输出新增最小头部布局回归测试。

已执行并通过：

```bash
g++ -std=c++14 -Isource -Isource/source_base/module_container -c source/source_io/test/write_hs_r_compat_test.cpp -o /tmp/write_hs_r_compat_test.o
cmake --build build-phase0 --target MODULE_IO_write_hs_r_compat_test -j8
ctest --test-dir build-phase0 -R MODULE_IO_write_hs_r_compat_test -V
cmake --build build-phase0 --target MODULE_IO_cal_pLpR_test -j8
ctest --test-dir build-phase0 -R MODULE_IO_cal_pLpR_test -V
cmake --build build-phase0 --target io_advanced -j8
cmake -B build-phase1-nompi -DBUILD_TESTING=OFF -DENABLE_MPI=OFF
cmake --build build-phase1-nompi --target io_advanced -j8
cmake --build build-phase0 --target MODULE_IO_single_R_test -j8
ctest --test-dir build-phase0 -R MODULE_IO_single_R_test -V
cmake --build build-phase0 --target MODULE_IO_sparse_matrix_test -j8
ctest --test-dir build-phase0 -R MODULE_IO_sparse_matrix_test -V
```

注意：`MODULE_IO_cal_pLpR_test` 和 `MODULE_IO_single_R_test` 在 sandbox 内会因 OpenMPI/PMIx socket 初始化失败而无法运行，已在非 sandbox 环境重跑并通过。非 MPI + `BUILD_TESTING=ON` 配置当前会失败在既有 `source_hsolver/test` 目标依赖缺失的 `MPI::MPI_CXX`，因此本阶段使用 `BUILD_TESTING=OFF` 覆盖生产代码非 MPI 编译。

## Phase 2：抽取 HContainer CSR 公共 writer

### 目的

把 H/S 和 DMR 已经高度重复的 CSR 写出逻辑统一起来。这个阶段收益最高，也最适合先做。

### 建议内部结构

新增一个内部 options 结构，例如：

```cpp
struct HContainerCsrWriteOptions
{
    std::string filename;
    std::string matrix_name;      // H, S, DM 等
    std::string matrix_symbol;    // H(R), S(R), DM(R) 等
    int precision = 8;
    int istep = 0;
    int ispin = 0;
    int nspin = 1;
    bool append_file = false;
    double sparse_threshold = 1e-10;
};
```

如果维护者希望少增文件，可以先把 options 和 writer 放在 `.cpp` 的匿名命名空间里。接口稳定后再拆成独立文件。

### 具体任务

1. 从以下函数中抽公共逻辑：
   - `write_hcontainer_csr()`；
   - `write_dmr_csr()`。
2. 保留现有函数作为 wrapper：
   - `write_hcontainer_csr()` 仍然可调用；
   - `write_dmr_csr()` 仍然可调用；
   - wrapper 输出必须完全一致。
3. 第一轮不要强行合并 `Overlap::output_SR_async_csr()`，因为它的 `IONIC_STEP:` 头是独立格式。
4. 增加打开文件失败检查：
   - 如果文件是必须输出，使用 `ModuleBase::WARNING_QUIT`；
   - 如果附近已有“警告后继续”的历史行为，则保持相同失败策略。

### 验证命令

```bash
ctest --test-dir build -R hcontainer -V
ctest --test-dir build -R dmr -V
ctest --test-dir build -R scf_out_hsr -V --timeout 1700
ctest --test-dir build -R scf_out_dmr_dmk -V --timeout 1700
```

### 合入门禁

H/S 和 DMR CSR 文件必须和重构前字节级一致。

## Phase 3：集中管理文件名、目录和 append 策略

### 目的

把“写到哪里、叫什么名、是否 append、离子步怎么进文件名”从矩阵写出函数里剥离出来。

### 建议内部结构

新增内部策略对象，例如：

```cpp
struct RealspaceMatrixFilePolicy
{
    std::string directory;
    std::string basename;
    bool append_to_existing_file;
};
```

### 具体任务

1. 抽取实空间矩阵文件名策略 helper。
2. 默认保持现有行为：
   - `write_hsr()` 的 H/S 仍写当前目录；
   - legacy sparse 输出仍写当前目录；
   - `syns_nao.csr` 不改名；
   - MD append/non-append 的现有文件名后缀不改。
3. 只替换重复拼接逻辑，不改变结果字符串。
4. 保持现有步数编号：
   - `--- Ionic Step` 显示 `istep + 1`；
   - legacy `STEP:` 保持当前值；
   - 文件名里的 `g` 后缀保持当前每类输出自己的 0-based 或 1-based 行为。
5. 如果发现行为不一致，先写进兼容性矩阵，不在本阶段修正。

### 验证命令

```bash
ctest --test-dir build -R md_out_syns -V --timeout 1700
ctest --test-dir build -R nscf_out_hsr_tr_rr -V --timeout 1700
```

### 合入门禁

输出路径和文件名不能变。

## Phase 4：在兼容层内重构 legacy sparse-map writer

### 目的

减少 `save_sparse()`、`save_dH_sparse()`、`output_single_R()` 中的重复逻辑，但保持旧格式输出完全不变。

### 建议内部结构

新增 legacy options：

```cpp
struct LegacySparseCsrOptions
{
    std::string filename;
    std::string label;
    int istep = 0;
    bool binary = false;
    bool reduce = true;
    double sparse_threshold = 1e-10;
    int precision = 16; // 保持 legacy 当前默认行为
};
```

### 具体任务

1. 给 `output_single_R()` 增加可选 precision 参数，但默认值保持 16。
2. 抽取非零元统计逻辑：
   - 单个 sparse map 的 nonzero count；
   - MPI reduce；
   - 输出 R 块数量统计；
   - R 块头写出。
3. 重构 dH/dS 的 x/y/z 三轴重复逻辑：
   - 用 axis descriptor 表达 x/y/z；
   - 避免三份文件打开、文件头写出、CSR 块写出代码。
4. 保持旧文件头、空行、row pointer、二进制写出类型不变。

### 验证命令

```bash
ctest --test-dir build -R scf_out_dh_t -V --timeout 1700
ctest --test-dir build -R get_s -V --timeout 1700
ctest --test-dir build -R nscf_out_hsr_tr_rr -V --timeout 1700
```

### 合入门禁

旧 sparse 输出参考文件必须字节级一致。

## Phase 5：优化 Output_HContainer 的 R 遍历

### 目的

`Output_HContainer::write()` 当前会先求 R 的 min/max，然后三重循环扫描整个 R 立方范围。对于稀疏 R 集合，这可能做很多无效 `find_R()`。可以优化，但必须保证 R 块顺序不变。

### 具体任务

1. 先增加 R 块顺序测试，覆盖代表性 H/S 文件。
2. 新增内部方法：收集真实存在的 R 向量，并排序成当前输出顺序。
3. 排序结果和当前三重循环顺序一致后，再替换扫描逻辑。
4. 保持 `write_empty` 的行为不变。

### 验证命令

```bash
ctest --test-dir build -R hcontainer -V
ctest --test-dir build -R scf_out_hsr -V --timeout 1700
ctest --test-dir build -R scf_out_hsr_spin4 -V --timeout 1700
```

### 合入门禁

CSR R 块顺序必须不变。

## Phase 6：可选行为改进，必须 opt-in

这些不属于纯重构，只能在前面阶段完成后独立做。

可选项：

1. 让旧 sparse 输出真正使用输入精度：
   - `out_mat_t 1 5`;
   - `out_mat_dh 1 5`;
   - `out_mat_ds 1 5`;
   - `out_mat_r 1 5`。
2. 为 `get_s` 增加新 HContainer CSR 格式输出选项。
3. 为 T/dH/dS/r 增加带 UnitCell 信息的新 CSR 格式输出选项。
4. 为 H/S 在 MD 且 `out_app_flag=false` 时提供写入 `OUT.*/matrix/` 的 opt-in 行为。
5. 增加统一 CSR reader，兼容读取：
   - `STEP:`;
   - `--- Ionic Step`;
   - `IONIC_STEP:`。

每个可选项都必须满足：

- 默认旧行为不变；
- 有新参数或明确 opt-in 路径；
- 有迁移文档；
- 有旧行为和新行为两套测试；
- 有 release note。

## 推荐 PR 拆分

### PR 1：行为冻结和测试护栏

范围：

- 只加测试；
- 只记录当前兼容性行为；
- 不改生产逻辑。

风险：低。

合入标准：

- 现有参考不变；
- CI 通过。

### PR 2：修复 P0/P1 缺陷和低风险清理

范围：

- 修复 `save_sparse()` binary header 空指针写入；
- 修复 `save_mat()` 非 MPI 文本分支 stream 遮蔽；
- 修复 `write_Vxc_R()` 资源泄漏；
- 修复 `AngularMomentumCalculator` 日志指针空值风险；
- 删除或隔离未使用 helper；
- 不改变输出契约和公开行为。

风险：中。虽然是缺陷修复，但触碰基础输出路径。

合入标准：

- 新增缺陷回归测试通过；
- Phase 0 characterization tests 通过；
- 默认文本输出和文件名不变；
- binary 修复有明确测试说明。

### PR 3：抽取 HContainer CSR writer

范围：

- 新增公共 HContainer CSR writer；
- `write_hcontainer_csr()` 和 `write_dmr_csr()` 保留为 wrapper；
- 增加打开文件失败检查。

风险：中。会触碰 H/S 和 DMR 输出。

合入标准：

- H/S 和 DMR 输出字节级一致；
- 当前调用方无需修改或只做等价迁移。

### PR 4：集中管理文件名和输出策略

范围：

- 抽取目录、文件名、append、step 后缀策略；
- 保持所有结果字符串不变。

风险：中，尤其是 MD 场景。

合入标准：

- MD append/non-append 测试通过；
- 文件路径和文件名不变。

### PR 5：legacy sparse writer 清理

范围：

- 重构 `save_sparse()`、`save_dH_sparse()`、`output_single_R()`；
- 引入 options 和 axis descriptor；
- 保持旧格式。

风险：中高。旧路径覆盖 T、dH、dS、r、get_s。

合入标准：

- legacy 输出参考文件不变。

### PR 6：Output_HContainer R 遍历优化

范围：

- 用排序后的真实 R 列表替换三重范围扫描；
- 保证顺序不变。

风险：中。

合入标准：

- R 块顺序测试通过；
- H/S 参考文件不变；
- 性能不回退。

### PR 7：opt-in 行为增强

范围：

- 新格式输出；
- 精度修正；
- 目录策略修正；
- reader 增强。

风险：视具体项而定。

合入标准：

- 默认旧行为仍有测试；
- 新行为必须显式开启；
- 文档齐全。

## 详细依赖关系

必须按以下依赖推进：

1. 没有 characterization tests，不做公共 writer 抽取。
2. 没有针对性回归测试，不修 P0/P1 输出缺陷。
3. 没有 P0/P1 缺陷修复，不做大规模 writer 抽取。
4. 没有 writer 抽取，不做目录策略统一。
5. 没有 legacy 输出测试，不动 `save_sparse()` 和 `output_single_R()` 的结构。
6. 没有 R 块顺序测试，不动 `Output_HContainer::write()` 的遍历方式。
7. 没有 opt-in 设计，不做用户可见行为变化。

## 建议的最终代码结构

可选拆分如下：

```text
source/source_io/module_hs/
  write_HS_R.cpp
  write_HS_R.h
  realspace_csr_writer.cpp       # 可选新增：HContainer CSR 公共 writer
  realspace_csr_writer.h         # 可选新增：内部 writer 声明
  legacy_sparse_csr_writer.cpp   # 可选新增：legacy sparse 兼容 writer
  legacy_sparse_csr_writer.h     # 可选新增：内部声明

source/source_io/module_dm/
  write_dmr.cpp                  # 保留 wrapper，复用公共 writer
```

如果维护者不想一次新增太多文件，可以先在 `.cpp` 匿名命名空间里放 helper。等结构稳定后再拆文件。

## 最小测试矩阵

每个生产逻辑 PR 至少运行：

```bash
cmake --build build -j8
ctest --test-dir build -R hcontainer -V
ctest --test-dir build -R source_io -V
ctest --test-dir build -R scf_out_hsr -V --timeout 1700
ctest --test-dir build -R scf_out_hsr_spin4 -V --timeout 1700
ctest --test-dir build -R nscf_out_hsr_tr_rr -V --timeout 1700
ctest --test-dir build -R scf_out_dh_t -V --timeout 1700
ctest --test-dir build -R get_s -V --timeout 1700
ctest --test-dir build -R scf_out_dmr_dmk -V --timeout 1700
ctest --test-dir build -R md_out_syns -V --timeout 1700
```

如果当前构建启用了 MPI，还应补跑 MPI 相关代表性 CTest。

## 代码风格与质量治理清单

以下事项不一定都要在同一个 PR 完成，但每次触碰对应文件时应顺手向目标状态推进。

### 所有权和资源管理

- 优先使用栈对象、`std::vector`、`std::unique_ptr`，减少 `new/delete` 和 `new[]/delete[]`。
- 对跨函数传入的裸指针明确所有权语义：
  - 非 owning 指针用注释或 `const` 表达；
  - 可空指针必须在入口检查；
  - 必须存在的对象优先使用引用。
- 文件流打开后必须检查状态：
  - 对必须生成的输出文件，失败时应明确报错；
  - 对历史上允许失败后继续的路径，保持旧策略但补注释说明。

### 输出格式和 I/O

- `STEP:`、`--- Ionic Step`、`IONIC_STEP:` 三类头部属于兼容格式，不混用、不合并。
- 文本和二进制 writer 应共享“打开文件、写 header、写 R 块头”的小 helper，但不要共享到改变空行、精度、类型宽度或顺序。
- 避免用固定临时文件名拼接中间结果；如果暂时保留，文件名必须包含 rank 和足够的调用上下文。
- `output_single_R()` 后续建议改为内存中缓存 CSR indices，替代 `temp_sparse_indices.dat`。

### API 和命名

- 新增参数优先收敛到 options struct，避免函数签名继续变长。
- 对 `const int&`、`const bool&`、`const double&` 这类基础类型参数，新增代码优先按值传递；旧签名可在兼容 wrapper 中保留。
- include guard 不使用 `__NAME__` 形式，新增头文件使用普通项目风格宏或 `#pragma once`，按项目约定执行。
- `output_mat_sparse()` 中的未使用参数要么恢复语义，要么通过 wrapper 逐步迁移，避免继续扩大无效 API。

### 复杂度控制

- `cal_r_overlap_R.cpp` 后续拆分目标：
  - orbital table 构造；
  - nonlocal projector 准备；
  - 单个 R 的 r(R) 三方向计算；
  - r(R) 文件输出。
- `write_vxc.hpp`、`write_vxc_lip.hpp`、`write_vxc_r.hpp` 后续拆分目标：
  - Potential/Veff 构造；
  - AO/MO 矩阵变换；
  - orbital energy 计算；
  - 文件输出。
- 避免把长业务逻辑继续放入 `.hpp`，模板必须留在头文件时，也应把非模板 helper 下沉到 `.cpp`。

### 测试策略

- 每修一个 review 风险点，至少补一个直接失败/通过的最小单测。
- 涉及 MPI rank0 写文件时，至少覆盖 `reduce=true` 和 `reduce=false` 的 writer 行为。
- 涉及 `ENABLE_MPI=OFF` 的分支时，必须有非 MPI 配置的编译或测试记录。
- 涉及 binary 输出时，必须检查字节布局，不只检查文件存在。

## Review Checklist

每个 PR 都要检查：

- 现有输出文件名不变；
- 现有输出目录不变；
- 现有文本文件头不变；
- 现有二进制布局不变；
- Ionic step 编号不变；
- spin index 和 spin 文件数量不变；
- CSR R 块顺序不变；
- CSR row pointer 的文本/二进制类型不变；
- 默认精度不变；
- 当前调用方仍能编译；
- 新 helper 有测试；
- binary writer 不从空指针写数据；
- file stream 打开后再写，打开失败策略明确；
- 新增或保留的 raw pointer 有清晰所有权；
- `new/delete` 优先替换为 RAII，无法替换时说明原因；
- 临时文件名不会跨 rank、跨并发测试冲突；
- 没有无理由更新旧参考文件；
- 新行为如果存在，必须 opt-in，且文档说明默认兼容。

## 回滚策略

每个 PR 应能独立回滚：

- PR 1 只有测试，回滚无生产影响。
- PR 2 只修明确缺陷和低风险清理，回滚后恢复旧实现。
- PR 3 保留旧 wrapper，回滚后恢复直接实现。
- PR 4 只集中路径策略，回滚后恢复原地拼接。
- PR 5 只改 legacy 内部实现，参考文件能捕获问题。
- PR 6 只改 `Output_HContainer` 遍历方式，可独立回滚。
- PR 7 是 opt-in 行为增强，应能通过关闭新选项或回滚 PR 恢复旧行为。

## 完成标准

满足以下条件才认为重构完成：

- `write_HS_R.cpp` 不再直接承担过多不相关职责；
- H/S 和 DMR 共享同一套 HContainer CSR writer；
- legacy sparse-map 输出有集中 helper，但旧格式保持不变；
- 文件名、目录、append、step 策略显式且有测试；
- 默认行为下集成参考文件不变；
- 如有新增行为，必须 opt-in 并有迁移文档。
