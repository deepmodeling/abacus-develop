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
2. 做不改变 API 和输出的低风险清理。
3. 抽取 HContainer CSR 公共 writer。
4. 抽取文件名、目录、append 策略。
5. 在兼容层内重构旧 sparse-map writer。
6. 优化 `Output_HContainer` 的 R 遍历内部实现。
7. 最后再考虑 opt-in 行为改进。

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

## Phase 1：低风险清理，不改 API 和输出

### 目的

先清掉明显无用代码和局部重复，为后续抽公共 writer 降低噪声。

### 具体任务

1. 删除 `write_HS_R.cpp` 匿名命名空间中未使用的 `sparse_map_to_hcontainer()`，或移到后续专门分支。
2. 保持现有函数签名不变：
   - 暂不删除未使用的 `kv` 参数；
   - 暂不删除 `binary` 参数；
   - 不改变默认参数。
3. 在安全前提下降低 include 耦合：
   - 能放到 `.cpp` 的重 include 不放在 `.h`；
   - 但不能破坏当前调用方编译。
4. 增加小型内部 helper，但必须保证输出字节级一致：
   - `open_text_or_append_file()`；
   - `should_append_matrix_step(istep, append)`；
   - `serial_output_rank()`。

### 验证命令

```bash
cmake --build build -j8
ctest --test-dir build -R hcontainer -V
ctest --test-dir build -R source_io -V
```

### 合入门禁

所有 characterization tests 通过，测试覆盖到的输出必须字节级一致。

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

### PR 2：抽取 HContainer CSR writer

范围：

- 新增公共 HContainer CSR writer；
- `write_hcontainer_csr()` 和 `write_dmr_csr()` 保留为 wrapper；
- 增加打开文件失败检查。

风险：中。会触碰 H/S 和 DMR 输出。

合入标准：

- H/S 和 DMR 输出字节级一致；
- 当前调用方无需修改或只做等价迁移。

### PR 3：集中管理文件名和输出策略

范围：

- 抽取目录、文件名、append、step 后缀策略；
- 保持所有结果字符串不变。

风险：中，尤其是 MD 场景。

合入标准：

- MD append/non-append 测试通过；
- 文件路径和文件名不变。

### PR 4：legacy sparse writer 清理

范围：

- 重构 `save_sparse()`、`save_dH_sparse()`、`output_single_R()`；
- 引入 options 和 axis descriptor；
- 保持旧格式。

风险：中高。旧路径覆盖 T、dH、dS、r、get_s。

合入标准：

- legacy 输出参考文件不变。

### PR 5：Output_HContainer R 遍历优化

范围：

- 用排序后的真实 R 列表替换三重范围扫描；
- 保证顺序不变。

风险：中。

合入标准：

- R 块顺序测试通过；
- H/S 参考文件不变；
- 性能不回退。

### PR 6：opt-in 行为增强

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
2. 没有 writer 抽取，不做目录策略统一。
3. 没有 legacy 输出测试，不动 `save_sparse()` 和 `output_single_R()`。
4. 没有 R 块顺序测试，不动 `Output_HContainer::write()` 的遍历方式。
5. 没有 opt-in 设计，不做用户可见行为变化。

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
- 没有无理由更新旧参考文件；
- 新行为如果存在，必须 opt-in，且文档说明默认兼容。

## 回滚策略

每个 PR 应能独立回滚：

- PR 1 只有测试，回滚无生产影响。
- PR 2 保留旧 wrapper，回滚后恢复直接实现。
- PR 3 只集中路径策略，回滚后恢复原地拼接。
- PR 4 只改 legacy 内部实现，参考文件能捕获问题。
- PR 5 只改 `Output_HContainer` 遍历方式，可独立回滚。
- PR 6 是 opt-in 行为增强，应能通过关闭新选项或回滚 PR 恢复旧行为。

## 完成标准

满足以下条件才认为重构完成：

- `write_HS_R.cpp` 不再直接承担过多不相关职责；
- H/S 和 DMR 共享同一套 HContainer CSR writer；
- legacy sparse-map 输出有集中 helper，但旧格式保持不变；
- 文件名、目录、append、step 策略显式且有测试；
- 默认行为下集成参考文件不变；
- 如有新增行为，必须 opt-in 并有迁移文档。
