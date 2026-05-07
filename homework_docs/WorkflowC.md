# Workflow C: SIMD 向量化与缓存复用实施清单

## 1. 目标与范围

本工作流聚焦 ABACUS 平面波基组模块中的两类优化：

1. **gather/scatter 重排循环的 SIMD 向量化**
   - 目标是提升 `PW_Basis::gatherp_scatters` 和 `PW_Basis::gathers_scatterp` 中局部重排段的单核吞吐。
   - 优先优化规则、连续、可预测的数组拷贝循环。
   - 尽量保持原有数据布局和接口不变，避免牵连 FFT、MPI 或 Gamma Only 路径。

2. **G 矢量与 `|G+K|^2` 等不可变数据的缓存复用**
   - 目标是减少 `PW_Basis_K` 在初始化和重复使用中的重复计算。
   - 重点对象包括 `kvec_c`、`gk2`、`gcar`、`igl2isz_k`、`igl2ig_k`、`ig2ixyz_k`。
   - 缓存必须定义清楚构建时机、失效条件和只读访问边界。

### 不在本工作流内处理的内容

- MPI 非阻塞通信替换
- Gamma Only 扩展与半谱 FFT
- `count_pw_st` 并行化
- FFT 通信-计算重叠

---

## 2. 预备阅读与代码定位

在动手修改之前，先把以下文件和函数看清楚：

### 2.1 gather/scatter 路径

- `source/source_basis/module_pw/pw_gatherscatter.h`
  - `PW_Basis::gatherp_scatters`
  - `PW_Basis::gathers_scatterp`

需要确认：
- 哪些循环只是在搬运数据
- 哪些循环在做索引转换
- 哪些循环与 MPI 通信绑定在一起
- 哪些循环在 `poolnproc == 1` 时可以独立优化

### 2.2 缓存路径

- `source/source_basis/module_pw/pw_basis_k.h`
- `source/source_basis/module_pw/pw_basis_k.cpp`
  - `PW_Basis_K::initparameters`
  - `PW_Basis_K::setuptransform`
  - `PW_Basis_K::setupIndGk`
  - `PW_Basis_K::collect_local_pw`
  - `PW_Basis_K::getgk2`
  - `PW_Basis_K::getgcar`
  - `PW_Basis_K::getgdirect`
  - `PW_Basis_K::getgpluskcar`
  - `PW_Basis_K::get_gk2_data`
  - `PW_Basis_K::get_gcar_data`

需要确认：
- 哪些成员变量是“初始化后不再改变”的
- 哪些成员变量会在重新设定网格、k 点、截断能、分发方式时失效
- 哪些接口当前是只读语义，哪些地方可能误用成可写共享状态

### 2.3 现有测试

- `source/source_basis/module_pw/test_serial/pw_basis_k_test.cpp`
- `source/source_basis/module_pw/test/test-other.cpp`
- `source/source_basis/module_pw/test/test1-4.cpp`
- `source/source_basis/module_pw/test/test4-4.cpp`

这些测试里已经覆盖了：
- `getgk2`、`getgcar`、`getgdirect` 的数值一致性
- FFT 往返的基础正确性
- 多 k 点场景下的基础构造与转换

---

## 3. 工作拆分原则

为了避免一次性改动太大，建议把工作流 C 拆成两个逻辑子任务：

### 子任务 C1: gather/scatter SIMD 化

目标：在不改变输出语义的前提下，把重排循环改成编译器更容易向量化的形式。

### 子任务 C2: 缓存语义整理

目标：明确 `PW_Basis_K` 中缓存的生命周期、复用边界和读取方式，避免重复计算和陈旧数据。

这两个子任务可以并行推进，但最终提交时建议拆成两个提交，便于回归定位。

---

## 4. 详细实施步骤

## 4.1 C1: gather/scatter 的 SIMD 向量化

### Step 1: 找出纯拷贝热区

先在 `pw_gatherscatter.h` 中区分以下几类循环：

- **纯拷贝循环**
  - 例如把 `inp[iz]` 逐元素拷贝到 `outp[iz]`
  - 这类循环最适合 `omp simd`

- **索引准备循环**
  - 例如通过 `istot2ixy`、`startg`、`startr`、`startz`、`numg`、`numr` 计算段起始地址
  - 这类循环一般不直接向量化，但可以作为拷贝循环的前置准备

- **通信相关循环**
  - `MPI_Alltoallv` 前后的打包/解包
  - 这部分不和 SIMD 强行绑定，只做本地重排段优化

### Step 2: 优先处理 `poolnproc == 1` 分支

这条路径最稳定，适合先做第一版优化：

- `gatherp_scatters` 中：
  - `for (int is = 0; is < this->nst; ++is)` 内的逐 `iz` 拷贝
- `gathers_scatterp` 中：
  - 先清零 `out`
  - 再把 `in` 写入 `out`

实施要求：
- 保持结果与原实现一致
- 内层循环尽量只保留线性访问
- 如果编译器识别失败，再考虑增加局部指针、临时变量、`#pragma omp simd` 或更明确的循环结构

### Step 3: 再处理 `__MPI` 分支的局部重排

`__MPI` 分支中有两类本地数组重排：

- 通信前打包到连续缓冲区
- 通信后从连续缓冲区解包回目标布局

这里的原则是：

1. **SIMD 只覆盖本地 contiguous copy**
2. **不把 MPI 调用本身与向量化混在同一个循环里**
3. **优先保持当前内存布局**，不轻易改通信协议

### Step 4: 处理可能阻碍向量化的因素

如果编译器无法生成向量化，通常是以下原因：

- 循环中存在不规则索引
- 指针别名不明确
- 目标数据和源数据有潜在重叠风险
- 循环体太复杂，混入了索引计算或条件判断

对应处理方式：

- 把索引提前计算到局部变量
- 将指针绑定到局部常量指针
- 把复杂逻辑拆到循环外部
- 保持循环体只有赋值语句

### Step 5: 保持接口不变

C1 阶段不改这些接口：

- `PW_Basis::gatherp_scatters`
- `PW_Basis::gathers_scatterp`

也不改它们的输入输出形状和调用顺序。优化必须是“同接口、同语义、低侵入”。

---

## 4.2 C2: 缓存语义整理

### Step 1: 明确缓存对象清单

`PW_Basis_K` 中需要整理语义的对象包括：

- `kvec_c`
- `gk2`
- `gcar`
- `igl2isz_k`
- `igl2ig_k`
- `ig2ixyz_k`
- 以及它们对应的 CPU/GPU 镜像指针

这里要特别区分：

- **逻辑缓存**：算法上可复用的数据
- **设备镜像**：为 CPU/GPU 访问而存在的复制版本

### Step 2: 定义缓存的构建时机

建议把以下时机视为缓存构建点：

- `initparameters(...)`
  - k 点、截断能、分发方式变化时，需要重建相关数据

- `setuptransform()`
  - FFT 网格或分布映射变化时，需要重建 `igl2isz_k`、`igl2ig_k`、`ig2ixyz_k`

- `collect_local_pw(...)`
  - 当需要更新局部平面波属性时，重建 `gk2`、`gcar`

### Step 3: 定义失效条件

缓存失效的典型条件：

- 网格变化：`initgrids` 被重新调用
- k 点变化：新的 `kvec_d` 输入
- 截断能变化：`gk_ecut` 变化
- 分发方式变化：`distribution_type` 或 `xprime` 变化
- FFT 布局变化：`setuptransform` 重新执行

如果某个对象在这些条件下会改变，就不能当作长期稳定缓存。

### Step 4: 把现有一次性生成逻辑整理成显式缓存语义

当前 `collect_local_pw` 已经在做一次性生成，但还缺“它是什么缓存、什么时候失效、谁能读谁不能写”的说明。

建议落实为：

- `getgk2`、`getgcar`、`getgdirect`、`getgpluskcar` 只读访问
- `get_gk2_data`、`get_gcar_data` 仅作为底层只读数据出口
- 不允许上层直接修改这些缓存
- 若未来确实需要重建，提供明确的 `invalidate` 或 `rebuild` 语义，而不是依赖外部手动清空指针

### Step 5: 处理设备侧镜像

如果当前代码需要维护 CPU/GPU 镜像，要明确：

- 哪些指针是 host ownership
- 哪些指针只是 view
- 哪些数据由 `resmem_*` / `syncmem_*` / `castmem_*` 管理

原则上，缓存语义应该先定义在逻辑层，再映射到设备镜像，不要反过来让设备内存结构决定算法边界。

---

## 5. 推荐的代码修改顺序

### 顺序 A: 先 SIMD，后缓存

适用于希望先快速看到性能收益的情况。

1. 修改 `pw_gatherscatter.h`
2. 通过现有测试验证结果一致
3. 再整理 `pw_basis_k.cpp` 的缓存语义
4. 最后补测试与文档

### 顺序 B: 先缓存，后 SIMD

适用于希望先把生命周期和边界定义清楚的情况。

1. 先整理 `pw_basis_k.h/.cpp`
2. 补 `collect_local_pw` 的复用说明
3. 再处理 `pw_gatherscatter.h`
4. 最后统一回归测试

### 推荐

如果团队里并行开发压力较大，建议采用**顺序 B**，因为缓存边界定义清楚后，后续 SIMD 改动更容易判断是否触碰到共享状态。

---

## 6. 具体改动点清单

### 6.1 `pw_gatherscatter.h`

建议重点检查并可能修改以下内容：

- `poolnproc == 1` 分支中的内层拷贝循环
- `__MPI` 分支中通信前后的一维连续拷贝循环
- `#pragma omp parallel for` 与 `#pragma omp simd` 的组合位置
- 循环体是否可以只保留简单赋值
- 是否需要将部分索引预计算为局部变量

### 6.2 `pw_basis_k.h`

建议整理：

- `gk2`、`gcar`、`kvec_c` 的语义说明
- 设备镜像和 host 缓存之间的关系
- `getgk2`、`getgcar` 等接口的只读约定
- 是否需要增加 `invalidate` / `rebuild` 类接口声明

### 6.3 `pw_basis_k.cpp`

建议重点关注：

- `initparameters(...)`
  - 何时重置相关成员
  - 何时清空旧缓存

- `setuptransform()`
  - 何时依赖网格变化重新生成映射

- `setupIndGk()`
  - `npwk`、`igl2isz_k`、`igl2ig_k` 的构建边界

- `collect_local_pw(...)`
  - `gk2`、`gcar` 的生成逻辑
  - 是否需要补充缓存状态标记

- `get_ig2ixyz_k()`
  - GPU 路径相关的索引缓存是否与当前设计一致

---

## 7. 测试计划

### 7.1 正确性测试

必须通过以下几类测试：

1. **FFT 往返一致性**
   - 优化前后结果一致
   - 误差保持在现有测试允许范围内

2. **缓存数值一致性**
   - `getgk2` 和 `getgcar` 返回结果与当前实现一致
   - `getgdirect` 和 `getgpluskcar` 结果与几何关系一致

3. **重复初始化一致性**
   - 重复调用 `initparameters`、`setuptransform`、`collect_local_pw` 后不出现陈旧数据

4. **多 k 点回归**
   - 确保不同 k 点下的缓存数组索引没有串位

### 7.2 重点测试文件

优先复用：

- `source/source_basis/module_pw/test_serial/pw_basis_k_test.cpp`
- `source/source_basis/module_pw/test/test-other.cpp`
- `source/source_basis/module_pw/test/test1-4.cpp`
- `source/source_basis/module_pw/test/test4-4.cpp`

### 7.3 建议补充的测试点

如果现有测试不够明确，建议补充：

- gather/scatter 前后数组完全一致的回归测试
- 缓存重建后指针语义和数值语义都正确的测试
- 对 `get_gk2_data` / `get_gcar_data` 的只读访问验证
- 对重复调用 `collect_local_pw` 后结果不变化的测试

---

## 8. 性能验证方法

### 8.1 SIMD 收益验证

测量对象：

- `gatherp_scatters`
- `gathers_scatterp`

建议指标：

- 单次调用耗时
- 大小规模扩展时的带宽表现
- OpenMP 线程数变化下的加速比

### 8.2 缓存收益验证

测量对象：

- `collect_local_pw`
- `setupIndGk`
- `getgk2` / `getgcar` 的重复访问场景

建议指标：

- 初始化耗时
- 重复初始化耗时
- 内存分配次数是否减少
- 是否避免重复计算 `G+K` 模长平方

### 8.3 建议对比方式

- 优化前 vs 优化后
- 单线程 vs 多线程
- 小体系 vs 大体系

如果有 benchmark 框架，建议单独记录：
- `distributeg`
- `real2recip`
- `recip2real`
- `gather/scatter`
- `collect_local_pw`

---

## 9. 验收标准

工作流 C 完成时，至少满足以下条件：

1. `gather/scatter` 重排结果与原实现一致
2. FFT 往返结果不引入新的数值偏差
3. `gk2`、`gcar`、`kvec_c` 等缓存对象的生命周期和失效条件清楚
4. 重复初始化后不会读取旧缓存
5. 新增或修改的测试全部通过
6. 如果有性能基准，能看到 gather/scatter 或缓存复用至少一项有实际收益

---

## 10. 风险与回退方案

### 10.1 风险点

- SIMD 优化可能被不规则索引限制，收益不明显
- 缓存语义若定义不清，容易与后续 Gamma Only 或 GPU 路径冲突
- 过度改动循环结构可能引入边界错误
- 设备侧镜像和 host 缓存如果没有分清，可能出现悬空指针或重复释放

### 10.2 回退原则

- 保留原始接口
- 先局部改造，再扩展语义
- 单次提交只解决一个主题
- 如果某个优化难以通过测试，优先回退到更保守的实现

---

## 11. 建议提交拆分

如果最终要做成可审查的代码提交，建议拆成两次：

### Commit 1: gather/scatter SIMD

内容：
- 仅修改 `pw_gatherscatter.h`
- 增加或调整相关测试
- 证明重排结果不变

### Commit 2: 缓存语义整理

内容：
- 修改 `pw_basis_k.h` / `pw_basis_k.cpp`
- 明确缓存生命周期和只读访问边界
- 补充重复初始化相关测试

---

## 12. 最终交付物

工作流 C 最终应交付：

- 一份 SIMD 向量化实现
- 一份缓存语义整理实现
- 对应的测试补充
- 性能和正确性验证记录
- 一份简短的改动说明，说明做了什么、收益在哪、风险如何控制

---

## 13. 执行顺序建议

如果要马上开始做，建议按下面顺序执行：

1. 读 `pw_gatherscatter.h`
2. 读 `pw_basis_k.cpp`
3. 先改本地纯拷贝循环
4. 跑测试确认不回归
5. 再整理缓存语义
6. 再补测试
7. 最后做性能统计和文档收尾
