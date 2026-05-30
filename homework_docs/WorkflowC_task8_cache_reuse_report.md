# WorkflowC 题目 8：缓存复用优化总结报告

## 1. 题目背景与优化目标

本报告对应 `homework_docs/01_plane_wave.md` 中 WorkflowC 的第 8 题“平面波预计算与缓存优化”。题目的核心目标是识别平面波模块中可重复利用的中间结果，把原本在自洽场迭代中反复发生的重复构造，改造成“按需构建、状态不变时重复复用、输入变化时重新失效重建”的缓存机制。

这一题适合做缓存优化，原因在于平面波模块中有一批数据具备共同特征：

- 它们由网格、晶格和截断参数决定，一旦对象状态不变，结果就是确定的。
- 它们在后续 SCF、FFT 和 `(k, G)` 相关流程中会被反复访问。
- 它们的构造往往涉及对大量平面波或 `(k, G)` 组合的遍历，重复计算成本明显。

因此，这次优化的重点不是引入新的物理模型，而是减少平面波相关重复构造，补齐缓存失效条件，并把缓存逻辑从原型推进到更完整、可验证、可维护的工程实现。

## 2. 优化前的主要问题

从 `PW_Basis` 和 `PW_Basis_K` 的职责看，优化前最适合缓存的几类数据包括：

- `PW_Basis` 中的 `gg`、`gdirect`、`gcar`
- `PW_Basis` 中的 `ig2igg`、`gg_uniq`
- `PW_Basis_K` 中的 `gcar`、`gk2`

这些数据在对象状态不变时并不会发生变化，但原始实现更接近“每次调用重新构造”的模式，因此会产生三类问题。

第一类问题是基础平面波信息的重复构造。`PW_Basis::collect_local_pw()` 每次都可能重新生成局域平面波对应的 `|G|^2`、整数格点坐标和笛卡尔坐标。对于后续会多次访问的模块来说，这意味着同一批 `G` 相关数据会被重复遍历和重复写入。

第二类问题是基于 `G` 向量继续加工出来的数据也会重复构造。`PW_Basis::collect_uniqgg()` 需要在 `gg` 的基础上建立去重后的 `gg_uniq` 和索引映射 `ig2igg`。如果底层 `gg` 没变，这部分排序、去重和映射其实也没有必要重复做。

第三类问题出现在 `PW_Basis_K`。这里不再只处理 `npw` 量级的 `G` 向量，而是要处理 `nks * npwk` 量级的 `(k, G)` 组合。`gcar` 与 `gk2` 的构造成本更高，而 `setupIndGk()` 还存在“先判断一遍是否满足截断条件，再在建索引时再判断一遍”的结构性重复遍历。

从题目 8 的角度看，以上问题说明这里确实存在一类典型的优化机会：不是算法公式变了，而是同样的构造工作在对象生命周期里做了太多次。

## 3. 缓存复用实现的演进

### 3.1 `906db2944 add cache-reuse`：缓存原型起点

提交 `906db2944 add cache-reuse` 是这次题目 8 工作的起点。根据提交记录，这一版已经把缓存复用的基本方向建立起来：

- `PW_Basis::collect_local_pw()` 开始在缓存有效时直接返回。
- `PW_Basis::collect_uniqgg()` 开始复用已经构造好的唯一 `|G|^2` 结构。
- `PW_Basis_K::collect_local_pw()` 开始尝试缓存 `gcar` 与 `gk2`。
- 相关单测开始补充“再次调用时不重新构造”的验证。

这一阶段的意义在于，它第一次把原本每次都重算的平面波中间结果，明确识别成“可缓存对象状态”。也就是说，题目 8 不再停留在概念层面，而是已经具备了可运行的缓存原型。

### 3.2 `82749621d Refine cache reuse...`：从原型走向可交付实现

相比原型版本，`82749621d Refine cache reuse, and refine test script...` 这次提交把缓存逻辑补到了更完整的工程形态。对应 `source/source_basis/module_pw` 下的实际代码改动，可以归纳为以下几个方面。

#### 3.2.1 `PW_Basis`：增加统计、RAII 和线程安全保护

在 `pw_basis.h` 中新增了 `CacheStats`，统计字段包括：

- `local_pw_hits`
- `local_pw_misses`
- `uniqgg_hits`
- `uniqgg_misses`
- `cache_bytes`

同时新增：

- `get_cache_stats()`
- `reset_cache_stats()`

这意味着缓存不再只是“做了”，而是可以被定量观察。尤其是 `cache_bytes` 提供了一个直接入口，用于估算为了复用而保留下来的主机端缓存空间。

在资源管理上，这一版把缓存数组交由 `std::unique_ptr` 管理，包括：

- `gg_cache_storage`
- `gdirect_cache_storage`
- `gcar_cache_storage`
- `ig2igg_cache_storage`
- `gg_uniq_cache_storage`

析构函数不再分别手工 `delete[]` 这些缓存数组，而是统一走 `clear_owned_cache()`。这使缓存所有权更清晰，也更符合题目要求里“遵循 RAII 原则”的方向。

在并发安全上，这一版引入：

- `std::atomic<bool> local_pw_cache_valid`
- `std::atomic<bool> uniqgg_cache_valid`
- `std::mutex cache_mutex`

`collect_local_pw()` 与 `collect_uniqgg()` 都采用了“先无锁检查、后加锁二次检查”的模式。这样做的结果是：

- 缓存已有效时，快速命中并直接返回。
- 首次构建时，通过互斥锁避免多个线程同时进入构建路径。
- 命中和未命中次数通过原子计数器记录。

这部分改动说明，题目 8 的缓存实现已经不只是“功能上能复用”，而是开始正面处理首次构建竞态和缓存可观测性问题。

#### 3.2.2 `PW_Basis`：补齐缓存失效链路

缓存可用的前提，是对象状态没有变化。因此这一轮很重要的一项工作，是把缓存失效条件显式补全。

`PW_Basis::invalidate_cache()` 统一负责把：

- `local_pw_cache_valid`
- `uniqgg_cache_valid`

标记为失效。

随后，在会改变平面波分布或网格状态的路径中，都补上了 `invalidate_cache()`，包括：

- `initmpi()`
- `initgrids()` 的两个重载路径
- `initparameters()`
- `setfullpw()`
- `get_ig2isz_is2fftixy()`

这部分的工程意义很直接：缓存机制只有在“状态变了就一定失效”的前提下才是正确的。补齐这些失效入口，解决的是缓存“对不对”的问题，而不只是“快不快”的问题。

#### 3.2.3 `collect_uniqgg()`：优先复用已有 `gg`

`collect_uniqgg()` 的细化也是这轮优化的重要点。新实现中，如果 `local_pw` 缓存已经有效且 `gg` 已存在，函数会优先直接读取已有的 `gg` 来构造排序和去重所需的临时数组，而不是再通过坐标回推一次 `|G|^2`。

这一点虽然不像“整段跳过重算”那样显眼，但它本质上仍然是在减少重复工作：在已有基础数据存在时，不再重复推导同样的能量模长。

#### 3.2.4 `PW_Basis_K`：把 `gcar` 与 `gk2` 拆成独立缓存

`PW_Basis_K` 的提升，是这轮实现中最能体现“缓存粒度细化”的部分。

在 `pw_basis_k.h` 中新增了 `KCacheStats`，在继承 `PW_Basis::CacheStats` 的基础上补充：

- `gcar_hits`
- `gcar_misses`
- `gk2_hits`
- `gk2_misses`

同时新增：

- `get_k_cache_stats()`
- `reset_k_cache_stats()`

更重要的是，`PW_Basis_K` 不再把 `gcar` 和 `gk2` 视作必须整体重建的一组数据，而是拆成两个独立有效位：

- `gcar_cache_valid`
- `gk_cache_valid`

对应地，缓存存储也拆成：

- `k_gcar_cache_storage`
- `k_gk2_cache_storage`

这种拆分的直接效果是：`gcar` 只依赖 `G` 向量和网格/晶格信息，而 `gk2` 还额外依赖 `erf_ecut`、`erf_height`、`erf_sigma`。因此当 `erf_*` 参数变化时，没有必要连 `gcar` 一起重算。

新实现的 `collect_local_pw()` 正是按这个逻辑工作的：

- 如果 `gcar` 有效，则直接命中。
- 如果 `gk2` 有效且 `erf_*` 参数未变，则直接命中。
- 如果只改了 `erf_*`，则保留 `gcar`，只重建 `gk2`。

这说明这轮优化已经不是简单的“有缓存/没缓存”二元判断，而是在 `PW_Basis_K` 中开始区分不同数据的依赖边界，避免不必要的连带重建。

#### 3.2.5 `setupIndGk()`：减少重复筛选

`setupIndGk()` 的结构也做了实际优化。新实现先为每个 `k` 点收集满足截断条件的 `selected_ig`，随后直接用这份结果去生成：

- `npwk`
- `igl2isz_k`
- `igl2ig_k`

这样做的意义在于，把原本“先遍历一次做计数，再遍历一次重新判断条件建索引”的流程，收敛成“第一次筛选时就把结果保存下来，后面直接消费这份筛选结果”。它没有改变数学结果，但减少了重复判断。

#### 3.2.6 设备侧同步拆分为更细粒度路径

在 `PW_Basis_K` 中，设备同步被拆成：

- `sync_gcar_device_cache()`
- `sync_gk2_device_cache()`

这使得当只有 `gk2` 重建、`gcar` 没变时，不必重复同步 `gcar` 的设备侧缓存；反之亦然。这部分改动属于缓存优化在 CPU/GPU 数据流层面的延伸，目标仍然是减少不必要的重复工作。

#### 3.2.7 单元测试补充了缓存复用行为验证

缓存机制如果没有测试，很难确认命中与失效是否符合预期。这一轮在两个测试文件中都补上了更具体的断言：

- `source/source_basis/module_pw/test/test1-1-1.cpp`
- `source/source_basis/module_pw/test_serial/pw_basis_k_test.cpp`

新增验证包括：

- 首次构建后 miss 计数是否为 1
- 再次调用时 hit 计数是否增加
- 二次调用时缓存指针是否保持不变
- 样本值是否保持一致
- `cache_bytes` 是否大于 0
- `PW_Basis_K` 中修改 `erf_*` 后是否只触发 `gk2` 的局部重建，而 `gcar` 继续复用

这些测试让题目 8 的缓存优化从“看代码推测成立”，进一步变成“有行为断言支持的实现”。

## 4. benchmark 方法与数据来源

本报告中的 benchmark 结论只基于仓库中已经存在的真实结果，不引用任何未提交或无法在当前仓库中复现的数据。

数据来源包括：

- `homework_docs/test_cases/README_task8_cache_benchmark.md`
- `homework_docs/test_cases/task8_cache_suite_runs/cache_20260530_221730`
- `homework_docs/test_cases/task8_cache_suite_runs/baseline_20260530_225026`
- `homework_docs/test_cases/task8_cache_suite_runs/combined_baseline_vs_cache/suite_summary.csv`
- `homework_docs/test_cases/task8_cache_suite_runs/combined_baseline_vs_cache/suite_summary.md`

其中，bench 验证链路的建立对应以下提交：

- `434b12cc4 bench: add bench scripts`
- `87e1304e4 bench: complete cache-reuse bench`
- `7c821b026 merge baseline bench`

需要强调的是，这几次提交的作用是补齐验证工作流，而不是题目 8 算法优化本身。也就是说，它们解决的是“如何稳定收集和比对结果”的问题，而不是“缓存逻辑怎么设计”的问题。

根据 `README_task8_cache_benchmark.md`，这套 benchmark 的固定规则是：

- 测试 case 为 `gaas_small` 和 `gaas_medium`
- MPI 进程数取 `1 / 2 / 4`
- OpenMP 线程数取 `1 / 2 / 4`
- 每组配置包含 `1` 次 warmup 与 `3` 次正式 repeat
- 汇总结果只统计 repeat，不把 warmup 计入最终性能数据
- baseline 与 cache 的比较以 median time 为准

加速比定义为：

```text
speedup_vs_baseline = baseline_median_time_s / cache_median_time_s
```

因此，本报告讨论的不是单次最优时间，而是按同一测试矩阵、同一统计口径得到的中位数结果。

## 5. benchmark 结果与如实结论

### 5.1 汇总结果

当前 `combined_baseline_vs_cache` 汇总表中，共有 18 组 baseline/cache 对比，对应：

- 2 个 case：`gaas_small`、`gaas_medium`
- 3 个 MPI 配置：`np = 1 / 2 / 4`
- 3 个 OpenMP 配置：`omp = 1 / 2 / 4`

按照汇总表中的 `speedup_vs_baseline` 统计：

- 1 组比 baseline 更快
- 13 组与 baseline 持平
- 4 组比 baseline 更慢

`speedup_vs_baseline` 的范围是：

- 最小值 `0.833`
- 最大值 `1.111`

全部 18 组的中位数为：

- `1.0`

也就是说，基于这批已经入库的 suite 数据，最保守且准确的结论是：缓存分支在当前测试矩阵下，大多数组合与 baseline 持平，少数组合略快或略慢，尚不能据此得出“总体性能稳定显著提升”的结论。

### 5.2 与 baseline 存在差异的 5 组配置

汇总表中只有 5 组配置的 `speedup_vs_baseline` 不等于 `1.0`，具体如下：

| Case | MPI | OMP | baseline median(s) | cache median(s) | speedup_vs_baseline |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gaas_medium` | 2 | 1 | 17.000 | 18.000 | 0.944 |
| `gaas_medium` | 4 | 1 | 9.000 | 10.000 | 0.900 |
| `gaas_medium` | 4 | 4 | 10.000 | 9.000 | 1.111 |
| `gaas_small` | 1 | 1 | 20.000 | 21.000 | 0.952 |
| `gaas_small` | 4 | 2 | 5.000 | 6.000 | 0.833 |

从这 5 组可以看到，差异本身是存在的，但方向并不统一：既有略快，也有略慢。尤其是在所有组合里，只有 `gaas_medium, np4, omp4` 这一组的中位数结果显示 cache 分支更快，其余 4 组都略慢于 baseline。

因此，这些数据更适合支持如下判断：

- 题目 8 的缓存机制已经进入真实工作负载测试；
- 当前 wall-clock 结果仍然比较接近，说明总时间受多种模块共同影响；
- 缓存优化的理论收益，尚未在这批总耗时指标上表现出稳定的一致优势。

### 5.3 细粒度 timer 的限制

根据 `suite_summary.md` 与 `README_task8_cache_benchmark.md` 的说明，当前汇总脚本会尝试解析：

- `pw_init`
- `plane_wave_generation`
- `gcar_or_gk2`
- `real2recip`
- `recip2real`

但在这次保存下来的题目 8 suite 结果中：

- `pw_init(s)`
- `pw_gen(s)`
- `gcar/gk2(s)`

三列全部为 `NA`。

这意味着当前日志并没有暴露出足够细的通用 timer，无法直接从这套结果里单独量化“平面波初始化阶段节省了多少时间”或者“`gcar/gk2` 复用具体减少了多少耗时”。

因此，现阶段可以直接引用的细粒度时间，主要仍然是 `real2recip` 和 `recip2real` 这类日志中已有的通用计时项，但它们并不能单独证明题目 8 缓存逻辑的收益大小。

## 6. 这次优化的实际价值

如果只看当前 suite 的总 wall-clock 结果，这次题目 8 工作并没有形成“性能已经稳定显著改善”的强证据。但如果结合代码改动本身来看，这轮工作依然有明确的工程价值。

第一，它把缓存复用从初步原型推进成了更完整的机制实现。`PW_Basis` 和 `PW_Basis_K` 中可复用的数据已经被显式建模，首次构建、命中返回、状态变化失效、析构释放等路径都更清晰了。

第二，它补齐了正确性和可维护性边界。缓存失效逻辑被集中化，缓存存储采用 RAII 管理，懒加载路径增加了线程安全保护，这些都让缓存机制更接近可长期维护的工程实现，而不是仅在局部场景下“碰巧能用”的优化补丁。

第三，它补齐了验证链路。命中/未命中统计接口、缓存字节数估算、单元测试行为断言，以及固定矩阵的 benchmark 脚本，使后续继续分析题目 8 的真实收益时，不必重新搭建基础设施。

换句话说，这次优化的成果更偏向“减少重复计算的机制建设”和“工程可靠性完善”，而不是已经完成了一次可以直接下结论的性能冲刺。

## 7. 结论与后续工作

综合当前分支中的真实代码、提交记录和 benchmark 结果，可以给出如下结论:

第一，题目 8 的缓存复用机制已经从 `906db2944 add cache-reuse` 的原型阶段，演进到 `82749621d Refine cache reuse...` 所代表的更完整工程实现。实际改动覆盖了缓存结构、缓存失效、RAII、线程安全、设备同步细化、统计接口和单元测试。

第二，bench 工作流已经通过 `434b12cc4`、`87e1304e4` 和 `7c821b026` 补齐，当前分支中已经能直接看到 baseline 与 cache 的成套对比结果。这为后续继续分析和复现实验提供了稳定基础。

第三，当前 `task8_cache_suite_runs` 中保存的 18 组对比结果，并不支持“总体性能稳定显著提升”的结论。更准确的说法是：大部分组合持平，少数组合略快或略慢，细粒度平面波初始化相关 timer 当前仍不可得。

因此，对这次题目 8 优化的最合适表述应当是：

- 缓存机制已经被实现、测试并工程化；
- 它的设计方向和代码质量都有明显进展；
- 但真实整机收益仍需更敏感的 timer、更多重负载样例。

