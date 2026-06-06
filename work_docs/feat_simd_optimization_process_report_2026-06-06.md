# `feat/simd` 优化过程报告

## 1. 总览

本报告结合 `feat/simd` 相对 `develop` 的提交历史与 **2026-06-06** 采集的性能数据，还原并总结这条优化分支的演进过程。

该分支的目标并非修改 `module_pw` 的数学结果，而是降低 transform 路径中复数缓冲区重复搬运的开销，重点覆盖以下阶段：

- gather planes / scatter sticks
- gather sticks / scatter planes
- `real2recip` 与 `recip2real` 周围的顶层缓冲区拷贝

最终性能呈”部分显著成功、部分仍需验证”的状态：

- `PW_Basis` 路径：约 **1.5x** 的明确提升
- `PW_Basis_K.medium`：功能正确，但未稳定跑出性能收益

这条分支的价值不仅在于最终加速比，更在于其演进过程本身——它展示了一条从”先试一个低层 copy 优化”起步，逐步推进到”可维护、可测试、可量化”的完整优化路径。

## 2. 按 commit 展开的优化过程

### `7c58a45ce` `have a try`

这条分支的首次优化尝试，直接瞄准 `pw_gatherscatter.h` 中高频重复的复数拷贝循环。

#### 背景与问题

原始实现按 `std::complex<T>` 逐元素复制：每轮循环只处理一个复数，且循环在多个 gather/scatter 路径里反复出现。对编译器而言，这种写法未能明确表达”底层其实是连续交错的实部/虚部标量流”，在 transform 多次调用下，拷贝开销容易成为热点。

#### 改动

在 `gatherp_scatters` 与 `gathers_scatterp` 的每个内部循环中（[pw_gatherscatter.h](../source/source_basis/module_pw/pw_gatherscatter.h)），将原始逐元素复数拷贝改写为：

- 复数缓冲区通过 `reinterpret_cast` 转为连续交错的标量数组
- 以 `2 * count` 个标量取代 `count` 个复数对象进行复制
- 引入 `__restrict__` 指针消除别名冲突
- 对 GCC 构建增加 `#pragma GCC ivdep` 提示向量化

```cpp
// 原始实现：每轮一个复数
for (int iz = 0; iz < nz_; ++iz) { outp[iz] = inp[iz]; }

// 优化后：连续标量流 + __restrict__ + ivdep
T* __restrict__ outp_r = reinterpret_cast<T*>(outp);
const T* __restrict__ inp_r = reinterpret_cast<const T*>(inp);
#pragma GCC ivdep
for (int iz = 0; iz < 2 * nz_; ++iz) { outp_r[iz] = inp_r[iz]; }
```

核心思路：将复数拷贝重写为更接近”连续标量流”的形式，降低编译器向量化门槛。

#### 局限

方向正确，但版本不成熟：

- 优化逻辑在多个位置重复展开
- 依赖 GCC 特定 pragma，可移植性一般
- 仅覆盖 gather/scatter 内部循环，未触及 `pw_transform.cpp` / `pw_transform_k.cpp` 的顶层 copy
- 缺乏针对性的正确性测试

该 commit 定位为”验证思路是否值得推进”的试探性步骤，而非最终形态。

### `c268969f9` `refine complex buffer copies in module_pw`

本次分支最核心的一次重构：将前一个提交的”局部尝试”抽象为可复用机制。

#### 背景与问题

第一轮尝试后，面临两个主要问题：

1. 低层 copy 优化代码在多个位置重复，维护成本高
2. 优化范围不完整，未覆盖整个 transform 拷贝链路

#### 改动

在 [pw_gatherscatter.h](../source/source_basis/module_pw/pw_gatherscatter.h) 中引入两个 helper（第 12–44 行）：

- `detail::copy_complex_buffer` — 基于 `std::copy_n` 的串行半精度复数拷贝
- `detail::copy_complex_buffer_parallel` — 按 1024 元素分块、配合 OpenMP `#pragma omp parallel for` 的并行版本

```cpp
// 串行版本（pw_gatherscatter.h:12-20）
template <typename T>
inline void copy_complex_buffer(const std::complex<T>* in,
                                std::complex<T>* out, const int count) {
    if (count <= 0) return;
    std::copy_n(in, count, out);
}

// 并行版本（pw_gatherscatter.h:25-44）：大于 chunk_size 时启动 OpenMP
template <typename T>
inline void copy_complex_buffer_parallel(const std::complex<T>* in,
                                         std::complex<T>* out, const int count) {
    constexpr int chunk_size = 1024;
    if (count <= chunk_size) { copy_complex_buffer(in, out, count); return; }
    #pragma omp parallel for schedule(static)
    for (int offset = 0; offset < count; offset += chunk_size) {
        int chunk_count = std::min(chunk_size, count - offset);
        std::copy_n(in + offset, chunk_count, out + offset);
    }
}
```

并将 helper 应用于以下路径：

| 方法 | 位置 | 调用 |
|---|---|---|
| `PW_Basis::gatherp_scatters` | [pw_gatherscatter.h:72,94,133](../source/source_basis/module_pw/pw_gatherscatter.h#L72) | `copy_complex_buffer` |
| `PW_Basis::gathers_scatterp` | [pw_gatherscatter.h:176,205,249](../source/source_basis/module_pw/pw_gatherscatter.h#L176) | `copy_complex_buffer` |
| `PW_Basis::real2recip` | [pw_transform.cpp:37](../source/source_basis/module_pw/pw_transform.cpp#L37) | `copy_complex_buffer_parallel` |
| `PW_Basis::recip2real` | [pw_transform.cpp:196](../source/source_basis/module_pw/pw_transform.cpp#L196) | `copy_complex_buffer_parallel` |
| `PW_Basis_K::real2recip` | [pw_transform_k.cpp:36](../source/source_basis/module_pw/pw_transform_k.cpp#L36) | `copy_complex_buffer_parallel` |
| `PW_Basis_K::recip2real` | [pw_transform_k.cpp:197](../source/source_basis/module_pw/pw_transform_k.cpp#L197) | `copy_complex_buffer_parallel` |

至此，优化目标不再局限于局部 gather/scatter 循环，而是扩展到 FFT transform 前后整条复数缓冲区搬运链的 **6 个调用点**。

#### 为何是关键转折

最终 benchmark 印证了这一设计的重要性：`PW_Basis` 端到端加速比明显高于单个 gather/scatter timer 的提升。这说明最终收益并非仅来自局部循环，而是来自 **gather/scatter 内部 copy** 与 **transform 顶层 staging copy** 两端的共同改善。

“沿整条 copy 链统一优化”的思路，正是从这个 commit 开始建立的。

#### 仍缺什么

性能路径已经收拢，但正确性保障不足。copy-heavy 优化极易出现”看起来能跑，但数值细节有偏差”的问题，需要补充有针对性的 round-trip 测试。

### `754fe85bb` `add module_pw complex transform round-trip tests`

补上优化过程中最缺的一块：针对 transform 正确性的专项测试。

#### 背景与问题

复数拷贝优化存在几个典型风险：

- 实部和虚部次序被破坏
- 拷贝长度正确但布局错误
- round-trip 表面上能工作，但数值发生微小漂移

若无专门测试，这些问题往往在后期才暴露。

#### 改动

为两条路径分别增加串行单元测试：

- `PWBasisTEST.ComplexTransformRoundTrip` — [test_serial/pw_basis_test.cpp:367-403](../source/source_basis/module_pw/test_serial/pw_basis_test.cpp#L367)
- `PWBasisKTEST.ComplexTransformRoundTrip` — [test_serial/pw_basis_k_test.cpp:195-236](../source/source_basis/module_pw/test_serial/pw_basis_k_test.cpp#L195)

测试流程：构造确定性的 reciprocal-space 输入 → 调用 `recip2real` → 再调用 `real2recip` → 对恢复的 reciprocal-space 数据逐元素比较。

```cpp
// pw_basis_test.cpp:367-403 — 简化示意
std::vector<std::complex<double>> recip_in(pwb.npw);
// 用确定性格子 (ig % 11 - 5) / 7.0 + i * (ig % 13 - 6) / 9.0 填充
pwb.recip2real(recip_in.data(), real_space.data());
pwb.real2recip(real_space.data(), recip_out.data());
for (int ig = 0; ig < pwb.npw; ++ig) {
    EXPECT_NEAR(recip_in[ig].real(), recip_out[ig].real(), 1e-10);
    EXPECT_NEAR(recip_in[ig].imag(), recip_out[ig].imag(), 1e-10);
}
```

选择 reciprocal-space 作为输入源是刻意的：若直接取任意 real-space 数据做 round-trip，plane-wave cutoff 投影本身即可能导致数据无法精确恢复，从而制造与优化无关的假失败（详见 [pw_basis_test.cpp:383-384](../source/source_basis/module_pw/test_serial/pw_basis_test.cpp#L383) 注释）。

#### 为何重要

这是整条优化链中”建立安全护栏”的一步。只有将 round-trip 正确性明确固定下来，后续的性能改写才不会演变成”为了更快，悄悄牺牲数值行为”。

### `25ebe2e30` `document module_pw copy helpers and tests`

改动量不大，但让分支从”能跑”走向”容易理解和维护”。

#### 背景与问题

helper 与测试补上之后，若设计意图缺乏说明，后续维护者仍可能误判：

- 为何 helper 要按交错标量流的思路来写
- 为何测试从 reciprocal-space 输入开始，而非任意 real-space 数组
- 为何同时保留串行 helper 和顶层 parallel helper

#### 改动

补充几类说明性注释（详见 [pw_gatherscatter.h:22-23](../source/source_basis/module_pw/pw_gatherscatter.h#L22)、[pw_basis_test.cpp:383-384](../source/source_basis/module_pw/test_serial/pw_basis_test.cpp#L383)）：

- helper 设计意图：让编译器更容易处理连续的实部/虚部数据
- round-trip 测试起点：reciprocal-space 输入才是合理选择
- 并行策略：顶层 transform copy 拥有自己的 OpenMP 区域；gather/scatter 内部循环则在已有并行区中调用非并行 helper

#### 价值

它本身不直接提升性能，但降低了后续维护中”误删优化””误改测试””看不懂所以回退”的风险。对于底层性能优化而言，可解释性不是附属品，而是长期稳定性的组成部分。

### `f3a0b6b4c` `Merge branch 'deepmodeling:develop' into feat/simd`

该 merge commit 并非 SIMD 优化的实现部分，但在演进过程中承担了”集成与收敛”的角色。

#### 意义

- 将分支与当时的 `develop` 同步
- 确保 SIMD 改动能与主线近期演化共存
- 降低后续合并时大规模冲突的风险

因此，它应被理解为”将性能优化保持在主线可集成状态”的必要步骤，而非单独的性能改进。

### `3245d2d31` `remove pragma GCC ivdep and use std::copy_n`

这是分支最后一个关键的”收口”与”工程化”提交。

#### 背景与问题

此前 helper 仍较依赖 GCC 风格：

- `reinterpret_cast` 到标量流
- `#pragma GCC ivdep`

这种写法虽可能有效，但存在明显问题：编译器耦合强、可读性弱、维护者理解成本高。

#### 改动

**标准库化** — 将 helper 调整为更可移植的实现（[pw_gatherscatter.h:12-44](../source/source_basis/module_pw/pw_gatherscatter.h#L12)）：

- `copy_complex_buffer` → 内层使用 `std::copy_n`，移除 `reinterpret_cast`、`__restrict__` 与 `#pragma GCC ivdep`
- `copy_complex_buffer_parallel` → 按 1024 元素分块，在大 buffer 上以 OpenMP `#pragma omp parallel for` 调度 `std::copy_n`
- 显式补充 `#include <algorithm>`

```cpp
// 最终版：std::copy_n 替代手工循环（pw_gatherscatter.h:19）
std::copy_n(in, count, out);
```

**可观测性** — 在 copy-sensitive 路径补充 `ModuleBase::timer` 时间戳（[pw_gatherscatter.h](../source/source_basis/module_pw/pw_gatherscatter.h) 的 6 处 start/end：[L63-L74](../source/source_basis/module_pw/pw_gatherscatter.h#L63)、[L85-L96](../source/source_basis/module_pw/pw_gatherscatter.h#L85)、[L120-L136](../source/source_basis/module_pw/pw_gatherscatter.h#L120)、[L157-L178](../source/source_basis/module_pw/pw_gatherscatter.h#L157)、[L192-L208](../source/source_basis/module_pw/pw_gatherscatter.h#L192)、[L239-L251](../source/source_basis/module_pw/pw_gatherscatter.h#L239)）；并在 `real2recip`/`recip2real` 顶层（[pw_transform.cpp:30-66](../source/source_basis/module_pw/pw_transform.cpp#L30)、[pw_transform_k.cpp:32-68](../source/source_basis/module_pw/pw_transform_k.cpp#L32)）也统一添加 timer。

**验证加强**：

- 在 [pw_basis_k_test.cpp:213-214](../source/source_basis/module_pw/test_serial/pw_basis_k_test.cpp#L213) 增加 `npwk` 合法性检查：`ASSERT_NE(basis_k.npwk, nullptr)` 与 `ASSERT_GT(basis_k.npwk[0], 0)`
- 在 [pw_basis_k_test.cpp:238-274](../source/source_basis/module_pw/test_serial/pw_basis_k_test.cpp#L238) 增加可选开启的 copy benchmark（`CopyComplexBufferTimerBenchmark`，通过环境变量 `ABACUS_PW_SIMD_TIMER_TEST` 控制）

#### 为何重要

这是分支从”手工调优尝试”走向”标准库驱动、便于测量、便于维护”的关键一步：

- 不再依赖 GCC 特有 pragma 表达优化意图
- 代码语义更清晰，更接近标准 C++
- 性能热点变得可测量，而非仅凭感觉判断

本次 `feat/simd` 与 `develop` 的性能对比之所以能顺利完成，正是因为这个提交将 timer 与 helper 结构整理到了适合 benchmark 的状态。

## 3. 性能结果到提交的映射

最终 benchmark 可反推各阶段的贡献分布。

### `PW_Basis` 路径

端到端约 **1.53x–1.55x** 的提升，表明以下提交的组合有效：

- `7c58a45ce` — 首次优化尝试
- `c268969f9` — 核心抽象与链路扩展
- `3245d2d31` — 工程化收口

尤其值得关注的是：gather/scatter copy timer 有提升，而 `real2recip` / `recip2real` 顶层 timer（[pw_transform.cpp:30-66](../source/source_basis/module_pw/pw_transform.cpp#L30)、[pw_transform_k.cpp:32-68](../source/source_basis/module_pw/pw_transform_k.cpp#L32)）提升更明显。这恰好印证了 `c268969f9` 的核心思路——收益来自整条 transform copy 链的统一优化，而非局部循环的孤立改进。

### `PW_Basis_K` 路径

`PW_Basis_K.medium` 的结果更为复杂：

- `754fe85bb` 保证了正确性
- `25ebe2e30` 和 `3245d2d31` 提升了可解释性与可测量性
- 但当前串行中等规模基准未给出净加速结果

这并不意味着分支没有价值，而是揭示了明显的证据分层：对 `PW_Basis` 收益明确，对本次测试的 `PW_Basis_K` case 尚不充分。从优化过程本身来看，这同样是有价值的——分支不仅做了优化，还通过测试和 timer 明确暴露了”哪些地方收益明显、哪些地方还没跑出来”。

## 4. 总结

`feat/simd` 展示了一条完整的性能优化路径：

1. **试探** — 将热点复数拷贝循环改写为利于向量化的低层形式
2. **抽象** — 将局部技巧封装为 helper，扩展至整条 transform 拷贝链路
3. **验证** — 增加 round-trip 测试，确保优化不破坏数值正确性
4. **解释** — 补充注释与设计说明，确保可维护性
5. **工程化** — 以 `std::copy_n` 替换编译器特定 pragma，提升可移植性
6. **可观测** — 补充 timer，使性能结论可量化验证

与 `develop` 对比后，该分支在 `PW_Basis` 路径上的成功是明确的：带来了可重复、可解释的性能提升。在 `PW_Basis_K.medium` 测试口径下，它保持了正确性，但尚未取得稳定加速的充分证据。

这条分支最值得肯定的，不仅是它让某些热点更快了，更在于它将优化本身做成了一个完整的工程过程——有实验、有抽象、有测试、有注释、有 timer，也有对收益边界的诚实呈现。
