# `feat/simd` 与 `develop` 在 `module_pw` 上的性能对比报告

## 1. 对比范围

本报告对比了 **2026-06-06** 当天 `feat/simd` 与 `develop` 两个分支在 `module_pw` 模块上的性能表现。

- 优化分支：`feat/simd`
- 基线分支：`develop`
- 对比重点：SIMD 优化涉及到的复数缓冲区拷贝路径
- 基准入口：`source/source_basis/module_pw/test_serial/pw_simd_bench.cpp`
- 构建选项：`-DENABLE_MPI=OFF -DUSE_OPENMP=OFF -DUSE_ELPA=OFF -DBUILD_TESTING=OFF`

为了保证对比公平，`develop` 分支是在独立临时 worktree `/home/aunixt/abacus-develop-develop` 中测得，临时只补充了两类内容：

1. 与 `feat/simd` 完全一致的 `MODULE_PW_simd_bench` 基准程序
2. 与优化分支同名的 `ModuleBase::timer` 时间戳，用于记录相同的 gather/scatter 拷贝阶段

对 `develop` 没有回移植任何 SIMD 优化逻辑，基线分支的算法行为保持不变。

## 2. 基准设计

本次基准使用 reciprocal-space round-trip transform 作为统一测试口径。

- `PW_Basis.medium`
  - 晶格：`Matrix3(1,0,1; 0,2,0; 0,0,2)`
  - `gridecut=30.0`，`pwecut=20.0`
  - `nrxx=320`，`npw=49`
  - 重复次数：`4096`
- `PW_Basis.large`
  - 晶格：`Matrix3(2,0,0; 0,2,0; 0,0,2)`
  - `gridecut=40.0`，`pwecut=25.0`
  - `nrxx=729`，`npw=147`
  - 重复次数：`2048`
- `PW_Basis_K.medium`
  - 单个 k 点：`{0,0,0}`
  - 几何参数与 `PW_Basis.medium` 相同
  - `nrxx=320`，`npwk=49`
  - 重复次数：`4096`

两个分支都在同一台机器上各运行 **3 次**，最终报告采用 **3 次结果的中位数** 作为比较依据。

## 3. 中位数结果

### 3.1 端到端 round-trip 耗时

| 用例 | 指标 | `develop` 中位数 | `feat/simd` 中位数 | 加速比 |
|---|---:|---:|---:|---:|
| `PW_Basis.medium` | ms/op | 0.001819651 | 0.001192664 | **1.526x** |
| `PW_Basis.large` | ms/op | 0.004269503 | 0.002761381 | **1.546x** |
| `PW_Basis_K.medium` | ms/op | 0.001135719 | 0.001236483 | **0.919x** |

### 3.2 与拷贝路径相关的 timer 分解结果

| 用例 | Timer 名称 | `develop` 中位数 (s) | `feat/simd` 中位数 (s) | 加速比 |
|---|---|---:|---:|---:|
| `PW_Basis.medium` | `real2recip` | 0.002785 | 0.001803 | **1.545x** |
| `PW_Basis.medium` | `recip2real` | 0.004378 | 0.002761 | **1.586x** |
| `PW_Basis.medium` | `gatherp_copy_serial` | 0.000396 | 0.000311 | **1.273x** |
| `PW_Basis.medium` | `gathers_copy_serial` | 0.000528 | 0.000341 | **1.548x** |
| `PW_Basis.large` | `real2recip` | 0.003565 | 0.002367 | **1.506x** |
| `PW_Basis.large` | `recip2real` | 0.005033 | 0.003151 | **1.597x** |
| `PW_Basis.large` | `gatherp_copy_serial` | 0.000361 | 0.000264 | **1.367x** |
| `PW_Basis.large` | `gathers_copy_serial` | 0.000333 | 0.000293 | **1.137x** |
| `PW_Basis_K.medium` | `real2recip` | 0.001949 | 0.002108 | 0.925x |
| `PW_Basis_K.medium` | `recip2real` | 0.002479 | 0.002634 | 0.941x |
| `PW_Basis_K.medium` | `gatherp_copy_serial` | 0.000342 | 0.000354 | 0.966x |
| `PW_Basis_K.medium` | `gathers_copy_serial` | 0.000342 | 0.000309 | 1.107x |

## 4. 结果解读

### 4.1 哪些部分获得了明显提升

在本次测试口径下，`feat/simd` 在 `PW_Basis` 路径上表现出比较稳定的收益，端到端 round-trip 性能大约提升 **1.5 倍**。这种提升不仅体现在总耗时上，也体现在顶层 transform timer 上：

- `PW_Basis.medium`：`real2recip` 与 `recip2real` 都提升了约 **1.55x**
- `PW_Basis.large`：`real2recip` 与 `recip2real` 提升约 **1.51x 到 1.60x**

这与分支中的实现修改是一致的。`pw_gatherscatter.h`、`pw_transform.cpp` 和 `pw_transform_k.cpp` 的改动并没有改变 FFT 的数学流程，而是降低了 transform 前后重复复数缓冲区搬运的成本。

### 4.2 性能收益主要来自哪里

从 copy-phase timer 可以直接看到，优化分支确实改善了串行 gather/scatter 阶段的拷贝开销：

- `gatherp_copy_serial`：提升约 **1.27x 到 1.37x**
- `gathers_copy_serial`：提升约 **1.14x 到 1.55x**

而端到端加速比比这些局部 timer 还高，说明收益不只来自 gather/scatter 内部循环，还来自顶层 transform 中的连续缓冲区拷贝优化，具体包括：

- `PW_Basis::real2recip`
- `PW_Basis::recip2real`
- `PW_Basis_K::real2recip`
- `PW_Basis_K::recip2real`

也正因为覆盖了更完整的 copy 链路，所以整体 transform 路径的收益会大于单个 gather/scatter 子阶段。

### 4.3 哪些部分在本次基准中没有体现提升

`PW_Basis_K.medium` 在这组串行微基准下没有表现出净提升，反而中位数上出现了轻微回退：

- `develop`：`0.001135719 ms/op`
- `feat/simd`：`0.001236483 ms/op`
- 比值：`0.919x`

内部 timer 也说明这个 case 当前更接近噪声区间：

- `gatherp_copy_serial` 基本持平
- `gathers_copy_serial` 有小幅提升
- 顶层 `real2recip` / `recip2real` 略慢于基线

因此，从现有证据出发，更准确的结论应当是：

- `PW_Basis`：SIMD/copy 重构收益明确
- `PW_Basis_K.medium`：在这组基准下，收益暂时没有被证明出来

## 5. 原始三轮数据

### 5.1 `feat/simd`

| 轮次 | `PW_Basis.medium` ms/op | `PW_Basis.large` ms/op | `PW_Basis_K.medium` ms/op |
|---|---:|---:|---:|
| 1 | 0.001192664 | 0.002739758 | 0.001266641 |
| 2 | 0.001468912 | 0.002972134 | 0.001236483 |
| 3 | 0.001141960 | 0.002761381 | 0.001045821 |

### 5.2 `develop`

| 轮次 | `PW_Basis.medium` ms/op | `PW_Basis.large` ms/op | `PW_Basis_K.medium` ms/op |
|---|---:|---:|---:|
| 1 | 0.001819651 | 0.004574592 | 0.001130419 |
| 2 | 0.001776332 | 0.004053817 | 0.001147463 |
| 3 | 0.002436848 | 0.004269503 | 0.001135719 |

## 6. 最终结论

从 `module_pw` 这次对比来看，`feat/simd` 在 `PW_Basis` 路径上带来了**明确且可重复的性能提升**，在本次串行 round-trip 微基准中，中位数加速比约为 **1.53x 到 1.55x**。

但在当前这组 `PW_Basis_K.medium` 基准中，并没有观察到同样明确的收益。因此更客观的总结应该是：

- `PW_Basis`：优化成功，收益已经被清楚证明
- `PW_Basis_K.medium`：功能正确性保持不变，但当前基准下性能收益尚未建立

如果后续还需要继续补充性能论证，下一步更值得做的是增加更大的 `PW_Basis_K` 用例，或者在 MPI / OpenMP 配置下继续测试，以便观察更大问题规模下 helper 开销被摊薄之后，是否能更充分体现 SIMD/copy 优化的价值。
