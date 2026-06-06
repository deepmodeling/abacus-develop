# module_pw 缓存复用性能对比

日期：2026-06-06  
对比分支：`develop` vs `feat/cache-reuse`  
对比范围：`source/source_basis/module_pw`

## 1. 本次补充的计时点

为了定位缓存复用收益，本次在 `module_pw` 内部补了最小必要的 timer：

- `PW_Basis::collect_local_pw`
- `PW_Basis::collect_local_pw_cache_hit`
- `PW_Basis::collect_local_pw_cache_build`
- `PW_Basis::collect_uniqgg`
- `PW_Basis::collect_uniqgg_cache_hit`
- `PW_Basis::collect_uniqgg_cache_build`
- `PW_Basis_K::collect_local_pw`
- `PW_Basis_K::collect_local_pw_cache_hit`
- `PW_Basis_K::collect_local_pw_build_gcar`
- `PW_Basis_K::collect_local_pw_build_gk2`

说明：

- `feat/cache-reuse` 上保留了命中/构建分支的细分 timer。
- `develop` 基线 worktree 只补了等价的“构建路径” timer，用来做公平基准，不包含缓存实现本身。

## 2. 基准方法

### 2.1 构建方式

为避免 MPI 环境对 `timer.cpp` 的影响，最终性能数据使用串行构建：

```bash
cmake -S /home/aunixt/abacus-develop -B /home/aunixt/abacus-develop/build-bench-serial \
  -DENABLE_MPI=OFF -DUSE_OPENMP=OFF -DUSE_ELPA=OFF -DBUILD_TESTING=OFF
cmake --build /home/aunixt/abacus-develop/build-bench-serial --target MODULE_PW_cache_bench
```

`develop` 基线在独立 worktree 中执行同样流程：

```bash
git -C /home/aunixt/abacus-develop worktree add /home/aunixt/abacus-develop-develop develop
cmake -S /home/aunixt/abacus-develop-develop -B /home/aunixt/abacus-develop-develop/build-bench-serial \
  -DENABLE_MPI=OFF -DUSE_OPENMP=OFF -DUSE_ELPA=OFF -DBUILD_TESTING=OFF
cmake --build /home/aunixt/abacus-develop-develop/build-bench-serial --target MODULE_PW_cache_bench
```

### 2.2 基准程序

新增基准程序：`MODULE_PW_cache_bench`

测试内容：

- `PW_Basis.setuptransform`
- `PW_Basis.collect_local_pw` 首次调用
- `PW_Basis.collect_local_pw` 重复 2000 次
- `PW_Basis.collect_uniqgg` 首次调用
- `PW_Basis.collect_uniqgg` 重复 2000 次
- `PW_Basis_K.setuptransform`
- `PW_Basis_K.collect_local_pw` 首次调用
- `PW_Basis_K.collect_local_pw` 重复 2000 次
- `PW_Basis_K.collect_local_pw(1.0, 0.5, 0.2)` 重复 2000 次

统计口径：

- 外层 wall time：基准程序直接测量
- 内层 timer：`ModuleBase::timer` 累积结果
- 每个分支各跑 3 次，结论使用中位数，避免单次抖动

## 3. 结果汇总

### 3.1 关键结论

1. `setuptransform` 基本无变化，说明优化没有破坏初始化主路径。
2. `PW_Basis.collect_local_pw` 的重复调用从持续重建，变成几乎纯命中路径，中位数加速约 `255.5x`。
3. `PW_Basis.collect_uniqgg` 的重复调用收益最大，中位数加速约 `2284.4x`。
4. `PW_Basis_K.collect_local_pw` 的重复调用中位数加速约 `342.7x`。
5. 即使传入新 `erf` 参数导致 `gk2` 需要重建，`PW_Basis_K` 仍然复用了 `gcar`，该路径中位数加速约 `463.0x`。

### 3.2 中位数对比表

单位：秒

| 场景 | develop | feat/cache-reuse | 提升 |
| --- | ---: | ---: | ---: |
| `PW_Basis.setuptransform.wall` | 0.001487681 | 0.001502684 | 0.99x |
| `PW_Basis.collect_local_pw.first.wall` | 0.000117131 | 0.000124510 | 0.94x |
| `PW_Basis.collect_local_pw.repeat.wall` | 0.088178164 | 0.000345102 | 255.52x |
| `PW_Basis.collect_uniqgg.first.wall` | 0.000440547 | 0.000401324 | 1.10x |
| `PW_Basis.collect_uniqgg.repeat.wall` | 0.754634649 | 0.000330313 | 2284.40x |
| `PW_Basis_K.setuptransform.wall` | 0.000259432 | 0.000217700 | 1.19x |
| `PW_Basis_K.collect_local_pw.first.wall` | 0.000134854 | 0.000121473 | 1.11x |
| `PW_Basis_K.collect_local_pw.repeat.wall` | 0.109138201 | 0.000318489 | 342.68x |
| `PW_Basis_K.collect_local_pw.gk2_rebuild.wall` | 0.193014060 | 0.000416850 | 463.04x |

## 4. 原始样本

### 4.1 feat/cache-reuse

#### Run 1

| 指标 | 数值 |
| --- | ---: |
| `PW_Basis.setuptransform.wall` | 0.001464156 |
| `PW_Basis.collect_local_pw.repeat.wall` | 0.000462526 |
| `PW_Basis.collect_uniqgg.repeat.wall` | 0.000362283 |
| `PW_Basis_K.collect_local_pw.repeat.wall` | 0.000318489 |
| `PW_Basis_K.collect_local_pw.gk2_rebuild.wall` | 0.000416850 |

#### Run 2

| 指标 | 数值 |
| --- | ---: |
| `PW_Basis.setuptransform.wall` | 0.001733625 |
| `PW_Basis.collect_local_pw.repeat.wall` | 0.000345102 |
| `PW_Basis.collect_uniqgg.repeat.wall` | 0.000330313 |
| `PW_Basis_K.collect_local_pw.repeat.wall` | 0.000373748 |
| `PW_Basis_K.collect_local_pw.gk2_rebuild.wall` | 0.000529150 |

#### Run 3

| 指标 | 数值 |
| --- | ---: |
| `PW_Basis.setuptransform.wall` | 0.001502684 |
| `PW_Basis.collect_local_pw.repeat.wall` | 0.000314556 |
| `PW_Basis.collect_uniqgg.repeat.wall` | 0.000287104 |
| `PW_Basis_K.collect_local_pw.repeat.wall` | 0.000317902 |
| `PW_Basis_K.collect_local_pw.gk2_rebuild.wall` | 0.000416157 |

### 4.2 develop

#### Run 1

| 指标 | 数值 |
| --- | ---: |
| `PW_Basis.setuptransform.wall` | 0.001473400 |
| `PW_Basis.collect_local_pw.repeat.wall` | 0.086802023 |
| `PW_Basis.collect_uniqgg.repeat.wall` | 0.754634649 |
| `PW_Basis_K.collect_local_pw.repeat.wall` | 0.109138201 |
| `PW_Basis_K.collect_local_pw.gk2_rebuild.wall` | 0.193014060 |

#### Run 2

| 指标 | 数值 |
| --- | ---: |
| `PW_Basis.setuptransform.wall` | 0.001487681 |
| `PW_Basis.collect_local_pw.repeat.wall` | 0.088178164 |
| `PW_Basis.collect_uniqgg.repeat.wall` | 0.745032054 |
| `PW_Basis_K.collect_local_pw.repeat.wall` | 0.104775062 |
| `PW_Basis_K.collect_local_pw.gk2_rebuild.wall` | 0.176847280 |

#### Run 3

| 指标 | 数值 |
| --- | ---: |
| `PW_Basis.setuptransform.wall` | 0.001664659 |
| `PW_Basis.collect_local_pw.repeat.wall` | 0.105299161 |
| `PW_Basis.collect_uniqgg.repeat.wall` | 0.965178448 |
| `PW_Basis_K.collect_local_pw.repeat.wall` | 0.134283688 |
| `PW_Basis_K.collect_local_pw.gk2_rebuild.wall` | 0.266846180 |

## 5. 内部 timer 观察

### 5.1 feat/cache-reuse

代表性现象：

- `timer.PW_Basis.collect_local_pw_cache_hit.calls = 2000`
- `timer.PW_Basis.collect_uniqgg_cache_hit.calls = 2000`
- `timer.PW_Basis_K.collect_local_pw_cache_hit.calls = 3999`
- `timer.PW_Basis_K.collect_local_pw_build_gcar.calls = 1`
- `timer.PW_Basis_K.collect_local_pw_build_gk2.calls = 2`

解释：

- `PW_Basis` 的两条缓存路径在首轮构建后，后续 2000 次全部命中。
- `PW_Basis_K` 在默认参数重复调用时，只有首轮需要构建。
- 改变 `erf` 参数后，`gk2` 需要重新构建，但 `gcar` 仍然保持复用。

### 5.2 develop

代表性现象：

- `timer.PW_Basis.collect_local_pw_cache_build.calls = 2001`
- `timer.PW_Basis.collect_uniqgg_cache_build.calls = 2001`
- `timer.PW_Basis_K.collect_local_pw_build_gcar.calls = 4001`
- `timer.PW_Basis_K.collect_local_pw_build_gk2.calls = 4001`

解释：

- 基线分支每次调用都走完整构建，没有任何命中路径。
- 这与重复调用 wall time 的数量级差异完全一致。

## 6. 验证记录

在 `feat/cache-reuse` 上新增并通过的回归测试：

- `PWBasisTEST.CacheCollectionRecordsTimers`
- `PWBasisKTEST.CollectLocalPWRecordsTimers`

测试命令：

```bash
cd /home/aunixt/abacus-develop/build-tests-mpi/source/source_basis/module_pw/test_serial
./MODULE_PW_basis_pw_serial --gtest_filter=PWBasisTEST.CacheCollectionRecordsTimers
./MODULE_PW_basis_pw_k_serial --gtest_filter=PWBasisKTEST.CollectLocalPWRecordsTimers
```

结果：两项均通过。

## 7. 最终结论

`feat/cache-reuse` 在 `module_pw` 中的缓存复用优化是有效且收益非常显著的，主要收益集中在重复调用路径：

- `PW_Basis.collect_local_pw`
- `PW_Basis.collect_uniqgg`
- `PW_Basis_K.collect_local_pw`

其中最关键的是：

- `develop` 重复调用仍然持续分配并重建数据。
- `feat/cache-reuse` 已经将这部分开销压缩到首轮构建，后续主要变成 cache hit。
- `PW_Basis_K` 在参数部分变化时还能保留 `gcar` 复用，说明缓存粒度设计是合理的。

如果后续需要，我建议直接把这份基准保留在仓库里，后面可以继续扩展成 CI 可重复的 micro-benchmark。
