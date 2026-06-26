# Phase 4 近邻搜索性能 Benchmark 报告

## 测试口径

- 基线：ABACUS `3.11.0-beta.1`，提交 `78d71c9134a47dd63f76e69cc1a7025f302cc09d`。
- 优化版：`final-project`，提交 `241c28a54fa194650ccf6a8f0044563ae31e464b`。
- 优化版工作树：包含待提交修改；neighbor diff SHA-256 为 `c2be3bda733fa8a43ddec20bff38b67718fb5267cf987995c80d2d79ea2f7b33`。
- 两个版本使用相同编译器、依赖、Release 配置、算例和 CPU 绑定。
- 当前版在解决上游 PR 冲突时包含少量 develop 接口适配、测试依赖和 mock 链接修复；这些修改未改变近邻搜索核心算法。
- 本报告为单节点测试，不据此声称多节点扩展效率。
- 时间采用正式重复运行的中位数；提交版原始数据位于 `docs/benchmark_results/phase4/`。
- beta.1 与当前 conda ScaLAPACK/MKL 组合存在 ABI 崩溃，因此整体版本对比统一使用单 rank LAPACK；MPI 由独立 neighbor probe 测量。

## 环境与方法

- CPU：Intel(R) Core(TM) i9-14900HX，16 核/32 线程；内存约 15Gi。
- 环境：WSL 单节点，GCC 14.3.0、OpenMPI 5.0.10、Conda `bxcx`。
- 构建：Release、MPI/LCAO/OpenMP 开启，ELPA 和 native optimization 关闭。
- kernel：4,096/32,768/131,072 原子，每个模式预热一次并正式运行 10 次。
- 三路径：beta.1 legacy、current Full27、current Half14 使用相同 PBC 晶胞、坐标和 cutoff；beta.1 因 O(N²) 仅运行到 32,768 原子，其中 4,096 原子为 10 次正式重复、32,768 原子为 3 次正式重复。
- MPI：固定 65,536 原子的强扩展，以及每 rank 4,096 原子的弱扩展；覆盖 1/2/4/8/16 ranks 和 orthogonal/triclinic，每组运行 5 次。
- 版本对比：Si16/Si64 LCAO SCF，预热一次、正式运行 3 次，使用中位数。
- Verlet：Si2 relax 5 步、Si8 MD 10 步和 Si64 MD 20 步，比较 `neighbor_skin=0.0/3.0`。

## 结论摘要

- 完整 SCF 用于验证端到端物理结果；其总耗时会被 Hamiltonian、求解器等主流程主导，不能单独代表近邻搜索热点收益。
- 三路径核心搜索显示 beta.1 legacy 在 32K 原子时约需 97-104 秒，而当前 Half14 为 83-181 毫秒；131K 下 beta.1 因 O(N²) 成本不再实际运行。
- 重构后的 Half14 在 4K/32K/131K 原子上分别比 Full27 快约 1.68/2.05/2.24 倍，pair hash 完全一致，峰值 RSS 也未增加。
- 取消默认 Full27 reference 继续提供明确收益：131,072 原子时避免约 81 MiB 峰值 RSS 和额外 reference 构建。
- Verlet 显著降低重建次数；Si8 小体系时间持平，Si64 20 步 MD 从 301 秒降至 297 秒，收益约 1.3%。
- 单节点 MPI strong/weak 在 1/2/4 ranks 表现较好，8/16 ranks 受到 ghost 增长、通信和 WSL 单节点资源竞争影响。

## 三路径核心搜索对比

| Cell | Radius | Atoms | beta.1 legacy total (ms) | Current Full27 total (ms) | Current Half14 total (ms) | beta.1/Half14 | Full27/Half14 | Pair hash |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| orthogonal | 2.2 | 4,096 | 843.703 | 10.851 | 10.336 | 81.624 | 1.050 | 一致 |
| orthogonal | 2.2 | 32,768 | 96748.051 | 98.977 | 83.758 | 1155.090 | 1.182 | 一致 |
| orthogonal | 2.2 | 131,072 | N/A | 445.159 | 334.744 | N/A | 1.330 | 一致（beta.1 大规模跳过） |
| orthogonal | 3.4 | 4,096 | 855.606 | 30.941 | 18.988 | 45.062 | 1.630 | 一致 |
| orthogonal | 3.4 | 32,768 | 102702.009 | 301.244 | 181.366 | 566.269 | 1.661 | 一致 |
| orthogonal | 3.4 | 131,072 | N/A | 1203.959 | 737.776 | N/A | 1.632 | 一致（beta.1 大规模跳过） |
| triclinic | 2.2 | 4,096 | 839.460 | 9.884 | 9.726 | 86.315 | 1.016 | 一致 |
| triclinic | 2.2 | 32,768 | 102872.726 | 87.648 | 82.612 | 1245.252 | 1.061 | 一致 |
| triclinic | 2.2 | 131,072 | N/A | 385.928 | 329.613 | N/A | 1.171 | 一致（beta.1 大规模跳过） |
| triclinic | 3.4 | 4,096 | 839.059 | 27.013 | 18.393 | 45.618 | 1.469 | 一致 |
| triclinic | 3.4 | 32,768 | 103614.053 | 273.945 | 166.281 | 623.126 | 1.647 | 一致 |
| triclinic | 3.4 | 131,072 | N/A | 1130.997 | 685.332 | N/A | 1.650 | 一致（beta.1 大规模跳过） |

## 版本整体对比

| Case | Ranks | beta.1 中位时间 (s) | 当前版中位时间 (s) | 加速比 | beta.1 RSS (MiB) | 当前版 RSS (MiB) | 能量差 (eV) |
|---|---:|---:|---:|---:|---:|---:|---:|
| si16 | 1 | 48.906 | 48.622 | 1.006 | 470.016 | 470.844 | 0.00000000 |
| si64 | 1 | 220.738 | 219.600 | 1.005 | 1551.422 | 1553.066 | 0.00000000 |

## 当前版本核心路径消融

| Atoms | Half14 search (ms) | Full27 search (ms) | Full27/Half14 | Half14+reference (ms) | reference/Half14 | Half14 RSS (MiB) | reference RSS (MiB) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 4,096 | 3.660 | 6.139 | 1.677 | 11.406 | 3.116 | 27.363 | 32.195 |
| 32,768 | 33.986 | 69.565 | 2.047 | 108.008 | 3.178 | 66.875 | 104.000 |
| 131,072 | 141.951 | 317.243 | 2.235 | 478.841 | 3.373 | 207.824 | 288.500 |

## MPI 单节点扩展

### orthogonal 强扩展

| Ranks | 总原子数 | Exchange max (ms) | Search max (ms) | 加速比/弱扩展比 | 效率 | Ghost 总数 | Payload sent 总数 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 65,536 | 0.003 | 64.615 | 1.000 | 1.000 | 0 | 0 |
| 2 | 65,536 | 0.882 | 31.971 | 1.967 | 0.983 | 4096 | 36864 |
| 4 | 65,536 | 1.650 | 17.987 | 3.291 | 0.823 | 12800 | 115200 |
| 8 | 65,536 | 4.776 | 13.928 | 3.455 | 0.432 | 22592 | 203328 |
| 16 | 65,536 | 3.115 | 12.630 | 4.104 | 0.257 | 27776 | 249984 |

### triclinic 强扩展

| Ranks | 总原子数 | Exchange max (ms) | Search max (ms) | 加速比/弱扩展比 | 效率 | Ghost 总数 | Payload sent 总数 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 65,536 | 0.004 | 51.489 | 1.000 | 1.000 | 0 | 0 |
| 2 | 65,536 | 0.912 | 26.172 | 1.901 | 0.951 | 4096 | 36864 |
| 4 | 65,536 | 1.786 | 14.868 | 3.092 | 0.773 | 12800 | 115200 |
| 8 | 65,536 | 5.679 | 12.605 | 2.816 | 0.352 | 22592 | 203328 |
| 16 | 65,536 | 3.330 | 9.901 | 3.892 | 0.243 | 27776 | 249984 |

### orthogonal 弱扩展

| Ranks | 总原子数 | Exchange max (ms) | Search max (ms) | 加速比/弱扩展比 | 效率 | Ghost 总数 | Payload sent 总数 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 4,096 | 0.002 | 3.434 | 1.000 | 1.000 | 0 | 0 |
| 2 | 8,192 | 0.167 | 4.101 | 0.805 | 0.805 | 1024 | 9216 |
| 4 | 16,384 | 0.536 | 4.593 | 0.670 | 0.670 | 4352 | 39168 |
| 8 | 32,768 | 2.674 | 8.966 | 0.295 | 0.295 | 13888 | 124992 |
| 16 | 65,536 | 3.152 | 12.882 | 0.214 | 0.214 | 27776 | 249984 |

### triclinic 弱扩展

| Ranks | 总原子数 | Exchange max (ms) | Search max (ms) | 加速比/弱扩展比 | 效率 | Ghost 总数 | Payload sent 总数 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 4,096 | 0.002 | 3.133 | 1.000 | 1.000 | 0 | 0 |
| 2 | 8,192 | 0.150 | 3.260 | 0.919 | 0.919 | 1024 | 9216 |
| 4 | 16,384 | 0.484 | 3.410 | 0.805 | 0.805 | 4352 | 39168 |
| 8 | 32,768 | 3.045 | 6.105 | 0.343 | 0.343 | 13888 | 124992 |
| 16 | 65,536 | 3.204 | 9.727 | 0.242 | 0.242 | 27776 | 249984 |


## Verlet 增量重建消融

| Case | skin=0 时间 (s) | skin=3 时间 (s) | 时间比 (skin0/skin3) | skin=0 重建/复用 | skin=3 重建/复用 | 能量差 (eV) |
|---|---:|---:|---:|---:|---:|---:|
| md_si64 | 301.000 | 297.000 | 1.013 | 20/0 | 1/19 | 0.00000000 |
| md_si8 | 17.000 | 17.000 | 1.000 | 10/0 | 1/9 | 0.00000000 |
| relax_si2 | 6.000 | 7.000 | 0.857 | 4/0 | 1/3 | 0.00000000 |

## 回归测试

- 完整 `MODULE_CELL` CTest 36/36 通过，包括 MPI wrapper。
- `MODULE_CELL_NEIGHBOR_mpi_grid` 在 1/2/4/8/16 ranks 下通过。
- `MODULE_LCAO_record_adj_mpi` 在 1/2/4 ranks 下通过。
- no-MPI 配置下 `cell`、`neighbor`、`hamilt_lcao` 和 `esolver` 构建通过。

## 正确性与限制

- Half14、Full27 和显式 Full27 reference 的 pair hash 必须逐规模一致。
- 版本对比仅在两个版本均正常结束且物理结果满足容差时计入。
- 峰值 RSS 是完整独立进程的最大常驻集，适合版本/模式相对比较，不等同于容器理论容量。
- WSL 单节点上的 MPI 数据主要用于观察算法和通信趋势，不能替代集群多节点 Benchmark。

## 原始结果

- `docs/benchmark_results/phase4/environment.json`：硬件、软件、命令参数和提交号。
- `docs/benchmark_results/phase4/summary.json`：各配置的中位数、最小值和最大值。
- `docs/benchmark_results/phase4/version.csv`：beta.1 与当前版完整 ABACUS 运行结果。
- `docs/benchmark_results/phase4/grid_paths.csv`：beta.1 legacy、current Full27 与 current Half14 的直接对比。
- `docs/benchmark_results/phase4/kernel.csv`：Half14、Full27 与 reference 消融结果。
- `docs/benchmark_results/phase4/mpi.csv`：MPI local/ghost、payload 和耗时结果。
- `docs/benchmark_results/phase4/verlet.csv`：skin=0/3 的真实 relax/MD 对照。
- 完整标准输出保存在本地 Benchmark 工作目录的 `logs/`，不纳入仓库。
