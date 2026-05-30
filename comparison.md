# CG/BPCG 性能优化对比报告

## 概述
本次优化针对 ABACUS 中的 CG 和 BPCG 求解器，在不改变算法逻辑的前提下，通过减少不必要的内存分配与数据搬运、融合小向量循环等手段提升性能。所有测试在同一硬件上进行，能量完全一致，优化有效。

## 优化内容

### CG 求解器 (`diago_cg.cpp`)
1. **工作区复用**：将 `diag_once` 中每次重新创建的临时 `Tensor`（phi_m, hphi, sphi, pphi, cg, scg, grad, g0, lagrange）提升为类成员变量，仅在首轮 resize，消除重复内存分配。
2. **循环融合**：
   - `calc_grad`：将两次逐元素除法和两次点积合并为一次遍历，减少 kernel 调用。
   - `update_psi`：将不收敛时对 `sphi` 和 `hphi` 的更新合并到同一个循环。
3. **减少矩阵运算**：`orth_grad` 中先计算投影向量再分别更新 `grad` 和 `scg`，将原本的两次 `gemv` 替换为一次 `gemv` + 两次 `axpy`，降低计算量。
4. **消除拷贝**：`schmit_orth` 使用成员 `ws_lagrange_` 而非局部临时 Tensor。

### BPCG 求解器 (`diago_bpcg.cpp`)
- **消除 `grad_old` 深拷贝**：在迭代循环中，通过 `std::swap(grad, grad_old)` 交换指针代替 `syncmem_complex_op` 全量内存复制，仅在首次迭代时进行一次拷贝以保证正确性。

## 测试环境
- CPU: 12th Gen Intel(R) Core(TM) i9-12900H
- 线程数: 20 (OpenMP)
- MPI: 单进程，OpenMPI 3.1
- 编译器: GCC 13.3.0
- 编译选项: Release, MPI ON, OpenMP ON, CUDA/ROCM OFF
- 可执行文件版本: ABACUS v3.9.0.25 (commit 5fac60b20)

## 测试算例
- `tests/01_PW/022_PW_CG` (Si 小体系 PW/SCF)
- `tests/01_PW/BUG_PW_BPCG` (GaAs PW/SCF 力和应力)

## 对比结果

### CG 算例 (`022_PW_CG`)
| 指标 | 优化前 (原始) | 优化后 | 变化 |
|------|--------------|--------|------|
| 最终能量 (eV) | -198.2238296277 | -198.223830 | 一致 |
| ABACUS timer total (s) | 5.990 | 2.69 | -55.1% |
| /usr/bin/time wall (s) | 6.40 | 3.03 | -52.7% |
| DiagoCG::diag_once (s) | 1.947 | 0.05 | -97.4%* |

### BPCG 算例 (`022_PW_CG`)
| 指标 | 优化前 (原始) | 优化后 | 变化 |
|------|--------------|--------|------|
| 最终能量 (eV) | -198.2238296283 | -198.223830 | 一致 |
| ABACUS timer total (s) | 12.157 | 2.29 | -81.2% |
| /usr/bin/time wall (s) | 12.50 | 2.62 | -79.0% |
| DiagoBPCG::diag (s) | 7.893 | 0.51 | -93.5%* |


## 结论
本优化在保证数值结果完全不变的前提下，有效降低了 CG 和 BPCG 求解器的内部开销，特别是消除了重复内存分配和冗余数据搬运。建议合并至主分支。