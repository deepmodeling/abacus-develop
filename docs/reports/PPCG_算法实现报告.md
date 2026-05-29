# ABACUS PPCG（Projected Preconditioned Conjugate Gradient）算法实现与单元测试报告

> 项目：abacus-develop（HSolver 子模块）
>
> 分支：PPCG
>
> 日期：2026-05-29

## 1. 背景与目标

在 ABACUS 的本征值求解模块（HSolver）中，已存在 CG/BPCG（Block Preconditioned Conjugate Gradient）等对角化求解器。为提升子空间迭代求解能力并丰富算法选型，本工作实现 PPCG（Projected Preconditioned Conjugate Gradient，投影预条件共轭梯度）求解器，并优先配套单元测试以验证正确性。

本阶段目标：

1. 参照现有 CG/BPCG 的工程结构与接口风格，实现 PPCG 求解器类。
2. 将 PPCG 接入 CMake/CTest，补充与 BPCG 类似风格的单元测试。
3. 优先跑通编译与测试框架（可运行），并逐步修正数值问题使测试通过。

## 2. 算法概述（实现采用的思路）

本实现采用 LOBPCG/PPCG 常见的“子空间投影 + 广义 Rayleigh-Ritz（RR）”框架。

### 2.1 基本符号

- 目标：求解 Hermitian 本征问题 $H x = \lambda x$（单元测试里采用稠密 Hermitian 矩阵）。
- $X \in \mathbb{C}^{n\times b}$：当前 block 近似本征向量（b = nband）。
- $HX = H X$。
- 残差：$R = HX - X\Lambda$（$\Lambda$ 为对角 Ritz 值）。
- 预条件方向：$W \approx -M^{-1}R$，其中 $M$ 为对角预条件器。
- 共轭方向：$P$（上一轮的搜索方向/子空间补充）。

### 2.2 子空间构造与投影 RR

每次外层迭代构造子空间：

- 首次迭代：$V = [X, W]$（列数 $2b$）
- 后续迭代：$V = [X, W, P]$（列数 $3b$）

并计算投影矩阵：

- $H_c = V^\dagger (H V) = V^\dagger HV$
- $S_c = V^\dagger V$

解广义本征值问题：

$$(H_c) c = (S_c) c \Lambda$$

取对应最小的 $b$ 个本征对，更新：

- $X \leftarrow V c_{1:b}$
- $HX \leftarrow HV c_{1:b}$

并按系数块更新搜索方向 $P$（来自 $W,P$ 部分）。

### 2.3 投影与正交化策略

为避免子空间病态与方向退化，实现中使用：

- 投影：将 $W$（以及更新后的 $P$）投影到 $X$ 与 $P$ 的补空间。
- 块正交化（Cholesky）：对 $P$、$W$ 做块正交化以改善条件数。

注意：若对 $W$ 做块正交化，则必须对 $HW$ 做一致变换，保持 $HW = H W$，否则投影矩阵 $V^\dagger HV$ 不再对应真实子空间。

## 3. 工程设计与文件结构

### 3.1 新增/修改的核心文件

- `source/source_hsolver/diago_ppcg.h`
  - 定义 `hsolver::DiagoPPCG` 类。
  - 对齐 BPCG 风格：`init_iter()` + `diag()`，并接收 `HPsiFunc` 形式的矩阵-向量（块）乘。

- `source/source_hsolver/diago_ppcg.cpp`
  - PPCG 主流程实现：
    - 初始 RR（仅在 $X$ 子空间上）
    - 外层迭代：残差/预条件、构造子空间、投影 RR、更新 $X/P$、收敛检查
  - 复用/对齐内核：
    - 使用 `hsolver::normalize_op / precondition_op / apply_eigenvalues_op`（来自 `source/source_hsolver/kernels/bpcg_kernel_op.*`）
    - 使用 `ModuleBase::gemm_op / axpy_op / dot_real_op` 等基础算子

- `source/source_hsolver/test/diago_ppcg_test.cpp`
  - PPCG 单元测试：
    1. `TwoByTwo`：2x2 Hermitian 矩阵（应快速正确）
    2. `readH`：读取数据文件 `H-KPoints-Si2.dat` 并与 LAPACK 对比
    3. `RandomHamilt`：随机 Hermitian（通过 LAPACK `zheev_` 得到参考本征值）

- `source/source_hsolver/test/CMakeLists.txt`
  - 新增 `MODULE_HSOLVER_ppcg` 测试 target，并通过 CTest 注册。

- `source/source_hsolver/CMakeLists.txt`
  - 将 `diago_ppcg.cpp` 加入 hsolver objects。

### 3.2 与 BPCG/CG 的接口一致性

`DiagoPPCG` 的外部接口与 `DiagoBPCG` 对齐：

- `init_iter(nband, nband_l, nbasis, ndim)`：初始化问题规模与 workspace
- `diag(hpsi_func, psi_in, eigenvalue_out, ethr_band)`：执行对角化/迭代

测试中的 `hpsi_func` 写法与 BPCG 单元测试保持一致，均通过 `ModuleBase::gemm_op` 完成稠密矩阵乘。

## 4. 单元测试设计与运行方式

### 4.1 测试判据

单元测试使用 LAPACK 输出作为参考，逐带比较：

- `EXPECT_NEAR(en[i], e_lapack[i], threshold)`

其中 `threshold` 随测试用例设置（例如 `TwoByTwo` 更严格，`RandomHamilt/readH` 较宽松）。

### 4.2 运行命令

在已 configure 的 build 目录下运行：

```bash
cmake --build build -j8 --target MODULE_HSOLVER_ppcg
ctest --test-dir build -V -R MODULE_HSOLVER_ppcg
```

## 5. 当前进度与结果（截至 2026-05-29）

### 5.1 已完成

- PPCG 求解器代码已完成“可编译、可链接、可运行”状态。
- `MODULE_HSOLVER_ppcg` 测试可以被 CTest 发现并执行。
- `TwoByTwo` 用例已通过。

### 5.2 当前问题（测试失败现象）

- `readH` 与 `RandomHamilt` 仍失败：计算得到的本征值与 LAPACK 参考值偏差较大。
- 在失败输出中，部分 `en[i]` 会出现接近 0 或极小值（如 `~1e-310`），表明当前迭代结果可能未正确收敛或某些更新步骤仍存在数值/布局错误。

### 5.3 已定位并修复过的关键工程性问题

- 内核接口签名：`normalize_op/precondition_op/apply_eigenvalues_op` 的调用方式与其真实接口不一致（已按 `bpcg_kernel_op.cpp` 真实签名修正）。
- `HW` 一致性：在对 $W$ 进行块正交化时同步对 $HW$ 施加同变换，保持 $HW=HW$ 的物理含义。
- 去除不必要依赖：移除 PPCG 中对 `DiagoBPCG` 的 fallback 依赖，避免测试 target 链接错误，并保证单测真正测试 PPCG 本身。

## 6. 根因分析（当前仍需继续攻关的数值点）

结合现有现象与实现流程，当前 PPCG 单测失败可能来自以下一个或多个原因（需进一步通过日志与断点验证）：

1. **投影/正交化策略是否与 RR 一致**：
   - `project_out()` 当前采用 `coeff = basis^H vecs`，默认 basis 列正交归一；若某一步 basis 未严格正交，投影会偏离。

2. **子空间系数块（vcc）的使用是否与 LAPACK 返回布局匹配**：
   - `hegvd_op` 输出 `vcc` 为列主序本征向量；在 `update_from_projected()` 中对系数块的行/列偏移必须严格正确。

3. **收敛与阈值设置**：
   - PPCG 外层迭代上限来自 `DiagoIterAssist::PW_DIAG_NMAX`；若算法参数或更新策略不当，可能需要更多迭代或更稳健的正交策略。

## 7. 后续计划

为尽快跑通单测（与 LAPACK 对齐），后续建议按以下顺序推进：

1. 在 `diag()` 每轮迭代打印/记录：`eval[0..b)`、`||R||` 与 `not_conv` 变化，确认迭代是否在正确下降。
2. 对 `project_out()` 改为严格投影（基于 $S = basis^H basis$ 解小线性系统），或确保 basis 在投影前块正交化。
3. 复核 `update_from_projected()` 中 `P/HP` 更新公式是否正确（系数块切片与 stride）。
4. 逐步调小测试规模并与 LAPACK 比对中间量（例如对 $H_c,S_c$ 做一致性检查）。

## 8. 附录：关键实现要点摘录

- PPCG 子空间：`V=[X,W,P]`（或首轮 `V=[X,W]`）
- RR 求解：通过 `hsolver::hegvd_op` 解 $(V^\dagger HV)c=(V^\dagger V)c\Lambda$
- 预条件：`precondition_op` 使用对角预条件向量与 Ritz 值近似构造

---

（本报告为阶段性实现与测试进度总结；算法数值正确性与鲁棒性仍在迭代完善中。）
