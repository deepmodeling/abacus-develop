# LCAO DeltaSpin Subspace 加速：修改记录与使用指南

## 一、Commit 历史与修改内容

### 1. `9c49503` feat(deltaspin): opt-in subspace acceleration for LCAO spin-constrained DFT

**核心改动**：将 subspace/first_order 加速改为 opt-in 模式，全对角化恢复为默认路径。

- 新增输入参数 `sc_acceleration_mode`（off/first_order/subspace）和 `sc_acceleration_rms_thr`
- 在 `cal_mw_from_lambda()` 中实现三条执行路径：
  1. `i_step == -2`：全对角化 + 构建子空间缓存（仅首次激活时调用）
  2. `accel_enabled + subspace_built`：加速路径（subspace/first_order）
  3. 其他：全 HSolverLCAO 对角化（默认）
- 修复关键 bug：subspace 对角化不能替代全求解器，只应在 RMS 低于阈值后用于小扰动近似
- 添加 `diagonalization_engine.cpp` 到 CMakeLists.txt

### 2. `b2c9d04` fix(deltaspin): fix LCAO subspace acceleration bugs and add diagnostic

**三个关键修复**：

**① 特征值偏移符号错误**
```cpp
// 错误：
ekb_new = ekb_old + spin_sign * delta_lambda * P_diag;
// 正确（与 calculate_delta_hcc_lcao 一致）：
ekb_new = ekb_old - spin_sign * delta_lambda * P_diag;
```

**② 跨 SCF 迭代的子空间缓存过期**
- 问题：一个 SCF 迭代构建的子空间被复用到后续迭代，此时 H/S/C 已变化，导致 Mi 冻结
- 修复：在 `run_lambda_loop()` 入口处重置所有加速状态标志并释放缓存

**③ 加速激活后缺少最终全对角化**
- 问题：`update_psi_charge()` 在加速激活时未做全对角化，psi 停留在参考 lambda 处
- 修复：如果加速曾经激活过，在收敛 lambda 处执行一次全对角化后再 psiToRho

**新增诊断工具**：`run_trace_vs_dmr_diagnostic()` 对比三种 Mi 计算方法随 Δλ 的变化。

### 3. `e61fc8f` fix(deltaspin): replace cal_PI_sub with calculate_PI_sub_from_hr

**核心问题**：旧的 `cal_PI_sub` 使用 $D_I^\dagger D_I$ 路径构建子空间投影矩阵，与 DMR 路径的数据源不一致。

- 旧路径：`B_I_data`（onsite 重叠积分缓存）→ $D_I$ → $P_{I,\text{sub}} = D_I^\dagger D_I$
- 问题：`B_I_data` 的邻域截断与 `pre_hr` 不完全一致
- 实测偏差：BCC Fe AFM 中原子 1 的 Mi 差 4.66 μB（符号反转）

**新方案**：`calculate_PI_sub_from_hr` 使用与 $H_\text{sub}$ 完全对称的计算流程：
1. `folding_HR(pre_hr[iat])` → $P_I(k)$（2D-block 分布）
2. `pzgemm('N','N')`：temp = $P_I(k) \times C(k)$
3. `pzgemm('C','N')`：$P_{I,\text{sub}} = C^\dagger(k) \times$ temp
4. `gather_sub_matrix_to_all`：汇总到所有进程

**关键保证**：$P_I(k)$ 来自 `pre_hr[iat]`，与 `cal_moment` 使用的数据源完全相同，因此迹公式与 DMR 路径精确一致。

**新增访问器**：
```cpp
dspin_op->get_pre_hr(iat)        // 获取实空间投影算符 HContainer
dspin_op->get_constraint_atom_list()  // 获取约束原子列表
```

### 4. `5735ea6` feat(deltaspin): simplify execution strategy with sc_strategy parameter

**简化用户接口**：用 `sc_strategy` 参数替代手动配置 `sc_acceleration_mode` + `sc_acceleration_rms_thr`。

| 策略 | sc_acceleration_mode | sc_acceleration_rms_thr | 行为 |
|------|---------------------|------------------------|------|
| `normal`（默认） | subspace | 1e-2 | 全对角化起步，RMS 降到阈值后自动切换 |
| `fast` | subspace | 1e10 | 第一步即激活 subspace（跳过收敛判断） |
| `accuracy` | off | -1.0 | 始终全对角化，不启用加速 |

**默认值变更**：
- `nsc`: 100 → 5
- `sc_scf_thr`: 1e-3 → 10
- `sc_scf_thr_mode`: threshold → immediate

用户仍可显式设置 `sc_acceleration_mode` 和 `sc_acceleration_rms_thr` 覆盖策略默认值。

---

## 二、算法流程

```
run_lambda_loop 入口
│
├── 重置加速状态（每个 SCF 迭代独立）
│
├── BFGS 内循环 (i_step = -1, 0, 1, ..., nsc-1)
│   │
│   ├── i_step = -1: 全对角化计算初始 Mi
│   │
│   ├── i_step >= 0:
│   │   ├── 更新 lambda = initial + delta_lambda
│   │   ├── cal_mw_from_lambda(i_step)
│   │   │   │
│   │   │   ├── [加速首次激活] RMS < sc_acceleration_rms_thr
│   │   │   │   ├── 全对角化 → C_ref, ε_ref
│   │   │   │   ├── calculate_lcao_sub_hs → H0_sub, S_sub
│   │   │   │   ├── calculate_PI_sub_from_hr → P_I_sub
│   │   │   │   │   ├── folding_HR(pre_hr[iat]) → P_I(k)
│   │   │   │   │   ├── pzgemm: temp = P_I(k) × C(k)
│   │   │   │   │   ├── pzgemm: P_I_sub = C† × temp
│   │   │   │   │   └─ gather_sub_matrix_to_all
│   │   │   │   └── 迹公式/DMR路径 计算 Mi(V=I)
│   │   │   │
│   │   │   ├── [子空间加速] acceleration_active && subspace_built
│   │   │   │   ├── H_sub = H0_sub + δH
│   │   │   │   ├── diag_hegvd → V, ε  (本地 zhegvd, O(n³))
│   │   │   │   ├── calculate_weights(ε)
│   │   │   │   └── DMR路径 Mi (rotate_psi → cal_dm_psi → cal_DMR → cal_mi_lcao)
│   │   │   │
│   │   │   └── [全对角化] 默认路径
│   │   │       └── HSolverLCAO::solve → calculate_weights → cal_mi_lcao
│   │   │
│   │   └── 收敛检查 / CG 方向更新 / 试探步
│   │
│   └── 收敛或达到最大步数
│
├── update_psi_charge
│   ├── 如果子空间激活过：最终全对角化(收敛 lambda) → psiToRho
│   └── 否则：直接 psiToRho
│
└── 返回 SCF 主循环
```

### 并行架构

| 阶段 | 复杂度 | 并行方式 | 通信 |
|------|--------|----------|------|
| 构建子空间 | $O(N^2 \cdot n)$ | ScaLAPACK pzgemm (2D-block) | k-pool 内 MPI_Allreduce |
| 子空间对角化 | $O(n^3)$ | 本地 LAPACK zhegvd（无通信） | 无 |
| 波函数旋转 | $O(N \cdot n^2)$ | ScaLAPACK pzgemm | k-pool 内 |
| Mi 计算 | 同 cal_mi_lcao | DMR 路径 | k-pool 内 |

**KPAR 支持**：每个 k-pool 独立处理自己的 k 点集，子空间构建和对角化天然并行。

---

## 三、使用方式

### 推荐：sc_strategy（三个预设）

```text
# INPUT
sc_mag_switch    1
sc_thr           1e-4
nsc              50
nsc_min          2
alpha_trial      0.01
sccut            3.0
sc_scf_thr       10

# 选择执行策略（三选一）
sc_strategy      fast       # 最快：subspace 从第一步开始
sc_strategy      normal     # 平衡：RMS < 1e-2 后自动切换（默认）
sc_strategy      accuracy   # 最精确：始终全对角化
```

### 高级：手动配置

```text
# 手动覆盖 sc_strategy 的默认值
sc_acceleration_mode    subspace      # off / first_order / subspace
sc_acceleration_rms_thr 0.5           # RMS 低于此值时激活加速
```

### 典型场景

| 场景 | 推荐配置 |
|------|----------|
| 大体系快速扫描 | `sc_strategy fast` |
| 生产计算 | `sc_strategy normal` |
| 调试/验证 | `sc_strategy accuracy` |
| 需要精细控制阈值 | `sc_strategy normal` + 手动设 `sc_acceleration_rms_thr` |
