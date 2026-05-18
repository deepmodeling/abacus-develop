# DeltaSpin nspin=2 (Collinear) 磁约束方案

## 概述

在 nspin=2（共线磁）模式下，DeltaSpin 的磁约束能力与其他参数组合有根本性差异。核心原因是：**共线磁只有 z 分量**，`constrain_` 向量被强制为 `(0, 0, Mz)`，导致 `sc_direction_only` 投影将 lambda 归零。本文档总结所有 nspin=2 磁约束方案的参数组合、行为和适用场景。

---

## 参数速查表

| 参数 | 默认值 | 选项 | 说明 |
|---|---|---|---|
| `sc_mag_switch` | 0 | 1 | 开启 DeltaSpin |
| `sc_scf_thr_mode` | `"threshold"` | `"threshold"` / `"immediate"` / `"off"` | 控制何时激活 lambda 循环 |
| `sc_scf_thr` | 1e-3 | 正实数 | 仅 `threshold` 模式下有效：drho 低于此值时激活 lambda 循环 |
| `sc_direction_only` | false | true/false | 仅约束磁矩方向，不约束大小 |
| `sc_dir_phase1_steps` | 5 | ≥2 整数 | `direction_only + nspin=2` 两阶段策略的 Phase 1 迭代步数 |
| `mixing_restart` | 0 (自动设置) | 实数 | 根据 `sc_scf_thr_mode` 自动设置（见下表） |

### mixing_restart 自动设置

| 条件 | mixing_restart 值 | 说明 |
|---|---|---|
| `direction_only=true` | 0（禁用） | 两阶段策略在 Phase 1→2 过渡时自行管理 mixing reset |
| `sc_scf_thr_mode="threshold"` | `= sc_scf_thr` | 在 lambda 循环激活时重启混合 |
| `sc_scf_thr_mode="immediate"` | `= scf_thr / 10` | 早期激活，轻微提前重启混合 |
| `sc_scf_thr_mode="off"` | 0（禁用） | 无 lambda 优化，无需重启混合 |

---

## 方案一：标准约束（约束磁矩大小和方向）

**参数组合：** `sc_direction_only=false` + `sc_scf_thr_mode` 选一种

### 1a. threshold 模式（默认，推荐）

```
sc_mag_switch    1
sc_scf_thr_mode  threshold
sc_scf_thr       1e-3        # 建议 10-100x scf_thr
```

**行为：**
- SCF 迭代中，当电荷密度误差 `drho < sc_scf_thr` 时激活 lambda 优化循环
- lambda 循环通过 BFGS 优化 lambda_z，使磁矩 Mz 逼近目标值
- 一旦激活（`mag_converged=true`），后续每步都继续运行 lambda 循环
- `mixing_restart` 自动设为 `sc_scf_thr`，在 lambda 循环激活时重启混合

**适用场景：** 一般用途，lambda 在 SCF 收敛后逐步优化磁矩

### 1b. immediate 模式（立即激活）

```
sc_mag_switch    1
sc_scf_thr_mode  immediate
# sc_scf_thr 不再需要设置，immediate 模式不检查 drho
```

**行为：**
- 从第 2 步 SCF 迭代起立即激活 lambda 循环（第 1 步无波函数，无法计算磁矩）
- 等效于旧版 `sc_scf_thr=10`（而 `drho < 10` 几乎总成立）
- `mixing_restart` 自动设为 `scf_thr / 10`

**适用场景：** PW 基组（第 1 步无法计算磁矩）、希望尽早引入约束的场景

### 1c. off 模式（关闭 lambda 优化）

```
sc_mag_switch    1
sc_scf_thr_mode  off
# sc_scf_thr 不再需要设置，off 模式不运行 lambda 循环
```

**行为：**
- 永不激活 lambda 循环
- lambda 值从 STRU 文件读入，作为常量约束加到哈密顿量中
- 正常 SCF 收敛即可，磁约束力固定不变
- `mixing_restart` 自动设为 0（禁用）

**适用场景：** 已知 lambda 值，只需求解固定约束下的电子结构；等效于旧版 `sc_scf_thr=1e-10`

---

## 方案二：仅约束方向（direction_only）+ nspin=2

**参数组合：** `sc_direction_only=true` + `nspin=2` + `sc_dir_phase1_steps ≥ 2`

```
sc_mag_switch         1
sc_direction_only     1
sc_dir_phase1_steps   5          # Phase 1 迭代步数，默认 5
# sc_scf_thr 和 sc_scf_thr_mode 在两阶段策略中不被检查
```

### 核心问题

对于 nspin=2，`direction_only` 投影会将 lambda **完全归零**：

- 目标方向 `dir = target / |target| = (0, 0, 1)`（共线只有 z 分量）
- 投影移除 lambda 的平行分量：`lambda_z -= lambda_z × 1 = 0`
- 结果：lambda = 0，约束力为零，BFGS 优化器无法工作

### 两阶段策略

Phase 1 + Phase 2 解决方案：

| 阶段 | 步骤 | 行为 | 效果 |
|---|---|---|---|
| **Phase 1** | iter ≤ `sc_dir_phase1_steps` | 临时禁用 `direction_only`，运行 BFGS 约束磁矩大小 Mz 逼近目标值 | 建立正确的反铁磁序和合理范围的 lambda (~0.3 eV/μB) |
| **Phase 1→2 过渡** | iter = `sc_dir_phase1_steps + 1` | 重置混合历史（`mix_reset()`），清除 Broyden 缓冲 | Phase 1 的 BFGS DM 更新与 Phase 2 的 SCF 电荷混合不兼容 |
| **Phase 2** | iter > `sc_dir_phase1_steps` | lambda 每步衰减 `λ *= 0.5^(1/3) ≈ 0.794`（约每 3 步减半） | 约束力逐渐释放，系统自然弛豫到磁基态 |

### ⚠️ 风险提示

> **`direction_only + nspin=2` 方案在复杂磁结构体系中可能无法保持磁矩方向不变，特别是磁矩大小 < 1.0 μB 的原子。**
>
> 原因：
> 1. Phase 2 中 lambda 衰减 → 约束力减弱 → 自旋可能翻转到其他方向
> 2. 磁矩小的原子，交换劈裂弱，更容易受微扰影响
> 3. 该方案本质上是一种近似：先锁定方向，再逐渐松开，希望系统自然选择正确方向
>
> **优势是计算速度快**——Phase 2 无需运行 BFGS 内循环，仅做简单标量乘法。

### PW 基组的限制

> **⚠️ PW 基组 + nspin=2 + `direction_only` 目前不支持两阶段策略。**
>
> PW 路径 `deltaspin_pw.cpp` 中没有 Phase 1/2 Logic，`direction_only` 投影会将 lambda 归零，约束完全失效。
> 如需在 PW 基组中使用 `direction_only`，请使用 **nspin=4**。

### 参数注意事项

| 参数 | 注意 |
|---|---|
| `sc_scf_thr` | 两阶段策略不检查此参数，设置与否均无影响 |
| `sc_scf_thr_mode` | 两阶段策略不检查此参数，设置与否均无影响 |
| `mixing_restart` | 被强制设为 0（禁用），Phase 1→2 过渡时代码自行 reset |
| `nsc` | Phase 1 步数由 `sc_dir_phase1_steps` 控制，不是 `nsc`；`nsc` 控制 BFGS 内循环最大迭代 |
| `sc_dir_phase1_steps` | 必须 ≥ 2；值越大 Phase 1 越充分，但 Phase 2 起始 lambda 也可能过大 |

---

## 方案三：线性扫描（lambda_strategy = linear_scan）

```
sc_mag_switch        1
sc_lambda_strategy   linear_scan
```

**行为：** lambda 从初始值开始按固定步长递增，每步扫描一个 lambda 值。

> 注意：此策略代码尚未完成（已被删除），仅供参考。

---

## nspin=2 参数组合决策树

```
是否需要 DeltaSpin？
├─ 否 → sc_mag_switch = 0
└─ 是 → 是否需要优化 lambda？
    ├─ 否（使用已知 lambda）→ sc_scf_thr_mode = "off"
    └─ 是 → 是否仅约束方向？
        ├─ 是 → sc_direction_only = 1
        │    └─ nspin = 2 → 两阶段策略（见方案二）
        │    └─ nspin = 4 → sc_scf_thr_mode = "immediate" 或 "threshold"
        └─ 否 → 约束大小和方向
             ├─ PW 基组 → sc_scf_thr_mode = "immediate"
             ├─ LCAO 基组 → sc_scf_thr_mode = "threshold"（默认）
             └─ 两者均可使用 "immediate"
```

---

## 测试覆盖

| 测试 | nspin | 基组 | 方案 | sc_scf_thr_mode |
|---|---|---|---|---|
| 12_PW_DS_S2_Z | 2 | PW | 标准约束 | immediate |
| 24_LCAO_DS_S2_Z | 2 | LCAO | direction_only | — (两阶段) |
| 36_PW_DS_S2_ReadLam_Z | 2 | PW | off (ReadLam) | immediate |
| 38_PW_DS_S2_Thr1e10_Z | 2 | PW | off (常量 lambda) | off |
| 40_PW_DS_S2_Thr10_Z | 2 | PW | 标准约束 | immediate |
| 44_PW_DFTU_DS_S2_Thr10_Z | 2 | PW | 标准+DFT+U | immediate |
| 65_LCAO_DS_S2_DirectionOnly_Z | 2 | LCAO | direction_only | — (两阶段) |

---

## 参数迁移对照表（旧 → 新）

| 旧参数 | 含义 | 新参数 |
|---|---|---|
| `sc_scf_thr = 10` | 立即激活 lambda 循环（magic number） | `sc_scf_thr_mode = "immediate"` |
| `sc_scf_thr = 1e-10` | 永不激活 lambda 循环（magic number） | `sc_scf_thr_mode = "off"` |
| `sc_scf_thr = 1e-3` | threshold 模式正常阈值 | `sc_scf_thr = 1e-3` + `sc_scf_thr_mode = "threshold"`（默认） |