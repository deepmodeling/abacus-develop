# DeltaSpin Lambda 优化方案：数值稳定性分析与改进

## 1. 失败算例特征

| 类别 | 算例 | 基组 | DFT+U | 约束方向 | 状态 |
|------|------|------|-------|---------|------|
| 纯 DS (PW) | 14, 17 | PW | 否 | XYZ [111] | **FAIL** |
| DFT+U + DS (PW) | 19, 22 | PW | 是 | XY [100] | **FAIL** |
| 纯 DS (LCAO) | 25, 26, 27, 28, 29 | LCAO | 否 | XY/XYZ/Z | **全部FAIL** |
| DFT+U + DS (LCAO) | 31, 32, 33, 34, 35 | LCAO | 是 | XY/XYZ/Z | **全部FAIL** |
| DirectionOnly | 58 | LCAO | 否 | XY | **FAIL** |

**共同特征**：全部使用 BCC Fe，2原子原胞，`lambda=4`，`sc_thr=1e-4`，非共线自旋，关闭对称性。

## 2. Case 14 根因分析

### 2.1 关键日志数据

| 步骤 | 磁矩偏差 (x,y,z) | 总磁矩 | 能量变化 | RMS |
|------|-----------------|--------|----------|-----|
| DS8 (SCF restart前) | `2.65e-02, 2.07e-02, 2.76e-02` | 8.23 | -0.0218 | 0.25 |
| DS9 (SCF restart后, 内循环初始) | `-2.93e-02, -3.53e-02, -2.82e-02` | 8.31 | -1.662 | 0.31 |

**内循环初始态 (DS9 step -1)**：
```
Atom 1 实际磁矩: (-1.80, -1.81, -1.80)  指向 [-1,-1,-1]
Atom 2 实际磁矩: ( 1.94,  1.95,  1.94)  指向 [ 1, 1, 1]

Atom 1 目标磁矩: ( 1.155,  1.155,  1.155)  指向 [ 1, 1, 1]
Atom 2 目标磁矩: (-1.155, -1.155, -1.155)  指向 [-1,-1,-1]

初始 RMS = 5.25 uB  -- 磁矩方向完全相反！
```

### 2.2 问题链条

```
SCF 接近收敛 (drho < sc_scf_thr=1e-3)
  |
  v
触发 SCF restart: 清空混合历史，重新开始
  |
  v
电荷密度落入相反的磁态 (能量简并的局部极小值)
  |-- [111] 非易磁化轴，能量景观比 [100]/[001] 更复杂
  |-- 自旋翻转需要跨越势垒 (~180° 旋转)
  |
  v
Lambda 内循环启动，初始 RMS = 5.25
  |
  v
内循环仅做子空间对角化 (H += becp†·λ·becp)
  |-- 电荷密度在整个内循环中保持不变
  |-- 用线性响应处理高度非线性问题
  |
  v
100步后 RMS 仅从 5.25 降到 1.22，无法收敛
```

### 2.3 根本原因

**子空间近似的局限**：`cal_mw_from_lambda` 在 PW 基组下只做子空间对角化，电荷密度固定。对于需要翻转磁矩（跨越势垒）的情况，子空间旋转只能在当前电子密度的线性响应范围内调整波函数，无法改变底层电荷分布，导致优化器被困在错误的磁态盆地中。

## 3. 优化方案

### 方案 2：初始态检测与预处理

**目标**：在进入 lambda 内循环前，检测磁矩方向是否严重偏离，若需要翻转则施加引导 lambda。

**修改文件**：`source/source_lcao/module_deltaspin/lambda_loop.cpp`

**算法流程**：
```
在 step -1 (初始化) 之后、step 0 (优化) 之前：

for each constrained atom ia:
    dot = spin[ia] · target_mag[ia]
    if dot < 0:  // 方向相反（夹角 > 90°）
        // 计算引导 lambda：指向目标方向
        residual = target_mag[ia] - spin[ia]
        guiding_lambda = guide_scale * |residual|
        lambda_[ia] += guiding_lambda * normalize(residual)
        
        // 可选：立即调用 cal_mw_from_lambda() 验证翻转效果
        cal_mw_from_lambda()
        if new_dot > 0:
            print "[DS-DIAG] Flip detection succeeded, pre-guiding applied"
```

**新增参数**：

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `sc_flip_detect` | bool | true | 是否启用翻转检测 |
| `sc_guide_lambda_scale` | double | 2.0 | 引导 lambda 的缩放系数 |

**参数输入方式**：通过 INPUT 文件或 STRU 的 `sc` 行扩展。

**优点**：
- 独立模块，不影响现有收敛路径
- 不增加额外 SCF 调用
- 对正常收敛情况无副作用（仅 dot < 0 时触发）

**限制**：
- 不保证一定能翻转成功（能量势垒过高时仍会失败）
- 需要配合方案 4 的优化器切换才能在大误差区有效推进

---

### 方案 4：自适应优化器切换

**目标**：根据当前 RMS 误差自动切换优化策略，避免大误差区 CG 方法的震荡。

**修改文件**：`source/source_lcao/module_deltaspin/lambda_loop.cpp`

**策略切换表**：

| RMS 范围 | 优化策略 | 说明 |
|----------|---------|------|
| RMS > `sc_rms_cg_switch_high` (3.0) | 纯最速下降 | `search = delta_spin`，无 CG 历史，无线性插值 |
| 中间区 (1.0 ~ 3.0) | PR-CG + 步长限制 | 使用 CG 搜索方向，禁用线性插值，仅用 `restrict_current_` 限制步长 |
| RMS < `sc_rms_cg_switch_mid` (1.0) | 当前完整算法 | PR-CG + 线性插值 + alpha 自适应 |

**代码修改位置**：`lambda_loop.cpp` 第 306-396 行（主循环搜索方向更新和 trial step 部分）

**伪代码**：
```cpp
// 阶段判断
bool in_high_rms = (rms_error > sc_rms_cg_switch_high_);
bool in_mid_rms  = (rms_error > sc_rms_cg_switch_mid_);

if (in_high_rms) {
    // 纯最速下降
    search = delta_spin;  // 已在第 204 行设置
    
    // 步长限制
    this->check_restriction(search, alpha_trial);
    
    // 直接推进，跳过线性插值
    dnu_last_step = dnu;
    add_scalar_multiply_2d(dnu, search, alpha_trial, dnu);
    delta_lambda = dnu;
    where_fill_scalar_else_2d(constrain_, 0, zero, delta_lambda, delta_lambda);
    add_scalar_multiply_2d(initial_lambda, delta_lambda, one, lambda_);
    
    cal_mw_from_lambda();
    spin_plus = Mi_;
    
    // 简单 alpha 自适应：如果 spin_plus 更接近目标，增大 alpha
    // 否则减小 alpha
    ...
}
else if (in_mid_rms) {
    // CG 搜索方向
    if (i_step >= 2) {
        beta = mean_error / mean_error_old;
        add_scalar_multiply_2d(search, search_old, beta, search);
    }
    this->check_restriction(search, alpha_trial);
    
    // 推进但不做线性插值
    dnu_last_step = dnu;
    add_scalar_multiply_2d(dnu, search, alpha_trial, dnu);
    ...
}
else {
    // 当前完整算法（保持不变）
    ...
}
```

**新增参数**：

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `sc_rms_cg_switch_high` | double | 3.0 | 切换到最速下降的 RMS 阈值 |
| `sc_rms_cg_switch_mid` | double | 1.0 | 切换到完整 CG 的 RMS 阈值 |

**优点**：
- 不增加额外计算量（最速下降阶段甚至减少计算）
- 改动集中在一个函数内部
- 对已收敛的算例无影响

**限制**：
- 最速下降区仍受限于固定电荷密度的子空间近似
- 需要合理设置阈值参数

---

### 方案 5：SCF-Lambda 交替策略

**目标**：在 SCF 未完全收敛时就开始 lambda 优化，避免等到电荷密度"冻结"后磁态难以翻转。

**修改文件**：
- `source/source_pw/module_pwdft/deltaspin_pw.cpp`
- `source/source_lcao/module_deltaspin/lambda_loop.cpp`（接口扩展）

**核心思想**：

```
当前流程：
  SCF 循环 -> drho < sc_scf_thr (1e-3) -> 触发 lambda 循环 -> 磁态已冻结

新流程：
  SCF 循环 -> drho < sc_scf_thr * sc_alt_factor (1e-2) -> 交替执行
    |-- 轻量 lambda: 子空间更新，不刷新电荷 (pw_solve=false)
    |-- 完整 lambda: 每 N 次调用一次 HSolverPW (pw_solve=true)
    |-- SCF 继续推进，逐步收紧 drho 要求
```

**修改 deltaspin_pw.cpp**：
```cpp
bool run_deltaspin_lambda_loop(const int iter,
                               const double drho,
                               const Input_para& inp)
{
    if (!inp.sc_mag_switch) return false;
    
    auto& sc = spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();
    
    // Case 0: linear_scan (保持不变)
    if (inp.sc_lambda_strategy == "linear_scan") {
        sc.run_lambda_linear_scan(iter);
        return true;
    }
    
    // [NEW] Case 0.5: 交替策略（放宽阈值，轻量更新）
    if (inp.sc_alt_enabled && !sc.mag_converged() && drho > 0 
        && drho < inp.sc_scf_thr * inp.sc_alt_factor) {
        
        // 每 sc_alt_interval 次执行一次完整更新
        static int alt_count = 0;
        bool full_update = (alt_count % inp.sc_alt_interval == 0);
        
        sc.run_lambda_loop(iter - 1, true, full_update);
        alt_count++;
        return true;
    }
    
    // Case 1: 标准 lambda 循环（保持不变）
    if (!sc.mag_converged() && drho > 0 && drho < inp.sc_scf_thr) {
        sc.run_lambda_loop(iter - 1);
        sc.set_mag_converged(true);
        return true;
    }
    
    // Case 2: 已收敛，继续优化（保持不变）
    if (sc.mag_converged()) {
        sc.run_lambda_loop(iter - 1);
        return true;
    }
    
    return false;
}
```

**修改 lambda_loop.cpp**：扩展 `run_lambda_loop` 签名以控制最终更新阶段：
```cpp
void run_lambda_loop(int outer_step, bool rerun = true, bool pw_solve = true);
// 新增 pw_solve 参数：控制最终 update_psi_charge 是否调用 HSolverPW
```

**新增参数**：

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `sc_alt_enabled` | bool | false | 是否启用交替策略 |
| `sc_alt_factor` | double | 10.0 | 放宽进入 lambda 循环的 drho 阈值倍数 |
| `sc_alt_interval` | int | 3 | 每几次轻量循环执行一次完整更新 |

**优点**：
- 在 SCF 早期就开始调整 lambda，磁态更灵活
- 完整更新定期刷新电荷，避免子空间近似失效

**限制**：
- 改动较大，涉及跨文件接口修改
- 完整更新（HSolverPW）计算开销较大
- 需要确保 `alt_count` 在 SCF 重启时正确清零
- 可能在 SCF 早期引入不稳定

---

## 4. 推荐实施顺序

| 优先级 | 方案 | 复杂度 | 计算开销 | 风险 |
|--------|------|--------|----------|------|
| 1 | 方案 4（自适应优化器） | 低 | 无增加（甚至减少） | 低 |
| 2 | 方案 2（初始态检测） | 低 | 1-2 次额外 cal_mw | 低 |
| 3 | 方案 5（SCF-lambda 交替） | 高 | 显著增加 | 中 |

**理由**：
- 方案 4 和 2 都是"局部优化"，改动范围小，易于测试和回滚
- 方案 4 解决的是大误差区 CG 方法不稳定的问题，是通用性改进
- 方案 2 针对的是具体的磁矩翻转场景
- 方案 5 涉及 SCF 流程层面的改动，需要更多验证
- 建议先实施方案 4+2，评估效果后再决定是否引入方案 5

## 5. 参数汇总

所有新增参数统一通过 INPUT 文件读取：

```
sc_flip_detect 1                # 启用翻转检测
sc_guide_lambda_scale 2.0       # 引导 lambda 缩放系数
sc_rms_cg_switch_high 3.0       # 最速下降阈值
sc_rms_cg_switch_mid 1.0        # 完整 CG 阈值
sc_alt_enabled 0                # 禁用交替策略（默认）
sc_alt_factor 10.0              # 交替策略放宽因子
sc_alt_interval 3               # 交替间隔
```
