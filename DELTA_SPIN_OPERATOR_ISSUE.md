# DeltaSpin 哈密顿量处理问题

## 问题：DeltaSpin HR 贡献在多 k 点循环中重复累积

### 背景

LCAO 基组下，每个 SCF 迭代需要遍历所有 k 点计算哈密顿量。`HamiltLCAO::updateHk(ik)` 在每个 k 点调用 `OperatorLCAO::init(ik)`，后者根据 `cal_type` 分发到不同 operator 的 `contributeHR()` 和 `contributeHk()`。

DeltaSpin 的 `contributeHR()` 使用 `+=` 累加模式向共享 `hR` 写入贡献。在多 k 点循环中，当 `hr_done` 在 k 点间传播为 `true` 时，DeltaSpin 的 `contributeHR()` 仍被无条件调用，导致贡献被累积多次。

## 当前修复方案：DeltaSpin 独立 `sc_hr_done` 标志

### 核心思路

DeltaSpin 不再依赖共享的 `OperatorLCAO::hr_done`，而是维护自己的 `sc_hr_done` 标志，独立控制每个 k 点/每个 spin 通道的 HR 计算。

### 修改文件

| 文件 | 修改内容 |
|------|----------|
| `dspin_lcao.h` | 1. 添加 `bool sc_hr_done = false;` 成员<br>2. `update_lambda()` 中重置 `sc_hr_done = false`<br>3. 新增 `set_current_spin()` 函数，spin 切换时重置 `sc_hr_done` |
| `dspin_lcao.cpp` | 1. `contributeHR()` 中用 `sc_hr_done` 替代 `this->hr_done` 控制计算逻辑<br>2. `contributeHR()` 末尾设置 `sc_hr_done = true` |
| `hamilt_lcao.cpp` | `updateHk()` 中添加缺失的 `set_current_spin()` 调用（与 tmp 对齐） |

### 核心逻辑

**dspin_lcao.h:**
```cpp
bool sc_hr_done = false;  // 独立于 OperatorLCAO::hr_done

void update_lambda() {
    for(int is=0; is<spin_num; is++)
        update_lambda_[is] = true;
    sc_hr_done = false;  // lambda 变化时重置
}

void set_current_spin(const int current_spin_in) {
    if (this->current_spin != current_spin_in)
        sc_hr_done = false;  // spin 切换时重置
    OperatorLCAO<TK, TR>::set_current_spin(current_spin_in);
}
```

**dspin_lcao.cpp - contributeHR():**
```cpp
if (!sc_hr_done) {
    // Case 1: 首次或 spin 切换或 lambda 变化，清零 lambda_save
    lambda_save.assign(nat * 3, 0.0);
}
else if (sc_hr_done && !update_lambda_[current_spin]) {
    // Case 2: HR 已计算且 lambda 未变化，跳过
    return;
}
// Case 3: 计算 HR 贡献 +=
// ...
sc_hr_done = true;
```

**hamilt_lcao.cpp - updateHk():**
```cpp
this->current_spin = this->kv->isk[ik];
// 补充 tmp 中有的调用，使 spin 切换能触发 set_current_spin()
dynamic_cast<hamilt::OperatorLCAO<TK, TR>*>(this->ops)->set_current_spin(this->kv->isk[ik]);
```

### 执行流程（nspin=2, lambda 循环）

```
Lambda loop entry:
  refresh_times=0 (无 HR reset)
  update_lambda() → sc_hr_done=false

k=0 (spin-up, isk=0):
  updateHk(0): 无 spin switch (current_spin=0)
  set_current_spin(0): current_spin 不变 → sc_hr_done 仍为 false
  init(0): first node, hr_done=true → hR 不清零 (但 operator_lcao.cpp 中
           仅 first node 且 !hr_done 时清零 hR，这里 hr_done=true 所以不清)
  DeltaSpin contributeHR(): !sc_hr_done → 计算 spin-up ΔHR → sc_hr_done=true

k=1 (spin-down, isk=1):
  updateHk(1): spin switch, hR 指针切换到 spin-down 半区
               refresh_times=0 → hr_done 保持 true（其他 operator 跳过）
  set_current_spin(1): current_spin 0→1 → sc_hr_done 重置为 false ✓
  init(1): hr_done=true → hR 不清零 (spin-down 半区未被清零)
           Veff/DFTU: !hr_done=false → 跳过
           DeltaSpin: !sc_hr_done → 计算 spin-down ΔHR → sc_hr_done=true
```

### 为什么需要 `set_current_spin()` 调用

当前版本 `hamilt_lcao.cpp:526` 仅更新了 `HamiltLCAO::current_spin`，但**没有**调用 operator 链的 `set_current_spin()`。tmp 版本有此行：

```cpp
// tmp: hamilt_lcao.cpp:437
dynamic_cast<hamilt::OperatorLCAO<TK, TR>*>(this->ops)->set_current_spin(this->kv->isk[ik]);
```

没有这行调用时，DeltaSpin 的 `current_spin` 始终为 0，lambda 循环中 spin switch 无法触发 `sc_hr_done` 重置，导致 spin-down 的 HR 计算被错误跳过。

### 当前状态

- ✅ 编译通过
- ❌ LCAO DeltaSpin 测试 lambda 内循环**仍无法收敛**
- 需进一步排查
