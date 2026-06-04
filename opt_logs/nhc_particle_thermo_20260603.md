# Nose_Hoover::particle_thermo() OpenMP 并行化

日期：2026-06-03

## 修改文件

- `source/source_md/nhchain.cpp`，第 556-561 行

## 修改前代码

```cpp
    for (int i = 0; i < ucell.nat; ++i)
    {
        vel[i] *= scale;
    }
```

## 修改后代码

```cpp
    const int nat = ucell.nat;
#pragma omp parallel for schedule(static) if (nat >= 256)
    for (int i = 0; i < nat; ++i)
    {
        vel[i] *= scale;
    }
```

## 改动说明

- 此循环位于 `Nose_Hoover::particle_thermo()` 末尾。`scale` 已由前段的 thermostat chain 积分循环（NHC 链变量 `eta`/`v_eta`/`g_eta` 递推）串行计算完毕。
- thermostat chain 自身的积分具有严格递推依赖（`m` 从 `md_tchain-1` 向 0 迭代 + 多步 Yoshida-Suzuki 分解），保持完全串行。
- 此循环仅将最终 `scale` 应用到每个原子速度上，每原子独立。
- 使用 `schedule(static)` + `if (nat >= 256)` 阈值。

## OpenMP 指令

- `#pragma omp parallel for schedule(static) if (nat >= 256)`

## 调用链

- `Nose_Hoover::first_half()` / `second_half()` → `particle_thermo()`

## 潜在风险

- NHC/NPT 物理量对数值误差敏感。并行不改变 `scale` 的计算过程，只改变乘法应用顺序（独立写入，无归约差异）。
- 多线程一致性需通过 `MODULE_MD_nhc` 单测验证。
