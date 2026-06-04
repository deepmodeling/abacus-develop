# Verlet::thermalize() OpenMP 并行化

日期：2026-06-03

## 修改文件

- `source/source_md/verlet.cpp`，第 122-127 行

## 修改前代码

```cpp
    for (int i = 0; i < ucell.nat; ++i)
    {
        vel[i] *= fac;
    }
```

## 修改后代码

```cpp
    const int nat = ucell.nat;
#pragma omp parallel for schedule(static) if (nat >= 256)
    for (int i = 0; i < nat; ++i)
    {
        vel[i] *= fac;
    }
```

## 改动说明

- 此循环位于 `Verlet::thermalize()` 中，`fac` 由上层温控策略（rescaling / rescale_v / berendsen）在循环外串行计算完毕。
- 每次迭代只写自己的 `vel[i]`，读取共享标量 `fac`，迭代之间无数据依赖。
- 使用 `schedule(static)` 保证可重复性和缓存局部性。
- 使用 `if (nat >= 256)` 避免小体系时线程启动开销。

## OpenMP 指令

- `#pragma omp parallel for schedule(static) if (nat >= 256)`

## 调用链

- `Verlet::apply_thermostat()` → `thermalize()`（rescaling / rescale_v / berendsen 分支）
- 不涉及 Anderson 分支（该分支使用 `std::rand()` 和 `gaussrand()`，保持串行）。

## 潜在风险

- 无。纯独立写入，不涉及归约或共享状态。
