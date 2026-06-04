# MSST::rescale() OpenMP 并行化

日期：2026-06-03

## 修改文件

- `source/source_md/msst.cpp`，第 265-270 行

## 修改前代码

```cpp
    for (int i = 0; i < ucell.nat; ++i)
    {
        vel[i][sd] *= dilation[sd];
    }
```

## 修改后代码

```cpp
    const int nat = ucell.nat;
#pragma omp parallel for schedule(static) if (nat >= 256)
    for (int i = 0; i < nat; ++i)
    {
        vel[i][sd] *= dilation[sd];
    }
```

## 改动说明

- 此循环位于 `MSST::rescale()` 末尾，在此之前已完成晶胞矩阵的 dilation 更新和 `unitcell::setup_cell_after_vc()` 调用（均保持串行）。
- 每次迭代只写 `vel[i][sd]`，读取共享标量 `dilation[sd]` 和常量 `sd`。
- `sd` 为 MSST 冲击方向（0/1/2），循环内所有迭代使用相同的 `sd` 和 `dilation[sd]`。
- 使用 `schedule(static)` 保证可重复性和缓存局部性。

## OpenMP 指令

- `#pragma omp parallel for schedule(static) if (nat >= 256)`

## 调用链

- `MSST::second_half()` → `rescale()`

## 潜在风险

- 无。不改变晶胞更新顺序，只并行原子速度的纯缩放写入。
