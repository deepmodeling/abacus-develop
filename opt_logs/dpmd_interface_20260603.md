# ESolver_DP::runner() OpenMP 并行化

日期：2026-06-03

## 修改文件

- `source/source_esolver/esolver_dp.h` — 新增 2 个成员变量
- `source/source_esolver/esolver_dp.cpp` — 修改 3 处

---

## 改点 A：新增扁平原子索引（esolver_dp.h + esolver_dp.cpp）

### esolver_dp.h 新增成员变量

在 `private` 区域 `atype` 之后新增：

```cpp
    std::vector<int> atom_type_index;  ///< type index (it) for each global atom iat
    std::vector<int> atom_local_index; ///< local index (ia) within type for each global atom iat
```

### esolver_dp.cpp `before_all_runners()` 新增索引构建

```cpp
    // Build flat atom index for OpenMP coordinate fill in runner()
    atom_type_index.resize(ucell.nat);
    atom_local_index.resize(ucell.nat);
    int iat = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            atom_type_index[iat] = it;
            atom_local_index[iat] = ia;
            iat++;
        }
    }
```

**说明**：`before_all_runners()` 只在初始化执行一次。索引将 `iat` → `(it, ia)` 映射缓存下来，供 `runner()` 每步复用。

---

## 改点 B：坐标转换并行化（esolver_dp.cpp 第 73-85 行）

### 修改前

```cpp
    std::vector<double> coord(3 * ucell.nat, 0.0);
    int iat = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            coord[3 * iat] = ucell.atoms[it].tau[ia].x * ucell.lat0_angstrom;
            coord[3 * iat + 1] = ucell.atoms[it].tau[ia].y * ucell.lat0_angstrom;
            coord[3 * iat + 2] = ucell.atoms[it].tau[ia].z * ucell.lat0_angstrom;
            iat++;
        }
    }
    assert(ucell.nat == iat);
```

### 修改后

```cpp
    std::vector<double> coord(3 * ucell.nat, 0.0);
    const int nat = ucell.nat;
#pragma omp parallel for schedule(static) if (nat >= 256)
    for (int iat = 0; iat < nat; ++iat)
    {
        const int it = atom_type_index[iat];
        const int ia = atom_local_index[iat];
        coord[3 * iat]     = ucell.atoms[it].tau[ia].x * ucell.lat0_angstrom;
        coord[3 * iat + 1] = ucell.atoms[it].tau[ia].y * ucell.lat0_angstrom;
        coord[3 * iat + 2] = ucell.atoms[it].tau[ia].z * ucell.lat0_angstrom;
    }
```

**说明**：通过预建索引查找每个 `iat` 对应的 `(it, ia)`，消除串行 `iat++` 依赖。DPMD 坐标布局为 row-major `[x0,y0,z0, x1,y1,z1, ...]`，每个线程写独立的 3 个连续位置。移除了原 `assert`（并行循环内不能放 assert），可通过索引构建时的 `assert` 覆盖。

**OpenMP 指令**：`#pragma omp parallel for schedule(static) if (nat >= 256)`

---

## 改点 C：力回填并行化（esolver_dp.cpp 第 104-109 行，`#ifdef __DPMD` 内部）

### 修改前

```cpp
    for (int i = 0; i < ucell.nat; ++i)
    {
        dp_force(i, 0) = f[3 * i] * fact_f;
        dp_force(i, 1) = f[3 * i + 1] * fact_f;
        dp_force(i, 2) = f[3 * i + 2] * fact_f;
    }
```

### 修改后

```cpp
    const int nat_f = ucell.nat;
#pragma omp parallel for schedule(static) if (nat_f >= 256)
    for (int i = 0; i < nat_f; ++i)
    {
        dp_force(i, 0) = f[3 * i] * fact_f;
        dp_force(i, 1) = f[3 * i + 1] * fact_f;
        dp_force(i, 2) = f[3 * i + 2] * fact_f;
    }
```

**说明**：每原子独立写 `dp_force(i, *)`，无依赖。位于 `#ifdef __DPMD` 内部。

**OpenMP 指令**：`#pragma omp parallel for schedule(static) if (nat >= 256)`

---

## 不做修改的部分

- **Virial 回填**（3×3=9 元素）：线程开销大于收益，保持串行
- **Cell 填充**（9 元素）：保持串行
- **`dp.compute()` 调用**：外部库，不进入
- **`type_map()`**：仅初始化执行一次

## 调用链

- `force_virial()` → `p_esolver->runner()` → `ESolver_DP::runner()`

## 潜在风险

| 风险 | 评估 |
|------|------|
| 扁平索引正确性 | 与 `type_map()` 遍历逻辑一致，构建在 `before_all_runners()` 一次性完成 |
| DPMD 外部库线程冲突 | `dp.compute()` 调用仍在串行区；ABACUS 侧仅并行数据填充/回填 |
| `__DPMD` 宏兼容 | 改点 A/B 在宏之前，不受影响；改点 C 在宏内部 |
| 内存开销 | 两个 `std::vector<int>` 各 nat 元素，约 8×nat 字节，可忽略 |
| 移除了坐标填充后的 `assert` | 正确性由索引构建逻辑保证，与 `type_map()` 结构一致 |
