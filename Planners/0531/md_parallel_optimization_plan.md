# ABACUS MD 模块并行优化方案

日期：2026-05-31
范围：以 `source/source_md` 为核心，覆盖其直接关联的 `source/source_esolver`、`source/source_cell`、构建与测试代码。
原则：优先选择实现容易、收益明确、正确性风险低的并行点；OpenMP 优先，MPI 只作为兼容与后续扩展方向；本阶段不安排 CUDA 编程。

---

## 一、相对既有项目计划的调整

本方案以 `Planners/06_md.md` 的总任务为目标约束，以 `Planners/Project_plan.md` 为基础，但根据当前代码实际做如下调整：

| 调整项 | `Project_plan.md` 中的表述 | 本方案调整 | 原因 |
| --- | --- | --- | --- |
| CUDA | 作为有余力时的扩展尝试 | 本阶段明确暂缓，不列入实现主线 | 用户要求暂时不考虑 CUDA；且当前易加速点主要是 CPU 侧数据转换和归约 |
| 势函数专项 | DPMD 或 NEP 二选一深入 | 建议优先 NEP 接口层，DPMD 作为同类同步优化对象 | NEP 当前有能量 `accumulate`、力回填、virial 求和三类典型循环，OpenMP 收益更直接；DPMD 主要是坐标与力回填，virial 仅 3x3 |
| Langevin/Andersen | 原计划把 Langevin 作为候选阅读对象 | 本方案暂不并行随机数生成循环，只保留阻尼力/后处理改造的预研 | `std::rand()` 和 `MD_func::gaussrand()` 依赖全局状态与调用顺序，直接并行会改变随机序列并引入数据竞争 |
| MPI | 课程说明中包含 MPI 层级 | 本方案不做 MD 域分解，只保证 OpenMP 与现有 MPI 广播兼容 | 当前 MD 积分器在 `my_rank == 0` 上更新再 `MPI_Bcast`，引入 MPI 原子分解会改变数据所有权，工作量和风险明显高于本阶段收益 |
| 测试 | 建立测试与 benchmark 框架 | 先复用已有 `MODULE_MD_*` 单测并增加线程一致性脚本，再做完整 benchmark | 当前仓库已有 MD 单测入口，先复用成本最低 |
| 新增候选 | 原计划未强调 LJ 势 | 将 LJ 作为可选 P2 优化和 benchmark 参考 | LJ 不依赖 DPMD/NEP 外部库，适合在没有模型文件时评估 MD 主循环和 OpenMP 正确性 |

---

## 二、当前代码事实与并行边界

MD 主循环在 `source/source_md/run_md.cpp` 中完成。`Run_MD::md_line()` 根据 `md_type` 创建 `FIRE`、`Nose_Hoover`、`Verlet`、`Langevin`、`MSST`，每步依次执行 `first_half()`、`MD_func::force_virial()`、`second_half()`、`compute_stress()`、`current_temp()`、输出与 restart。关键依赖顺序在 `run_md.cpp:57-136`，不建议改变。

积分器公共状态在 `MD_base`。`MD_base::update_pos()` 的位置增量循环位于 `md_base.cpp:95-120`，`MD_base::update_vel()` 的速度半步循环位于 `md_base.cpp:126-144`。两者当前只在 `my_rank == 0` 上串行更新，然后广播给所有 MPI rank。OpenMP 应只包住 rank 0 内部循环，广播和 `unitcell::update_pos_taud()` 保持原顺序。

统计与力转换集中在 `MD_func`。`kinetic_energy()` 位于 `md_func.cpp:43-52`，`temp_vector()` 位于 `md_func.cpp:499-518`，`force_virial()` 中 `force_temp` 到 `Vector3` 的回填位于 `md_func.cpp:275-313`。这些都是 OpenMP 的低风险入口。

DPMD 接口在 `source/source_esolver/esolver_dp.cpp`。每步包括 cell 写入、坐标转换、外部 `dp.compute()`、能量/力/virial 单位转换。坐标转换在 `esolver_dp.cpp:76-88`，力回填在 `esolver_dp.cpp:108-113`。

NEP 接口在 `source/source_esolver/esolver_nep.cpp`。坐标转换在 `esolver_nep.cpp:73-87`，能量求和在 `esolver_nep.cpp:102-104`，力回填在 `esolver_nep.cpp:107-113`，virial 求和在 `esolver_nep.cpp:115-133`。这是机器学习势函数接口层最值得优先处理的区域。

随机数相关循环需要保守处理。`Verlet::apply_thermostat()` 的 Anderson 分支在 `verlet.cpp:74-97` 使用 `std::rand()` 和 `MD_func::gaussrand()`；`Langevin::post_force()` 在 `langevin.cpp:86-107` 中每原子调用 `std::rand()`；`MD_func::rand_vel()` 在 `md_func.cpp:157-211` 中依赖串行高斯随机数。这些循环本阶段不直接 OpenMP 化。

---

## 三、优化目标

1. 在不改变 MD 主循环和物理算法的前提下，提高大原子数体系中积分器、统计归约、势函数接口层数据转换的 CPU 利用率。
2. 保持 OpenMP 关闭时行为完全兼容；OpenMP 开启时与现有 MPI 广播逻辑兼容。
3. 对确定性计算路径做到多线程结果在数值容差内一致；对随机数路径不改变随机序列。
4. 优先完成能被已有单元测试覆盖的改动，再推进依赖外部 DPMD/NEP 库的测试。

---

## 四、优先级总览

| 优先级 | 优化点 | 文件/函数 | OpenMP 策略 | 收益判断 | 风险 |
| --- | --- | --- | --- | --- | --- |
| P0 | 基础积分器位置更新 | `MD_base::update_pos()` | `parallel for schedule(static)`，仅 rank 0 | 每步必经，所有 MD 类型共享 | 低 |
| P0 | 基础积分器速度更新 | `MD_base::update_vel()` | `parallel for schedule(static)`，仅 rank 0 | 每步至少两次，所有 MD 类型共享 | 低 |
| P0 | 动能与温度归约 | `MD_func::kinetic_energy()` | `reduction(+:ke)` | `current_temp()` 多处调用 | 低，浮点归约有微小顺序差 |
| P0 | 温度张量归约 | `MD_func::temp_vector()` | 9 个局部标量 reduction 或线程局部 3x3 后合并 | NPT/MSST stress 计算依赖 | 中低 |
| P0 | 力矩阵回填 | `MD_func::force_virial()` | `parallel for` 回填 `force[i][j]` | 每步势函数后必经 | 低 |
| P1 | Verlet 温控缩放 | `Verlet::thermalize()` | `parallel for` 缩放 `vel[i]` | NVT rescaling/berendsen 低成本提速 | 低 |
| P1 | FIRE 归约和速度更新 | `FIRE::check_force()`、`FIRE::check_fire()` | `max`、`+` reduction + `parallel for` | 弛豫任务常用；循环清晰 | 中低 |
| P1 | NEP 接口层 | `ESolver_NEP::runner()` | 坐标/力回填并行，能量/virial reduction | 机器学习势主线，收益明确 | 中 |
| P1 | DPMD 接口层 | `ESolver_DP::runner()` | 坐标/力回填并行 | 与 NEP 同类，可同步完成 | 中低 |
| P2 | NHC 原子速度缩放 | `particle_thermo()`、`vel_baro()` | 串行计算 scale/factor，原子循环并行 | NVT/NPT 常用 | 中 |
| P2 | MSST 原子循环 | `vel_sum()`、`rescale()`、`propagate_vel()` | `reduction` + `parallel for` | 特定场景提速 | 中 |
| P2 | LJ 势参考优化 | `ESolver_LJ::runner()` | 原子邻居循环并行，能量/virial reduction | 无外部库 benchmark | 中高 |
| 暂缓 | 随机数循环 | `rand_vel()`、Anderson、Langevin | 不直接并行 | 避免改变随机序列和线程安全问题 | 高 |

---

## 五、详细优化方案

### 5.1 OpenMP 接入与通用写法

仓库顶层 `CMakeLists.txt` 已有 `USE_OPENMP` 选项，默认开启，并在 `CMakeLists.txt:398-400` 查找和链接 `OpenMP::OpenMP_CXX`。测试框架 `cmake/Testing.cmake` 也会给单测目标链接 OpenMP。因此不需要新增全局构建选项，只需在代码中使用标准 pragma。

建议写法：

```cpp
#pragma omp parallel for schedule(static) if (natom >= 256)
for (int i = 0; i < natom; ++i)
{
    ...
}
```

实现要求：

1. 不新增必须依赖 OpenMP API 的逻辑；只使用 pragma 时，未开启 OpenMP 也可编译为串行。
2. 对短循环增加 `if (natom >= 256)` 或同类阈值，避免小体系线程启动开销。
3. 默认使用 `schedule(static)`，保持可重复性和缓存局部性。
4. 所有 reduction 接受浮点求和顺序差异，测试使用绝对/相对容差，不要求 bitwise 相同。
5. 不在已经可能被外部库 OpenMP 并行的区域包大范围 `parallel`，避免嵌套并行过度订阅。

### 5.2 P0：`MD_base::update_pos()` 并行化

位置：`source/source_md/md_base.cpp:95-120`

当前逻辑：

1. rank 0 遍历所有原子；
2. 根据 `ionmbl[i][k]` 计算速度导致的位置增量；
3. 用 `ucell.GT` 转为直接坐标增量；
4. MPI 广播 `pos`；
5. 调用 `unitcell::update_pos_taud()` 更新 `UnitCell` 中的原子位置。

并行方案：

```cpp
if (my_rank == 0)
{
    const int nat = ucell.nat;
    #pragma omp parallel for schedule(static) if (nat >= 256)
    for (int i = 0; i < nat; ++i)
    {
        ...
    }
}
```

正确性依据：

- 每次迭代只写 `pos[i]`，读取 `vel[i]`、`ionmbl[i]`、`md_dt`、`ucell.lat0`、`ucell.GT`。
- 不改变 MPI 广播位置；所有 rank 仍收到 rank 0 的完整 `pos`。
- 不并行 `unitcell::update_pos_taud()`，因为它修改 `ucell.atoms`，应先保持原串行语义。

预期收益：所有 MD 类型每步至少调用一次，原子数大时可获得稳定但有限的加速；对机器学习势主导的总耗时贡献不一定最大，但风险最低。

### 5.3 P0：`MD_base::update_vel()` 并行化

位置：`source/source_md/md_base.cpp:126-144`

并行方案与 `update_pos()` 类似：

- rank 0 内部对 `i` 并行；
- 每个线程只更新 `vel[i]`；
- `force[i]`、`allmass[i]`、`ionmbl[i]` 只读；
- 保持 `MPI_Bcast(vel, ...)` 在并行区之后。

注意点：

- `update_vel()` 每个常规时间步通常调用两次，是最值得先做的共享积分器优化。
- 固定自由度 `ionmbl[i][k] == 0` 保持不更新。
- 不要把 `MPI_Bcast` 放入 OpenMP 并行区。

### 5.4 P0：MD 统计归约并行化

#### 5.4.1 `kinetic_energy()`

位置：`source/source_md/md_func.cpp:43-52`

并行方案：

```cpp
double ke = 0.0;
#pragma omp parallel for reduction(+:ke) schedule(static) if (natom >= 256)
for (int ion = 0; ion < natom; ++ion)
{
    ke += 0.5 * allmass[ion] * vel[ion].norm2();
}
```

调用链：

- `calc_kinetic_state()`；
- `current_temp()`；
- `Verlet::apply_thermostat()`；
- `Nose_Hoover::first_half()`、`second_half()`；
- `MSST::second_half()`；
- `MD_base::print_md()`。

#### 5.4.2 `temp_vector()`

位置：`source/source_md/md_func.cpp:499-518`

当前 `t_vector(i,j) += ...` 会对同一个 3x3 矩阵累加，不能直接 `parallel for` 写共享矩阵。

建议方案：使用 9 个局部 reduction 标量，最后写回矩阵。

```cpp
double t00 = 0.0, t01 = 0.0, ... , t22 = 0.0;
#pragma omp parallel for reduction(+:t00,t01,t02,t10,t11,t12,t20,t21,t22) schedule(static) if (natom >= 256)
for (int ion = 0; ion < natom; ++ion)
{
    const double m = allmass[ion];
    ...
}
t_vector(0,0) = t00;
...
```

正确性注意：

- `t_vector.create(3, 3)` 后不要依赖其是否清零，最终应显式赋值 9 个元素。
- 这个改动会影响 `calc_stress_state()` 和 NPT/MSST 压力相关路径，测试应覆盖 `MODULE_MD_nhc`、`MODULE_MD_msst`。

### 5.5 P0：`force_virial()` 力回填并行化

位置：`source/source_md/md_func.cpp:275-313`

`p_esolver->cal_force()` 返回 `ModuleBase::matrix force_temp` 后，当前串行拷贝到 `ModuleBase::Vector3<double>* force`。可并行化 `md_func.cpp:303-309` 的回填循环：

```cpp
const int nat = unit_in.nat;
#pragma omp parallel for schedule(static) if (nat >= 256)
for (int i = 0; i < nat; ++i)
{
    force[i][0] = force_temp(i, 0);
    force[i][1] = force_temp(i, 1);
    force[i][2] = force_temp(i, 2);
}
```

保持不变的部分：

- `p_esolver->runner()`、`cal_energy()`、`cal_force()`、`cal_stress()` 的调用顺序；
- `potential *= 0.5`、`force_temp *= 0.5`、`virial *= 0.5` 的单位转换顺序。

可选后续：若 `ModuleBase::matrix::operator*=(double)` 内部仍是串行，可单独评估是否替换为并行缩放。但这属于共享基础类型优化，先不扩散范围。

### 5.6 P1：Verlet 温控缩放并行化

位置：`source/source_md/verlet.cpp:110-126`

`Verlet::thermalize()` 在 rescaling、rescale_v、berendsen 温控中对所有原子速度乘同一 `fac`。这是独立写 `vel[i]` 的典型循环，可直接 `parallel for`。

不并行的部分：

- `Verlet::apply_thermostat()` 的 Anderson 分支，位置 `verlet.cpp:74-97`，因为它调用 `std::rand()` 和 `MD_func::gaussrand()`。

### 5.7 P1：FIRE 归约与更新并行化

位置：

- `source/source_md/fire.cpp:156-176`，`check_force()`
- `source/source_md/fire.cpp:180-236`，`check_fire()`

建议拆分：

1. `check_force()` 使用 `reduction(max:max_force)` 获取最大力，再赋给成员 `max`。
2. `check_fire()` 第一段对 `P`、`sumforce`、`normvel` 使用 `reduction(+:P,sumforce,normvel)`。
3. `sumforce` 和 `normvel` 开方后，第二段速度更新使用 `parallel for`。
4. `P <= 0` 分支中的速度清零也使用 `parallel for`。

注意点：

- `md_dt`、`alpha`、`negative_count` 的更新必须保持串行，在归约之后执行。
- 若 `sumforce == 0`，原代码已有潜在除零风险，本方案不顺手改算法，只建议在实现时加测试覆盖或保留原行为。

### 5.8 P1：NEP 接口层 OpenMP 优化

位置：`source/source_esolver/esolver_nep.cpp:56-137`

NEP 是本方案建议优先做的机器学习势函数方向。优化点包括坐标转换、能量求和、力回填和 virial 求和。

#### 5.8.1 坐标转换

当前坐标转换使用 `iat++` 串行编号，不能直接 `parallel for`。建议在 `before_all_runners()` 中构建一次扁平索引：

- `atom_type_of_iat[iat] = it`
- `atom_local_index_of_iat[iat] = ia`

之后每步：

```cpp
const int nat = ucell.nat;
#pragma omp parallel for schedule(static) if (nat >= 256)
for (int iat = 0; iat < nat; ++iat)
{
    const int it = atom_type_of_iat[iat];
    const int ia = atom_local_index_of_iat[iat];
    nep_coord[iat] = ucell.atoms[it].tau[ia].x * ucell.lat0_angstrom;
    nep_coord[iat + nat] = ucell.atoms[it].tau[ia].y * ucell.lat0_angstrom;
    nep_coord[iat + 2 * nat] = ucell.atoms[it].tau[ia].z * ucell.lat0_angstrom;
}
```

这个小索引也可用于 DPMD，避免多个接口重复手写 `iat++` 并行变体。

#### 5.8.2 能量求和

当前 `std::accumulate(_e.begin(), _e.end(), 0.0)` 位于 `esolver_nep.cpp:102-104`。改为 OpenMP reduction：

```cpp
double e_sum = 0.0;
#pragma omp parallel for reduction(+:e_sum) schedule(static) if (nat >= 256)
for (int i = 0; i < nat; ++i)
{
    e_sum += _e[i];
}
nep_potential = fact_e * e_sum;
```

注意：多线程浮点归约会与串行 `std::accumulate` 有末位差异，测试使用容差。

#### 5.8.3 力回填

位置：`esolver_nep.cpp:107-113`

每个原子只写 `nep_force(i, 0..2)`，可直接 `parallel for`。

#### 5.8.4 Virial 求和

位置：`esolver_nep.cpp:115-133`

当前外层 9 个分量、内层 nat 累加。建议写成 9 个分量各自 reduction，或者对 `j` 串行、`i` 并行 reduction：

```cpp
for (int j = 0; j < 9; ++j)
{
    double sum = 0.0;
    #pragma omp parallel for reduction(+:sum) schedule(static) if (nat >= 256)
    for (int i = 0; i < nat; ++i)
    {
        sum += _v[j * nat + i];
    }
    nep_virial_sum[j] = sum;
}
```

如果担心 9 次 parallel region 开销，可使用线程局部 `double local[9]` 后合并，但实现复杂度略高。第一版建议用简单写法，之后根据 profile 决定是否合并并行区。

### 5.9 P1：DPMD 接口层 OpenMP 优化

位置：`source/source_esolver/esolver_dp.cpp:61-125`

可复用 NEP 的扁平原子索引思路。

优化点：

1. `esolver_dp.cpp:76-88` 坐标转换并行化；
2. `esolver_dp.cpp:108-113` `dp_model_force` 到 `dp_force` 的单位转换和回填并行化；
3. `dp_virial` 只有 3x3，不需要 OpenMP；
4. `type_map()` 位于 `esolver_dp.cpp:162-208`，只在初始化执行，不作为并行重点。

注意：

- 不进入 DeePMD-kit 的 `dp.compute()` 内部；外部库可能已有自己的线程策略。
- 如果 DeePMD-kit 使用 OpenMP，需要 benchmark `OMP_NUM_THREADS` 对 ABACUS 外层和库内层的影响，避免线程过度订阅。

### 5.10 P2：Nose-Hoover 原子循环并行化

位置：

- `Nose_Hoover::particle_thermo()` 的速度缩放，`nhchain.cpp:555-559`
- `Nose_Hoover::vel_baro()` 的速度缩放与剪切修正，`nhchain.cpp:690-720`
- `Nose_Hoover::update_baro()` 间接依赖 `MD_func::temp_vector()`，`nhchain.cpp:643-688`

方案：

1. thermostat chain 自身的 `eta`、`v_eta`、`g_eta` 积分循环保持串行，因为链变量有递推依赖。
2. 计算出最终 `scale` 后，对所有 `vel[i] *= scale` 使用 `parallel for`。
3. `vel_baro()` 中 `factor[3]` 先串行计算，然后对 `i` 并行；每个原子只写自己的速度。
4. 不改 `update_volume()`，因为它修改晶格元素和调用 `unitcell::setup_cell_after_vc()`，顺序依赖强。

风险控制：

- NHC/NPT 物理量对数值误差敏感，必须运行 `MODULE_MD_nhc` 多线程一致性测试。
- 先合入 P0/P1 后再做本项，便于定位误差来源。

### 5.11 P2：MSST 原子循环并行化

位置：

- `MSST::vel_sum()`，`msst.cpp:239-249`
- `MSST::rescale()` 中速度缩放，`msst.cpp:251-269`
- `MSST::propagate_vel()`，`msst.cpp:272-309`

方案：

1. `vel_sum()` 使用 `reduction(+:vsum)`。
2. `rescale()` 中晶胞修改和 `setup_cell_after_vc()` 保持串行，只并行其后的 `vel[i][sd] *= dilation[sd]`。
3. `propagate_vel()` 在 rank 0 内对原子并行，`const_C`、`const_D`、`expd` 都是线程局部变量，结束后保留原 MPI 广播。

风险：

- MSST 依赖 `vsum`、体积、应力和速度传播顺序；只改每原子独立部分。
- 先跑 `MODULE_MD_msst`，再跑短 MD 案例比较 `MD_dump`。

### 5.12 P2：LJ 势作为参考 benchmark 的可选并行化

位置：`source/source_esolver/esolver_lj.cpp:58-159`

LJ 不属于机器学习势，但适合作为不依赖外部模型的 benchmark。当前 `runner()` 中按原子和邻居遍历，累加 `lj_potential` 和 `lj_virial`，并写 `lj_force(index, *)`。

可行方案：

1. 建立扁平原子索引，避免 `index++` 串行依赖。
2. 每个原子 `index` 的 `lj_force(index, *)` 只由该原子邻居循环写入，可并行。
3. `lj_potential` 用 reduction。
4. `lj_virial` 使用 9 个局部 reduction 标量，替代共享 `LJ_virial()` 直接累加。

暂不作为主线原因：

- `NeighborSearch` 构建本身可能是主要耗时之一；
- 力是否依赖完整双计数邻居表需要测试确认；
- 该项与机器学习势函数目标相关性弱于 NEP/DPMD 接口层。

---

## 六、暂缓或不建议并行的点

1. `MD_func::gaussrand()`：使用静态变量 `v1`、`v2`、`S`、`phase`，不是线程安全函数。
2. `MD_func::rand_vel()`：初始化速度需要可重复随机序列、总动量归零和温度重标定，直接并行会改变结果。
3. `Verlet::apply_thermostat()` Anderson 分支：随机碰撞选择和高斯随机数均依赖串行 RNG。
4. `Langevin::post_force()`：每原子随机力调用 `std::rand()`，直接并行有数据竞争和序列改变。
5. `Run_MD::md_line()` 时间步循环：MD 时间推进具有严格前后依赖，不能并行时间步。
6. `unitcell::setup_cell_after_vc()`、`unitcell::update_pos_taud()`：涉及 `UnitCell` 结构更新，先保持串行。
7. DPMD/NEP 核心推理：位于外部库，本阶段只优化 ABACUS 接口层。

---

## 七、实现顺序与分支划分

### PR 1：MD 公共循环与统计归约

建议分支：`opt/md-openmp-base`

修改文件：

- `source/source_md/md_base.cpp`
- `source/source_md/md_func.cpp`
- `source/source_md/verlet.cpp`
- `source/source_md/fire.cpp`

内容：

1. 并行 `update_pos()`、`update_vel()`；
2. 并行 `kinetic_energy()`、`temp_vector()`；
3. 并行 `force_virial()` 的力回填；
4. 并行 `Verlet::thermalize()`；
5. 并行 FIRE 的 max/sum reduction 和速度更新。

验证：

```bash
cmake --build build -j
OMP_NUM_THREADS=1 ctest --test-dir build -R 'MODULE_MD_(func|verlet|fire)' --output-on-failure
OMP_NUM_THREADS=4 ctest --test-dir build -R 'MODULE_MD_(func|verlet|fire)' --output-on-failure
```

### PR 2：NEP/DPMD 接口层 OpenMP

建议分支：`opt/md-ml-interface-openmp`

修改文件：

- `source/source_esolver/esolver_nep.h`
- `source/source_esolver/esolver_nep.cpp`
- `source/source_esolver/esolver_dp.h`
- `source/source_esolver/esolver_dp.cpp`

内容：

1. 增加扁平原子索引缓存；
2. 并行 NEP/DPMD 坐标转换；
3. 并行 NEP 能量和 virial 求和；
4. 并行 NEP/DPMD 力回填。

验证：

- 编译不带 `__NEP`/`__DPMD` 时仍通过；
- 有模型文件时运行 1、2、4、8 线程短 MD；
- 比较总能、最大力差、virial/stress 差。

### PR 3：NHC/MSST 与可选 LJ benchmark

建议分支：`opt/md-openmp-extended`

修改文件：

- `source/source_md/nhchain.cpp`
- `source/source_md/msst.cpp`
- 可选：`source/source_esolver/esolver_lj.cpp`

内容：

1. NHC velocity scale 和 barostat velocity update；
2. MSST velocity sum、rescale、propagate_vel；
3. 可选 LJ 原子邻居循环并行化。

验证：

```bash
OMP_NUM_THREADS=1 ctest --test-dir build -R 'MODULE_MD_(nhc|msst|LJ_pot)' --output-on-failure
OMP_NUM_THREADS=4 ctest --test-dir build -R 'MODULE_MD_(nhc|msst|LJ_pot)' --output-on-failure
```

---

## 八、测试与 benchmark 方案

### 8.1 单元测试

优先使用现有 MD 单元测试：

- `MODULE_MD_func`
- `MODULE_MD_verlet`
- `MODULE_MD_fire`
- `MODULE_MD_nhc`
- `MODULE_MD_msst`
- `MODULE_MD_lgv`
- `MODULE_MD_LJ_pot`

每个 PR 至少运行：

```bash
for t in 1 2 4 8; do
  OMP_NUM_THREADS=$t OMP_PROC_BIND=close OMP_PLACES=cores \
  ctest --test-dir build -R 'MODULE_MD' --output-on-failure
done
```

`MODULE_MD_lgv` 虽然本阶段不并行 Langevin 随机数循环，但仍要跑，确保 P0 公共循环没有破坏 Langevin 流程。

### 8.2 线程一致性检查

建议新增一个轻量脚本，例如 `tools/benchmark/md_openmp_check.sh` 或放在课程报告目录中，自动执行：

1. 使用 `OMP_NUM_THREADS=1` 生成基准输出；
2. 使用 `OMP_NUM_THREADS=2/4/8` 生成对比输出；
3. 提取 `#TOTAL ENERGY#`、`Energy (Ry)`、`Temperature (K)`、最大力；
4. 用容差比较：
   - 能量相对误差 `<= 1e-10`；
   - 力最大绝对误差 `<= 1e-10` 到 `1e-8`，按势函数精度调整；
   - 温度相对误差 `<= 1e-10`；
5. 对随机数相关测试只检查程序稳定和物理量有限，不要求跨线程 bitwise 相同。

### 8.3 性能 benchmark

建议至少三类规模：

| 规模 | 原子数 | 用途 |
| --- | --- | --- |
| small | 100-500 | 检查线程开销，验证阈值不会拖慢小体系 |
| medium | 2,000-10,000 | 主要性能对比 |
| large | 20,000+ | 观察内存带宽和并行效率 |

线程数：

- `OMP_NUM_THREADS=1`
- `2`
- `4`
- `8`
- 机器允许时 `16`

统计指标：

- `MD_base::update_pos`
- `MD_base::update_vel`
- `MD_func::force_virial`
- `ESolver_NEP::runner` 或 `ESolver_DP::runner`
- 完整 MD 每步平均耗时
- 加速比 `S_p = T_1 / T_p`
- 并行效率 `E_p = S_p / p`

ABACUS 已有 `ModuleBase::timer`，第一版可直接利用已有 timer 输出；如果粒度不够，再给 NEP/DPMD 接口层局部循环增加临时 timer，最终提交前移除过细的临时输出。

---

## 九、正确性风险与处理

| 风险 | 影响 | 处理 |
| --- | --- | --- |
| 浮点归约顺序变化 | 能量、温度、stress 末位差异 | 测试用容差；报告说明非 bitwise 差异来源 |
| 小体系线程开销 | 性能变慢 | `if (nat >= 256)` 阈值；benchmark 后调整 |
| 外部库也使用 OpenMP | 线程过度订阅 | benchmark 记录外部库线程设置；必要时文档建议控制库内线程 |
| MPI + OpenMP 混合 | rank 0 线程计算后广播 | 不改变广播位置；运行 MPI 单测 |
| 随机数线程安全 | 数据竞争、随机序列改变 | 本阶段不并行随机数循环 |
| 可变胞更新顺序 | NPT/MSST 结果漂移 | 只并行原子独立循环，不并行晶胞更新函数 |

---

## 十、预期收益

P0 的积分器和统计优化属于“每步必经”的低风险收益。对于 ML 势主导的大体系，总耗时加速可能不如势函数核心明显，但可以减少 ABACUS 接口层和后处理的串行尾部，提升多线程下的整体伸缩。

P1 的 NEP/DPMD 接口层优化更贴合机器学习分子动力学目标。NEP 中能量、力、virial 后处理均与原子数线性相关，且每步调用一次；大体系下预期比单纯积分器循环更有可见收益。

P2 的 NHC/MSST/LJ 优化用于扩展覆盖面和 benchmark 支撑。建议在 P0/P1 稳定后再做，避免一次性改动过多导致数值误差难定位。

---

## 十一、最终交付物建议

1. 代码改动 PR：
   - `opt/md-openmp-base`
   - `opt/md-ml-interface-openmp`
   - `opt/md-openmp-extended`
2. 测试脚本：
   - 多线程单测脚本；
   - 短 MD 输出一致性检查脚本；
   - benchmark 数据收集脚本。
3. 报告数据：
   - 每个优化点的改动说明；
   - 1/2/4/8 线程耗时表；
   - 加速比和并行效率；
   - 正确性误差表；
   - 暂缓随机数和 CUDA 的理由。

本方案的主线应先完成 P0 和 NEP/DPMD P1。若时间有限，至少交付 `MD_base`、`MD_func`、`Verlet/FIRE`、`ESolver_NEP` 四块；这些点覆盖 MD 主循环、统计归约和机器学习势接口层，且最符合“容易实现、提效显著、OpenMP 优先”的要求。
