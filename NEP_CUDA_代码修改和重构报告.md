# NEP CUDA 后处理优化代码修改和重构报告

## 1. 修改目标

本次大作业实现选择两种机器学习势函数中的 **NEP** 作为第一阶段修改对象。选择 NEP 的原因是：

- 当前 ABACUS 对 NEP 的接入代码集中在 `source/source_esolver/esolver_nep.cpp`，调用链比 DPMD 更短，适合先做一个可控原型。
- NEP 外部库返回的是每原子能量 `_e`、每原子力 `_f`、每原子 virial `_v`，ABACUS 原代码需要再做一次能量求和、力单位换算、virial 汇总，这部分是规则的线性数组处理，适合 GPU 并行。
- 项目文档和前期算法文档已经确认，当前 NEP 核心计算位于外部 NEP 库中。本次修改不重写外部 NEP 势函数核心，而是在 ABACUS 接入层先完成后处理加速和结构重构。

本次修改不是完整解决 NEP 势函数 GPU 化问题，而是完成一个边界清楚、可继续扩展的第一阶段版本：

```text
NEP 外部库 compute()
  -> 返回 _e, _f, _v
  -> ABACUS 后处理
      CPU 路径: 原等价逻辑
      CUDA 路径: 并行能量求和、力转换、virial 汇总
```

## 2. 修改文件

本次新增和修改的文件如下：

```text
source/source_esolver/CMakeLists.txt
source/source_esolver/esolver_nep.h
source/source_esolver/esolver_nep.cpp
source/source_esolver/esolver_nep_postprocess.h
source/source_esolver/esolver_nep_postprocess.cpp
source/source_esolver/esolver_nep_postprocess.cu
```

## 3. 原始代码问题

原始 `ESolver_NEP::runner()` 中同时承担了四类职责：

```text
1. 构造 NEP cell 数组
2. 构造 NEP coord 数组
3. 调用 nep.compute()
4. 对 _e/_f/_v 做后处理并写回 nep_potential/nep_force/nep_virial
```

这种写法的问题是：

- 每个 MD step 都临时创建 `std::vector<double> cell(9)` 和 `std::vector<double> coord(3 * nat)`。
- 后处理逻辑和 NEP 外部库调用混在同一个函数里，后续很难加 CPU/GPU 双路径。
- 能量、力、virial 的处理都是线性数组操作，但原代码只在 CPU 串行执行。
- 计时粒度只有整个 `runner`，不方便分析输入准备、外部库计算、后处理各自耗时。

## 4. 核心重构

### 4.1 持久化输入缓冲区

在 `ESolver_NEP` 中新增成员：

```cpp
std::vector<double> cell;
std::vector<double> coord;
```

并在 `before_all_runners()` 中按体系大小初始化：

```cpp
cell.resize(9);
coord.resize(3 * ucell.nat);
```

这样每个 MD step 只更新数组内容，不再反复构造临时 vector。

### 4.2 拆分输入准备函数

新增：

```cpp
void ESolver_NEP::prepare_input_buffers(const UnitCell& ucell);
```

该函数专门负责把 ABACUS 的 `UnitCell` 转换成 NEP 需要的数据布局：

```text
cell:
  column-major 3x3 matrix

coord:
  [x0, x1, ..., xN-1,
   y0, y1, ..., yN-1,
   z0, z1, ..., zN-1]
```

同时增加 timer：

```cpp
ModuleBase::timer::start("ESolver_NEP", "prepare_input");
ModuleBase::timer::end("ESolver_NEP", "prepare_input");
```

### 4.3 拆分后处理函数

新增：

```cpp
void ESolver_NEP::postprocess_outputs(const UnitCell& ucell);
```

该函数负责单位换算和 CPU/GPU 路径选择：

```text
if compiled with CUDA and INPUT has device gpu:
    postprocess_nep_cuda(...)
else:
    postprocess_nep_cpu(...)
```

同时增加 timer：

```cpp
ModuleBase::timer::start("ESolver_NEP", "postprocess");
ModuleBase::timer::end("ESolver_NEP", "postprocess");
```

## 5. CPU 后处理路径

新增文件：

```text
source/source_esolver/esolver_nep_postprocess.cpp
```

提供函数：

```cpp
void postprocess_nep_cpu(...);
```

CPU 路径保持原始逻辑等价：

```text
nep_potential = sum(_e) * fact_e
nep_force(i, 0) = _f[i] * fact_f
nep_force(i, 1) = _f[i + nat] * fact_f
nep_force(i, 2) = _f[i + 2 * nat] * fact_f
nep_virial(i, j) = sum(_v[(3*i+j)*nat : ...]) * fact_v
```

这样即使不启用 CUDA，代码行为也应与原版保持一致。

## 6. CUDA 后处理路径

新增文件：

```text
source/source_esolver/esolver_nep_postprocess.cu
```

提供函数：

```cpp
void postprocess_nep_cuda(...);
```

CUDA kernel 的并行粒度为“每个线程处理一个或多个原子”。每个原子线程完成：

```text
1. atomicAdd 到总能量
2. 写出该原子的三维力
3. atomicAdd 到 9 个 virial 分量
```

核心映射：

```text
thread i:
  potential += _e[i] * fact_e
  force[3*i + 0] = _f[i] * fact_f
  force[3*i + 1] = _f[i + nat] * fact_f
  force[3*i + 2] = _f[i + 2*nat] * fact_f
  virial[j] += _v[j*nat + i] * fact_v
```

当前 CUDA 版本是一个教学/原型实现，特点是：

- 代码结构简单，容易检查正确性。
- 能体现能量求和、力转换、virial 汇总的 GPU 并行方式。
- 暂时每次后处理都会申请和释放显存，后续可进一步改成持久化 device buffer。

## 7. 构建系统修改

修改：

```text
source/source_esolver/CMakeLists.txt
```

新增 CPU 文件：

```cmake
esolver_nep_postprocess.cpp
```

当 `USE_CUDA` 开启时，额外编译：

```cmake
esolver_nep_postprocess.cu
```

因此默认 CPU 构建不会依赖 CUDA 文件。

## 8. 当前调用链

修改后的 NEP `runner()` 逻辑变为：

```text
ESolver_NEP::runner()
  -> prepare_input_buffers()
  -> nep.compute(atype, cell, coord, _e, _f, _v)
  -> postprocess_outputs()
      -> postprocess_nep_cpu()
      -> or postprocess_nep_cuda()
```

相比原始版本，`runner()` 现在更像调度函数，具体的数据转换和后处理被拆分出去，后续扩展更容易。

## 9. 正确性与验证

已完成的验证：

```text
g++ -std=c++11 -I.../source -c source/source_esolver/esolver_nep_postprocess.cpp
g++ -std=c++11 -I.../source -c source/source_esolver/esolver_nep.cpp
```

这两个 CPU 侧文件均已通过编译语法检查。

未完成完整 CMake 构建，原因是当前环境默认 CMake 版本过低：

```text
项目要求: CMake >= 3.16
当前 /usr/local/bin/cmake: 3.14.5
当前 /usr/bin/cmake: 3.10.2
```

CUDA 文件也未实际编译，原因是当前环境未发现 `nvcc`。

## 10. 性能收益分析

本次修改的收益分两部分：

1. CPU 路径收益
   减少每步临时分配 `cell/coord`，同时将后处理逻辑拆出，便于进一步优化和测试。

2. GPU 路径收益
   当使用 CUDA 构建并设置 `device gpu` 时，NEP 后处理中的能量求和、力写回、virial 汇总会转到 GPU 并行执行。

需要注意的是，当前 NEP 外部库 `nep.compute()` 仍然是主要瓶颈。如果外部 NEP 库本体仍是 CPU 实现，那么本次 GPU 加速只能覆盖 ABACUS 接入层的后处理部分，不能代表完整 NEP 势函数核心已经 GPU 化。

## 11. 后续改进方向

后续可以继续做三类增强：

```text
1. 将 CUDA 后处理中的 cudaMalloc/cudaFree 改成 ESolver_NEP 成员级持久化 device buffer。
2. 为 postprocess_nep_cpu/postprocess_nep_cuda 增加单元测试，对比能量、力、virial 输出。
3. 如果课程要求更深入的 GPU 加速，需要修改或替换外部 NEP 库，让 nep.compute() 内部的邻域、描述符和神经网络推理也运行在 GPU 上。
```

## 12. 小结

本次修改选择 NEP 作为第一阶段目标，完成了 ABACUS 接入层的结构重构和 CUDA 后处理原型。它保留原有 `ESolver` 接口，不改变 MD 主循环，也不改变外部输入方式。

这是一版适合课程大作业继续推进的中间成果：代码改动集中、风险较低、能解释清楚加速边界，并为后续把 NEP 外部库核心迁移到 GPU 留出了接口位置。
