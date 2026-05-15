# 迭代法求解特征值的并行优化

## 大作业说明

---

## 一、背景介绍

### 0.1 特征值问题基础

#### 0.1.1 什么是特征值问题？

**特征值问题**是线性代数中的核心问题，在科学计算和工程应用中具有广泛的应用。对于一个 $n \times n$ 的矩阵 $A$，特征值 $\lambda$ 和对应的特征向量 $v$ 满足：

$$A v = \lambda v$$"

**在ABACUS中的应用**：
- **电子结构计算**：求解哈密顿量的本征值和本征函数
- **分子动力学**：计算振动频率
- **结构优化**：确定分子和晶体的稳定结构
- **光谱计算**：模拟材料的光学性质

#### 0.1.2 特征值求解方法

**传统方法**：
- **直接法**：如QR算法、特征值分解，计算复杂度 $O(n^3)$
- **迭代法**：如幂法、Lanczos算法、适合大规模稀疏矩阵

**ABACUS中的特征值求解器**：
- **DiagoCG**：基于共轭梯度的求解器
- **DiagoDavidson**：Davidson迭代法

#### 0.1.3 迭代法的优势

**迭代法特别适合**：
- **大规模稀疏矩阵**：如LCAO基组下的哈密顿量
- **只需要部分特征值**：如费米面附近的能级
- **分布式内存环境**：易于并行化
- **内存受限系统**：内存使用与矩阵大小线性相关

**主要迭代方法**：

| 方法 | 适用场景 | 优势 | 计算复杂度 |
|------|---------|------|-----------|
| **幂法** | 求最大特征值 | 简单高效 | $O(n^2)$ per iteration |
| **Davidson** | 大规模稀疏矩阵 | 收敛快 | $O(n^2)$ per iteration |

---

### 1.1 问题由来

在ABACUS的电子结构计算中，特征值求解是计算瓶颈之一。随着体系规模的增大，传统的直接求解方法面临以下挑战：

1. **计算复杂度高**：直接法的 $O(n^3)$ 复杂度限制了可处理的体系大小
2. **内存需求大**：存储完整矩阵和特征向量需要大量内存
3. **并行效率低**：直接法的并行扩展性有限
4. **收敛困难**：金属体系的费米面附近能级密集，传统方法收敛慢

迭代法为解决这些问题提供了有效途径，但现有实现仍有优化空间：

- **并行性能**：MPI和OpenMP并行效率有待提高
- **异构计算**：GPU加速尚未充分利用
- **精度控制**：混合精度计算潜力未发挥
- **算法选择**：缺乏自适应的算法选择机制
- **代码结构**：需要更模块化、可测试的设计

### 1.2 现有代码结构

#### 1.2.1 特征值求解器架构

ABACUS的特征值求解器采用插件式架构：

```
source/source_hsolver/
├── hsolver.h/cpp          # 哈密顿量求解器基类
├── hsolver_lcao.cpp       # LCAO基组求解器
├── hsolver_pw.cpp         # 平面波基组求解器
├── diago_*.cpp            # 各种特征值求解器实现
│   ├── diago_cg.cpp       # 共轭梯度求解器
│   ├── diago_davidson.cpp # Davidson迭代法
│   ├── diago_elpa.cpp     # ELPA求解器
│   └── diago_pexsi.cpp    # PEXSI求解器
└── module_diag/           # 特征值求解相关模块
```

#### 1.2.2 核心接口

```cpp
// source/source_hsolver/hsolver.h
class HSolver
{
public:
    virtual ~HSolver() = default;
    
    // 求解哈密顿量
    virtual void solve(hamilt::Hamilt<T>* phamilt, psi::Psi<T>& psi, double* eigenvalue) = 0;
    
    // 设置求解参数
    virtual void set_parameters(const int& npw, const int& nev) = 0;
};

// 特征值求解器接口
class Diago
{
public:
    virtual ~Diago() = default;
    
    // 对角化求解
    virtual void diag(hamilt::Hamilt<T>* phamilt, psi::Psi<T>& psi, double* eigenvalue) = 0;
    
    // 设置迭代参数
    virtual void set_iterations(int max_iter, double tol) = 0;
};
```

#### 1.2.3 现有迭代法实现

**Davidson迭代法**：
```cpp
// source/source_hsolver/diago_davidson.cpp
void DiagoDavidson<T>::diag(hamilt::Hamilt<T>* phamilt, psi::Psi<T>& psi, double* eigenvalue)
{
    // 初始化 Davidson 子空间
    // 迭代求解
    for (int iter = 0; iter < max_iter; ++iter)
    {
        // 计算残差
        // 扩展子空间
        // 求解小型特征值问题
        // 收敛判断
    }
}
```

**共轭梯度法**：
```cpp
// source/source_hsolver/diago_cg.cpp
void DiagoCG<T>::diag(hamilt::Hamilt<T>* phamilt, psi::Psi<T>& psi, double* eigenvalue)
{
    // 初始化
    // CG 迭代
    for (int iter = 0; iter < max_iter; ++iter)
    {
        // 矩阵-向量乘积
        // 计算残差
        // 更新搜索方向
        // 线搜索
        // 收敛判断
    }
}
```

### 1.3 性能瓶颈分析

#### 1.3.1 计算瓶颈

| 瓶颈 | 位置 | 原因 |
|------|------|------|
| **矩阵-向量乘积** | `hamilt_*.cpp` | 计算量最大，占总时间的60-80% |
| **子空间求解** | `diago_*.cpp` | 小型矩阵对角化，占10-20% |
| **残差计算** | `diago_*.cpp` | 向量操作，占5-10% |
| **收敛判断** | `diago_*.cpp` | 向量范数计算，占1-5% |

#### 1.3.2 并行瓶颈

| 瓶颈 | 原因 | 影响 |
|------|------|------|
| **MPI通信** | 进程间数据传输 | 随着进程数增加，通信开销增大 |
| **内存访问** | 非连续内存访问 | 缓存命中率低，影响计算效率 |
| **负载均衡** | 工作分配不均 | 部分进程空闲，并行效率下降 |
| **同步开销** | 进程间同步 | 等待时间增加，特别是在异构环境 |

---

## 二、建议可以做的事情（共 8 题）

### 题目 1：PPCG 方法实现

**难度**：⭐⭐⭐

#### 题目描述

实现 PPCG（Projected Preconditioned Conjugate Gradient）方法求解特征值问题，这是一种高效的预条件共轭梯度法。

#### 现有代码位置

- `source/source_hsolver/diago_bpcg.h` - BPCG方法实现
- `source/source_hsolver/diago_bpcg.cpp` - BPCG方法实现
- `source/source_hsolver/diago_cg.cpp` - 共轭梯度法实现

#### 具体要求

1. **算法实现**
   - 实现 PPCG 方法，包括预条件器设计
   - 确保算法的数值稳定性
   - 优化收敛策略和预条件器

2. **接口设计**
   - 遵循现有特征值求解器接口
   - 支持不同基组（LCAO和平面波）
   - 提供合理的参数配置

3. **性能测试**
   - 测试不同体系规模的收敛速度
   - 对比与现有方法（如CG、Davidson）的性能
   - 分析计算复杂度和加速比

4. **正确性验证**
   - 与传统方法对比结果
   - 测试不同类型的矩阵
   - 验证收敛性和精度

5. **单元测试要求**
   - 编写单元测试验证 PPCG 算法正确性
   - 测试边界情况和特殊矩阵
   - 验证与现有求解器的结果一致性

6. **代码重构（加分项）**
   - 将 PPCG 方法抽象为可插拔的策略类
   - 实现预条件器的自动选择
   - 设计统一的迭代法接口

### 题目 2：混合精度求解器

**难度**：⭐⭐⭐

#### 题目描述

实现混合精度的特征值求解器，利用单精度计算提高性能，双精度保证精度。

#### 现有代码位置

- `source/source_hsolver/hsolver.h` - 求解器基类
- `source/source_hsolver/diago_*.cpp` - 现有求解器实现

#### 具体要求

1. **精度分析**
   - 分析不同计算步骤的精度需求
   - 确定哪些步骤可以使用单精度
   - 评估混合精度的精度损失

2. **实现方案**
   - 实现float/double混合精度计算
   - 优化精度切换策略
   - 确保最终结果的精度

3. **性能测试**
   - 对比单精度、双精度和混合精度的性能
   - 测试不同体系规模的加速比
   - 分析内存带宽节省

4. **正确性验证**
   - 确保混合精度结果与双精度一致（误差 < 1e-6）
   - 测试不同类型的矩阵
   - 验证收敛性

5. **单元测试要求**
   - 编写单元测试验证混合精度的正确性
   - 测试不同精度组合的效果
   - 验证精度切换的边界情况

6. **代码重构（加分项）**
   - 使用模板实现精度无关的代码
   - 设计精度选择策略
   - 支持运行时精度配置

### 题目 3：MPI并行优化

**难度**：⭐⭐⭐

#### 题目描述

优化特征值求解器的MPI并行实现，提高并行效率和扩展性。

#### 现有代码位置

- `source/source_hsolver/diago_*.cpp` - 特征值求解器
- `source/source_hsolver/module_diag/` - 相关模块

#### 具体要求

1. **并行分析**
   - 分析现有MPI并行实现的瓶颈
   - 识别通信密集型操作
   - 评估负载均衡情况

2. **优化实现**
   - 使用非阻塞通信减少等待
   - 实现计算与通信重叠
   - 优化数据分布和负载均衡

3. **性能测试**
   - 测试不同进程数的加速比
   - 分析并行效率和扩展性
   - 对比优化前后的性能

4. **正确性验证**
   - 确保并行结果与串行一致
   - 测试不同进程数的正确性
   - 验证边界情况

5. **单元测试要求**
   - 编写单元测试验证MPI并行的正确性
   - 测试不同进程数的结果一致性
   - 验证通信错误处理

6. **代码重构（加分项）**
   - 将MPI通信抽象为独立接口
   - 实现通信策略的可配置性
   - 设计自适应的并行策略

### 题目 4：OpenMP多线程加速

**难度**：⭐⭐

#### 题目描述

实现特征值求解器的OpenMP多线程并行，提高共享内存系统的性能。

#### 现有代码位置

- `source/source_hsolver/diago_*.cpp` - 特征值求解器
- `source/source_hsolver/module_diag/` - 相关模块

#### 具体要求

1. **并行化分析**
   - 分析计算密集型操作的并行潜力
   - 识别可并行的循环和操作
   - 评估数据依赖关系

2. **OpenMP实现**
   - 使用`#pragma omp parallel for`实现并行计算
   - 优化线程分配和负载均衡
   - 处理线程私有变量和归约操作

3. **性能测试**
   - 测试不同线程数的加速比
   - 分析并行效率
   - 对比优化前后的性能

4. **正确性验证**
   - 确保并行结果与串行一致
   - 测试不同线程数的正确性
   - 验证线程安全

5. **单元测试要求**
   - 编写单元测试验证OpenMP并行的正确性
   - 测试不同线程数的结果一致性
   - 验证线程同步的正确性

6. **代码重构（加分项）**
   - 将并行计算逻辑抽象为独立模块
   - 实现线程池管理
   - 支持动态线程数调整

### 题目 5：GPU异构加速

**难度**：⭐⭐⭐⭐

#### 题目描述

实现特征值求解器的GPU加速，利用CUDA提高计算性能。

#### 现有代码位置

- `source/source_hsolver/diago_*.cpp` - 特征值求解器
- `source/source_hsolver/module_diag/` - 相关模块

#### 具体要求

1. **GPU加速分析**
   - 分析适合GPU加速的计算部分
   - 评估内存传输开销
   - 设计GPU计算方案

2. **CUDA实现**
   - 实现GPU版本的核心计算
   - 优化内存访问模式
   - 使用CUDA流实现计算与数据传输重叠

3. **性能测试**
   - 对比CPU和GPU版本的性能
   - 测试不同体系规模的加速比
   - 分析内存传输开销

4. **兼容性**
   - 保持与现有代码的接口兼容
   - 支持CPU/GPU自动切换
   - 处理GPU不可用的情况

5. **单元测试要求**
   - 编写单元测试验证GPU计算的正确性
   - 对比CPU和GPU版本的结果一致性
   - 测试不同GPU设备的兼容性

6. **代码重构（加分项）**
   - 将计算设备抽象为独立接口
   - 实现设备选择策略
   - 支持多GPU并行

### 题目 6：代码重构与模块化

**难度**：⭐⭐⭐

#### 题目描述

重构特征值求解器的代码结构，提高模块化程度和可维护性。

#### 现有代码位置

- `source/source_hsolver/` - 求解器相关代码

#### 具体要求

1. **代码分析**
   - 分析现有代码的结构和依赖关系
   - 识别重复代码和设计问题
   - 设计模块化架构

2. **重构实现**
   - 将公共功能提取为独立模块
   - 实现依赖反转和接口抽象
   - 优化代码结构和命名

3. **模块设计**
   - 设计清晰的模块边界
   - 定义明确的接口
   - 减少模块间依赖

4. **测试验证**
   - 确保重构后功能与原代码一致
   - 测试边界情况
   - 验证性能不劣化

5. **单元测试要求**
   - 编写单元测试验证重构后的模块
   - 测试模块间接口的正确性
   - 验证依赖注入的有效性

6. **代码质量**
   - 遵循项目代码规范
   - 添加详细的文档和注释
   - 确保代码可读性

### 题目 7：单元测试框架

**难度**：⭐⭐

#### 题目描述

设计并实现特征值求解器的单元测试框架，确保代码质量和功能正确性。

#### 题目背景

现有特征值求解器缺乏全面的单元测试，这使得代码修改和优化存在风险。建立一个完善的单元测试框架对于保证代码质量至关重要。

#### 具体要求

1. **测试框架设计**
   - 设计适合特征值求解器的单元测试框架
   - 定义测试用例和测试方法
   - 实现测试结果的自动验证

2. **测试用例实现**
   - 编写迭代法求解的测试用例
   - 编写并行计算的测试用例
   - 编写混合精度的测试用例

3. **测试覆盖**
   - 确保关键功能的测试覆盖
   - 测试边界情况和异常处理
   - 验证不同并行配置的正确性

4. **性能测试**
   - 实现性能基准测试
   - 监控优化效果
   - 提供性能分析工具

5. **集成与自动化**
   - 集成到CI/CD流程
   - 实现测试的自动化运行
   - 提供测试报告生成

6. **代码重构（加分项）**
   - 将测试框架抽象为独立的模块
   - 实现测试数据的自动生成
   - 支持测试结果的可视化

### 题目 8：效率提升与算法优化

**难度**：⭐⭐⭐

#### 题目描述

优化特征值求解器的算法和实现，提高计算效率和收敛速度。

#### 现有代码位置

- `source/source_hsolver/diago_*.cpp` - 特征值求解器

#### 具体要求

1. **算法分析**
   - 分析现有迭代法的收敛特性
   - 识别计算瓶颈
   - 评估优化潜力

2. **优化实现**
   - 改进收敛加速策略
   - 优化预条件器
   - 实现自适应算法参数

3. **性能测试**
   - 测试不同优化策略的效果
   - 分析收敛速度和计算时间
   - 对比优化前后的性能

4. **正确性验证**
   - 确保优化后结果与原代码一致
   - 测试不同类型的矩阵
   - 验证收敛性和稳定性

5. **单元测试要求**
   - 编写单元测试验证优化后的算法
   - 测试不同优化策略的正确性
   - 验证边界情况

6. **代码重构（加分项）**
   - 实现算法参数的自动调优
   - 设计自适应的收敛策略
   - 支持多种预条件器

---

## 三、测试环境与基准数据

### 3.1 推荐测试体系

| 体系 | 原子数 | 基组 | 矩阵大小 | 推荐测试规模 |
|------|--------|------|----------|-------------|
| H₂O 分子 | 3 | LCAO | ~100 | 初级测试 |
| Si 晶体 | 64 | LCAO | ~1000 | 基准测试 |
| Al 金属 | 128 | LCAO | ~2000 | 性能测试 |
| TiO₂ | 192 | LCAO | ~3000 | 大规模测试 |

### 3.2 性能基准

| 优化项 | 当前时间 | 目标时间 | 最低加速比 |
|--------|---------|---------|-----------|
| PPCG方法 | T₁ | T₁/2 | 2x |
| 混合精度 | T₂ | T₂/1.5 | 1.5x |
| MPI 并行 | T₃ | T₃/4 | 4x (4进程) |
| OpenMP 并行 | T₄ | T₄/4 | 4x (4线程) |
| GPU 加速 | T₅ | T₅/10 | 10x |
| 算法优化 | T₆ | T₆/2 | 2x |

### 3.3 测试脚本参考

```bash
#!/bin/bash
# benchmark_diago.sh - 特征值求解性能测试

export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8

for nproc in 1 2 4 8 16; do
    for nthread in 1 2 4 8; do
        echo "Testing: nproc=$nproc, nthread=$nthread"
        export OMP_NUM_THREADS=$nthread
        mpirun -np $nproc ./abacus INPUT > log_p${nproc}_t${nthread}.out 2>&1
        grep "eigenvalue calculation" log_p${nproc}_t${nthread}.out | tail -1
    done
done

# GPU测试
if [ -n "$CUDA_VISIBLE_DEVICES" ]; then
    echo "Testing with GPU"
    mpirun -np 1 ./abacus INPUT_gpu > log_gpu.out 2>&1
    grep "eigenvalue calculation" log_gpu.out | tail -1
fi
```

---

## 四、代码规范与提交流程

### 4.1 代码规范

1. **命名规范**
   - 遵循项目现有的命名风格
   - 新增函数需添加文档注释

2. **模块化设计**
   - 独立功能封装为独立函数/类
   - 便于单元测试

3. **错误处理**
   - 检查所有 MPI 调用返回值
   - 妥善处理异常情况

4. **并行代码规范**
   - 明确并行区域和同步点
   - 避免死锁和竞争条件
   - 注释并行策略和通信模式

### 4.2 提交流程

#### 4.2.1 推荐方式：GitHub Pull Request ⭐

为了更好地模拟真实软件开发流程，我们**强烈推荐**使用 GitHub 进行代码提交和协作。具体方式如下：

1. **Fork 仓库**
   - Fork ABACUS deepmodeling仓库到你自己的 GitHub 账户
   - 地址：`https://github.com/deepmodeling/abacus-develop`

2. **创建分支**
   ```bash
   git checkout -b feature/eigen-solver-optimization
   ```

3. **少量多次提交**
   ```bash
   # 每次完成一个小功能就提交
   git add source/source_hsolver/
   git commit -m "Add Jacobi-Davidson solver implementation"
   git push origin feature/eigen-solver-optimization
   ```

4. **提交 Pull Request**
   - 在 GitHub 上创建 Pull Request
   - 描述你做了哪些优化
   - 请求代码 Review

#### 4.2.2 提交策略

| 原则 | 说明 |
|------|------|
| **少量多次** | 每完成一个小功能就提交，不要等到最后一次性提交 |
| **问题导向** | 每个 PR 解决一个具体问题 |
| **文档完善** | PR 描述中说明解决了什么瓶颈、预期性能提升 |
| **可验证** | 提交时附带测试结果或性能数据 |

#### 4.2.3 代码接受标准

**你的代码被官方仓库接受将获得额外加分**：

| 🌟 代码被 merged | PR 被接受并合并到主分支 |
| 🌟 代码可运行 | 通过基本编译和测试 |

#### 4.2.4 评分原则

> **核心原则：以实际解决问题的质量和数量作为评价标准**

- 代码不被接受也可以获得分数，取决于工作量和完成质量
- 重点关注：是否真正解决了实际问题、是否有创新性、代码是否健壮
- 不以"是否被接受"作为唯一标准

---

### 4.3 报告格式要求

```latex
\documentclass[12pt,a4paper]{article}

\title{迭代法求解特征值的并行优化}
\author{姓名}
\date{\today}

\begin{document}
\maketitle

\section{引言}
% 描述问题背景和优化目标

\section{现有代码分析}
% 分析当前实现的瓶颈

\section{优化方案}
% 描述实现的优化方法

\section{性能测试}
% 包含测试结果和图表

\section{结论}
% 总结优化效果和心得

\end{document}
```

---

## 五、参考资料

### 5.1 代码位置索引

| 文件 | 路径 | 说明 |
|------|------|------|
| 求解器基类 | `source/source_hsolver/hsolver.h` | 哈密顿量求解器基类 |
| Davidson求解器 | `source/source_hsolver/diago_davidson.cpp` | Davidson迭代法 |
| CG求解器 | `source/source_hsolver/diago_cg.cpp` | 共轭梯度法 |

### 5.2 推荐阅读

1. **迭代法**：《Iterative Methods for Sparse Linear Systems》- Y. Saad
2. **特征值算法**：《Numerical Linear Algebra》- T. G. Kolda et al.
3. **并行计算**：《Parallel Programming with MPI》- P. S. Pacheco
4. **CUDA编程**：《Professional CUDA C Programming》- J. Cheng et al.
5. **Davidson方法**："Davidson's method for eigenvalue problems" - E. R. Davidson
6. **Jacobi-Davidson方法**："Jacobi-Davidson style QR and QZ algorithms for the reduction of matrix pencils" - G. L. G. Sleijpen et al.

---

## 六、致谢

本大作业题目设计参考了以下资源：

1. ABACUS 软件源代码 (https://github.com/abacusmodeling/abacus-develop)
2. 特征值求解算法相关文献
3. 并行计算最佳实践
4. 高性能科学计算经验

---

**最后更新**：2026-04-21

**版本**：v1.0
