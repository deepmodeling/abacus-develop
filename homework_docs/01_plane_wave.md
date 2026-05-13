# 平面波生成与应用的并行优化

## 大作业说明

---

## 一、背景介绍

### 0.1 平面波基组基础

#### 0.1.1 什么是平面波基组？

在密度泛函理论（DFT）中，平面波基组是一种常用的基函数展开方式。电子波函数可以展开为平面波的线性组合：

$$\psi_k(\mathbf{r}) = \sum_{\mathbf{G}} c_{\mathbf{k,G}} e^{i(\mathbf{k+G})\cdot\mathbf{r}}$$

其中：
- $\mathbf{k}$ 是k点矢量
- $\mathbf{G}$ 是倒格子矢量
- $c_{\mathbf{k,G}}$ 是展开系数

**平面波基组的特点**：
- 形式简单，易于实现快速傅里叶变换（FFT）
- 收敛性好，增加截断能量即可提高精度
- 适合周期性体系的计算
- 计算量与平面波数量成正比

#### 0.1.2 平面波与FFT

平面波方法的核心是利用FFT实现实空间和倒易空间的快速变换：

```
实空间波函数 ψ(r)
       ↓ FFT
倒易空间波函数 ψ(G)
       ↓ 对易哈密顿量操作
       ↓ IFFT
更新后的实空间波函数 ψ(r)
```

在ABACUS中，平面波变换主要由以下模块处理：
- `PW_Basis`: 单k点平面波基组
- `PW_Basis_K`: 多k点平面波基组
- `pw_transform.cpp`: 实空间与倒易空间变换

#### 0.1.3 布里渊区与k点

##### 什么是布里渊区？

布里渊区（Brillouin Zone）是倒易空间中定义的一个重要概念，它是倒格子原胞在动量空间中的对应区域。

**第一布里渊区**定义为：从原点出发，到最近倒格子点的垂直平分面所围成的最小区域。

```
正格子空间                    倒易空间（布里渊区）
    ↑                              ↑
    |    ● ← 原子                  |    ■ ← 倒格子点
    |                              |   / \
----●----→                       ----■----→ (边界)
    |    ● ← 原子                  |    /   \
    |                              |   /     \
```

对于不同晶格类型，第一布里渊区的形状不同：
- **简单立方**: 立方体
- **面心立方(FCC)**: 截角八面体（十四面体）
- **体心立方(BCC)**: 菱形十二面体

##### k点的重要性

在周期性体系的量子力学计算中，我们不能直接处理连续的k点，而只能在布里渊区内选取离散的k点进行计算。这基于**Bloch定理**：

$$\psi_{n\mathbf{k}}(\mathbf{r}) = e^{i\mathbf{k}\cdot\mathbf{r}} u_{n\mathbf{k}}(\mathbf{r})$$

其中：
- $\mathbf{k}$ 是晶格动量（Bloch动量）
- $u_{n\mathbf{k}}(\mathbf{r})$ 具有与晶格相同的周期性
- $n$ 是能带索引

**为什么需要多个k点？**

1. **能带结构**：不同k点对应不同的电子动量态
2. **积分计算**：布里渊区积分需要多k点采样
3. **精度要求**：k点越密，结果越精确

##### k点选取方法

常见的k点选取方法包括：

| 方法 | 说明 | 适用场景 |
|------|------|---------|
| Monkhorst-Pack | 等间距网格采样 | 一般计算 |
| Gamma-centered | Gamma点为中心的网格 | 高对称性体系 |
| Line-path | 沿着高对称性路径 | 能带计算 |

```bash
# ABACUS中k点设置示例（KPT文件）
K_POINTS
0                         # 0表示 Monkhorst-Pack 格式
Gamma                     # 可选：Gamma 或 Monkhorst-Pack
4 4 4 0 0 0               # n1 n2 n3 s1 s2 s3
                          # n1,n2,n3: 各方向k点数
                          # s1,s2,s3: 位移参数
```

或使用Gamma Only：

```bash
K_POINTS
0
Gamma
1 1 1 0 0 0               # 只有一个Gamma点
```

##### Gamma点（k=0）的特殊地位

Gamma点是指布里渊区的中心点，即 $\mathbf{k} = (0, 0, 0)$。

**Gamma点的特点**：
1. **最高对称性**：Gamma点具有完整的点群对称性
2. **时间反演对称**：$T\psi(\mathbf{r}) = \psi^*(\mathbf{r})$
3. **计算效率高**：可利用对称性减少计算量
4. **内存节省**：可以采用gamma_only模式节省50%存储

**Gamma点适用条件**：
- 非金属体系（费米面附近无需精细采样）
- 高对称性大晶胞
- 初步结构优化

**什么时候不能用Gamma点**：
- 金属体系（需要细密k点采样费米面）
- 需要精确能带（沿高对称性路径）
- 自旋轨道耦合（破坏时间反演对称）

##### 自学资源

同学们可以通过以下方式深入理解布里渊区和k点：

1. **教材参考**
   - 《Solid State Physics》- Ashcroft & Mermin，第2章
   - 《Introduction to Solid State Physics》- Kittel，第5章
   - 《Electronic Structure》- Martin，第3章

2. **在线资源**
   - Wikipedia: Brillouin zone
   - AFLOW: http://aflow.org/aflowWiki/
   - Materials Project: https://materialsproject.org/

3. **代码探索**
   - 阅读ABACUS中k点相关代码：`source/source_cell/klist*`
   - 查看KPT文件格式：`source/source_io/module_io/write_kpoints.cpp`

4. **可视化工具**
   - XCrysDen: 可视化布里渊区和k点路径
   - VESTA: 查看晶体结构
   - SeeK-path: http://seekpath.materialsvirtuallab.org/

### 0.2 Gamma Only 详解

#### 0.2.1 什么是 Gamma Only？

当k点只取Gamma点（$\mathbf{k}=0$）时，由于体系哈密顿量的时间反演对称性，波函数满足：

$$\psi(\mathbf{r}) = \psi^*(\mathbf{r})$$

这意味着波函数是**实函数**，对应的傅里叶变换满足：

$$F(-\mathbf{G}) = F^*(\mathbf{G})$$

**核心结论**：对于实函数，只需要存储一半的傅里叶分量即可完整描述原函数。

#### 0.2.2 Gamma Only 的物理意义

在实空间中，实函数的特点是：
- 实函数 $f(\mathbf{r})$ 的FFT满足 $f(-\mathbf{G}) = f^*(\mathbf{G})$
- 当 $f(\mathbf{G})$ 是复数时，其共轭 $f^*(-\mathbf{G})$ 包含了冗余信息
- 因此可以只存储 $\mathbf{G} \geq 0$ 的分量

```
普通平面波（k≠0）:
存储: [G=0, G=1, G=2, ..., G=N] 全部 N 个平面波
每个平面波系数是复数 c[G] = a + bi

Gamma Only (k=0):
存储: [G=0, G=1, G=2, ..., G=N/2] 只存储 N/2 个平面波
每个平面波系数是复数，但由于共轭对称性，
另一半可以从已存储的系数推导
```

#### 0.2.3 Gamma Only 在代码中的实现

ABACUS代码中，`gamma_only` 标志控制FFT网格的分配：

```cpp
// source/source_basis/module_pw/pw_init.cpp:221-230
this->gamma_only = gamma_only_in;
// if use gamma point only, when convert real function f(r) to F(k) = FFT(f),
// we have F(-k) = F(k)*, so that only half of planewaves are needed.
this->fftny = this->ny;
this->fftnx = this->nx;
if (this->gamma_only)
{
    if(this->xprime) this->fftnx = int(this->nx / 2) + 1;
    else            this->fftny = int(this->ny / 2) + 1;
}
```

当 `gamma_only=true` 时：
- 如果 `xprime=true`：在x方向只扫描正半轴，fftnx = nx/2 + 1
- 如果 `xprime=false`：在y方向只扫描正半轴，fftny = ny/2 + 1

#### 0.2.4 Gamma Only 节省的资源

| 数据类型 | 普通模式 | Gamma Only | 节省比例 |
|---------|---------|------------|---------|
| 波函数系数 | N个复数 | N/2个复数 | 50% |
| 电荷密度(倒空间) | N个复数 | N/2个复数 | 50% |
| 势函数(倒空间) | N个复数 | N/2个复数 | 50% |
| 内存带宽 | 2N | N | 50% |
| FFT计算量 | N log N | (N/2) log(N/2) | ~60% |

#### 0.2.5 Gamma Only 的限制

Gamma Only 只能在以下条件满足时使用：
1. **k点只取Gamma点**：布里渊区中心
2. **没有磁场**：时间反演对称性成立
3. **非自旋极化或自旋简并**：没有破坏时间反演对称性
4. **非相对论近似**：自旋轨道耦合会破坏时间反演对称性

---

### 1.1 问题由来

在ABACUS中，平面波的生成和变换是计算的核心部分，主要涉及：

1. **平面波分布**（`pw_distributeg.cpp`）：将平面波分配到不同进程
2. **FFT变换**（`pw_transform.cpp`）：实空间与倒易空间变换
3. **Gather/Scatter操作**（`pw_gatherscatter.h`）：进程间数据交换
4. **波函数操作**：哈密顿量作用、能量计算等

这些操作中存在大量的并行优化空间，包括：

| 操作类型 | 现有并行情况 | 优化方向 |
|---------|------------|---------|
| 平面波分布 | 部分MPI并行 | MPI+OpenMP混合 |
| FFT变换 | OpenMP并行 | SIMD/异构加速 |
| Gather/Scatter | MPI Alltoall | 非阻塞通信 |
| 数据重排 | 串行 | OpenMP并行 |

### 1.2 现有代码结构

#### 1.2.1 平面波基类结构

```cpp
// source/source_basis/module_pw/pw_basis.h
class PW_Basis
{
public:
    int npw = 0;                    // 当前进程的平面波数量
    int npwtot = 0;                 // 总平面波数量
    int nstot = 0;                  // 总stick数量
    int nst = 0;                    // 当前进程的stick数量

    bool gamma_only = false;        // 是否使用gamma only模式
    double ggecut = 0.0;           // 能量截断

    int *ig2isz = nullptr;          // 平面波到FFT网格的映射
    int *is2fftixy = nullptr;       // stick到FFT坐标的映射
    int *fftixy2ip = nullptr;       // FFT坐标到进程的映射

    ModuleBase::Vector3<double> *gcar = nullptr;  // G矢量笛卡尔坐标
    double *gg = nullptr;          // |G|^2
};
```

#### 1.2.2 平面波分布流程

```
count_pw_st()      - 统计每个stick的平面波数量
       ↓
distribute_g()     - 将平面波分配到不同进程
       ↓
get_ig2isz_is2fftixy() - 构建索引映射
```

关键代码在 `pw_distributeg.cpp:44-143`：

```cpp
// count_pw_st() 函数中的三重循环
for (int ix = ix_start; ix <= ix_end; ++ix)
{
    for (int iy = iy_start; iy <= iy_end; ++iy)
    {
        for (int iz = iz_start; iz <= iz_end; ++iz)
        {
            // 计算 |G|^2 并判断是否在截断球内
            double modulus = f * (this->GGT * f);
            if (modulus <= this->ggecut || this->full_pw)
            {
                // 计数平面波
            }
        }
    }
}
```

**性能瓶颈**：
1. 三重嵌套循环，串行执行
2. 每次迭代都要做矩阵乘法 $f \cdot (GGT \cdot f)$
3. 没有利用OpenMP并行化

#### 1.2.3 FFT变换流程

```cpp
// pw_transform.cpp:real2recip()
real2recip() {
    // 1. 数据拷贝（OpenMP并行）
    #pragma omp parallel for
    for (int ir = 0; ir < this->nrxx; ++ir) { ... }

    // 2. XY方向的FFT（库函数，内部可能无并行）
    this->fft_bundle.fftxyfor(...);

    // 3. Gather操作（MPI Alltoallv，阻塞通信）
    this->gatherp_scatters(...);

    // 4. Z方向的FFT（库函数）
    this->fft_bundle.fftzfor(...);

    // 5. 提取平面波系数（OpenMP并行）
    #pragma omp parallel for
    for (int ig = 0; ig < this->npw; ++ig) { ... }
}
```

**性能瓶颈**：
1. `gatherp_scatters` 使用阻塞的 `MPI_Alltoallv`
2. 所有进程必须等待通信完成才能继续
3. 没有利用非阻塞通信重叠计算与通信

#### 1.2.4 Gather/Scatter 操作详解

```cpp
// pw_gatherscatter.h:gatherp_scatters()
// 单进程情况
if(this->poolnproc == 1)
{
    // 简单的内存拷贝
    for(int is = 0; is < this->nst; ++is) { ... }
}

// 多进程情况
#ifdef __MPI
// 第一步：数据重排（OpenMP并行）
#pragma omp parallel for
for (int istot = 0; istot < nstot; ++istot) { ... }

// 第二步：MPI Alltoallv 通信（阻塞）
MPI_Alltoallv(out, numr, startr, MPI_DOUBLE_COMPLEX,
              in, numg, startg, MPI_DOUBLE_COMPLEX,
              this->pool_world);

// 第三步：再次数据重排（OpenMP并行 collapse）
#pragma omp parallel for collapse(2)
for (int ip = 0; ip < this->poolnproc; ++ip)
    for (int is = 0; is < this->nst; ++is) { ... }
#endif
```

**性能瓶颈**：
1. 两处循环可进一步优化
2. `MPI_Alltoallv` 是阻塞调用
3. 没有实现计算与通信的重叠

---

### 1.3 Gamma Only 的当前实现

#### 1.3.1 Gamma Only 的FFT处理

在 `pw_transform.cpp` 中，gamma_only 模式有特殊处理：

```cpp
// pw_transform.cpp:79-136
void PW_Basis::real2recip(const FPTYPE* in, std::complex<FPTYPE>* out, ...)
{
    if (this->gamma_only)
    {
        // 使用 r2c 变换（实数到复数）
        const int npy = this->ny * this->nplane;
        #pragma omp parallel for collapse(2)
        for (int ix = 0; ix < this->nx; ++ix)
        {
            for (int ipy = 0; ipy < npy; ++ipy) { ... }
        }
        this->fft_bundle.fftxyr2c(...);  // 实数FFT
    }
    else
    {
        // 使用 c2c 变换（复数到复数）
        #pragma omp parallel for
        for (int ir = 0; ir < this->nrxx; ++ir) { ... }
        this->fft_bundle.fftxyfor(...);   // 复数FFT
    }
}
```

#### 1.3.2 PW_Basis_K 中 Gamma Only 的处理

```cpp
// pw_basis_k.cpp:81-99
this->gamma_only = gamma_only_in;
if (kmaxmod > 0)
{
    this->gamma_only = false; // 非gamma点不使用gamma_only
}
```

这意味着：
- 只有当所有k点都是Gamma点时才启用gamma_only
- 否则自动关闭gamma_only优化

#### 1.3.3 Gamma Only 的限制

从代码中可以看到：

```cpp
// pw_transform_k.cpp:34
assert(this->gamma_only == false);
```

`PW_Basis_K::real2recip` 对gamma_only有断言检查，这说明：
- 多k点情况下不支持gamma_only
- Gamma_only只在单k点（PW_Basis）情况下生效

---

### 1.4 并行架构分析

#### 1.4.1 两级并行架构

ABACUS使用MPI+OpenMP两级并行：

```
MPI_COMM_WORLD
    ├── POOL_WORLD (k点并行)
    │       ├── Pool 0: 进程 0,1,2... (平面波并行)
    │       ├── Pool 1: 进程 ...
    │       └── ...
    └── 其他进程组
```

关键变量：
| 变量 | 说明 |
|------|------|
| `GlobalV::NPROC_IN_POOL` | 每个pool的进程数 |
| `GlobalV::MY_POOL` | 当前pool编号 |
| `GlobalV::RANK_IN_POOL` | pool内rank |
| `pool_world` | pool通信子 |

#### 1.4.2 平面波并行策略

每个进程的平面波数量计算：

```cpp
// 截断球内平面波数量
npw = number of G vectors satisfying |G+k|^2 < Ecut
```

不同进程的平面波分配通过 `distribute_g()` 实现。

---

## 二、建议可以做的事情（共 8 题）

### 题目 1：平面波分布的 OpenMP 并行化

**难度**：⭐⭐

#### 题目描述

优化 `pw_distributeg.cpp` 中的 `count_pw_st` 函数，实现其 OpenMP 并行版本。

#### 现有代码位置

- `pw_distributeg.cpp:44-143` - `count_pw_st` 函数
- `pw_distributeg.cpp:13-32` - `distribute_g` 函数

#### 具体要求

1. **循环依赖分析**
   - 分析三重循环的数据依赖
   - 确定可并行的循环层级
   - 验证结果正确性

2. **OpenMP 实现**
   - 使用 `#pragma omp parallel for collapse(2)` 实现外层循环并行
   - 处理线程私有变量和归约
   - 确保线程安全

3. **性能测试**
   - 测试不同线程数（1, 2, 4, 8）下的加速比
   - 分析 Amdahl 定律
   - 对比优化前后性能

4. **代码质量**
   - 遵循项目代码规范
   - 添加必要注释

5. **单元测试要求**
   - 编写单元测试验证count_pw_st输出的正确性
   - 测试不同网格大小和截断能量下的边界情况
   - 验证并行结果与串行结果的一致性

6. **代码重构（加分项）**
   - 将平面波判断逻辑抽象为独立的函数对象
   - 使用依赖注入便于测试不同GGT矩阵实现
   - 设计可配置的截断条件判断器

#### 示例：平面波判断逻辑抽象

```cpp
// 抽象截断条件判断接口
class IPWCriterion {
public:
    virtual ~IPWCriterion() = default;
    virtual bool is_in_sphere(const ModuleBase::Vector3<double>& g) const = 0;
};

// 基于能量截断的实现
class EnergyCutoffCriterion : public IPWCriterion {
private:
    double ggecut_;
    ModuleBase::Matrix3 GGT_;
public:
    EnergyCutoffCriterion(double ggecut, const ModuleBase::Matrix3& GGT)
        : ggecut_(ggecut), GGT_(GGT) {}

    bool is_in_sphere(const ModuleBase::Vector3<double>& g) const override {
        return g * (GGT_ * g) <= ggecut_;
    }
};

// 测试中使用Mock对象
class MockPWCriterion : public IPWCriterion {
public:
    mutable std::vector<bool> results;
    bool is_in_sphere(const ModuleBase::Vector3<double>& g) const override {
        return results[g.x() * 100 + g.y() * 10 + g.z()];
    }
};
```

---

### 题目 2：MPI_Alltoallv 的非阻塞优化

**难度**：⭐⭐⭐

#### 题目描述

优化 `pw_gatherscatter.h` 中的 `gatherp_scatters` 函数，使用非阻塞通信重叠计算与通信。

#### 现有代码位置

- `pw_gatherscatter.h:14-89` - `gatherp_scatters` 函数

#### 具体要求

1. **当前实现分析**
   - 分析阻塞通信的瓶颈
   - 识别可重叠的计算部分

2. **非阻塞实现**
   - 使用 `MPI_Isend`/`MPI_Irecv` 替代阻塞通信
   - 设计计算与通信重叠的方案
   - 确保通信结果正确性

3. **性能测试**
   - 测试不同进程数下的性能提升
   - 分析通信开销占比
   - 对比阻塞与非阻塞版本

4. **兼容性**
   - 保持与其他模块的接口兼容
   - 处理边界情况

#### 参考提示

```cpp
// 非阻塞通信示例
MPI_Request requests[2];
MPI_Status statuses[2];

// 异步发送
MPI_Isend(send_buf, count, MPI_DOUBLE_COMPLEX, dest,
          tag, pool_world, &requests[0]);

// 在等待时执行其他计算
do_local_computation();

// 等待通信完成
MPI_Waitall(2, requests, statuses);
```

---

### 题目 3：FFT 变换的 OpenMP 优化

**难度**：⭐⭐

#### 题目描述

优化 `pw_transform.cpp` 中 FFT 变换相关的数据拷贝和重排部分，提升 OpenMP 并行效率。

#### 现有代码位置

- `pw_transform.cpp:24-136` - `real2recip` 函数
- `pw_transform.cpp:145-286` - `recip2real` 函数

#### 具体要求

1. **性能热点分析**
   - 使用 timer 识别热点函数
   - 分析 OpenMP 的并行效率

2. **优化实现**
   - 优化循环分块以提高 cache 命中率
   - 使用 SIMD 指令加速（可选）
   - 优化内存访问模式

3. **性能测试**
   - 对比不同线程数的加速比
   - 分析并行效率

4. **正确性验证**
   - 确保浮点结果一致
   - 处理边界情况

#### 参考提示

```cpp
// 当前实现
#pragma omp parallel for schedule(static)
for (int ir = 0; ir < this->nrxx; ++ir)
{
    auxr[ir] = in[ir];
}

// 优化方向
// 1. 考虑 cache line 对齐
// 2. 调整 schedule 策略
// 3. 使用向量化
```

---

### 题目 4：Gamma Only 模式的完整实现

**难度**：⭐⭐⭐⭐

#### 题目描述

在 `PW_Basis_K` 中实现完整的 gamma_only 模式，使得多k点计算时也能利用gamma_only节省一半存储。

#### 题目背景

当前实现中，`PW_Basis_K::real2recip` 对 gamma_only 有断言检查，说明多k点情况下 gamma_only 尚未实现。但理论上，当所有k点都位于高对称点（如Gamma点）时，可以利用时间反演对称性节省存储。

#### 具体要求

1. **理论分析**
   - 分析多k点情况下时间反演对称性的适用条件
   - 确定哪些操作可以共享 gamma_only 优化
   - 设计数据结构支持混合k点（gamma + non-gamma）

2. **代码实现**
   - 修改 `PW_Basis_K` 支持 gamma_only
   - 实现对应的 FFT 变换
   - 确保与现有代码兼容

3. **测试验证**
   - 对比 gamma_only 与非 gamma_only 的计算结果
   - 测试不同体系
   - 验证内存节省效果

4. **性能评估**
   - 测量内存节省比例
   - 测量计算效率变化

#### 参考提示

```cpp
// 关键修改点
// 1. pw_basis_k.h - 添加 gamma_only 支持
void initparameters(..., const bool gamma_only_in, ...);

// 2. pw_basis_k.cpp - 实现 gamma_only FFT
template <typename FPTYPE>
void PW_Basis_K::real2recip_gamma(const FPTYPE* in,
                                   std::complex<FPTYPE>* out,
                                   const int ik) const;

// 3. 利用 gamma_only 时：
//    - xprime=true: fftnx = nx/2 + 1
//    - xprime=false: fftny = ny/2 + 1
```

---

### 题目 5：平面波 gather 操作的向量化

**难度**：⭐⭐⭐

#### 题目描述

优化 `pw_gatherscatter.h` 中的数据重排循环，使用 SIMD 指令实现向量化加速。

#### 现有代码位置

- `pw_gatherscatter.h:44-54` - 第一处数据重排
- `pw_gatherscatter.h:71-85` - 第二处数据重排

#### 具体要求

1. **向量化分析**
   - 分析数据访问模式
   - 确定向量化可行性

2. **SIMD 实现**
   - 使用 `#pragma omp simd` 指导编译器向量化
   - 或使用 Intrinsics 直接编写向量化代码
   - 确保数据对齐

3. **性能测试**
   - 对比向量化前后的性能
   - 分析加速比

4. **跨平台支持**
   - 支持 AVX2/AVX512
   - 考虑 ARM NEON（可选）

#### 参考提示

```cpp
// 原始代码
for (int iz = 0; iz < nplane; ++iz)
{
    outp[iz] = inp[iz];
}

// 向量化版本
#pragma omp simd aligned(outp, inp: 32)
for (int iz = 0; iz < nplane; ++iz)
{
    outp[iz] = inp[iz];
}
```

---

### 题目 6：电荷密度与波函数的 Gamma Only 存储优化

**难度**：⭐⭐⭐

#### 题目描述

实现 gamma_only 模式下电荷密度、势函数和波函数的紧凑存储格式，进一步减少内存占用。

#### 题目背景

当使用 gamma_only 时，倒空间数据存在共轭对称性：
- $\rho(-\mathbf{G}) = \rho^*(\mathbf{G})$
- $V(-\mathbf{G}) = V^*(\mathbf{G})$
- $\psi(-\mathbf{G}) = \psi^*(\mathbf{G})$

可以只存储一半的数据来节省内存。

#### 具体要求

1. **存储格式设计**
   - 设计紧凑的存储格式
   - 保持与现有接口兼容
   - 处理边界情况（G=0）

2. **实现方案**
   - 实现紧凑型数据结构
   - 修改相关函数支持新格式
   - 添加数据压缩/解压函数

3. **测试验证**
   - 验证紧凑存储的正确性
   - 测试内存节省效果
   - 确保计算结果一致

4. **性能评估**
   - 测量内存节省
   - 测量压缩/解压开销

#### 参考提示

```cpp
// 紧凑存储示例
class CompactGammaData {
private:
    int ngk;                      // 存储的平面波数量 (ngk = npw/2)
    std::complex<double>* data;   // 紧凑存储的数据
public:
    // G>=0 的数据直接存储
    // G<0 的数据通过共轭关系推导
    std::complex<double> get_psi(int ig) const {
        if (ig <= ngk) return data[ig];
        else return std::conj(data[ngk - ig]);
    }
};
```

---

### 题目 7：FFT 通信与计算的Overlap优化

**难度**：⭐⭐⭐

#### 题目描述

实现 FFT 变换中 MPI 通信与计算的Overlap，使通信和计算可以同时进行。

#### 题目背景

当前 FFT 变换流程：
```
1. 本地数据拷贝
2. XY-FFT
3. MPI_Alltoallv (阻塞)
4. Z-FFT
5. 提取平面波系数
```

步骤3是阻塞的，所有进程必须等待通信完成。可以将本地计算与通信重叠来隐藏延迟。

#### 具体要求

1. **Overlap 设计**
   - 识别可与通信重叠的计算
   - 设计双缓冲机制
   - 实现流水线处理

2. **实现方案**
   - 使用 MPI_Isend/Irecv 实现非阻塞通信
   - 实现 Double Buffering
   - 协调进程间的数据流

3. **性能测试**
   - 测量通信隐藏比例
   - 测试端到端性能提升

4. **正确性验证**
   - 确保结果与原实现一致
   - 处理边界情况

#### 参考提示

```cpp
// Double Buffering 方案
// Buffer A: 用于当前通信
// Buffer B: 用于同时计算

// Step 1: 通信Buffer A, 计算Buffer B
MPI_Irecv(bufA, ..., &req_recv);
compute(bufB);

// Step 2: 等待通信完成, 通信Buffer B
MPI_Wait(&req_recv, ...);
MPI_Irecv(bufB, ..., &req_recv);
compute(bufA);

// Step 3: 重复...
```

---

### 题目 8：平面波预计算与缓存优化

**难度**：⭐⭐

#### 题目描述

优化平面波相关的重复计算，如G矢量、gk2模长等，实现计算结果缓存。

#### 题目背景

在自洽场迭代中，某些平面波相关的数据在多个迭代步骤中保持不变，可以一次性计算并缓存复用。

#### 具体要求

1. **缓存分析**
   - 识别可缓存的数据
   - 确定缓存生命周期

2. **实现方案**
   - 在 PW_Basis 类中添加缓存成员
   - 实现懒加载机制
   - 处理缓存失效

3. **性能测试**
   - 测量缓存命中率
   - 测试内存开销

4. **代码质量**
   - 确保线程安全
   - 遵循RAII原则

#### 参考提示

```cpp
// 可缓存的数据
class PW_Basis {
private:
    bool gcar_cached = false;
    ModuleBase::Vector3<double>* gcar_cache = nullptr;

public:
    // 懒加载获取 gcar
    const ModuleBase::Vector3<double>* get_gcar() {
        if (!gcar_cached) {
            compute_gcar();
            gcar_cached = true;
        }
        return gcar_cache;
    }
};
```

---

## 三、附加题：Gamma Only 算法实现

### 附加题：完整的 Gamma Only 模式

**难度**：⭐⭐⭐⭐⭐

#### 背景知识

##### G点只取Gamma时的优化原理

当k点固定为Gamma点（$\mathbf{k}=0$）时，哈密顿量具有时间反演对称性：

$$H(\mathbf{r}) = H^*(\mathbf{r})$$

对应的本征函数满足：

$$\psi(\mathbf{r}) = \psi^*(\mathbf{r})$$

这意味着实空间波函数的FFT满足共轭对称性：

$$c_{-\mathbf{G}} = c_{\mathbf{G}}^*$$

**数学推导**：

对于实函数 $f(\mathbf{r})$，其傅里叶变换为：

$$F(\mathbf{G}) = \int f(\mathbf{r}) e^{-i\mathbf{G}\cdot\mathbf{r}} d\mathbf{r}$$

取共轭：

$$F^*(\mathbf{G}) = \int f^*(\mathbf{r}) e^{i\mathbf{G}\cdot\mathbf{r}} d\mathbf{r} = \int f(\mathbf{r}) e^{i\mathbf{G}\cdot\mathbf{r}} d\mathbf{r} = F(-\mathbf{G})$$

因此：

$$F(-\mathbf{G}) = F^*(\mathbf{G})$$

##### 内存节省分析

| 数据类型 | 普通模式存储 | Gamma Only存储 | 节省比例 |
|---------|------------|---------------|---------|
| 波函数系数 | $\psi[0...N-1]$ | $\psi[0...N/2-1]$ | 50% |
| 电荷密度(倒空间) | $\rho[0...N-1]$ | $\rho[0...N/2-1]$ | 50% |
| 势能(倒空间) | $V[0...N-1]$ | $V[0...N/2-1]$ | 50% |

##### FFT优化

对于实数输入的FFT，存在专门的优化算法：

- **Real-to-Complex FFT (r2c)**：将 $N$ 点实数序列转换为 $N/2+1$ 点复数序列
- **Complex-to-Real FFT (c2r)**：反向变换

在FFTW库中，可以使用 `FFTW_R2HC` 格式直接获取半复谱。

##### 实现要点

1. **网格设置**
```cpp
// Gamma only 时的FFT网格
if (gamma_only) {
    if (xprime) {
        fftnx = nx / 2 + 1;  // x方向只存储正半部分
    } else {
        fftny = ny / 2 + 1;  // y方向只存储正半部分
    }
}
```

2. **数据索引映射**
```cpp
// 普通模式：G向量索引 ig 范围 [0, npw)
// Gamma only：只需存储 ig < npw/2 的系数
//另一半通过共轭关系获得
```

3. **逆变换处理**
```cpp
// 从紧凑存储恢复完整数据
for (int ig = 0; ig < npw; ++ig) {
    if (ig < npw_gamma) {
        psi_full[ig] = psi_gamma[ig];
    } else {
        psi_full[ig] = conj(psi_gamma[npw - ig]);
    }
}
```

#### 作业要求

##### 必做部分

1. **理论分析**（20分）
   - 推导 gamma_only 节省存储的数学原理
   - 分析 FFT r2c/c2r 变换的优化空间
   - 评估对计算精度的影响

2. **代码实现**（40分）
   - 在 `PW_Basis` 类中实现 gamma_only 的紧凑存储
   - 实现从紧凑格式到完整格式的转换函数
   - 修改 FFT 变换支持紧凑格式
   - 确保与现有接口兼容

3. **测试验证**（20分）
   - 对比 gamma_only 与普通模式的计算结果
   - 验证正确性（能量误差 < 1e-6 Ry）
   - 测试内存节省效果

##### 选做部分（额外加分）

4. **性能优化**（20分）
   - 实现 SIMD 加速的数据转换
   - 测量并对比性能
   - 提交性能报告

#### 参考代码结构

```cpp
// 紧凑存储的平面波系数
template <typename FPTYPE>
class GammaOnlyPW {
private:
    int npw_gamma;  // 紧凑存储的平面波数 = npw/2
    std::complex<FPTYPE>* psi_gamma;

public:
    // 从完整格式构造
    GammaOnlyPW(const std::complex<FPTYPE>* psi_full, int npw);

    // 恢复到完整格式
    void expand(std::complex<FPTYPE>* psi_full) const;

    // FFT支持
    void fft_r2c(const FPTYPE* rdata);
    void fft_c2r(FPTYPE* rdata) const;
};
```

#### 评分标准

| 完成度 | 分数 | 要求 |
|-------|------|------|
| 理论分析 | 20 | 正确推导数学原理 |
| 基础实现 | 40 | 实现紧凑存储和转换 |
| 正确性验证 | 20 | 结果误差 < 1e-6 |
| 性能优化（选做） | 20 | 性能提升明显 |

---

## 四、测试环境与基准数据

### 4.1 推荐测试体系

| 体系 | 原子数 | 平面波数量 | 网格大小 | 备注 |
|------|--------|-----------|---------|------|
| Al fcc | 1 | ~1000 | 36×36×36 | Gamma only 适用 |
| Si diamond | 2 | ~2000 | 40×40×40 | Gamma only 适用 |
| Fe bcc | 1 | ~1500 | 38×38×38 | 自旋极化，需验证 |
| NaCl rocksalt | 2 | ~1800 | 42×42×42 | Gamma only 适用 |

### 4.2 性能基准

| 优化项 | 当前时间 | 目标时间 | 最低加速比 |
|--------|---------|---------|-----------|
| count_pw_st | T₁ | T₁/4 | 4x (4线程) |
| MPI_Alltoallv | T₂ | T₂/2 (Overlap) | 2x |
| FFT 数据拷贝 | T₃ | T₃/4 | 4x (4线程) |
| Gamma Only 存储 | - | -50% | 50% 内存节省 |

### 4.3 测试脚本参考

```bash
#!/bin/bash
# benchmark_pw.sh - 平面波性能测试

export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8

for nproc in 1 2 4 8; do
    for nthread in 1 2 4 8; do
        echo "Testing: nproc=$nproc, nthread=$nthread"
        export OMP_NUM_THREADS=$nthread
        mpirun -np $nproc ./abacus INPUT > log_p${nproc}_t${nthread}.out 2>&1
        grep "plane wave" log_p${nproc}_t${nthread}.out | tail -1
    done
done
```

---

## 五、代码规范与提交流流程

### 5.1 代码规范

1. **命名规范**
   - 遵循项目现有的命名风格
   - 新增函数需添加文档注释

2. **模块化设计**
   - 独立功能封装为独立函数/类
   - 便于单元测试

3. **错误处理**
   - 检查所有 MPI 调用返回值
   - 妥善处理异常情况

### 5.2 提交流程

#### 5.2.1 推荐方式：GitHub Pull Request

1. **Fork 仓库**
   - Fork ABACUS 仓库到你自己的 GitHub 账户

2. **创建分支**
   ```bash
   git checkout -b feature/plane-wave-optimization
   ```

3. **少量多次提交**
   ```bash
   git add source/source_basis/module_pw/
   git commit -m "Add OpenMP parallelization for count_pw_st"
   git push origin feature/plane-wave-optimization
   ```

4. **提交 Pull Request**
   - 在 GitHub 上创建 Pull Request
   - 描述你的优化内容

#### 5.2.2 提交策略

| 原则 | 说明 |
|------|------|
| **少量多次** | 每完成一个小功能就提交 |
| **问题导向** | 每个 PR 解决一个具体问题 |
| **文档完善** | PR 描述中说明优化内容 |

---

## 六、参考资料

### 6.1 代码位置索引

| 文件 | 路径 | 说明 |
|------|------|------|
| PW_Basis 类 | `source/source_basis/module_pw/pw_basis.h` | 平面波基类定义 |
| PW_Basis_K 类 | `source/source_basis/module_pw/pw_basis_k.h` | 多k点平面波基类 |
| 平面波分布 | `source/source_basis/module_pw/pw_distributeg.cpp` | 平面波分布 |
| FFT 变换 | `source/source_basis/module_pw/pw_transform.cpp` | 实-倒空间变换 |
| Gather/Scatter | `source/source_basis/module_pw/pw_gatherscatter.h` | 进程间数据交换 |

### 6.2 推荐阅读

1. **FFT**: FFTW 文档 - Real-to-Complex FFT 优化
2. **MPI**: MPI_Alltoallv 非阻塞版本 MPI_Ialltoallv
3. **OpenMP**: SIMD 向量化编程指南
4. **Gamma Only**: ABACUS 理论手册 - 平面波基组章节

---

## 七、致谢

本大作业题目设计参考了以下资源：

1. ABACUS 软件源代码 (https://github.com/abacusmodeling/abacus-develop)
2. FFTW 库文档 - Real-to-Complex FFT 变换
3. MPI 标准文档 - 非阻塞通信
4. 《Computational Materials Science》- 平面波方法章节

---

**最后更新**：2026-04-21

**版本**：v1.0
