# PPCG 特征值求解算法阶段性文档

## 一、任务背景

本阶段选择的问题是实现 PPCG（Projected Preconditioned Conjugate Gradient）方法，用于优化 ABACUS 中特征值问题的迭代求解过程。特征值求解是电子结构计算中的核心步骤，尤其在平面波基组下，Hamiltonian 与波函数的乘法、残差计算和正交化会占用大量计算时间。因此，在已有 CG、BPCG 和 Davidson 方法的基础上理解原算法，是设计 PPCG 方法的前提。

目前我主要阅读了 `source_hsolver` 目录下与迭代对角化相关的代码，包括：

- `hsolver_pw.cpp`
- `diago_cg.h / diago_cg.cpp`
- `diago_bpcg.h / diago_bpcg.cpp`
- `diago_david.h / diago_david.cpp`
- `diago_dav_subspace.h / diago_dav_subspace.cpp`
- `diago_iter_assist.h / diago_iter_assist.cpp`
- `kernels/bpcg_kernel_op.cpp`

其中，`diago_bpcg.cpp` 与本题最相关，因为它已经实现了 block 形式的预条件共轭梯度方法，可以作为 PPCG 的主要参考。同时，Davidson 相关代码对理解“投影子空间”也很重要。

## 二、现有代码结构理解

在平面波基组下，特征值求解的入口主要在 `hsolver_pw.cpp` 中。程序会根据输入参数选择不同的对角化方法，例如：

```cpp
cg
bpcg
dav
dav_subspace
```

这些方法共享两个重要操作：

```text
hpsi_func : 计算 H * psi
spsi_func : 计算 S * psi
```

其中 `hpsi_func` 是最核心的计算步骤，因为它对应 Hamiltonian 与波函数的乘法，也是迭代法中最耗时的部分。`spsi_func` 用来处理广义特征值问题中的重叠矩阵 `S`。

预条件器由 `HSolverPW::update_precondition` 生成，主要和动能项 `g2kin` 有关。对于 CG 和 BPCG 方法，预条件器的形式大致为：

```text
M = 1 + g2kin + sqrt(1 + (g2kin - 1)^2)
```

后续求解过程中会通过除以这个对角预条件器来改善收敛速度。

## 三、CG 方法原理

`DiagoCG` 是当前代码中的逐能带预条件共轭梯度方法。它一次只处理一条 band，因此逻辑比较清晰，但并行性和矩阵块操作效率有限。

它的基本流程可以概括为：

1. 对初始波函数做子空间对角化，得到较好的初始猜测。
2. 对每一条 band 单独进行迭代。
3. 计算当前波函数的 `H psi` 和 `S psi`。
4. 根据残差构造预条件梯度。
5. 将梯度与已经求出的低能态正交。
6. 更新共轭方向。
7. 在当前波函数和共轭方向张成的二维空间内做线搜索。
8. 判断本征值变化是否小于阈值。

从数学上看，CG 方法求解的是：

```text
H x = lambda S x
```

残差可以理解为：

```text
r = Hx - lambda Sx
```

预条件的作用是近似求解：

```text
M^{-1} r
```

这样可以让搜索方向更接近误差方向，从而加快收敛。

CG 方法的优点是内存占用较低，算法比较稳定；缺点是逐 band 处理，无法充分利用 block BLAS 和多能带之间的整体信息。

## 四、BPCG 方法原理

`DiagoBPCG` 可以看作 CG 方法的 block 版本。它不再逐条 band 单独处理，而是把多个 band 组成一个波函数块一起迭代。

在代码中，BPCG 主要维护以下数据：

```text
psi       当前波函数
hpsi      H * psi
grad      当前梯度或搜索方向
grad_old  上一步搜索方向
hgrad     H * grad
hsub      子空间 Hamiltonian 小矩阵
eigen     当前本征值
err_st    每条 band 的误差
```

它的主要流程是：

1. 首先计算 `hpsi = H psi`。
2. 构造小矩阵：

```text
hsub = psi^H H psi
```

3. 对 `hsub` 做一次小规模对角化，并旋转波函数，改善初始波函数。
4. 计算每条 band 的残差：

```text
r_i = H psi_i - epsilon_i psi_i
```

5. 使用预条件器得到梯度方向：

```text
grad_i = - r_i / M
```

6. 加入上一轮方向，形成类似共轭梯度的更新：

```text
grad_i = - r_i / M + beta_i grad_old_i
```

7. 将 `grad` 对当前 `psi` 做正交投影。
8. 计算 `hgrad = H grad`。
9. 在 `psi_i` 和 `grad_i` 张成的二维空间内做线搜索。
10. 对整个 `psi` block 重新正交化。
11. 重复迭代直到误差满足阈值。

相比 `DiagoCG`，BPCG 的主要优势是 block 化，可以一次处理多条 band，更适合并行计算和矩阵乘法优化。

不过当前 BPCG 仍然存在一个限制：虽然数据结构是 block 的，但每条 band 的更新仍然主要是在二维空间 `span{psi_i, grad_i}` 内完成的，还没有真正构造更大的投影子空间。

## 五、Davidson 方法原理

ABACUS 中和 Davidson 有关的实现主要有两个：普通 Davidson，即 `DiagoDavid`；以及 `Diago_DavSubspace`，对应输入方法中的 `dav_subspace`。二者都属于投影子空间方法，基本思想是不断扩展一个较小的子空间，在这个子空间中求解小规模特征值问题。

### 5.1 普通 Davidson

普通 Davidson 的实现位于 `diago_david.cpp`。它求解的问题形式是：

```text
H X = S X Lambda
```

其核心思想可以概括为：

1. 先对初始波函数做 Schmidt 正交化，得到初始子空间基 `basis`。
2. 计算：

```text
H basis
S basis
```

3. 在当前子空间中构造小矩阵，并求解小规模特征值问题。
4. 根据本征值变化判断哪些 band 尚未收敛。
5. 对未收敛的 band 构造残差：

```text
r = (H - lambda S) x
```

6. 对残差做预条件，得到新的修正方向。
7. 将新的方向正交化后加入子空间。
8. 子空间过大时进行 refresh，用当前 Ritz 向量重启子空间。

普通 Davidson 的特点是子空间会逐步增长。每次迭代只对未收敛的 band 增加新的方向，因此在收敛过程中可以避免处理已经收敛的部分。它的关键步骤是残差修正：

```text
w = M^{-1} (H - lambda S) x
```

这里的 `M` 是对 Hamiltonian 的近似对角预条件器。这个思想和 PPCG 中的预条件残差 `W` 非常接近。

普通 Davidson 的优势是收敛通常比较稳健，尤其适合求解少量低能本征态；缺点是子空间维度会增长，需要定期重启，并且小矩阵对角化和正交化的开销会随子空间大小增加。

### 5.2 DavSubspace 方法

`Diago_DavSubspace` 是另一套 Davidson 子空间实现，代码位于 `diago_dav_subspace.cpp`。它和普通 `DiagoDavid` 的主要思想相同，但在子空间矩阵构造和小矩阵求解上更强调统一的子空间处理。

在 `dav_subspace` 中，程序显式维护：

```text
psi_iter  子空间基
hpsi      H * psi_iter
spsi      S * psi_iter
hcc       子空间 Hamiltonian 矩阵
scc       子空间 overlap 矩阵
vcc       子空间特征向量
```

每一轮迭代中，先在当前子空间中构造：

```text
H_c = V^H H V
S_c = V^H S V
```

然后求解小规模广义特征值问题：

```text
H_c c = lambda S_c c
```

得到 Ritz 值和 Ritz 向量后，再根据未收敛的 band 构造残差和修正方向。与普通 Davidson 相比，`dav_subspace` 更明确地把 `H_c` 和 `S_c` 都作为子空间矩阵维护，因此更适合处理广义特征值问题。

另外，`dav_subspace` 的小矩阵对角化后端可以选择不同实现：

```text
diag_subspace = 0 : LAPACK
diag_subspace = 1 : Gen-ELPA
diag_subspace = 2 : ScaLAPACK
```

这说明 `dav_subspace` 主要考虑的是当子空间矩阵较大或并行规模较大时，小矩阵对角化本身也可能成为性能瓶颈，需要使用并行对角化库。

从 PPCG 的角度看，`dav_subspace` 的参考价值在于：它展示了如何构造和维护投影子空间中的 `H_c`、`S_c`，以及如何在小空间中求解广义特征值问题。PPCG 也需要类似的小空间 Rayleigh-Ritz 过程，只是 PPCG 的子空间通常固定为：

```text
span{X, W, P}
```

而 Davidson 的子空间则会随迭代不断扩展。

## 六、PPCG 算法设计

根据对 CG、BPCG 和 Davidson 的理解，PPCG 可以设计为当前 BPCG 方法的进一步改进。它的核心区别是：不再只对每条 band 做二维线搜索，而是在由 `X`、`W`、`P` 构成的投影子空间中进行 Rayleigh-Ritz 对角化。

设当前波函数块为：

```text
X = [x_1, x_2, ..., x_n]
```

对应的本征值为：

```text
Lambda = diag(lambda_1, lambda_2, ..., lambda_n)
```

首先计算残差：

```text
R = H X - S X Lambda
```

然后对残差做预条件：

```text
W = - M^{-1} R
```

其中 `M` 可以先复用当前代码中的对角预条件器。

如果已有上一轮搜索方向 `P`，则构造投影子空间：

```text
K = [X, W, P]
```

第一轮没有 `P` 时，可以使用：

```text
K = [X, W]
```

接下来在该子空间内构造小矩阵：

```text
H_k = K^H H K
S_k = K^H S K
```

并求解小规模广义特征值问题：

```text
H_k C = S_k C Lambda
```

求得系数矩阵 `C` 后，用它更新波函数：

```text
X_new = K C
```

同时更新搜索方向 `P`，用于下一轮迭代。

因此，PPCG 每次迭代不是只在单条 band 的二维空间里寻找更优方向，而是在所有 band 共同构成的投影空间中统一优化。这也是它相比 BPCG 更有潜力的地方。

## 七、与现有算法的关系

当前 BPCG 的更新方式可以简化理解为：

```text
psi_i 在 span{psi_i, grad_i} 中更新
```

而 PPCG 的更新方式是：

```text
X 在 span{X, W, P} 中更新
```

普通 Davidson 的更新方式可以理解为：

```text
不断扩展 basis，并在 basis 中求解投影特征值问题
```

所以 PPCG 处在 CG/BPCG 和 Davidson 之间：它保留了预条件共轭梯度中的搜索方向 `P`，同时也使用 Davidson 类似的投影子空间思想。但它不像 Davidson 那样让子空间持续增长，而是每轮主要使用 `X`、`W`、`P` 组成的小空间。

这样做的好处是：

1. 比逐 band 线搜索能利用更多 block 内信息。
2. 对近简并本征值问题可能更稳定。
3. Rayleigh-Ritz 投影更新比单独二维线搜索更系统。
4. 子空间大小相对固定，内存开销比 Davidson 的增长型子空间更容易控制。

## 八、性能瓶颈分析

从现有代码和算法流程看，特征值迭代求解中的主要瓶颈集中在以下几个方面。

第一，`H * psi` 是最主要的计算开销。无论 CG、BPCG、Davidson 还是 PPCG，每轮迭代都需要多次调用 `hpsi_func`。在平面波基组下，这一步通常包含 FFT、局域势、非局域赝势等操作，因此是整体耗时的核心。

第二，正交化和子空间矩阵构造会带来较多全局归约。比如计算：

```text
psi^H H psi
K^H H K
K^H S K
```

都需要内积和矩阵乘法。在 MPI 并行下，这些操作往往伴随 `reduce` 或通信同步。进程数增加后，通信开销会逐渐明显。

第三，小矩阵对角化也可能成为瓶颈。对于 CG 和 BPCG，这个开销相对较小；但 Davidson 和 PPCG 都需要在投影子空间中求解小规模特征值问题。特别是 `dav_subspace` 中已经提供 LAPACK、Gen-ELPA、ScaLAPACK 等不同后端，说明当子空间维度较大时，小矩阵对角化需要并行库支持。

第四，内存访问和临时数组也会影响性能。BPCG、Davidson 和 PPCG 都需要保存 `psi`、`hpsi`、残差、搜索方向以及小空间矩阵。如果频繁复制或重排这些数组，会增加额外开销。GPU 情况下还要考虑 host/device 数据同步。

第五，收敛性本身也会影响总耗时。单次迭代快并不一定总时间最短，如果迭代步数很多，总体仍然较慢。PPCG 的目标就是通过更大的投影空间减少迭代次数，但它每轮的小空间构造和对角化又比 BPCG 更贵。因此 PPCG 的性能关键在于平衡“单步开销”和“收敛步数”。

综合来看，PPCG 的优化重点应该是减少不必要的 `H * psi` 调用、提高 block 矩阵操作效率、控制投影子空间大小，并尽量降低正交化和小矩阵对角化带来的通信开销。

## 九、阶段性总结

通过阅读现有代码，我认为 PPCG 最适合在 `DiagoBPCG` 的基础上理解和设计。当前 BPCG 已经具备 block 波函数、预条件残差、正交化和并行矩阵操作等基础，但它的核心更新仍然偏向逐 band 的二维线搜索。

Davidson 和 `dav_subspace` 则提供了投影子空间方法的参考：通过构造小空间矩阵并进行 Rayleigh-Ritz 对角化，可以在较小维度内获得更好的 Ritz 向量。PPCG 的主要思想正是把 BPCG 的预条件共轭梯度方向和 Davidson 的投影子空间更新结合起来。

因此，PPCG 的关键是引入 `span{X, W, P}` 投影子空间，并在该子空间中进行 Rayleigh-Ritz 对角化。这样可以更充分地利用 block 方法的优势，也更符合本题“Projected Preconditioned Conjugate Gradient”的算法思想。
