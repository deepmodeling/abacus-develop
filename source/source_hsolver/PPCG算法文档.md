# PPCG 算法文档

按照原论文，分为一个基础版本和在此基础上的若干改进，可以先实现基础版本，再逐步实现改进版本和并行版本.

## 基础版本

1. 算法输入：厄密特矩阵 $A\in\mathbb{C}^{n\times n}$，一个预条件器 $T$ 是对 $A^{-1}$ 的近似，想求的最小特征值个数 $k$.

2. 算法初始化：生成 $X\in\mathbb{C}^{n\times k}$ 作为特征向量的初始近似，其中 $X$ 还满足正交性 ${X}^{H}X=I$.[1]

3. 算法迭代：在未收敛的情况下，不断迭代：
    1. 计算 $W=T(AX-X(X^HAX))$
    2. 计算 $W=(I-XX^H)W$
    3. 计算 $P=(I-XX^H)W$
    4. 对 $j\in\{1, \ldots ,k\}$，计算：
        1. $S=[x_j,w_j,p_j]$
        2. 通过求解 $3\times 3$ 的特征值问题，得到 $\alpha_j,\beta_j,\gamma_j$. [2]
        3. $p_j=\beta_jw_j+\gamma_jp_j$
        4. $\bar{x}_j=\alpha_jx_j+p_j$
    5. 对 $\bar{X}$ 进行正交化，得到新的估计值 $X$. [3]

### 算法细节
[1] 这里的正交性如何保证？先生成随机的，再用正交化算法？直接用前 $k$ 个标准正交基可以吗？
[2] 这里具体是怎么求解？
- $\alpha_j,\beta_j,\gamma_j=\arg\min\limits_{||\bar{x}_j||=1}\bar{x}_j^H A \bar{x}_j$
令 $c=(\alpha_j,\beta_j,\gamma_j)^T$，则 $\bar{x}_j=Sc$，根据 Lagrange 乘子法，考虑 $f(c,\lambda)=c^HS^HASc-\lambda c^HS^HSc$，则 $\dfrac{\mathrm{d} f}{\mathrm{d} c}=2(S^HASc-\lambda S^HSc)$. 相当于求解广义的特征值问题 $S^HASc=\lambda S^HSc$，由于 $S$ 的列数为 3，所以是一个 $3\times 3$ 的特征值问题。调用 LAPACK 的函数进行求解.

[3] 这里使用对 $\bar{X}$ 进行 QR 分解，分解得到的 $Q$ 作为新的 $X$.

## 改进版本
### 改进一：使用分块对角阵加速 3. iv. 步
具体地，设分块对角阵 $C_X=\operatorname{diag}\{C_{X_1}, \ldots ,C_{X_s}\}$，$C_W=\operatorname{diag}\{C_{W_1}, \ldots ,C_{W_s}\}$，$C_P=\operatorname{diag}\{C_{P_1}, \ldots ,C_{P_s}\}$，设第 $i$ 个块大小为 $k_i$，用同样的块大小划分 $X,W,P$，3. iv. 步骤改为：
- 对 $j\in\{1, \ldots ,s\}$，计算：
    a. 令 $S=[X_j,W_j,P_j]$，$C=\begin{pmatrix}C_{X_j}\\C_{W_j}\\C_{P_j}\end{pmatrix}$
    b. 求前 $k_i$ 个广义特征值 $S^HASC=\Lambda S^HSC$
    c. 令 $P_j=W_jC_{W_j}+P_jC_{P_j}$
    d. 令 $X_j=X_jC_{X_j}+P_j$

大体上转化为求解 $s$ 个 $3k_i\times 3k_i$ 的前 $k_i$ 个广义特征值问题。**最需要讨论的点：如何优化 $k_i$ 的选取？** 单就一轮而言，肯定是 $k_i=1$ 达到最好的效果，回到了基础版本的情况。但是精心选取的 $k_i$ 可以减少迭代次数，从而提高效率。

### 改进二：引入额外特征向量
具体地，如果 $k^{\text{th}}$ 特征值和 $(k+1)^{\text{th}}$ 特征值之间的间隔较小，算法收敛会比较慢，因此可以考虑求解 $k'=k+l$ 个特征值，但是只关注前 $k$ 个特征值的收敛情况。一般取 $\frac{l}{k}=1\%\sim 5\%$.

### 改进三：正交化的再考虑

在 $\bar{X}$ 的正交性较差时，直接使用基于 Cholesky 分解的 QR 算法即可：求单位上三角阵 $R$ 使得 $\bar{X}^H\bar{X}=R^HR$，再迭代 $\bar{X}\leftarrow \bar{X}R^{-1}$

如果 $\bar{X}$ 的正交性已经较好，可以考虑基于 Taylor 展开的正交化算法：令 $\bar{X}=X(X^HX)^{-0.5}$，其中 $X^HX=I+Y$，$Y$ 的范数较小，根据 Taylor 展开就有
$$
\bar{X}\leftarrow \bar{X}(I-\frac{Y}{2}+\frac{3Y^2}{8}-\frac{5Y^3}{16}+\cdots),Y=\bar{X}^H\bar{X}-I
$$

文章还发现，其实每次跑到 3.v. 时 $\bar{X}$ 的正交性已经比较好，因此可以采取周期性正交化的方法，每 $l$ 次才执行一次正交化算法，其余时候直接用 $\bar{X}$ 来代替 $X$.

**额外的改进方法：开发一套快速判断 $\bar{X}$ 正交性的方法，如果判断出来正交性还不错，就不做正交化了**

### 改进四：引入周期性 Rayleigh-Ritz 步骤
定期对整个矩阵做 RR 步骤，来加速收敛。

### 改进五：锁定已收敛的特征向量
当某个特征向量已经收敛时，可以将其锁定。同时在迭代空间中去掉这个特征向量对应的子空间（通过投影算子 $I-X_{\text{lock}}^HX_{\text{lock}}$）。

### 改进后的伪代码
```
输入：厄密特阵 A，要求解的特征值个数 k，预条件器 T
超参：分块方案 k_i，额外特征值个数 l，RR 方法周期 rr_period
初始化：W:=AX-X(X^HAX),X_{lock}={},J_{lock}={}
while not converged do:
    W:=TW\
    W:=(I-XX^H)W; W:=(I-X_{lock}X_{lock}^H)W
    P:=(I-XX^H)W; P:=(I-X_{lock}X_{lock}^H)P
    for j in {1,...,s} do:
        S:=[X_j,W_j,P_j],C=(C_X \\ C_W \\ C_P)
        求解前 k_i 个广义特征值问题 S^HASC=\Lambda S^HSC
        P_j:=W_jC_W+P_jC_P
        X_j:=X_jC_X+P_j
    if iter mod rr_period == 0 do: #周期性 RR 步骤
        S:=[X,X_{lock}]
        求解前 k 个广义特征值问题 S^HASC=\Lambda S^HSC
        X:=SC
        W:=AX-X\Lambda
        根据 W 的范数，判断哪些已经收敛了，更新 X,X_{lock},J_{lock},W,P
        更新分块方案 k_i
    else do:
        对 X 进行正交化*
        W:=AX-X(X^HAX)
最后再做一次 RR，得到最后的特征值和特征向量.
```
