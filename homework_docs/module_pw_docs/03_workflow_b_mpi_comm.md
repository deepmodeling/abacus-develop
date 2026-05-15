# Workflow B: 数据分布式整合与 MPI 通信瓶颈分析

为了让规模庞大的平面波能在多计算节点间正确拆解并交由单节点执行 FFT，在跨越实空间（Real space Z-planes）与倒空间（Reciprocal space sticks）的过程中，产生了密集的全向通信（All-to-All）开销。本工作流致力于解除其中的阻塞瓶颈。

## 1. 痛点：阻塞式全局路由（`MPI_Alltoallv`）的性能灾难

### 旧算法表现与实现形式
在 `pw_gatherscatter.h` 里的 `PW_Basis::gatherp_scatters` 与 `PW_Basis::gathers_scatterp` 是平面波变换的核心通信路径：

1. **局部封包（Pack）**：将 `(nplane,fftxy)` 或 `(nz,nst)` 重新布局为连续的 `(nplane,nstot)` 或 `(numz[ip],nst)` 缓冲区。
2. **全局同步发收**：调用 `MPI_Alltoallv` 完成 stick/plane 的全局转置。
3. **局部解包（Unpack）**：将接收缓冲按 `startg/startr` 规则拆回 `(nz,nst)` 或 `(nplane,fftxy)`。

**核心问题在于“阻塞与同步”**。`MPI_Alltoallv` 会强制所有进程进入同步等待，必须等最慢链路或最慢进程完成后整体继续，导致算力空转。

## 2. 旧流水线：计算与通信的完全割裂

当前数据流是严格串行的：
```
[FFT 或局部拷贝] -> (MPI 全局同步) -> [拷贝/重排] -> [FFT]
```
随着核心数增加，这段 AlltoAllv 同步区会显著拖慢强弱扩展效率。

## 3. 面向计算/通信重叠（Overlap）的改进思路

本工作流旨在利用**非阻塞通信引擎**彻底推翻该行为：

**利用 `MPI_Isend` / `MPI_Irecv` 或 `MPI_Ialltoallv` 技术**：
- Pack 完成某个 `ip` 的数据块后立即发起异步通信。
- 在等待期间执行剩余 Pack 或本地重排。
- 在需要访问接收数据前统一 `MPI_Waitall`。

**补充细节（现有实现的关键参数）**：
- `numg/numr` 与 `startg/startr` 控制每个进程在 AlltoAllv 中的发送/接收 counts 与 offsets。
- 单进程路径 (`poolnproc == 1`) 直接走本地拷贝，避免通信，但拷贝仍是热区。
- 当前实现根据 `typeid(T)` 选择 `MPI_DOUBLE_COMPLEX` 或 `MPI_COMPLEX`，这对非阻塞路径也必须保留一致性。