# `source_base/parallel*` 模块关系文档

本文梳理 `source/source_base/` 下 `parallel*` 系列头文件/源文件的整体架构、相互依赖、
上层调用和已知技术债务，作为后续并行层更大重构的基线材料。

- 适用的文件族：`parallel_comm`、`parallel_reduce`、`parallel_common`、`parallel_global`、
  `parallel_grid`、`parallel_2d`、`parallel_device`、`parallel_cell`。
- 阅读对象：需要改动并行分发/通信、把全局 MPI 状态改成显式传参、或迁移到
  `source_base/module_parallel/` 新抽象（`Parallel::ParaWorld`）的开发者。
- 行号/引用数为撰写时快照，随代码演进可能变化；以仓库实际内容为准。

---

## 1. 总览

### 1.1 分层与依赖图

这些文件并非同一层，而是从“最底层通信子/归约原语”到“面向具体分布语义的封装”逐层叠加：

```
   ┌────────────────────────────────────────────────────────────────────────────┐
   │  上层业务模块: esolver / pw / lcao / io / hsolver / estate / hamilt / cell   │
   └────────────────────────────────────────────────────────────────────────────┘
                         │  include（见 1.3 节）
   ┌─────────────────────┼───────────────────────┬───────────────────┬──────────┐
   ▼                     ▼                       ▼                   ▼          ▼
┌──────────────┐  ┌──────────────┐   ┌──────────────────┐  ┌─────────────┐ ┌──────────────┐
│parallel_2d   │  │parallel_grid │   │parallel_global   │  │parallel_    │ │parallel_     │
│2D 块循环分布 │  │实空间 z 分布 │   │进程网格/域名划分 │  │reduce 归约  │ │device 设备通信│
└──────┬───────┘  └──────┬───────┘   └────────┬─────────┘  └──────┬──────┘ └──────┬───────┘
       │                 │                    │ include           │               │
       │                 │                    ▼                   │               ▼
       │                 │           ┌──────────────────┐         │      module_device/types.h
       │                 │           │ parallel_comm.h  │◄────────┘ (parallel_reduce.cpp)
       │                 │           │ 全局通信子       │
       │                 │           │ + MPICommGroup   │
       │                 │           └──────────────────┘
       │                 │                    ▲
       │                 │                    │ include (only .cpp)
       │                 └────────────────────┘
       │ include (.cpp)
       ▼
module_external/{blacs,scalapack}_connector.h

并行层内部还挂着两个“工具层”文件：
  · parallel_common  —— bcast_* 系列，.cpp 复用 parallel_reduce.h 中的 MPI_Type<T>
  · parallel_cell    —— CommunicationDomain，串行/并行统一的通信子封装（较新）
```

依赖方向（`A ──▶ B` 表示 A include / 依赖 B）：

```
parallel_global.h  ──▶ parallel_comm.h ──▶ mpi.h
parallel_comm.cpp  ──▶ parallel_global.h
parallel_reduce.h  ──▶ mpi.h
parallel_reduce.cpp──▶ parallel_comm.h
parallel_common.cpp──▶ parallel_reduce.h        (仅为复用 MPI_Type<T>)
parallel_grid.cpp  ──▶ parallel_comm.h , global_variable.h
parallel_2d.cpp    ──▶ module_external/blacs_connector.h , scalapack_connector.h
parallel_device.h  ──▶ module_device/types.h
```

要点：

- `parallel_comm.h` 是并行层的“公共底座”：它声明了 6 个全局 MPI 通信子和
  `MPICommGroup`，被 `parallel_global`、`parallel_reduce`、`parallel_grid` 等共同依赖。
- `parallel_reduce.h` 的模板 `MPI_Type<T>` 被 `parallel_common.cpp` 复用，是二者之间
  唯一的耦合点。
- `parallel_2d` 相对独立：它依赖 `module_external` 的 BLACS/ScaLAPACK 适配层，而不是
  `parallel_comm.h`（它自带 BLACS 上下文，不与 `POOL_WORLD` 等共享）。
- `parallel_device` 依赖 `module_device`，与其余文件只有“同属 `Parallel_Common` 命名空间”
  的名义关系（`parallel_common.h` 与 `parallel_device.h` 都定义 `namespace Parallel_Common`，
  但声明不同符号）。

### 1.2 `__MPI` 编译边界

`source/source_base/CMakeLists.txt` **无条件**把全部 `parallel*.cpp` 编进 `base` 目标，
MPI 的有无完全靠预处理宏区分（仅 `test_parallel/` 子目录由 `ENABLE_MPI` 控制）。因此
“串行构建”不是不编译文件，而是走宏保护下的空实现：

| 文件 | 串行（无 `__MPI`）行为 | 受 `__MPI` 保护的部分 |
| --- | --- | --- |
| `parallel_comm.h/.cpp` | 头文件为空（仅 include guard）；`.cpp` 整个文件被 `#if defined __MPI` 包裹，**不产生符号** | 全部内容（6 个通信子 + `MPICommGroup`） |
| `parallel_device.h/.cpp` | 头文件整体在 `#ifdef __MPI` 内，串行不产生符号 | 全部内容 |
| `parallel_2d.h/.cpp` | `class Parallel_2D` 始终存在，串行走 `set_serial()`；BLACS 描述符等为空 | `init()` / `set()` / `comm()` / `desc` / BLACS 私有成员 |
| `parallel_common.h/.cpp` | 声明始终存在，函数体为空 | 每个函数体的实现 |
| `parallel_global.h/.cpp` | 命名空间与函数声明存在；`mpi_number`/`omp_number` 全局量**仅在 `__MPI` 下定义** | `divide_pools`/`init_pools`/`divide_mpi_groups`/`finalize_mpi` 及两个全局量 |
| `parallel_grid.h/.cpp` | 类始终存在，`init()` 在设置完几何量后早返回 | `bcast()`/`reduce()`/`zpiece_distribute()`/`reduce_across_pools` 的 MPI 分支 |
| `parallel_reduce.h/.cpp` | 模板声明与显式实例化存在，函数体为空 | `MPI_Type<T>` 特化与全部函数体 |

> 注意：`parallel_global.cpp` 中 `int mpi_number; int omp_number;` 的定义被
> `#if defined __MPI` 包裹，而头文件里的 `extern` 声明**没有**包裹。这是当前一个
> 不对称点：串行构建下这两个符号只有声明、没有定义（尚未被引用，暂不报错）。

### 1.3 上层模块使用概览

引用文件数（含测试文件，粗略快照）与主要用途：

| 头文件 | 引用文件数 | 主要使用模块 | 典型用途 |
| --- | --- | --- | --- |
| `parallel_reduce.h` | ≈165 | pw·lcao·io·estate·hamilt·hsolver·cell | 跨进程/池归约（`reduce_all/pool/min/max`、`reduce_double_grid/diag`、`gather_int_all`） |
| `parallel_global.h` | ≈95 | main·cell·io·hsolver·psi·basis | MPI/OpenMP 初始化、进程域划分（`read_pal_param`、`init_pools`、`split_*_world`） |
| `parallel_comm.h` | ≈63 | pw·io·hsolver·estate·lcao·esolver | 直接使用 `POOL_WORLD`/`KP_WORLD`/`BP_WORLD`/`INT_BGROUP`/`GRID_WORLD`/`DIAG_WORLD` |
| `parallel_common.h` | ≈50 | cell·io·basis·estate·relax·pw | 基础类型广播（`bcast_int/double/complex/string`） |
| `parallel_device.h` | 24 | pw·hsolver·psi | 设备侧 `send/recv/bcast/reduce/gatherv`（含 CPU/GPU 适配） |
| `parallel_2d.h` | 23 | hsolver·lcao(rdmft/lr/bse)·io(module_hs)·basis | 2D block-cyclic 矩阵分布与 ScaLAPACK 描述符 |
| `parallel_grid.h` | 23 | io(chgpot/wf/output)·estate(charge)·esolver_fp | 实空间网格的跨 k-pool 广播/归约 |

按模块归类的高层画像：

- **esolver**：`esolver_fp` 用 `parallel_grid`；`esolver_sdft_pw` 用 `BP_WORLD`；
  `esolver_ks_pw_tddft` 用 `parallel_reduce`；`esolver_dp/nep` 用 `parallel_common`。
- **pw**：`parallel_reduce`（受力/应力/EXX/非局域，~50 处）、`parallel_device`（~12 处）、
  `parallel_comm`（EXX、stodft）；`stru_fac` 用 `parallel_grid`。
- **lcao**：`parallel_2d`（rdmft/lr/bse）、`parallel_reduce`（deepks/dftu/operator/rt）、
  `parallel_comm`（`DIAG_WORLD`、`POOL_WORLD`）。
- **io**：`parallel_grid`（电荷/波函数/cube）、`parallel_2d`（module_hs 稠密 HS）、
  `parallel_comm`、`parallel_reduce` 均大量使用。
- **hsolver**：`parallel_2d`（`parallel_k2d`、`diag_hs_para`）、`parallel_comm`
  （`diago_dav*` 用 `POOL_WORLD`）、`parallel_device`（迭代/正交化的设备通信）。
- **cell**：`parallel_common`（`unitcell`、`atom_*`、`klist`）、`parallel_global`
  （`parallel_kpoints`）、`parallel_reduce`（`magnetism`）。

---

## 2. 各模块职责

### parallel_2d（`parallel_2d.h` 206 行 / `parallel_2d.cpp` 228 行）

- **核心类**：`Parallel_2D`。输入全局维度 `(mg, ng)`、块大小 `nb` 和一个 `MPI_Comm`，
  产出一个 ScaLAPACK 风格的 2D block-cyclic 分布。
- **功能**：
  - `init()`：按 `n = p*q`（最接近且 `p<=q`）分解进程网格，建 BLACS 网格（`Csys2blacs_handle`
    + `Cblacs_gridinit`），并用 `numroc_`/`descinit_` 生成局部尺寸与描述符；`mode=true` 时交换
    `dim0/dim1`。
  - `set()`：复用外部给定的 BLACS 上下文（借用，不拥有）。
  - `set_serial()`：串行布局（`nb=1, dim=1`，映射为恒等）。
  - 索引映射：`global2local_row/col`、`local2global_row/col`；查询 `in_this_processor`、
    `owner_processor`、`blacs_in_this_processor`。
  - RAII：析构与移动赋值通过 `release_blacs_grid()` 释放自有的 BLACS 上下文
    （`owns_blacs_ctxt_` 区分“拥有”与“借用”）。
- **依赖**：`source_base/module_external/blacs_connector.h`、`scalapack_connector.h`（仅 `.cpp`）。
- **被谁用**：hsolver（`parallel_k2d`、`diag_hs_para`）、lcao（rdmft / lr / bse）、
  io（`module_hs` 稠密 HS 读写、`write_wfc_nao`）、basis（`module_ao/parallel_orbitals`）。
- **已知细节**：`nrow/ncol/nloc/nb/dim0/dim1` 为 public（源码中有 `FIXME` 注释），
  目前被广泛直接访问，暂未收紧。

### parallel_common（`parallel_common.h` 30 行 / `parallel_common.cpp` 111 行）

- **核心函数**：`namespace Parallel_Common` 下的 `bcast_*` 系列。
  - 数组版：`bcast_complex_double/double/int(ptr, n)`、`bcast_string(ptr, n)`、`bcast_char(ptr, n)`。
  - 标量版：`bcast_complex_double/string/double/int(bool&)`。
- **功能**：把基础类型从 rank 0 广播到所有进程（**固定 `MPI_COMM_WORLD`**，固定 root=0）。
  `bcast_string` 先广播长度再广播内容；`bcast_bool` 走 `int` 中转。
- **依赖**：`parallel_reduce.h`（仅 `.cpp`，复用模板 `Parallel_Reduce::MPI_Type<T>`）。
  **注意**：虽与 `parallel_comm.h` 同属并行层，但二者无直接 include 关系。
- **被谁用**：全局变量初始化、输入参数广播。cell（`unitcell`、`atom_*`、`klist`）、
  io（`module_qo`、`module_wannier`、波函数读写）、basis（`module_nao`）、
  estate（`module_charge`）、relax、esolver（dp/nep）等。
- **接口限制**：所有 `bcast_*` 都硬编码 `MPI_COMM_WORLD`，无法指定通信子/root；需要
  在子域内广播时必须改用 `parallel_device.h` 的 `bcast_data(..., comm, root)`。

### parallel_global（`parallel_global.h` 91 行 / `parallel_global.cpp` 336 行）

- **核心**：`namespace Parallel_Global`，是 ABACUS 并行的“启动与拓扑”入口。
- **功能**：
  - `read_pal_param()`：`MPI_Init`/`MPI_Init_thread`、取 `NPROC`/`MY_RANK`、按节点内进程数
    自动设定 OpenMP 线程数并做超订警告。参数（`NPROC/NTHREAD_PER_PROC/MY_RANK`）通过
    引用返回，**不直接写全局量**（用户后续再赋给 `GlobalV`）。
  - `divide_mpi_groups()`：把 `procs` 个进程均分为 `num_groups` 组，返回组内进程数与
    本进程的 `my_group/rank_in_group`；`even=true` 时强制整除。
  - `divide_pools()` / `init_pools()`：按 `KPAR`、`BNDPAR` 划分 k 池与带组，产出
    `POOL_WORLD`/`KP_WORLD`/`INT_BGROUP`/`BP_WORLD` 及对应的 `NPROC_IN_POOL`/`RANK_IN_POOL` 等
    兼容变量（借助 `parallel_comm.h` 的 `MPICommGroup`）。
  - `split_diag_world()` / `split_grid_world()`：由 `DIAGO_PROC` 切出 `DIAG_WORLD`/`GRID_WORLD`。
  - `finalize_mpi()`：释放上述通信子并 `MPI_Finalize()`。
  - 全局量 `mpi_number`（节点内进程数）、`omp_number`（线程数）。
- **依赖**：`parallel_comm.h`；`.cpp` 另用 `global_function.h`/`global_variable.h`/`tool_quit.h`。
- **被谁用**：`source_main/driver.cpp`（启动序列）、`cell`（k 点并行）、
  `hsolver/parallel_k2d`、io、psi、以及大量单测的 `main`。
- **函数调用点**：
  - `init_pools` → `driver.cpp`、多个测试；
  - `split_diag_world`/`split_grid_world` → `driver.cpp`；
  - `divide_mpi_groups` → `parallel_comm.cpp`（`MPICommGroup` 内部）、`parallel_k2d.cpp` 及测试。

### parallel_grid（`parallel_grid.h` 71 行 / `parallel_grid.cpp` 430 行）

- **核心类**：`Parallel_Grid`。管理实空间 FFT 网格 `(ncx, ncy, ncz)` 在 **z 方向** 按
  k-pool 与池内进程的二维切分。
- **功能**：
  - `init()`：记录网格几何量，按 `nprocgroup` 与 `KPAR` 计算每个 pool 的进程数、
    每个进程拥有的 z 层数（`numz`）、起始层（`startz`）、以及“第 iz 层归哪个全局/池内进程”
    （`whichpro`/`whichpro_loc`）。
  - `z_distribution()`：生成上述分布表（round-robin 按 `iz % nproc` 累加 `bz`）。
  - `bcast()` + `zpiece_distribute()`：从 root 把一个又一个 xy 平面（z-slice）分发到
    各 pool 的属主进程（两种变体：`is_sdft` 决定用 `INT_BGROUP` 还是 `MPI_COMM_WORLD`）。
  - `reduce()`：把各进程的局部 slab 通过 `POOL_WORLD` `MPI_Gatherv` 汇总到池 root，
    并重排为 Cube 输出期望的 `[xy][global_z]` 顺序。
  - `reduce_across_pools()`：跨 k-pool 求和。等大 pool 直接 `MPI_Allreduce(KP_WORLD)`；
    不齐 pool 则先 `POOL_WORLD` allgather 重建公共布局，再 `INT_BGROUP` 求和。
- **依赖**：`.cpp` 用 `global_variable.h`（`GlobalV::` 系列）和 `parallel_comm.h`。
- **被谁用**：estate（`module_charge`）、io（电荷/波函数/cube 输出）、`esolver_fp`、
  pw（`stru_fac`）、lcao（`xc_kernel`）。

### parallel_reduce（`parallel_reduce.h` 87 行 / `parallel_reduce.cpp` 180 行）

- **核心**：`namespace Parallel_Reduce`，全仓库最广使用的并行封装。
- **功能**（均为 `MPI_Allreduce` 封装）：
  - `reduce_all`：`MPI_COMM_WORLD` 上求和（标量/数组，模板 + 显式实例化）。
  - `reduce_pool`：`POOL_WORLD` 上求和。
  - `reduce_min`/`reduce_max`：全局极值；`reduce_max(double*, n)`、`reduce_min_pool`/`reduce_max_pool`。
  - `reduce_double_grid` → `GRID_WORLD`；`reduce_double_diag` → `DIAG_WORLD`。
  - `reduce_double_allpool`：跨 pool 的“先除池内进程数再求和”归约。
  - `gather_int_all`：`MPI_Allgather`。
  - `MPI_Type<T>`：`int/double/float/complex<double>/complex<float>/long long` → MPI 类型的映射。
- **依赖**：`.cpp` 用 `parallel_comm.h`（`POOL_WORLD`/`GRID_WORLD`/`DIAG_WORLD`）。
- **被谁用**：几乎全部物理模块（pw/lcao/io/estate/hamilt/hsolver/cell/psi/…）。由于模板
  需要显式实例化，新增支持类型时必须同步改 `.h`（特化）与 `.cpp`（实例化）。

### parallel_device（`parallel_device.h` 217 行 / `parallel_device.cpp` 439 行）

- **核心**：`namespace Parallel_Common` 下的设备侧通信原语（与 `parallel_common` 同名命名空间
  但符号独立）。
- **功能**：
  - 基础 MPI 包装：`isend_data` / `send_data` / `recv_data` / `bcast_data` / `reduce_data` /
    `gatherv_data`，按 `double/std::complex<double>/float/std::complex<float>` 重载。
  - 设备模板：`send_dev`/`isend_dev`/`recv_dev`/`bcast_dev`/`reduce_dev`/`gatherv_dev`，通过
    `object_cpu_point<T,Device>` 在 CPU（直通）与 GPU（d2h/h2d 拷贝）之间适配；
    `#ifdef __CUDA_MPI` 时直接走 MPI，否则经 CPU 中转。
  - NCCL 路径：`#if defined(__NCCL_PARALLEL_DEVICE)` 时提供 `nccl_bcast_data`/`nccl_reduce_data`/
    `nccl_gatherv_data`，用 `NcclCommRegistry`（`MPI_Comm → ncclComm_t` 的懒创建缓存，带互斥锁）
    管理 NCCL 通信子，`gatherv` 以 `ncclAllGather` 分块模拟。
- **依赖**：`module_device/types.h`（设备标签 `DEVICE_CPU`/`DEVICE_GPU`）、`.cpp` 另用
  `memory_op.h`；NCCL 分支用 `device_check.h` + CUDA runtime。
- **被谁用**：pw（受力/应力/非局域/EXX）、hsolver（迭代与正交化）、psi（`psi_prepare`）、
  `para_gemm`。
- **结构问题**：NCCL 实现（`.cpp` 约 21–254 行）与纯 MPI 包装（257–352 行）及
  `object_cpu_point` 特化（354–436 行）混在同一文件，见 §4。

### 相关但不在核心 7 文件内

- **`parallel_comm.h/.cpp`**：并行层的公共底座，声明/定义 6 个全局通信子与 `MPICommGroup`
  （`divide_group_comm` 内部调用 `Parallel_Global::divide_mpi_groups`）。
- **`parallel_cell.h/.cpp`**：`ModuleBase::CommunicationDomain`（+ `world_comm_domain()`），
  一个串行/并行统一的通信子值对象（`rank/size/max`），是 `rhog_io` 等迁移到“显式传通信域”
  时引入的轻量封装，可视为 `module_parallel` 思路的过渡产物。

---

## 3. 关键全局变量

并行层通过两类“全局状态”耦合上层：`GlobalV` 里的进程/池索引，以及 `parallel_comm.h`
里的裸 `MPI_Comm`。

### 3.1 `GlobalV` 并行变量

| 变量 | 引用文件数 | 主要分布（模块:文件数） | 语义 |
| --- | --- | --- | --- |
| `GlobalV::MY_RANK` | 143 | io:54, cell:25, basis:17, lcao:16, pw:8, base:8, esolver:6 | 全局 MPI rank |
| `GlobalV::NPROC` | 55 | io:25, cell:12, esolver:5, basis:5 | 全局进程数 |
| `GlobalV::NPROC_IN_POOL` | 36 | pw:10, io:9, estate:7, lcao:4 | 每个 k 池进程数 |
| `GlobalV::RANK_IN_POOL` | 32 | io:10, pw:5, lcao:5, estate:3 | 池内 rank |
| `GlobalV::KPAR` | 31 | io:11, estate:10, pw:4, cell:2, base:2 | k 点并行组数 |
| `GlobalV::MY_POOL` | 21 | io:11, lcao:2, estate:2, cell:2 | 本进程所属 k 池 |
| `GlobalV::MY_BNDGROUP` | 13 | io:5, pw:2, cell:2 | 本进程所属带组 |

- `parallel_grid.cpp` 内部直接读取 `GlobalV::KPAR / MY_POOL / RANK_IN_POOL / RANK_IN_BPGROUP /
  MY_RANK / ofs_warning`；`parallel_global` 的划分函数则通过**引用参数**返回结果，由调用方
  （`driver.cpp`）写入 `GlobalV`。
- 这些变量是跨层控制的典型来源：`parallel_reduce`、`parallel_grid` 等都隐式依赖调用方
  已经正确初始化了 `GlobalV`，函数签名无法体现这一前提。

### 3.2 `parallel_comm.h` 的全局通信子

| 通信子 | 引用文件数 | 语义 |
| --- | --- | --- |
| `POOL_WORLD` | ≈139 | 同 k-point、同 band、不同平面波（池内） |
| `BP_WORLD` | 28 | 同 k、同平面波、不同 band（stodft/bpcg） |
| `DIAG_WORLD` | 22 | 各 grid 组的首个进程组成的对角化域 |
| `KP_WORLD` | 17 | 不同 k、同 band、同平面波（跨池，pool 不齐时为 `MPI_COMM_NULL`） |
| `INT_BGROUP` | 17 | 同 band、不同 k、不同平面波（带组内部） |
| `GRID_WORLD` | 12 | 由 `DIAGO_PROC` 切出的实空间网格域 |

这些变量在 `parallel_comm.cpp` 中定义，模块通过 `#include "source_base/parallel_comm.h"`
直接访问，属于“应显式传参”的已知债务（见 §4）。

### 3.3 应当改为传参的候选

- `parallel_grid.cpp` 中的全部 `GlobalV::` 读取（`KPAR`/`MY_POOL`/`RANK_IN_POOL`/…）。
- `parallel_common.cpp` 中硬编码的 `MPI_COMM_WORLD`（应允许传入 `comm`）。
- 裸通信子 `POOL_WORLD` 等（迁移目标为 `Parallel::ParaWorld`，见 §5）。

---

## 4. 已知技术债务

1. **`parallel_grid.cpp::init` 的 `#ifndef __MPI` 分支** —— 已解决（任务 1）。
   原 `init()` 在串行时于中途 `return`，把“几何量设置”与“MPI 分布”耦合在一个函数里；
   现已拆为 `init_serial()`（通用几何量）+ `init_parallel()`（`__MPI` 专属分布），
   `init()` 只做调用编排。
   分支 `refactor-parallel-grid-init`，提交 `94c7126f2`（*Refactor Parallel_Grid::init into
   serial/parallel subroutines*）。
   > 说明：该修复当前在独立分支上，未合入本文所基于的 `develop`；本仓库主线仍是旧写法。
2. **`z_distribution` 依赖 `GlobalV::KPAR`** —— 已解决（任务 2）。
   原 `z_distribution()` 直接读 `GlobalV::KPAR`，行为依赖全局状态；现改为显式形参 `kpar`，
   由 `init()` 传入 `GlobalV::KPAR`，函数体不再读任何全局变量（新增 `assert(kpar > 0)`）。
   分支 `refactor-parallel-grid-kpar`，提交 `a5da988f7`。
3. **`parallel_global.cpp` 偏大（336 行）**，混合了四类职责：
   ①进程/线程初始化（`read_pal_param`）、②通用分组算法（`divide_mpi_groups`）、
   ③k/带池划分（`divide_pools`/`init_pools`）、④对角化/网格域切分（`split_*_world`）
   与⑤收尾（`finalize_mpi`）。可按上述职责拆成多个 `.cpp`，降低单文件复杂度。
4. **`parallel_device.cpp` 中 NCCL 与 MPI 代码混在一起**：NCCL 注册表与实现、纯 MPI
   重载包装、CPU 回退特化三段风格迥异，建议按后端拆分（如 `parallel_device_nccl.cpp`）。
5. **`parallel_common.cpp` 硬编码 `MPI_COMM_WORLD`**：`bcast_*` 无法指定通信子/root，
   在子域内广播只能绕行 `parallel_device.h`，接口不对称。
6. **裸通信子全局量**：`POOL_WORLD` 等被 60+ 文件直接引用，是跨层耦合的最大来源；
   `module_parallel` 正在用 `ParaWorld` 逐步替换（§5）。
7. **`parallel_global.cpp` 全局量的宏不对称**：`mpi_number`/`omp_number` 的定义在
   `#if defined __MPI` 内，而头文件 `extern` 声明在宏外（见 §1.2 注）。
8. **`parallel_grid.cpp::init` 仍读取 `GlobalV`**：任务 1 只解决了“串行分支耦合”，
   `init_parallel()`/`reduce_across_pools()` 等仍直接依赖 `GlobalV::KPAR/MY_POOL/
   RANK_IN_POOL/RANK_IN_BPGROUP/MY_RANK`，是下一步参数化的候选。

---

## 5. 演进方向：`source_base/module_parallel/`

仓库正在引入一套新的并行抽象，用来替代 `parallel_global` 的函数与 `parallel_comm` 的裸通信子：

- **`Parallel::ParaWorld`**（`para_world.h`）：把“域名标签 + 通信子 + rank + size”封成一个值对象，
  串行构建下 `rank()==0`、`size()==1`，使调用点无需 `#ifdef __MPI`。
- **`Parallel::ParaCollection`**：域的集合/查找容器。
- **`Parallel::setup_para_worlds()`**（`para_setup.h`）：一次构建完整的域层次树，显式替换
  `divide_pools` + `split_diag_world` + `split_grid_world`；`divide_mpi_groups`、
  `split_images`/`split_pools`/`split_diag_world`/`split_grid_world` 都有对应的
  `ParaWorld` 版本。
- **`para_bridge.h`**：临时桥接层（`make_pw_world`/`make_band_world`），从旧的
  `POOL_WORLD`/`BP_WORLD` 构造 `ParaWorld`，待 driver 初始化接管后删除。
- 各域已拆分为 `para_kmesh_world`/`para_pw_world`/`para_bgroup_world`/`para_diag_world`/
  `para_rgrid_world`/`para_matrix_world` 等，并有对应单测。

因此，本文档记录的“旧并行层”目标是逐步被 `ParaWorld` 收编：新代码应优先接收上层传入的
`ParaWorld&`，而不是新增对 `GlobalV::*`/`POOL_WORLD` 的依赖（对应 AGENTS.md 第 1 条）。

---

## 6. 测试与验证入口

- 单元测试目录：`source/source_base/test_parallel/`
  （`parallel_2d_test.cpp`、`parallel_common_test.cpp`、`parallel_device_test.cpp`、
  `parallel_domain_grid_test.cpp`、`parallel_global_test.cpp`、`parallel_reduce_test.cpp`）。
- 新抽象测试：`source/source_base/module_parallel/test/`。
- 这些测试目录在 `source/source_base/CMakeLists.txt` 中由 `ENABLE_MPI` 控制；
  运行示例：`ctest --test-dir build -V -R MODULE_BASE`（并按需设 `OMP_NUM_THREADS=1`）。
- 改动 `GlobalV` 依赖或通信子时，建议至少覆盖：串行构建（无 `__MPI`）编译、MPI 构建下单测、
  以及一个多进程集成算例。

---

## 附录：文件清单与规模（快照）

| 文件 | 行数 | 角色 |
| --- | --- | --- |
| `parallel_comm.h/.cpp` | 44 / 54 | 全局通信子 + `MPICommGroup`（整个 `.cpp` 仅 `__MPI`） |
| `parallel_reduce.h/.cpp` | 87 / 180 | reduce/gather 封装 + `MPI_Type<T>` |
| `parallel_common.h/.cpp` | 30 / 111 | `bcast_*`（`MPI_COMM_WORLD`） |
| `parallel_global.h/.cpp` | 91 / 336 | 启动、分组、域切分、收尾 |
| `parallel_grid.h/.cpp` | 71 / 430 | 实空间网格 z 方向分布与跨池归约 |
| `parallel_2d.h/.cpp` | 206 / 228 | 2D block-cyclic / BLACS / ScaLAPACK |
| `parallel_device.h/.cpp` | 217 / 439 | 设备侧通信 + 可选 NCCL |
| `parallel_cell.h/.cpp` | 33 / 69 | `CommunicationDomain`（串行/并行统一通信域） |
