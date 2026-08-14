# 抽象倒空间基类（K/Q 点统一）重构计划

- 状态：已批准（2026-08-14）
- 关联模块：`source/source_cell`、`source/source_pw/module_dfpt`
- 背景：DFPT 需要 q 点管理。q 点与 k 点大量共用倒空间网格/坐标/归约逻辑，
  但有两点关键差异：(1) q 点不涉及自旋；(2) q 点需要小群不可约表示分解。

## 目标

提取一个抽象的倒空间基类 `ModuleCell::ReciprocalGrid`，让 `K_Vectors` 与
`QList` 分别继承其共用功能，再各自加入不共用的功能：

- 基类：网格生成（Monkhorst-Pack）、坐标/权重、star 归约原语，spin-free。
- `K_Vectors`（继承基类）：自旋展开（isk、nspin 翻倍、SOC/磁群）、kstars、
  k 特有 MPI 分发。
- `QList`（继承基类）：q 网格生成、小群不可约表示（接口占位→完整实现）。

## 关键差异与决策（已与作者确认）

1. q 不含自旋：第一阶密度是标量，q 作为 `nspin=1` 的纯列表；自旋相关逻辑
   全部下沉到 `K_Vectors`，基类不触碰。
2. q 需要两级归约：
   - **star 归约**：与 k 完全相同（含 `q≡-q` 的 time-reversal 加倍，`kvec_ibz_kpoint`
     已自带）。
   - **小群不可约表示分解**：每个 q 的小群 `{R|t}: Rq≡q+G`，将 3N 原子位移分解
     到各 irrep 的代表模，只解代表模。仓库无现成实现，需新增。
3. 对称性模块已提供所需全部输入：`symm.kgmatrix[48]`（K 空间旋转）+
   `symm.gtrans[48]`（分数平移，供 `e^{-iq·t}` 相因子）。

## 类体系设计

```
ModuleCell::ReciprocalGrid（新抽象基类，spin-free，放 source_cell/）
 ├─ 数据(protected+getter): kvec_c/kvec_d/kvec_c_full, wk, ngk, nmp[3],
 │    kl_segids, kc_done/kd_done, nks/nkstot/nkstot_full, nspin=1, is_mp
 ├─ 方法: renew, Monkhorst_Pack(+formula), set_both_kvec, normalize_wk,
 │        print, reduce_ibz（star 归约原语）
 └─ 纯虚钩子: virtual void reduce_by_symmetry(...) = 0
        ↑                        ↑
 K_Vectors（public 继承）    QList（public 继承）
 ├─ isk, nspin, set_kup_and_kdw    ├─ generate_mesh / read_from_file
 ├─ kstars, ibz_index              ├─ reduce_by_symmetry = star 归约（复用原语）
 ├─ kvec_mpi_k（含 spin 广播）     └─ get_irreps = 小群 irrep 分解（新增）
 └─ reduce_by_symmetry 覆写        └─ nirr_ / irrep_modes_

ModuleSymmetry::LittleGroup（新独立组件，放 module_symmetry/）
 └─ 输入 kgmatrix+gtrans+q，输出小群操作/irrep 模式；QList 聚合它
```

要点：
- 基类不含自旋、不含 irrep；irrep 用独立组件类，避免把 q 特有复杂度带进基类。
- 基类显式传参、不新增 `GlobalV` 依赖（AGENTS.md 规则 1）。
- `K_Vectors` public 继承且不改外部 API（`kv.kvec_d` 等按名访问全部兼容）。

## 实施阶段

### Phase 1 — 提取抽象基类 + K_Vectors 迁移（行为不变）
1. 记录基线：`ctest` 现有 `klist_test`/`klist_test_para`/`parallel_kpoints_test`。
2. 新增 `source/source_cell/reciprocal_grid.h/.cpp`（命名空间 `ModuleCell`）：
   - 迁移 K_Vectors/KVectorUtils 中 spin-free 的网格生成、坐标转换、
     权重归一化、打印逻辑（逐行一致）。
   - `reduce_ibz(...)`：`kvec_ibz_kpoint` 的通用内核（restrict + MP k-lattice
     转换 + 旋转等价判断 + 权重计数 + `-q` 加倍）。
   - 纯虚 `reduce_by_symmetry(...)`。
3. `K_Vectors` 改为 `public ModuleCell::ReciprocalGrid`：
   - 私有保留：isk、koffset、k_kword、k_nkstot、kstars、ibz_index、para_k。
   - 覆写 `reduce_by_symmetry()`：构造 kgmatrix（含 nspin=4 磁群分支 +
     include_inv 加倍 + kstars），调基类 `reduce_ibz`，再执行 K 特有
     `update_use_ibz`（nspin 扩容）。
   - `set()` 改为调用覆写方法。
   - `KVectorUtils` 自由函数先保留为薄封装委托基类（保住
     `esolver_fp.cpp:178` 的 `set_after_vc` 与现有测试编译），随后测试迁移、
     删封装。
4. `source/source_cell/CMakeLists.txt` 接入新文件。
5. 回归：重构后重跑基线测试，输出必须一致；不一致立即回退对应文件。

### Phase 2 — QList 接入基类
- `class QList : public ModuleCell::ReciprocalGrid`。
- `generate_mesh`：基类 `Monkhorst_Pack` 建 q 网格 → `reduce_by_symmetry()`
  （恒加 `-q`，无磁群）→ 归约后补 `kvec_c`（笛卡尔坐标）→ 权重归一化
  → `use_irreps` 开关控制是否填充 `nirr_`/`irrep_modes_`。
- 保持 `get_nq/get_q/get_nirr/get_irrep_modes` 接口不变；删除 design-phase 桩注释。
- q 点补充功能（已完成）：
  - `read_from_file`：读 q 点文件（Gamma/Monkhorst-Pack 网格、Direct/Cartesian
    列表、Line_Direct/Line_Cartesian 插值路径），不做对称归约（接口无 symm）。
  - `print_qlists`：打印 q 点笛卡尔/直接坐标表（`Q-POINTS` 标签）。
  - `get_nirr/get_irrep_modes` 越界安全；`use_irreps=false` 时 irrep 数据为空。
  - `nkstot_full` = 归约前网格规模；`wk` 求和 = 1（网格/列表）或逐点 1（路径）。

### Phase 3 — 不可约表示接口（module_symmetry，先留接口）
- 新增 `source/source_cell/module_symmetry/little_group.h/.cpp`
  （命名空间 `ModuleSymmetry`），QList 聚合。
- 首版仅接口：`set_q(q,symm)`、`get_nirr()`（返回 1，全对称 A1）、
  `get_mode_basis(irrep)`、`get_little_group_ops()`。
- 完整 irrep 表/投影算符下一轮实现，配金刚石/闪锌矿 q=Γ/X/L 已知 irrep 表单测。

### Phase 4 — DFPT 接线（后续迭代）
- `DFPT_PW::init` 真正走 `generate_mesh`；`DFPT_PW_Data` 用
  `get_nirr`/`get_irrep_modes` 驱动逐 irrep SCF。

## 测试策略

- 新增 `source/source_cell/test/reciprocal_grid_test.cpp`：MP 生成、d/c 转换、
  权重归一化、`reduce_ibz` 在已知小群（fcc、金刚石）上的归约结果。
- 新增 `source/source_cell/test/qlist_test.cpp`：q 网格生成、star 归约（含 `-q`）、
  `get_nirr` 接口。
- 回归基准：现有 `klist_test`、`klist_test_para`、`parallel_kpoints_test`
  输出逐字节一致。

## 风险与边界

1. **K_Vectors 行为回归**：靠基线测试锁定；提取中任何逻辑漂移立即回退。
2. **命名空间**：基类放 `ModuleCell`，K_Vectors 保持全局命名空间继承
   （跨命名空间继承合法）；彻底统一命名空间列为后续清理项。
3. **`ngk` 归属**：PW_Basis 填充的每点平面波数，k/q 都需要，放基类。
4. **C++11 基线、LF 换行、新文件进 CMakeLists**：全程遵守 AGENTS.md。
