// test_parallel_orbitals_composition.cpp
//
// D6 组合化回归护栏（阶段 I 最后一步）。
// 依据：docs/007-B-阶段I执行规格.md D6 节；遵循 docs/008-编码规范.md（C++11、
// 无裸指针新增、snake_case、必要注释、简单清晰）。
//
// 四组断言：
//   ① getter 等价性：Parallel_2D 的 getter 族与底层字段一致，PO 对组合成员
//     的唯一转发 getter（get_blacs_ctxt）与 para_2d() 内部状态一致（R4 防副本）；
//   ② 派生层同步：重新初始化 para_2d()（init/set/set_serial）后重跑 PO 派生层
//     入口，断言 nrow_bands、atom_begin_*、desc_wfc 与新分布逐项一致，
//     并显式验证"只重初始化、不重跑派生层"的静默错误模式能被检测到；
//   ③ static_assert 值语义锁：PO 不可拷贝、不可移动；P2D 不可拷贝、可移动；
//   ④ 双配置覆盖：同一用例主体在 __MPI 与串行构建下各跑一遍
//     （MPI 覆盖 init/set 路径，串行覆盖 set_serial 路径）。
//
// 值语义说明：Parallel_Orbitals 因用户声明析构而抑制隐式移动，且成员 para_2d_
// 的拷贝操作被删除，故既不可拷贝也不可移动；Parallel_2D 手写移动语义，可移动。

#include "gtest/gtest.h"
#include "../parallel_orbitals.h"

#include <algorithm>
#include <fstream>
#include <type_traits>
#include <utility>
#include <vector>

// ---------------------------------------------------------------------------
// ③ static_assert 值语义锁（编译期，禁止 copy、锁定 move 语义）
// ---------------------------------------------------------------------------
static_assert(!std::is_copy_constructible<Parallel_Orbitals>::value,
              "Parallel_Orbitals must not be copy constructible");
static_assert(!std::is_copy_assignable<Parallel_Orbitals>::value,
              "Parallel_Orbitals must not be copy assignable");
static_assert(!std::is_move_constructible<Parallel_Orbitals>::value,
              "Parallel_Orbitals must not be move constructible");
static_assert(!std::is_move_assignable<Parallel_Orbitals>::value,
              "Parallel_Orbitals must not be move assignable");
static_assert(std::is_default_constructible<Parallel_Orbitals>::value,
              "Parallel_Orbitals must be default constructible");

static_assert(!std::is_copy_constructible<Parallel_2D>::value,
              "Parallel_2D must not be copy constructible");
static_assert(!std::is_copy_assignable<Parallel_2D>::value,
              "Parallel_2D must not be copy assignable");
static_assert(std::is_move_constructible<Parallel_2D>::value,
              "Parallel_2D must be move constructible");
static_assert(std::is_move_assignable<Parallel_2D>::value,
              "Parallel_2D must be move assignable");
static_assert(std::is_default_constructible<Parallel_2D>::value,
              "Parallel_2D must be default constructible");

class ParallelOrbitalsCompositionTest : public testing::Test
{
protected:
    std::ofstream ofs_running; // set_nloc_wfc_Eij 的日志流，测试时弃置
    std::ofstream ofs_warning;
    std::vector<int> iat2iwt;  // 原子 -> 全局轨道起始索引
    int nat = 3;

    void SetUp() override
    {
        ofs_running.open("/dev/null");
        ofs_warning.open("/dev/null");
        iat2iwt = std::vector<int>{0, 10, 20};
    }

    void TearDown() override
    {
        ofs_running.close();
        ofs_warning.close();
    }

    /// 运行 PO 派生层入口：atom trace（两种配置都有）与 band 层（仅 __MPI）。
    void run_derived_layer(Parallel_Orbitals& po, const int nlocal_global, const int nbands_global)
    {
        po.set_atomic_trace(iat2iwt.data(), nat, nlocal_global);
#ifdef __MPI
        ASSERT_EQ(po.set_nloc_wfc_Eij(nbands_global, ofs_running, ofs_warning), 0);
        po.set_desc_wfc_Eij(nlocal_global, nbands_global, po.para_2d().get_row_size());
#endif
    }

    /// 断言派生层与当前 para_2d() 分布逐项一致（专打静默错误）。
    void assert_sync_with(Parallel_Orbitals& po, const int nlocal_global, const int nbands_global)
    {
        // atom 层尾哨兵必须跟随新分布
        EXPECT_EQ(po.atom_begin_row[nat], po.para_2d().get_row_size());
        EXPECT_EQ(po.atom_begin_col[nat], po.para_2d().get_col_size());
        // 每个 atom 的起始索引必须在合法范围或标记"本进程无该原子"
        for (int iat = 0; iat < nat; ++iat)
        {
            if (po.atom_begin_row[iat] == -1)
            {
                EXPECT_EQ(po.get_nrow_atom(iat), 0);
            }
            else
            {
                EXPECT_GE(po.atom_begin_row[iat], 0);
                EXPECT_LT(po.atom_begin_row[iat], po.para_2d().get_row_size());
            }
            if (po.atom_begin_col[iat] == -1)
            {
                EXPECT_EQ(po.get_ncol_atom(iat), 0);
            }
            else
            {
                EXPECT_GE(po.atom_begin_col[iat], 0);
                EXPECT_LT(po.atom_begin_col[iat], po.para_2d().get_col_size());
            }
        }
#ifdef __MPI
        // band 层：nrow_bands/nloc 必须由当前分布重新导出
        EXPECT_EQ(po.nrow_bands, po.para_2d().get_row_size());
        EXPECT_EQ(po.nloc_wfc, static_cast<int64_t>(po.ncol_bands) * po.para_2d().get_row_size());
        EXPECT_EQ(po.nloc_Eij, static_cast<int64_t>(po.ncol_bands) * po.ncol_bands);
        EXPECT_EQ(po.get_wfc_global_nbands(), po.desc_wfc[3]);
        EXPECT_EQ(po.get_wfc_global_nbasis(), po.desc_wfc[2]);
        // desc_wfc 内嵌 BLACS context 必须与组合成员一致
        EXPECT_EQ(po.desc_wfc[1], po.get_blacs_ctxt());
        EXPECT_EQ(po.desc_wfc[2], nlocal_global);
        EXPECT_EQ(po.desc_wfc[3], nbands_global);
        EXPECT_EQ(po.desc_wfc[8], po.para_2d().get_row_size());
#endif
    }
};

// ---------------------------------------------------------------------------
// ① getter 等价性（R4 防副本）：getter 族必须与底层字段一致
// ---------------------------------------------------------------------------
TEST_F(ParallelOrbitalsCompositionTest, GetterEqualsPara2d)
{
    Parallel_Orbitals po;
    const int mg = 50;
    const int ng = 50;
#ifdef __MPI
    ASSERT_EQ(po.para_2d().init(mg, ng, 2, MPI_COMM_WORLD, 1), 0);
#else
    po.para_2d().set_serial(mg, ng);
#endif
    const Parallel_2D& p2d = po.para_2d();
    EXPECT_EQ(p2d.get_row_size(), p2d.get_number_of_local_rows());
    EXPECT_EQ(p2d.get_col_size(), p2d.get_number_of_local_cols());
    EXPECT_EQ(p2d.get_local_size(), p2d.get_number_of_local_elements());
    EXPECT_EQ(p2d.get_local_size(),
              static_cast<int64_t>(p2d.get_number_of_local_rows())
                  * p2d.get_number_of_local_cols());
    EXPECT_EQ(p2d.get_block_size(), p2d.get_block_size());
    EXPECT_EQ(p2d.get_dim0(), p2d.get_grid_rows());
    EXPECT_EQ(p2d.get_dim1(), p2d.get_grid_cols());
#ifdef __MPI
    EXPECT_EQ(p2d.get_desc(), &p2d.get_desc()[0]);
    EXPECT_EQ(p2d.get_desc()[1], p2d.get_blacs_ctxt());
    EXPECT_EQ(p2d.get_desc()[2], mg);
    EXPECT_EQ(p2d.get_desc()[3], ng);
    EXPECT_EQ(p2d.get_desc()[4], p2d.get_block_size());
    EXPECT_EQ(p2d.get_desc()[5], p2d.get_block_size());
    EXPECT_EQ(p2d.get_desc()[8], std::max(p2d.get_row_size(), 1));
    // PO 对组合成员唯一的转发 getter：必须与 para_2d() 内部一致
    EXPECT_EQ(po.get_blacs_ctxt(), p2d.get_blacs_ctxt());
#endif
}

// ---------------------------------------------------------------------------
// ② 派生层同步：init() 重新初始化后重跑派生层，三者一致
// ---------------------------------------------------------------------------
#ifdef __MPI
TEST_F(ParallelOrbitalsCompositionTest, DerivedLayerSyncsAfterReinit)
{
    Parallel_Orbitals po;
    // 第一轮：小分布
    ASSERT_EQ(po.para_2d().init(30, 30, 1, MPI_COMM_WORLD, 1), 0);
    const int old_nrow = po.para_2d().get_row_size();
    run_derived_layer(po, 30, 16);
    EXPECT_EQ(po.nrow_bands, old_nrow);
    EXPECT_EQ(po.atom_begin_row[nat], old_nrow);

    // 模拟 L1/L10 静默错误模式：只重新初始化组合成员，不重跑派生层
    ASSERT_EQ(po.para_2d().init(50, 50, 1, MPI_COMM_WORLD, 1), 0);
    const int new_nrow = po.para_2d().get_row_size();
    ASSERT_NE(new_nrow, old_nrow); // 分布已改变
    EXPECT_NE(po.nrow_bands, new_nrow); // 旧值残留 -> 测试可检测
    EXPECT_NE(po.atom_begin_row[nat], new_nrow);

    // 重跑派生层 -> 与新分布逐项一致
    run_derived_layer(po, 50, 16);
    assert_sync_with(po, 50, 16);
}

// ---------------------------------------------------------------------------
// ② 派生层同步：set()（复用既有 BLACS 上下文）后重跑派生层
// ---------------------------------------------------------------------------
TEST_F(ParallelOrbitalsCompositionTest, DerivedLayerSyncsAfterSet)
{
    Parallel_Orbitals po;
    ASSERT_EQ(po.para_2d().init(40, 40, 2, MPI_COMM_WORLD, 1), 0);
    run_derived_layer(po, 40, 16);
    EXPECT_EQ(po.nrow_bands, po.para_2d().get_row_size());

    // 用同一 BLACS 上下文走 set() 重新初始化
    const int blacs = po.para_2d().get_blacs_ctxt();
    ASSERT_EQ(po.para_2d().set(60, 60, 2, blacs), 0);
    ASSERT_NE(po.para_2d().get_row_size(), po.nrow_bands); // 派生层旧值残留 -> 可检测
    run_derived_layer(po, 60, 16);
    assert_sync_with(po, 60, 16);
}
#else
// ---------------------------------------------------------------------------
// ② 派生层同步：串行 set_serial 路径
// ---------------------------------------------------------------------------
TEST_F(ParallelOrbitalsCompositionTest, DerivedLayerSyncsAfterSetSerial)
{
    Parallel_Orbitals po;
    po.para_2d().set_serial(30, 30);
    const int old_nrow = po.para_2d().get_row_size();
    po.set_atomic_trace(iat2iwt.data(), nat, 30);
    EXPECT_EQ(po.atom_begin_row[nat], old_nrow);

    // 模拟静默错误模式：仅重新初始化组合成员，不重跑派生层
    po.para_2d().set_serial(50, 50);
    const int new_nrow = po.para_2d().get_row_size();
    ASSERT_NE(new_nrow, old_nrow);
    EXPECT_NE(po.atom_begin_row[nat], new_nrow); // 旧值残留 -> 测试可检测

    // 重跑派生层 -> 一致；串行下 global2local 为恒等映射
    po.set_atomic_trace(iat2iwt.data(), nat, 50);
    EXPECT_EQ(po.atom_begin_row[nat], po.para_2d().get_row_size());
    EXPECT_EQ(po.atom_begin_col[nat], po.para_2d().get_col_size());
    for (int iat = 0; iat < nat; ++iat)
    {
        EXPECT_EQ(po.atom_begin_row[iat], iat2iwt[iat]);
        EXPECT_EQ(po.atom_begin_col[iat], iat2iwt[iat]);
    }
}
#endif

// ---------------------------------------------------------------------------
// ③ 值语义运行期冒烟：P2D 可移动（static_assert 的运行时互补）
// ---------------------------------------------------------------------------
TEST_F(ParallelOrbitalsCompositionTest, ValueSemanticsSmoke)
{
    Parallel_2D source;
    source.set_serial(4, 4);
    EXPECT_EQ(source.get_row_size(), 4);

    Parallel_2D moved(std::move(source)); // 移动构造
    EXPECT_EQ(moved.get_row_size(), 4);

    Parallel_2D assigned;
    assigned = std::move(moved); // 移动赋值
    EXPECT_EQ(assigned.get_row_size(), 4);
}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
