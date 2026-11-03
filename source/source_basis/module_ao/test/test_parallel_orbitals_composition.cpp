// test_parallel_orbitals_composition.cpp
//
// Regression coverage for Parallel_Orbitals composition:
// - the composed BlockCyclicDistributionInformation is exposed only through para_2d();
// - derived atom and band metadata can be resynchronized after layout changes;
// - copy and move semantics remain explicit;
// - MPI and serial builds exercise their respective initialization paths.

#include "gtest/gtest.h"
#include "../parallel_orbitals.h"

#include <algorithm>
#include <fstream>
#include <type_traits>
#include <utility>
#include <vector>

static_assert(!std::is_convertible<const Parallel_Orbitals&, const BlockCyclicDistributionInformation&>::value,
              "Parallel_Orbitals must expose the distribution only through para_2d()");
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

static_assert(!std::is_copy_constructible<BlockCyclicDistributionInformation>::value,
              "BlockCyclicDistributionInformation must not be copy constructible");
static_assert(!std::is_copy_assignable<BlockCyclicDistributionInformation>::value,
              "BlockCyclicDistributionInformation must not be copy assignable");
static_assert(std::is_move_constructible<BlockCyclicDistributionInformation>::value,
              "BlockCyclicDistributionInformation must be move constructible");
static_assert(std::is_move_assignable<BlockCyclicDistributionInformation>::value,
              "BlockCyclicDistributionInformation must be move assignable");
static_assert(std::is_default_constructible<BlockCyclicDistributionInformation>::value,
              "BlockCyclicDistributionInformation must be default constructible");

class ParallelOrbitalsCompositionTest : public testing::Test
{
protected:
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    std::vector<int> iat2iwt{0, 10, 20};
    const int nat = 3;

    /// Refresh state derived from the current para_2d() layout.
    void run_derived_layer(Parallel_Orbitals& po, const int nlocal_global, const int nbands_global)
    {
        po.set_atomic_trace(iat2iwt.data(), nat, nlocal_global);
#ifdef __MPI
        ASSERT_EQ(po.set_nloc_wfc_Eij(nbands_global, ofs_running, ofs_warning), 0);
        po.set_desc_wfc_Eij(nlocal_global, nbands_global, po.para_2d().get_row_size());
#endif
    }

    /// Assert that derived state matches the current para_2d() layout.
    void assert_sync_with(const Parallel_Orbitals& po, const int nlocal_global, const int nbands_global)
    {
        // The trailing sentinels must follow the new distribution.
        EXPECT_EQ(po.atom_begin_row[nat], po.para_2d().get_row_size());
        EXPECT_EQ(po.atom_begin_col[nat], po.para_2d().get_col_size());
        // Each atom must either start locally or be marked absent on this rank.
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
        // Band metadata must be re-derived from the current distribution.
        EXPECT_EQ(po.nrow_bands, po.para_2d().get_row_size());
        EXPECT_EQ(po.nloc_wfc, static_cast<int64_t>(po.ncol_bands) * po.para_2d().get_row_size());
        EXPECT_EQ(po.nloc_Eij, static_cast<int64_t>(po.ncol_bands) * po.ncol_bands);
        EXPECT_EQ(po.get_wfc_global_nbands(), po.desc_wfc[3]);
        EXPECT_EQ(po.get_wfc_global_nbasis(), po.desc_wfc[2]);
        // The embedded descriptor context must match the composed layout.
        EXPECT_EQ(po.desc_wfc[1], po.get_blacs_ctxt());
        EXPECT_EQ(po.desc_wfc[2], nlocal_global);
        EXPECT_EQ(po.desc_wfc[3], nbands_global);
        EXPECT_EQ(po.desc_wfc[8], po.para_2d().get_row_size());
#endif
    }
};

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
    const BlockCyclicDistributionInformation& p2d = po.para_2d();
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

#ifdef __MPI
TEST_F(ParallelOrbitalsCompositionTest, DerivedLayerSyncsAfterReinit)
{
    Parallel_Orbitals po;
    // First use a small distribution.
    ASSERT_EQ(po.para_2d().init(30, 30, 1, MPI_COMM_WORLD, 1), 0);
    const int old_nrow = po.para_2d().get_row_size();
    run_derived_layer(po, 30, 16);
    EXPECT_EQ(po.nrow_bands, old_nrow);
    EXPECT_EQ(po.atom_begin_row[nat], old_nrow);

    // Reinitialize only the composed layout and leave derived state stale.
    ASSERT_EQ(po.para_2d().init(50, 50, 1, MPI_COMM_WORLD, 1), 0);
    const int new_nrow = po.para_2d().get_row_size();
    ASSERT_NE(new_nrow, old_nrow);
    EXPECT_NE(po.nrow_bands, new_nrow);
    EXPECT_NE(po.atom_begin_row[nat], new_nrow);

    // Refreshing derived state restores the invariant.
    run_derived_layer(po, 50, 16);
    assert_sync_with(po, 50, 16);
}

TEST_F(ParallelOrbitalsCompositionTest, DerivedLayerSyncsAfterSet)
{
    Parallel_Orbitals po;
    ASSERT_EQ(po.para_2d().init(40, 40, 2, MPI_COMM_WORLD, 1), 0);
    run_derived_layer(po, 40, 16);
    EXPECT_EQ(po.nrow_bands, po.para_2d().get_row_size());

    // Reinitialize through set() while reusing the existing BLACS context.
    const int blacs = po.para_2d().get_blacs_ctxt();
    ASSERT_EQ(po.para_2d().set(60, 60, 2, blacs), 0);
    ASSERT_NE(po.para_2d().get_row_size(), po.nrow_bands);
    run_derived_layer(po, 60, 16);
    assert_sync_with(po, 60, 16);
}
#else
TEST_F(ParallelOrbitalsCompositionTest, DerivedLayerSyncsAfterSetSerial)
{
    Parallel_Orbitals po;
    po.para_2d().set_serial(30, 30);
    const int old_nrow = po.para_2d().get_row_size();
    po.set_atomic_trace(iat2iwt.data(), nat, 30);
    EXPECT_EQ(po.atom_begin_row[nat], old_nrow);

    // Reinitialize only the composed layout and leave derived state stale.
    po.para_2d().set_serial(50, 50);
    const int new_nrow = po.para_2d().get_row_size();
    ASSERT_NE(new_nrow, old_nrow);
    EXPECT_NE(po.atom_begin_row[nat], new_nrow);

    // Refresh derived state; serial global-to-local maps are identities.
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

TEST_F(ParallelOrbitalsCompositionTest, ValueSemanticsSmoke)
{
    BlockCyclicDistributionInformation source;
    source.set_serial(4, 4);
    EXPECT_EQ(source.get_row_size(), 4);

    BlockCyclicDistributionInformation moved(std::move(source));
    EXPECT_EQ(moved.get_row_size(), 4);

    BlockCyclicDistributionInformation assigned;
    assigned = std::move(moved);
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
