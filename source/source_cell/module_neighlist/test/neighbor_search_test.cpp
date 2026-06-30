#include <gtest/gtest.h>
#include "../neighbor_search.h"
#include "../unitcell_lite.h"
#include <algorithm>

// Helper function to create a simple UnitCellLite for testing
static UnitCellLite make_test_ucell(double lat0, double omega,
                                     const ModuleBase::Matrix3& latvec,
                                     int ntype, const std::vector<int>& na,
                                     const std::vector<ModuleBase::Vector3<double>>& tau) {
    UnitCellLite ucell;
    ucell.set_lattice(lat0, omega, latvec);
    ucell.set_atoms(ntype, na, tau);
    return ucell;
}

TEST(NeighborSearchTest, TwoAtomsNeighbor)
{
    ModuleBase::Matrix3 latvec;
    latvec.e11 = 1; latvec.e12 = 0; latvec.e13 = 0;
    latvec.e21 = 0; latvec.e22 = 1; latvec.e23 = 0;
    latvec.e31 = 0; latvec.e32 = 0; latvec.e33 = 1;

    UnitCellLite ucell = make_test_ucell(
        1.0, 1.0, latvec, 1, {2},
        {{0.0, 0.0, 0.0}, {0.5, 0.0, 0.0}}
    );

    NeighborSearch ns;
    double cutoff = 1.0;

    ns.init(ucell, cutoff, 0);
    ns.build_neighbors();

    auto &list = ns.get_neighbor_list();

    ASSERT_EQ(list.get_nlocal(), 2);

    EXPECT_EQ(list.get_numneigh(0), 8);
    EXPECT_EQ(list.get_numneigh(1), 8);
}

TEST(NeighborSearchTest, NoNeighbor)
{
    ModuleBase::Matrix3 latvec;
    latvec.e11 = 1; latvec.e12 = 0; latvec.e13 = 0;
    latvec.e21 = 0; latvec.e22 = 1; latvec.e23 = 0;
    latvec.e31 = 0; latvec.e32 = 0; latvec.e33 = 1;

    UnitCellLite ucell = make_test_ucell(
        1.0, 1.0, latvec, 1, {2},
        {{0.0, 0.0, 0.0}, {5.0, 0.0, 0.0}}
    );

    NeighborSearch ns;

    // use a smaller search radius to avoid counting periodic-image neighbors
    ns.init(ucell, 0.1, 0);
    ns.build_neighbors();

    auto &list = ns.get_neighbor_list();

    EXPECT_EQ(list.get_numneigh(0), 0);
    EXPECT_EQ(list.get_numneigh(1), 0);
}

TEST(NeighborSearchTest, DistributedInputUsesOwnedCentersAndGhostNeighbors)
{
    std::vector<LocalAtom> owned_atoms;
    std::vector<LocalAtom> ghost_atoms;
    owned_atoms.push_back(LocalAtom(ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                    ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                    0,
                                    0,
                                    0,
                                    0,
                                    false));
    ghost_atoms.push_back(LocalAtom(ModuleBase::Vector3<double>(0.5, 0.0, 0.0),
                                    ModuleBase::Vector3<double>(0.5, 0.0, 0.0),
                                    0,
                                    1,
                                    1,
                                    1,
                                    true));

    NeighborSearch ns;
    ns.init_distributed(owned_atoms, ghost_atoms, 1.0, 1.0);
    ns.build_neighbors();

    const NeighborList& list = ns.get_neighbor_list();
    ASSERT_EQ(list.get_nlocal(), 1);
    ASSERT_EQ(list.get_numneigh(0), 1);
    const int neighbor_id = list.get_firstneigh(0)[0];
    ASSERT_GE(neighbor_id, 0);
    ASSERT_LT(neighbor_id, static_cast<int>(ns.get_all_atoms().size()));
    EXPECT_EQ(ns.get_all_atoms()[neighbor_id].global_id, 1);
    EXPECT_EQ(ns.get_all_atoms()[neighbor_id].owner_rank, 1);
}

TEST(NeighborSearchUnit, DistanceBox)
{
    NeighborSearch ns;
    // set a single cell region at x=0..1,y=0..1,z=0..1
    ns.set_position(0, 0, 0);
    ns.set_width(1.0, 1.0, 1.0);

    double inside = ns.distance(0.2, 0.5, 0.5, 0.0, 0.0, 0.0);
    EXPECT_DOUBLE_EQ(inside, 0.0);

    double outside = ns.distance(2.0, 0.5, 0.5, 0.0, 0.0, 0.0);
    // squared distance should be (2-1)^2 = 1
    EXPECT_DOUBLE_EQ(outside, 1.0);
}

TEST(NeighborSearchUnit, DecomposeCases)
{
    NeighborSearch ns;
    int nx, ny, nz;

    ns.decompose(8, nx, ny, nz);
    EXPECT_EQ(nx * ny * nz, 8);
    // expect somewhat balanced cube factors for 8
    EXPECT_EQ(nx, 2);
    EXPECT_EQ(ny, 2);
    EXPECT_EQ(nz, 2);

    ns.decompose(7, nx, ny, nz);
    EXPECT_EQ(nx * ny * nz, 7);
    EXPECT_EQ(nx, 1);
    EXPECT_EQ(ny, 1);
    EXPECT_EQ(nz, 7);
}

TEST(NeighborSearchUnit, DecomposePrimeNumber)
{
    NeighborSearch ns;
    int nx, ny, nz;
    ns.decompose(13, nx, ny, nz);
    EXPECT_EQ(nx * ny * nz, 13);
    EXPECT_EQ(nx, 1);
    EXPECT_EQ(ny, 1);
    EXPECT_EQ(nz, 13);
}

TEST(NeighborSearchUnit, DecomposeSkipsZeroSpanDirections)
{
    NeighborSearch ns;
    int nx, ny, nz;

    ns.decompose(8, 1.0, 1.0, 0.0, nx, ny, nz);
    EXPECT_EQ(nx * ny * nz, 8);
    EXPECT_EQ(nz, 1);

    ns.decompose(4, 4.0, 0.0, 0.0, nx, ny, nz);
    EXPECT_EQ(nx, 4);
    EXPECT_EQ(ny, 1);
    EXPECT_EQ(nz, 1);

    ns.decompose(4, 0.0, 0.0, 0.0, nx, ny, nz);
    EXPECT_EQ(nx, 1);
    EXPECT_EQ(ny, 1);
    EXPECT_EQ(nz, 1);
}

TEST(NeighborSearchUnit, NonOrthogonalLatticeExpand)
{
    ModuleBase::Matrix3 latvec;
    // skewed lattice
    latvec.e11 = 1; latvec.e12 = 0.3; latvec.e13 = 0.0;
    latvec.e21 = 0.1; latvec.e22 = 1.0; latvec.e23 = 0.0;
    latvec.e31 = 0.0; latvec.e32 = 0.0; latvec.e33 = 1.0;

    UnitCellLite ucell = make_test_ucell(
        1.0, 1.0, latvec, 1, {1},
        {{0.0, 0.0, 0.0}}
    );

    NeighborSearch ns;
    ns.init(ucell, 2.5, 0);
    // for skewed lattice, expansion layers should be >= 1
    EXPECT_GE(ns.get_glayerX(), 1);
    EXPECT_GE(ns.get_glayerY(), 1);
    EXPECT_GE(ns.get_glayerZ(), 1);
}

TEST(NeighborSearchInit_WideZero_CentralInside, SingleAtomCell)
{
    ModuleBase::Matrix3 latvec;
    latvec.e11 = 1; latvec.e12 = 0; latvec.e13 = 0;
    latvec.e21 = 0; latvec.e22 = 1; latvec.e23 = 0;
    latvec.e31 = 0; latvec.e32 = 0; latvec.e33 = 1;

    UnitCellLite ucell = make_test_ucell(
        1.0, 1.0, latvec, 1, {1},
        {{0.0, 0.0, 0.0}}
    );

    NeighborSearch ns;
    // choose sr small enough; with mpi_size fixed to 1 in init, wide_* become 0
    ns.init(ucell, 0.1, 0);
    // central cell atom should be counted as inside
    EXPECT_EQ(ns.get_inside_atoms().size(), 1);
    EXPECT_EQ(ns.get_neighbor_list().get_nlocal(), static_cast<int>(ns.get_inside_atoms().size()));
}

TEST(NeighborSearchInit_MpiRankIndexing, RankValues)
{
    ModuleBase::Matrix3 latvec;
    latvec.e11 = 1; latvec.e12 = 0; latvec.e13 = 0;
    latvec.e21 = 0; latvec.e22 = 1; latvec.e23 = 0;
    latvec.e31 = 0; latvec.e32 = 0; latvec.e33 = 1;

    UnitCellLite ucell = make_test_ucell(
        1.0, 1.0, latvec, 1, {1},
        {{0.0, 0.0, 0.0}}
    );

    NeighborSearch ns0;
    ns0.init(ucell, 0.5, 0);
    // with mpi_size fixed to 1 in init, nx=ny=nz=1; for rank 0 expect x=y=0,z=0
    EXPECT_EQ(ns0.get_x(), 0);
    EXPECT_EQ(ns0.get_y(), 0);
    EXPECT_EQ(ns0.get_z(), 0);
}

TEST(NeighborSearchInit_MpiOwnership, SingleAtomZeroSpanIsOwnedOnce)
{
    ModuleBase::Matrix3 latvec;
    latvec.e11 = 1; latvec.e12 = 0; latvec.e13 = 0;
    latvec.e21 = 0; latvec.e22 = 1; latvec.e23 = 0;
    latvec.e31 = 0; latvec.e32 = 0; latvec.e33 = 1;

    UnitCellLite ucell = make_test_ucell(
        1.0, 1.0, latvec, 1, {1},
        {{0.0, 0.0, 0.0}}
    );

    size_t total_inside = 0;
    for (int rank = 0; rank < 4; ++rank)
    {
        NeighborSearch ns;
        ns.init(ucell, 0.1, rank, 4);
        total_inside += ns.get_inside_atoms().size();
        EXPECT_EQ(ns.get_neighbor_list().get_nlocal(), static_cast<int>(ns.get_inside_atoms().size()));
        EXPECT_EQ(ns.get_inside_atoms().size(), rank == 0 ? 1U : 0U);
    }
    EXPECT_EQ(total_inside, 1U);
}

TEST(NeighborSearchInit_MpiOwnership, SplitsOnlyNonzeroSpanDirection)
{
    ModuleBase::Matrix3 latvec;
    latvec.e11 = 4; latvec.e12 = 0; latvec.e13 = 0;
    latvec.e21 = 0; latvec.e22 = 1; latvec.e23 = 0;
    latvec.e31 = 0; latvec.e32 = 0; latvec.e33 = 1;

    UnitCellLite ucell = make_test_ucell(
        1.0, 4.0, latvec, 1, {4},
        {{0.0, 0.0, 0.0},
         {1.0, 0.0, 0.0},
         {2.0, 0.0, 0.0},
         {3.0, 0.0, 0.0}}
    );

    size_t total_inside = 0;
    for (int rank = 0; rank < 4; ++rank)
    {
        NeighborSearch ns;
        ns.init(ucell, 0.1, rank, 4);
        total_inside += ns.get_inside_atoms().size();
        EXPECT_EQ(ns.get_y(), 0);
        EXPECT_EQ(ns.get_z(), 0);
        EXPECT_EQ(ns.get_inside_atoms().size(), 1U);
    }
    EXPECT_EQ(total_inside, 4U);
}

TEST(NeighborSearchInit_MpiLocalAtoms, LocalIdsAreValidAndAllAtomsShrink)
{
    ModuleBase::Matrix3 latvec;
    latvec.e11 = 4; latvec.e12 = 0; latvec.e13 = 0;
    latvec.e21 = 0; latvec.e22 = 4; latvec.e23 = 0;
    latvec.e31 = 0; latvec.e32 = 0; latvec.e33 = 4;

    std::vector<ModuleBase::Vector3<double>> tau;
    for (int ix = 0; ix < 4; ++ix)
    {
        for (int iy = 0; iy < 4; ++iy)
        {
            for (int iz = 0; iz < 4; ++iz)
            {
                tau.emplace_back(ix, iy, iz);
            }
        }
    }

    UnitCellLite ucell = make_test_ucell(1.0, 64.0, latvec, 1, {static_cast<int>(tau.size())}, tau);

    NeighborSearch serial;
    serial.init(ucell, 1.1, 0);
    serial.build_neighbors();
    const size_t serial_all_atoms = serial.get_all_atoms().size();
    size_t serial_neighbor_pairs = 0;
    for (int local_i = 0; local_i < serial.get_neighbor_list().get_nlocal(); ++local_i)
    {
        serial_neighbor_pairs += serial.get_neighbor_list().get_numneigh(local_i);
    }

    size_t total_inside = 0;
    size_t parallel_all_atoms_sum = 0;
    size_t parallel_all_atoms_max = 0;
    size_t parallel_neighbor_pairs = 0;
    for (int rank = 0; rank < 4; ++rank)
    {
        NeighborSearch ns;
        ns.init(ucell, 1.1, rank, 4);
        ns.build_neighbors();

        const auto& all_atoms = ns.get_all_atoms();
        const auto& list = ns.get_neighbor_list();
        total_inside += ns.get_inside_atoms().size();
        parallel_all_atoms_sum += all_atoms.size();
        parallel_all_atoms_max = std::max(parallel_all_atoms_max, all_atoms.size());
        for (int local_i = 0; local_i < list.get_nlocal(); ++local_i)
        {
            parallel_neighbor_pairs += list.get_numneigh(local_i);
        }

        for (size_t i = 0; i < all_atoms.size(); ++i)
        {
            EXPECT_EQ(all_atoms[i].atom_id, static_cast<int>(i));
        }

        for (int local_i = 0; local_i < list.get_nlocal(); ++local_i)
        {
            for (int ad = 0; ad < list.get_numneigh(local_i); ++ad)
            {
                const int neighbor_id = list.get_firstneigh(local_i)[ad];
                EXPECT_GE(neighbor_id, 0);
                EXPECT_LT(neighbor_id, static_cast<int>(all_atoms.size()));
            }
        }
    }

    EXPECT_EQ(total_inside, tau.size());
    EXPECT_EQ(parallel_neighbor_pairs, serial_neighbor_pairs);
    EXPECT_LT(parallel_all_atoms_max, serial_all_atoms);
    EXPECT_LT(parallel_all_atoms_sum, serial_all_atoms * 4);
}

TEST(NeighborSearchDistance_OutsideCases, VariousAxes)
{
    NeighborSearch ns;
    ns.set_position(0, 0, 0);
    ns.set_width(2.0, 3.0, 4.0);

    // position inside box along x (no dx), but outside along y by above high bound
    double d = ns.distance(0.5, 4.5, 1.0, 0.0, 0.0, 0.0);
    // dy = position_y - (y_low + (y+1)*wide_y) = 4.5 - 3.0 = 1.5 -> squared 2.25
    // dx = 0, dz = 0 -> total 2.25
    EXPECT_DOUBLE_EQ(d, 2.25);

    // position left of low bound on x
    double d2 = ns.distance(-1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
    // dx = x_low - position_x = 0 - (-1) = 1 -> squared 1
    EXPECT_DOUBLE_EQ(d2, 1.0);
}

TEST(NeighborSearchDecompose_SmallSizes, TwoAndOne)
{
    NeighborSearch ns;
    int nx, ny, nz;
    ns.decompose(2, nx, ny, nz);
    EXPECT_EQ(nx * ny * nz, 2);
    // possible decomposition is nx=1, ny=1, nz=2 (or nx=1, ny=2, nz=1 depending on algorithm)
    EXPECT_EQ(nx, 1);

    ns.decompose(1, nx, ny, nz);
    EXPECT_EQ(nx, 1);
    EXPECT_EQ(ny, 1);
    EXPECT_EQ(nz, 1);
}

TEST(NeighborSearchUnit, ExpansionLayersAndAtomCount)
{
    ModuleBase::Matrix3 latvec;
    latvec.e11 = 1; latvec.e12 = 0; latvec.e13 = 0;
    latvec.e21 = 0; latvec.e22 = 1; latvec.e23 = 0;
    latvec.e31 = 0; latvec.e32 = 0; latvec.e33 = 1;

    UnitCellLite ucell = make_test_ucell(
        1.0, 1.0, latvec, 1, {2},
        {{0.0, 0.0, 0.0}, {0.5, 0.0, 0.0}}
    );

    NeighborSearch ns;
    ns.init(ucell, 1.0, 0);

    // For identity lattice with search_radius=1 expected ceil produce values
    EXPECT_EQ(ns.get_glayerX(), 2);
    EXPECT_EQ(ns.get_glayerY(), 2);
    EXPECT_EQ(ns.get_glayerZ(), 2);
    EXPECT_EQ(ns.get_glayerX_minus(), 1);

    // Check atom count
    int images_x = ns.get_glayerX() + ns.get_glayerX_minus();
    int images_y = ns.get_glayerY() + ns.get_glayerY_minus();
    int images_z = ns.get_glayerZ() + ns.get_glayerZ_minus();
    int expected = images_x * images_y * images_z * 2; // 2 atoms per cell
    EXPECT_EQ(static_cast<int>(ns.get_all_atoms().size()), expected);
}

// end of additional tests
