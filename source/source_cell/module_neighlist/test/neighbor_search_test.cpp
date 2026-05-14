#include <gtest/gtest.h>
#include "../neighbor_search.h"
#include "../unitcell_plus.h"

TEST(NeighborSearchTest, TwoAtomsNeighbor)
{
    UnitCellPlus ucell;

    ucell.lat0 = 1.0;
    ucell.omega = 1.0;

    ucell.latvec.e11 = 1; ucell.latvec.e12 = 0; ucell.latvec.e13 = 0;
    ucell.latvec.e21 = 0; ucell.latvec.e22 = 1; ucell.latvec.e23 = 0;
    ucell.latvec.e31 = 0; ucell.latvec.e32 = 0; ucell.latvec.e33 = 1;

    ucell.ntype = 1;
    ucell.na = {2};
    ucell.nat = 2;

    ucell.tau = {
        {0.0, 0.0, 0.0},
        {0.5, 0.0, 0.0}
    };
    ucell.compute_naa();

    NeighborSearch ns;

    double cutoff = 1.0;

    ns.init(ucell, cutoff, 0);
    ns.build_neighbors();

    auto &list = ns.get_neighbor_list();

    ASSERT_EQ(list.numneigh.size(), 2);

    EXPECT_EQ(list.numneigh[0], 8);
    EXPECT_EQ(list.numneigh[1], 8);
}

TEST(NeighborSearchTest, NoNeighbor)
{
    UnitCellPlus ucell;

    ucell.lat0 = 1.0;
    ucell.omega = 1.0;

    ucell.latvec.e11 = 1; ucell.latvec.e12 = 0; ucell.latvec.e13 = 0;
    ucell.latvec.e21 = 0; ucell.latvec.e22 = 1; ucell.latvec.e23 = 0;
    ucell.latvec.e31 = 0; ucell.latvec.e32 = 0; ucell.latvec.e33 = 1;

    ucell.ntype = 1;
    ucell.na = {2};
    ucell.nat = 2;

    ucell.tau = {
        {0.0, 0.0, 0.0},
        {5.0, 0.0, 0.0}
    };
    ucell.compute_naa();

    NeighborSearch ns;

    // use a smaller search radius to avoid counting periodic-image neighbors
    ns.init(ucell, 0.1, 0);
    ns.build_neighbors();

    auto &list = ns.get_neighbor_list();

    EXPECT_EQ(list.numneigh[0], 0);
    EXPECT_EQ(list.numneigh[1], 0);
}

TEST(NeighborSearchUnit, UCellToInputAtoms)
{
    UnitCellPlus ucell;
    ucell.lat0 = 1.0;
    ucell.omega = 1.0;
    ucell.latvec.e11 = 1; ucell.latvec.e12 = 0; ucell.latvec.e13 = 0;
    ucell.latvec.e21 = 0; ucell.latvec.e22 = 1; ucell.latvec.e23 = 0;
    ucell.latvec.e31 = 0; ucell.latvec.e32 = 0; ucell.latvec.e33 = 1;

    ucell.ntype = 1;
    ucell.na = {2};
    ucell.nat = 2;
    ucell.tau = {
        {0.0, 0.0, 0.0},
        {0.5, 0.0, 0.0}
    };
    ucell.compute_naa();
    NeighborSearch ns;
    auto inputs = ns.ucell_to_input_atoms(ucell);

    EXPECT_EQ(inputs.n_atoms, 2);
    ASSERT_EQ(inputs.InputAtom.size(), 2);
    EXPECT_DOUBLE_EQ(inputs.InputAtom[0].position_x, 0.0);
    EXPECT_DOUBLE_EQ(inputs.InputAtom[1].position_x, 0.5);
    EXPECT_DOUBLE_EQ(inputs.x_low, 0.0);
    EXPECT_DOUBLE_EQ(inputs.x_high, 0.5);
}

TEST(NeighborSearchUnit, CheckExpandAndSetMembers)
{
    UnitCellPlus ucell;
    ucell.lat0 = 1.0;
    ucell.omega = 1.0;
    ucell.latvec.e11 = 1; ucell.latvec.e12 = 0; ucell.latvec.e13 = 0;
    ucell.latvec.e21 = 0; ucell.latvec.e22 = 1; ucell.latvec.e23 = 0;
    ucell.latvec.e31 = 0; ucell.latvec.e32 = 0; ucell.latvec.e33 = 1;

    ucell.ntype = 1;
    ucell.na = {2};
    ucell.nat = 2;
    ucell.tau = {
        {0.0, 0.0, 0.0},
        {0.5, 0.0, 0.0}
    };
    ucell.compute_naa();
    NeighborSearch ns;
    ns.search_radius = 1.0; // use search radius = 1 for Check_Expand_Condition
    ns.Check_Expand_Condition(ucell);

    // For identity lattice with search_radius=1 expected ceil produce values
    EXPECT_EQ(ns.glayerX, 2);
    EXPECT_EQ(ns.glayerY, 2);
    EXPECT_EQ(ns.glayerZ, 2);
    EXPECT_EQ(ns.glayerX_minus, 1);

    // Now populate all_atoms and check count
    ns.setMemberVariables(ucell);
    int images_x = ns.glayerX + ns.glayerX_minus; // iterations in x
    int images_y = ns.glayerY + ns.glayerY_minus;
    int images_z = ns.glayerZ + ns.glayerZ_minus;
    int expected = images_x * images_y * images_z * 2; // 2 atoms per cell
    EXPECT_EQ(static_cast<int>(ns.all_atoms.size()), expected);
}

TEST(NeighborSearchUnit, DistanceBox)
{
    NeighborSearch ns;
    // set a single cell region at x=0..1,y=0..1,z=0..1
    ns.x = 0; ns.y = 0; ns.z = 0;
    ns.wide_x = 1.0; ns.wide_y = 1.0; ns.wide_z = 1.0;

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

TEST(NeighborSearchUnit, UCellToInputAtomsMultipleTypes)
{
    UnitCellPlus ucell;
    ucell.lat0 = 1.0;
    ucell.omega = 1.0;
    ucell.latvec.e11 = 1; ucell.latvec.e12 = 0; ucell.latvec.e13 = 0;
    ucell.latvec.e21 = 0; ucell.latvec.e22 = 1; ucell.latvec.e23 = 0;
    ucell.latvec.e31 = 0; ucell.latvec.e32 = 0; ucell.latvec.e33 = 1;

    ucell.ntype = 2;
    ucell.na = {1, 2};
    ucell.nat = 3;
    ucell.tau = {
        {0.0, 0.0, 0.0},
        {0.5, 0.0, 0.0},
        {0.0, 0.5, 0.0}
    };
    ucell.compute_naa();
    NeighborSearch ns;
    auto inputs = ns.ucell_to_input_atoms(ucell);

    EXPECT_EQ(inputs.n_atoms, 3);
    ASSERT_EQ(inputs.InputAtom.size(), 3);
    EXPECT_DOUBLE_EQ(inputs.InputAtom[2].position_y, 0.5);
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

TEST(NeighborSearchUnit, NonOrthogonalLatticeExpand)
{
    UnitCellPlus ucell;
    ucell.lat0 = 1.0;
    ucell.omega = 1.0;
    // skewed lattice
    ucell.latvec.e11 = 1; ucell.latvec.e12 = 0.3; ucell.latvec.e13 = 0.0;
    ucell.latvec.e21 = 0.1; ucell.latvec.e22 = 1.0; ucell.latvec.e23 = 0.0;
    ucell.latvec.e31 = 0.0; ucell.latvec.e32 = 0.0; ucell.latvec.e33 = 1.0;

    ucell.ntype = 1;
    ucell.na = {1};
    ucell.nat = 1;
    ucell.tau = {{0.0, 0.0, 0.0}};
    ucell.compute_naa();

    NeighborSearch ns;
    ns.search_radius = 2.5;
    ns.Check_Expand_Condition(ucell);
    // for skewed lattice, expansion layers should be >= 1
    EXPECT_GE(ns.glayerX, 1);
    EXPECT_GE(ns.glayerY, 1);
    EXPECT_GE(ns.glayerZ, 1);
}

// --- additional tests to cover remaining branches in neighbor_search.cpp ---

TEST(NeighborSearchInit_WideZero_CentralInside, SingleAtomCell)
{
    UnitCellPlus ucell;
    ucell.lat0 = 1.0;
    ucell.omega = 1.0;
    ucell.latvec.e11 = 1; ucell.latvec.e12 = 0; ucell.latvec.e13 = 0;
    ucell.latvec.e21 = 0; ucell.latvec.e22 = 1; ucell.latvec.e23 = 0;
    ucell.latvec.e31 = 0; ucell.latvec.e32 = 0; ucell.latvec.e33 = 1;

    ucell.ntype = 1;
    ucell.na = {1};
    ucell.nat = 1;
    ucell.tau = {{0.0, 0.0, 0.0}};
    ucell.compute_naa();

    NeighborSearch ns;
    // choose sr small enough; with mpi_size fixed to 1 in init, wide_* become 0
    ns.init(ucell, 0.1, 0);
    // central cell atom should be counted as inside
    EXPECT_EQ(ns.inside_atoms.size(), 1);
    EXPECT_EQ(ns.neighbor_list.nlocal, static_cast<int>(ns.inside_atoms.size()));
}

TEST(NeighborSearchInit_MpiRankIndexing, RankValues)
{
    UnitCellPlus ucell;
    ucell.lat0 = 1.0;
    ucell.omega = 1.0;
    ucell.latvec.e11 = 1; ucell.latvec.e12 = 0; ucell.latvec.e13 = 0;
    ucell.latvec.e21 = 0; ucell.latvec.e22 = 1; ucell.latvec.e23 = 0;
    ucell.latvec.e31 = 0; ucell.latvec.e32 = 0; ucell.latvec.e33 = 1;

    ucell.ntype = 1;
    ucell.na = {1};
    ucell.nat = 1;
    ucell.tau = {{0.0, 0.0, 0.0}};
    ucell.compute_naa();

    NeighborSearch ns0;
    ns0.init(ucell, 0.5, 0);
    // with mpi_size fixed to 1 in init, nx=ny=nz=1; for rank 0 expect x=y=0,z=0
    EXPECT_EQ(ns0.x, 0);
    EXPECT_EQ(ns0.y, 0);
    EXPECT_EQ(ns0.z, 0);

}

TEST(NeighborSearchDistance_OutsideCases, VariousAxes)
{
    NeighborSearch ns;
    ns.x = 0; ns.y = 0; ns.z = 0;
    ns.wide_x = 2.0; ns.wide_y = 3.0; ns.wide_z = 4.0;

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

// end of additional tests
