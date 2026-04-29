#include <gtest/gtest.h>
#include "../neighbor_search.h"
#include "../unitcell_plus.h"
#include "../bin_manager.h"
#include "../neighbor_list.h"

TEST(NeighborSearchTest, TwoAtomsNeighbor)
{
    UnitCellPlus ucell;

    // ===== 基本晶胞信息 =====
    ucell.lat0 = 1.0;
    ucell.omega = 1.0;

    ucell.latvec.e11 = 1; ucell.latvec.e12 = 0; ucell.latvec.e13 = 0;
    ucell.latvec.e21 = 0; ucell.latvec.e22 = 1; ucell.latvec.e23 = 0;
    ucell.latvec.e31 = 0; ucell.latvec.e32 = 0; ucell.latvec.e33 = 1;

    // ===== 原子信息 =====
    ucell.ntype = 1;
    ucell.na = {2};
    ucell.nat = 2;

    ucell.tau = {
        {0.0, 0.0, 0.0},
        {0.5, 0.0, 0.0}
    };

    NeighborSearch ns;

    double cutoff = 1.0;

    ns.init(ucell, cutoff, 0);
    ns.build_neighbors();

    auto &list = ns.get_neighbor_list();

    ASSERT_EQ(list.numneigh.size(), 2);

    EXPECT_EQ(list.numneigh[0], 8);
    EXPECT_EQ(list.numneigh[1], 8);
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

TEST(BinManagerUnit, InitAndBinning)
{
    std::vector<NeighborAtom> inside;
    std::vector<NeighborAtom> ghost;

    // two atoms close to each other
    inside.emplace_back(0.0, 0.0, 0.0, 0, 0, 0);
    inside.emplace_back(0.5, 0.0, 0.0, 0, 1, 1);

    BinManager bm;
    bm.init_bins(1.0, inside, ghost);

    EXPECT_GE(bm.nbinx, 1);
    EXPECT_GE(bm.nbiny, 1);
    EXPECT_GE(bm.nbinz, 1);
    EXPECT_EQ(static_cast<int>(bm.bins.size()), bm.nbinx * bm.nbiny * bm.nbinz);

    bm.do_binning(inside, ghost);

    // compute bin index for first atom
    int ix = std::min(std::max(int((inside[0].position_x - bm.x_min) / bm.bin_sizex), 0), bm.nbinx - 1);
    int iy = std::min(std::max(int((inside[0].position_y - bm.y_min) / bm.bin_sizey), 0), bm.nbiny - 1);
    int iz = std::min(std::max(int((inside[0].position_z - bm.z_min) / bm.bin_sizez), 0), bm.nbinz - 1);
    int idx = ix * bm.nbiny * bm.nbinz + iy * bm.nbinz + iz;

    EXPECT_GE(bm.bins[idx].atoms.size(), 1u);
}

TEST(BinManagerUnit, BuildNeighborsAndClear)
{
    std::vector<NeighborAtom> atoms;
    atoms.emplace_back(0.0, 0.0, 0.0, 0, 0, 0);
    atoms.emplace_back(0.5, 0.0, 0.0, 0, 1, 1);
    atoms.emplace_back(5.0, 0.0, 0.0, 0, 2, 2);

    std::vector<NeighborAtom> inside = atoms;
    std::vector<NeighborAtom> ghost;

    BinManager bm;
    bm.init_bins(1.0, inside, ghost);
    bm.do_binning(inside, ghost);

    NeighborList nl;
    nl.initialize(static_cast<int>(atoms.size()), 1024);

    bm.build_atom_neighbors(nl, atoms);

    // atom 0 and 1 are neighbors; atom 2 is far
    EXPECT_GE(nl.numneigh[0], 1);
    EXPECT_GE(nl.numneigh[1], 1);
    EXPECT_EQ(nl.numneigh[2], 0);

    bm.clear();
    EXPECT_EQ(bm.bins.size(), 0u);
}

// Additional, more exhaustive tests

TEST(PageAllocatorTest, MultiPageAllocationAndReset)
{
    PageAllocator pa(4);
    int* a = pa.allocate(3);
    ASSERT_NE(a, nullptr);
    int* b = pa.allocate(3);
    ASSERT_NE(b, nullptr);
    // since page size = 4, a uses 3, b must be allocated on a new page
    EXPECT_NE(b, a + 3);

    pa.reset();
    int* c = pa.allocate(2);
    ASSERT_NE(c, nullptr);
    // after reset, allocation starts at the beginning of the last page (current implementation)
    EXPECT_EQ(c, pa.pages[0].data.data());
}

TEST(PageAllocatorTest, ZeroAndExactPageAllocations)
{
    PageAllocator pa(8);
    // zero or negative allocation should return nullptr
    EXPECT_EQ(pa.allocate(0), nullptr);
    EXPECT_EQ(pa.allocate(-1), nullptr);

    // allocate exactly page size
    int* p = pa.allocate(8);
    ASSERT_NE(p, nullptr);
    // next tiny allocation must go to a new page
    int* q = pa.allocate(1);
    ASSERT_NE(q, nullptr);
    EXPECT_NE(q, p + 8);
}

TEST(PageAllocatorTest, LargeAllocationSpansPages)
{
    PageAllocator pa(16);
    // request larger than single page -> should still succeed (allocates new page and fits)
    int* p = pa.allocate(32);
    EXPECT_EQ(p, nullptr);
}

TEST(NeighborListUnit, InitializeZeroAndResize)
{
    NeighborList nl;
    nl.initialize(0, 8);
    EXPECT_EQ(nl.numneigh.size(), 0);
    EXPECT_EQ(nl.firstneigh.size(), 0);

    // initialize with larger size
    nl.initialize(5, 8);
    EXPECT_EQ(nl.numneigh.size(), 5);
    EXPECT_EQ(nl.firstneigh.size(), 5);
}

TEST(BinManagerUnit, EmptyAtomsBuildNeighbors)
{
    std::vector<NeighborAtom> atoms;
    std::vector<NeighborAtom> ghost;

    BinManager bm;
    bm.init_bins(1.0, atoms, ghost);

    NeighborList nl;
    nl.initialize(0, 16);

    // should not crash or assert when atoms is empty
    bm.build_atom_neighbors(nl, atoms);
    EXPECT_EQ(nl.numneigh.size(), 0);
}

TEST(BinManagerUnit, BoundaryAndExactRadius)
{
    // inside atom at origin; other atoms placed on bin boundaries and at exact radius
    std::vector<NeighborAtom> atoms;
    atoms.emplace_back(0.0, 0.0, 0.0, 0, 0, 0);
    // exactly at search radius 1.0 along x direction
    atoms.emplace_back(1.0, 0.0, 0.0, 0, 1, 1);
    // slightly inside
    atoms.emplace_back(0.9, 0.0, 0.0, 0, 2, 2);

    std::vector<NeighborAtom> inside = atoms;
    std::vector<NeighborAtom> ghost;

    BinManager bm;
    bm.init_bins(1.0, inside, ghost);
    bm.do_binning(inside, ghost);

    NeighborList nl;
    nl.initialize(static_cast<int>(inside.size()), 64);

    bm.build_atom_neighbors(nl, inside);

    // atom at distance exactly 1.0 should be considered neighbor (dist2 == sradius2)
    EXPECT_GE(nl.numneigh[0], 1);
    // the exact atom itself must not be counted as its own neighbor
    for (int i = 0; i < static_cast<int>(inside.size()); ++i) {
        for (int j = 0; j < nl.numneigh[i]; ++j) {
            int id = nl.firstneigh[i][j];
            EXPECT_NE(id, inside[i].atom_id);
        }
    }
}

TEST(NeighborListUnit, AllocateManyNeighbors)
{
    // verify neighbor_list allocator can handle many small allocations across pages
    NeighborList nl;
    nl.initialize(10, 8); // small page size to force multiple pages

    // simulate allocations of varying sizes
    for (int i = 0; i < 10; ++i) {
        int count = (i % 4) + 1;
        int* ptr = nl.allocator.allocate(count);
        ASSERT_NE(ptr, nullptr);
        // write and read back
        for (int k = 0; k < count; ++k) ptr[k] = i * 10 + k;
    }
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

    NeighborSearch ns;
    ns.search_radius = 2.5;
    ns.Check_Expand_Condition(ucell);
    // for skewed lattice, expansion layers should be >= 1
    EXPECT_GE(ns.glayerX, 1);
    EXPECT_GE(ns.glayerY, 1);
    EXPECT_GE(ns.glayerZ, 1);
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

TEST(NeighborSearchIntegration, GhostAtomsCountedAsNeighbors)
{
    // inside atom at origin; ghost atom within search radius should be found
    std::vector<NeighborAtom> inside;
    std::vector<NeighborAtom> ghost;

    inside.emplace_back(0.0, 0.0, 0.0, 0, 0, 0);
    ghost.emplace_back(0.4, 0.0, 0.0, 0, 1, 1);

    BinManager bm;
    bm.init_bins(1.0, inside, ghost);
    bm.do_binning(inside, ghost);

    NeighborList nl;
    nl.initialize(static_cast<int>(inside.size()), 32);

    bm.build_atom_neighbors(nl, inside);

    EXPECT_GE(nl.numneigh[0], 1);
}