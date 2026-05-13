#include <gtest/gtest.h>
#include "../bin_manager.h"
#include "../neighbor_list.h"

TEST(BinManagerUnit, InitAndBinning)
{
    std::vector<NeighborAtom> inside;
    std::vector<NeighborAtom> ghost;

    // two atoms close to each other
    inside.emplace_back(0.0, 0.0, 0.0, 0, 0, 0);
    inside.emplace_back(0.5, 0.0, 0.0, 0, 1, 1);

    BinManager bm;
    bm.init_bins(1.0, inside, ghost);

    EXPECT_EQ(bm.nbinx, 1);
    EXPECT_EQ(bm.nbiny, 1);
    EXPECT_EQ(bm.nbinz, 1);
    EXPECT_EQ(static_cast<int>(bm.bins.size()), bm.nbinx * bm.nbiny * bm.nbinz);

    bm.do_binning(inside, ghost);

    // compute bin index for first atom
    int ix = std::min(std::max(int((inside[0].position_x - bm.x_min) / bm.bin_sizex), 0), bm.nbinx - 1);
    int iy = std::min(std::max(int((inside[0].position_y - bm.y_min) / bm.bin_sizey), 0), bm.nbiny - 1);
    int iz = std::min(std::max(int((inside[0].position_z - bm.z_min) / bm.bin_sizez), 0), bm.nbinz - 1);
    int idx = ix * bm.nbiny * bm.nbinz + iy * bm.nbinz + iz;

    EXPECT_GE(bm.bins[idx].atoms.size(), 1u);
}

TEST(BinManagerUnit, InitBins)
{
    std::vector<NeighborAtom> atoms;
    atoms.emplace_back(0.0, 0.0, 0.0, 0, 0, 0);
    atoms.emplace_back(0.5, 0.0, 0.0, 0, 1, 1);
    atoms.emplace_back(4.9, 0.0, 0.0, 0, 2, 2);

    std::vector<NeighborAtom> inside = atoms;
    std::vector<NeighborAtom> ghost;

    BinManager bm;
    bm.init_bins(1.0, inside, ghost);
    EXPECT_EQ(bm.nbinx, 5);
    EXPECT_EQ(bm.nbiny, 1);
    EXPECT_EQ(bm.nbinz, 1);
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
    EXPECT_EQ(bm.nbinx, 5);
    EXPECT_EQ(bm.nbiny, 1);
    EXPECT_EQ(bm.nbinz, 1);
    EXPECT_EQ(static_cast<int>(bm.bins.size()), bm.nbinx * bm.nbiny * bm.nbinz);

    bm.do_binning(inside, ghost);

    NeighborList nl;
    nl.initialize(static_cast<int>(atoms.size()), 1024);

    bm.build_atom_neighbors(nl, atoms);

    // atom 0 and 1 are neighbors; atom 2 is far
    EXPECT_EQ(nl.numneigh[0], 1);
    EXPECT_EQ(nl.numneigh[1], 1);
    EXPECT_EQ(nl.numneigh[2], 0);

    bm.clear();
    EXPECT_EQ(bm.bins.size(), 0u);
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
    EXPECT_EQ(nl.numneigh[0], 2);
    // the exact atom itself must not be counted as its own neighbor
    for (int i = 0; i < static_cast<int>(inside.size()); ++i) {
        for (int j = 0; j < nl.numneigh[i]; ++j) {
            int id = nl.firstneigh[i][j];
            EXPECT_NE(id, inside[i].atom_id);
        }
    }
}

TEST(BinManagerUnit, InitWithGhostOnly)
{
    // inside empty, ghosts present -> init_bins should still compute bounds from ghosts
    std::vector<NeighborAtom> inside;
    std::vector<NeighborAtom> ghost;

    ghost.emplace_back(-1.0, -1.0, -1.0, 0, 0, 0);
    ghost.emplace_back(2.0, 0.0, 0.0, 0, 1, 1);

    BinManager bm;
    bm.init_bins(1.0, inside, ghost);

    EXPECT_EQ(bm.nbinx, 3);
    EXPECT_EQ(bm.nbiny, 1);
    EXPECT_EQ(bm.nbinz, 1);
}

TEST(BinManagerUnit, BuildNeighborsNoNeighborsFirstneighNull)
{
    // atoms all far apart => no neighbors; firstneigh entries should be nullptr
    std::vector<NeighborAtom> atoms;
    atoms.emplace_back(0.0, 0.0, 0.0, 0, 0, 0);
    atoms.emplace_back(100.0, 100.0, 100.0, 0, 1, 1);

    std::vector<NeighborAtom> inside = atoms;
    std::vector<NeighborAtom> ghost;

    BinManager bm;
    bm.init_bins(1.0, inside, ghost);
    bm.do_binning(inside, ghost);

    NeighborList nl;
    nl.initialize(static_cast<int>(inside.size()), 8);

    bm.build_atom_neighbors(nl, inside);

    EXPECT_EQ(nl.numneigh[0], 0);
    EXPECT_EQ(nl.numneigh[1], 0);
    EXPECT_EQ(nl.firstneigh[0], nullptr);
    EXPECT_EQ(nl.firstneigh[1], nullptr);
}

TEST(BinManagerUnit, GhostAtomsAreCounted)
{
    // inside atom at origin; ghost atom within search radius should be found
    std::vector<NeighborAtom> inside;
    std::vector<NeighborAtom> ghost;

    inside.emplace_back(0.0, 0.0, 0.0, 0, 0, 0);
    ghost.emplace_back(0.4, 0.0, 0.0, 0, 1, 3);

    BinManager bm;
    bm.init_bins(1.0, inside, ghost);
    bm.do_binning(inside, ghost);

    NeighborList nl;
    nl.initialize(static_cast<int>(inside.size()), 32);

    bm.build_atom_neighbors(nl, inside);

    EXPECT_EQ(nl.numneigh.size(), 1);
    EXPECT_EQ(nl.numneigh[0], 1);
    // ensure neighbor id matches ghost atom id
    bool found = false;
    if (nl.numneigh[0] > 0 && nl.firstneigh[0] != nullptr) {
        for (int k = 0; k < nl.numneigh[0]; ++k) {
            if (nl.firstneigh[0][k] == 3) found = true;
        }
    }
    EXPECT_TRUE(found);
}

TEST(BinManagerUnit, MultipleBinsNeighborSearch)
{
    // Create a 3x3x3 grid of atoms spaced by 1.0 so multiple bins exist
    std::vector<NeighborAtom> atoms;
    int id = 0;
    for (int x = 0; x < 3; ++x)
        for (int y = 0; y < 3; ++y)
            for (int z = 0; z < 3; ++z)
                atoms.emplace_back(x * 1.0, y * 1.0, z * 1.0, 0, 0, id++);

    std::vector<NeighborAtom> inside = atoms;
    std::vector<NeighborAtom> ghost;

    BinManager bm;
    bm.init_bins(1.0, inside, ghost);
    bm.do_binning(inside, ghost);

    NeighborList nl;
    nl.initialize(static_cast<int>(inside.size()), 16);

    bm.build_atom_neighbors(nl, inside);

    // center atom (1,1,1) should have multiple neighbors
    int center_index = 13; // 1*9 + 1*3 + 1
    EXPECT_EQ(nl.numneigh[center_index], 6);
}
