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

    EXPECT_EQ(bm.get_nbinx(), 1);
    EXPECT_EQ(bm.get_nbiny(), 1);
    EXPECT_EQ(bm.get_nbinz(), 1);
    EXPECT_EQ(bm.get_total_bins(), bm.get_nbinx() * bm.get_nbiny() * bm.get_nbinz());

    bm.do_binning(inside, ghost);

    // After binning, at least one atom should be in some bin
    int total_atoms_in_bins = 0;
    for (int i = 0; i < bm.get_total_bins(); ++i) {
        total_atoms_in_bins += bm.get_bin_atom_count(i);
    }
    EXPECT_GE(total_atoms_in_bins, 2);
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
    EXPECT_EQ(bm.get_nbinx(), 5);
    EXPECT_EQ(bm.get_nbiny(), 1);
    EXPECT_EQ(bm.get_nbinz(), 1);
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
    EXPECT_EQ(bm.get_nbinx(), 5);
    EXPECT_EQ(bm.get_nbiny(), 1);
    EXPECT_EQ(bm.get_nbinz(), 1);
    EXPECT_EQ(bm.get_total_bins(), bm.get_nbinx() * bm.get_nbiny() * bm.get_nbinz());

    bm.do_binning(inside, ghost);

    NeighborList nl;
    nl.initialize(static_cast<int>(atoms.size()), 1024);

    bm.build_atom_neighbors(nl, atoms);

    // atom 0 and 1 are neighbors; atom 2 is far
    EXPECT_EQ(nl.numneigh[0], 1);
    EXPECT_EQ(nl.numneigh[1], 1);
    EXPECT_EQ(nl.numneigh[2], 0);

    bm.clear();
    EXPECT_EQ(bm.get_total_bins(), 0);
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

    EXPECT_EQ(bm.get_nbinx(), 3);
    EXPECT_EQ(bm.get_nbiny(), 1);
    EXPECT_EQ(bm.get_nbinz(), 1);
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