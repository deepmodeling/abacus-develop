#include <gtest/gtest.h>
#include "../bin_manager.h"
#include "../neighbor_list.h"

TEST(BinManagerUnit, InitAndBinning)
{
    std::vector<NeighborAtom> inside;
    std::vector<NeighborAtom> ghost;

    inside.emplace_back(0.0, 0.0, 0.0, 0, 0, 0);
    inside.emplace_back(0.5, 0.0, 0.0, 0, 1, 1);

    BinManager bm;
    bm.init_bins(1.0, inside);

    EXPECT_EQ(bm.get_nbinx(), 1);
    EXPECT_EQ(bm.get_nbiny(), 1);
    EXPECT_EQ(bm.get_nbinz(), 1);
    EXPECT_EQ(bm.get_total_bins(), bm.get_nbinx() * bm.get_nbiny() * bm.get_nbinz());

    bm.do_binning(inside);

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

    BinManager bm;
    bm.init_bins(1.0, atoms);
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

    BinManager bm;
    bm.init_bins(1.0, atoms);
    EXPECT_EQ(bm.get_nbinx(), 5);
    EXPECT_EQ(bm.get_nbiny(), 1);
    EXPECT_EQ(bm.get_nbinz(), 1);
    EXPECT_EQ(bm.get_total_bins(), bm.get_nbinx() * bm.get_nbiny() * bm.get_nbinz());

    bm.do_binning(atoms);

    NeighborList nl;
    nl.initialize(static_cast<int>(atoms.size()), 1024);

    bm.build_atom_neighbors(nl, atoms, atoms);

    EXPECT_EQ(nl.get_numneigh(0), 1);
    EXPECT_EQ(nl.get_numneigh(1), 1);
    EXPECT_EQ(nl.get_numneigh(2), 0);

    bm.clear();
    EXPECT_EQ(bm.get_total_bins(), 0);
}

TEST(BinManagerUnit, EmptyAtomsBuildNeighbors)
{
    std::vector<NeighborAtom> atoms;
    std::vector<NeighborAtom> ghost;

    BinManager bm;
    bm.init_bins(1.0, atoms);

    NeighborList nl;
    nl.initialize(0, 16);

    bm.build_atom_neighbors(nl, atoms, atoms);
    EXPECT_EQ(nl.get_nlocal(), 0);
}

TEST(BinManagerUnit, BoundaryAndExactRadius)
{
    std::vector<NeighborAtom> atoms;
    atoms.emplace_back(0.0, 0.0, 0.0, 0, 0, 0);
    atoms.emplace_back(1.0, 0.0, 0.0, 0, 1, 1);
    atoms.emplace_back(0.9, 0.0, 0.0, 0, 2, 2);

    BinManager bm;
    bm.init_bins(1.0, atoms);
    bm.do_binning(atoms);

    NeighborList nl;
    nl.initialize(atoms.size(), 64);

    bm.build_atom_neighbors(nl, atoms, atoms);

    EXPECT_EQ(nl.get_numneigh(0), 2);
    for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
        for (int j = 0; j < nl.get_numneigh(i); ++j) {
            int id = nl.get_firstneigh(i)[j];
            EXPECT_NE(id, atoms[i].atom_id);
        }
    }
}

TEST(BinManagerUnit, InitWithGhostOnly)
{
    std::vector<NeighborAtom> inside;
    std::vector<NeighborAtom> ghost;

    ghost.emplace_back(-1.0, -1.0, -1.0, 0, 0, 0);
    ghost.emplace_back(2.0, 0.0, 0.0, 0, 1, 1);

    BinManager bm;
    bm.init_bins(1.0, ghost);

    EXPECT_EQ(bm.get_nbinx(), 3);
    EXPECT_EQ(bm.get_nbiny(), 1);
    EXPECT_EQ(bm.get_nbinz(), 1);
}

TEST(BinManagerUnit, BuildNeighborsNoNeighborsFirstneighNull)
{
    std::vector<NeighborAtom> atoms;
    atoms.emplace_back(0.0, 0.0, 0.0, 0, 0, 0);
    atoms.emplace_back(100.0, 100.0, 100.0, 0, 1, 1);

    BinManager bm;
    bm.init_bins(1.0, atoms);
    bm.do_binning(atoms);

    NeighborList nl;
    nl.initialize(atoms.size(), 8);

    bm.build_atom_neighbors(nl, atoms, atoms);

    EXPECT_EQ(nl.get_numneigh(0), 0);
    EXPECT_EQ(nl.get_numneigh(1), 0);
    EXPECT_EQ(nl.get_firstneigh(0), nullptr);
    EXPECT_EQ(nl.get_firstneigh(1), nullptr);
}

TEST(BinManagerUnit, GhostAtomsAreCounted)
{
    std::vector<NeighborAtom> inside;
    std::vector<NeighborAtom> ghost;

    inside.emplace_back(0.0, 0.0, 0.0, 0, 0, 0);
    ghost.emplace_back(0.4, 0.0, 0.0, 0, 1, 1, 3, 1);

    BinManager bm;
    std::vector<NeighborAtom> all_atoms = inside;
    all_atoms.insert(all_atoms.end(), ghost.begin(), ghost.end());
    bm.init_bins(1.0, all_atoms);
    bm.do_binning(all_atoms);

    NeighborList nl;
    nl.initialize(static_cast<int>(inside.size()), 32);

    bm.build_atom_neighbors(nl, inside, all_atoms);

    EXPECT_EQ(nl.get_nlocal(), 1);
    EXPECT_EQ(nl.get_numneigh(0), 1);
    bool found = false;
    if (nl.get_numneigh(0) > 0 && nl.get_firstneigh(0) != nullptr) {
        for (int k = 0; k < nl.get_numneigh(0); ++k) {
            if (nl.get_firstneigh(0)[k] == 1) found = true;
        }
    }
    EXPECT_TRUE(found);
}

TEST(BinManagerUnit, SamePositionDifferentAtomsAreNeighbors)
{
    std::vector<NeighborAtom> atoms;
    atoms.emplace_back(0.0, 0.0, 0.0, 0, 0, 0);
    atoms.emplace_back(0.0, 0.0, 0.0, 0, 1, 1);

    BinManager bm;
    bm.init_bins(1.0, atoms);
    bm.do_binning(atoms);

    NeighborList nl;
    nl.initialize(atoms.size(), 16);

    bm.build_atom_neighbors(nl, atoms, atoms);

    EXPECT_EQ(nl.get_numneigh(0), 1);
    EXPECT_EQ(nl.get_numneigh(1), 1);
    ASSERT_NE(nl.get_firstneigh(0), nullptr);
    ASSERT_NE(nl.get_firstneigh(1), nullptr);
    EXPECT_EQ(nl.get_firstneigh(0)[0], 1);
    EXPECT_EQ(nl.get_firstneigh(1)[0], 0);
}

TEST(BinManagerUnit, MultipleBinsNeighborSearch)
{
    std::vector<NeighborAtom> atoms;
    int id = 0;
    for (int x = 0; x < 3; ++x)
        for (int y = 0; y < 3; ++y)
            for (int z = 0; z < 3; ++z)
                atoms.emplace_back(x * 1.0, y * 1.0, z * 1.0, 0, 0, id++);

    BinManager bm;
    bm.init_bins(1.0, atoms);
    bm.do_binning(atoms);

    NeighborList nl;
    nl.initialize(atoms.size(), 16);

    bm.build_atom_neighbors(nl, atoms, atoms);

    int center_index = 13;
    EXPECT_EQ(nl.get_numneigh(center_index), 6);
}
