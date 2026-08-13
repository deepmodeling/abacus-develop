#include <gtest/gtest.h>
#include "source_cell/module_neighlist/bin_manager.h"
#include "source_cell/module_neighlist/neighbor_list.h"

#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{
std::vector<std::vector<int>> snapshot_neighbors(const NeighborList& list)
{
    std::vector<std::vector<int>> result(static_cast<std::size_t>(list.get_nlocal()));
    for (int i = 0; i < list.get_nlocal(); ++i)
    {
        const int count = list.get_numneigh(i);
        const int* first = list.get_firstneigh(i);
        if (count > 0)
        {
            result[static_cast<std::size_t>(i)].assign(first, first + count);
        }
    }
    return result;
}
} // namespace

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

#ifdef _OPENMP
TEST(BinManagerUnit, ParallelBuildPreservesSerialNeighborOrder)
{
    std::vector<NeighborAtom> centers;
    std::vector<NeighborAtom> binned_atoms;
    centers.reserve(300);
    binned_atoms.reserve(600);
    for (int i = 0; i < 300; ++i)
    {
        const double x = 0.2 * static_cast<double>(i % 10);
        const double y = 0.2 * static_cast<double>((i / 10) % 10);
        const double z = 0.2 * static_cast<double>(i / 100);
        centers.emplace_back(x, y, z, 0, i, i);
        binned_atoms.push_back(centers.back());
    }
    for (int i = 0; i < 300; ++i)
    {
        const NeighborAtom& center = centers[static_cast<std::size_t>(i)];
        binned_atoms.emplace_back(center.position_x + 0.05,
                                  center.position_y,
                                  center.position_z,
                                  0,
                                  i,
                                  300 + i,
                                  1000 + i,
                                  1);
    }

    BinManager bm;
    bm.init_bins(0.31, binned_atoms);
    bm.do_binning(binned_atoms);

    const int previous_dynamic = omp_get_dynamic();
    const int previous_threads = omp_get_max_threads();
    omp_set_dynamic(0);

    NeighborList serial_list;
    serial_list.initialize(centers.size(), centers.size() * 128);
    omp_set_num_threads(1);
    bm.build_atom_neighbors(serial_list, centers, binned_atoms);
    const std::vector<std::vector<int>> serial = snapshot_neighbors(serial_list);

    NeighborList parallel_list;
    parallel_list.initialize(centers.size(), centers.size() * 128);
    omp_set_num_threads(4);
    bm.build_atom_neighbors(parallel_list, centers, binned_atoms);
    const std::vector<std::vector<int>> parallel = snapshot_neighbors(parallel_list);

    omp_set_num_threads(previous_threads);
    omp_set_dynamic(previous_dynamic);

    EXPECT_EQ(parallel, serial);
    EXPECT_NE(std::find_if(serial[0].begin(), serial[0].end(),
                           [](int atom_id) { return atom_id >= 300; }),
              serial[0].end());
}

TEST(BinManagerUnit, ParallelBuildKeepsZeroNeighborPointersNull)
{
    std::vector<NeighborAtom> atoms;
    atoms.reserve(256);
    for (int i = 0; i < 256; ++i)
    {
        atoms.emplace_back(static_cast<double>(i), 0.0, 0.0, 0, i, i);
    }

    BinManager bm;
    bm.init_bins(0.4, atoms);
    bm.do_binning(atoms);

    const int previous_dynamic = omp_get_dynamic();
    const int previous_threads = omp_get_max_threads();
    omp_set_dynamic(0);
    omp_set_num_threads(4);

    NeighborList list;
    list.initialize(atoms.size(), 1024);
    bm.build_atom_neighbors(list, atoms, atoms);

    omp_set_num_threads(previous_threads);
    omp_set_dynamic(previous_dynamic);

    for (int i = 0; i < list.get_nlocal(); ++i)
    {
        EXPECT_EQ(list.get_numneigh(i), 0);
        EXPECT_EQ(list.get_firstneigh(i), nullptr);
    }
}
#endif
