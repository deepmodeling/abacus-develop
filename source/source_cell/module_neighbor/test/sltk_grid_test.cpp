#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <chrono>
#include <functional>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#define private public
#include "../sltk_grid.h"
#include "prepare_unitcell.h"
#include "source_io/module_parameter/parameter.h"
#include "synthetic_neighbor_unitcell.h"
#undef private
#include "source_cell/read_stru.h"
#ifdef __LCAO
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
LCAO_Orbitals::LCAO_Orbitals()
{
}
LCAO_Orbitals::~LCAO_Orbitals()
{
}
#endif
Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
    this->start_mag = nullptr;
}
Magnetism::~Magnetism()
{
    delete[] this->start_mag;
}

/************************************************
 *  unit test of sltk_grid
 ***********************************************/

/**
 * - Tested Functions:
 *   - Init: Grid::init()
 *     - setMemberVariables: really set member variables
 *       (like dx, dy, dz and d_minX, d_minY, d_minZ) by
 *       reading from getters of Atom_input, and construct the
 *       member Cell as a 3D array of CellSet
 */

void SetGlobalV()
{
    PARAM.input.test_grid = 0;
}

class SltkGridTest : public testing::Test
{
  protected:
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    std::ofstream ofs;
    std::ifstream ifs;
    bool pbc = true;
    double radius = ((8 + 5.01) * 2.0 + 0.01) / 10.2;
    int test_atom_in = 0;
    std::string output;
    void SetUp()
    {
        SetGlobalV();
        ucell = utp.SetUcellInfo();
    }
    void TearDown()
    {
        delete ucell;
    }
};

using SltkGridDeathTest = SltkGridTest;

namespace
{
using NeighborPairKey = std::tuple<int, int, int, int, int, int, int>;

std::vector<NeighborPairKey> collect_neighbor_pair_keys(const Grid& grid)
{
    std::vector<NeighborPairKey> keys;
    for (int center_type = 0; center_type < static_cast<int>(grid.all_adj_info.size()); ++center_type)
    {
        const auto& type_neighbors = grid.all_adj_info[center_type];
        for (int center_natom = 0; center_natom < static_cast<int>(type_neighbors.size()); ++center_natom)
        {
            const auto& atom_neighbors = type_neighbors[center_natom];
            for (const FAtom* atom: atom_neighbors)
            {
                keys.push_back(NeighborPairKey(center_type,
                                               center_natom,
                                               atom->type,
                                               atom->natom,
                                               atom->cell_x,
                                               atom->cell_y,
                                               atom->cell_z));
            }
        }
    }
    std::sort(keys.begin(), keys.end());
    return keys;
}

std::vector<ModuleNeighbor::NeighborPair> collect_legacy_neighbor_pairs(const Grid& grid)
{
    std::vector<ModuleNeighbor::NeighborPair> pairs;
    for (int center_type = 0; center_type < static_cast<int>(grid.all_adj_info.size()); ++center_type)
    {
        const auto& type_neighbors = grid.all_adj_info[center_type];
        for (int center_natom = 0; center_natom < static_cast<int>(type_neighbors.size()); ++center_natom)
        {
            const auto& atom_neighbors = type_neighbors[center_natom];
            for (const FAtom* atom: atom_neighbors)
            {
                ModuleNeighbor::NeighborPair pair;
                pair.center_type = center_type;
                pair.center_natom = center_natom;
                pair.neighbor_type = atom->type;
                pair.neighbor_natom = atom->natom;
                pair.cell_x = atom->cell_x;
                pair.cell_y = atom->cell_y;
                pair.cell_z = atom->cell_z;
                pairs.push_back(pair);
            }
        }
    }
    std::sort(pairs.begin(), pairs.end());
    return pairs;
}

long long elapsed_microseconds(const std::function<void()>& work)
{
    const auto begin = std::chrono::steady_clock::now();
    work();
    const auto end = std::chrono::steady_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
}

void expect_half14_full27_legacy_consistent(const UnitCell& ucell,
                                            const double radius,
                                            const bool pbc,
                                            const std::string& case_name)
{
    std::ofstream ofs_half("synthetic_half14.out");
    Grid grid_half(PARAM.input.test_grid);
    grid_half.init(ofs_half,
                   ucell,
                   radius,
                   pbc,
                   Grid::NeighborBuildMode::AtomPackAndLegacy,
                   Grid::NeighborSearchMode::Half14,
                   Grid::NeighborReferenceMode::Full27);
    ofs_half.close();

    std::ofstream ofs_full("synthetic_full27.out");
    Grid grid_full(PARAM.input.test_grid);
    grid_full.init(ofs_full,
                   ucell,
                   radius,
                   pbc,
                   Grid::NeighborBuildMode::AtomPackAndLegacy,
                   Grid::NeighborSearchMode::Full27,
                   Grid::NeighborReferenceMode::Full27);
    ofs_full.close();

    SCOPED_TRACE(case_name + ", pbc=" + std::to_string(pbc) + ", radius=" + std::to_string(radius));
    EXPECT_EQ(grid_half.neighbor_pairs, grid_half.neighbor_pairs_27);
    EXPECT_EQ(grid_full.neighbor_pairs, grid_full.neighbor_pairs_27);
    EXPECT_EQ(grid_half.neighbor_pairs, grid_full.neighbor_pairs);
    EXPECT_EQ(grid_half.neighbor_pairs, collect_legacy_neighbor_pairs(grid_half));
    EXPECT_EQ(grid_full.neighbor_pairs, collect_legacy_neighbor_pairs(grid_full));

    remove("synthetic_half14.out");
    remove("synthetic_full27.out");
}
} // namespace

TEST_F(SltkGridTest, Init)
{
    ofs.open("test.out");
    unitcell::check_dtau(ucell->atoms,ucell->ntype, ucell->lat0, ucell->latvec);
    test_atom_in = 2;
    PARAM.input.test_grid = 1;
    Grid LatGrid(PARAM.input.test_grid);
    LatGrid.init(ofs, *ucell, radius, pbc);
    EXPECT_FALSE(LatGrid.neighbor_pairs.empty());
    EXPECT_TRUE(LatGrid.neighbor_pairs_27.empty());
    EXPECT_EQ(LatGrid.getGlayerX(), 6);
    EXPECT_EQ(LatGrid.getGlayerY(), 6);
    EXPECT_EQ(LatGrid.getGlayerZ(), 6);
    EXPECT_EQ(LatGrid.getGlayerX_minus(), 5);
    EXPECT_EQ(LatGrid.getGlayerY_minus(), 5);
    EXPECT_EQ(LatGrid.getGlayerZ_minus(), 5);
    ofs.close();
    remove("test.out");
}

TEST_F(SltkGridTest, InitSmall)
{
    ofs.open("test.out");
    unitcell::check_dtau(ucell->atoms,ucell->ntype, ucell->lat0, ucell->latvec);
    test_atom_in = 2;
    PARAM.input.test_grid = 1;
    radius = 0.5;
    Grid LatGrid(PARAM.input.test_grid);
    LatGrid.init(ofs, *ucell, radius, pbc);
    EXPECT_TRUE(LatGrid.neighbor_pairs_27.empty());
    LatGrid.setMemberVariables(ofs,  *ucell);
    EXPECT_EQ(LatGrid.pbc, true);
    EXPECT_TRUE(LatGrid.pbc);
    EXPECT_DOUBLE_EQ(LatGrid.sradius2, radius * radius);
    EXPECT_DOUBLE_EQ(LatGrid.sradius2, 0.5 * 0.5);
    EXPECT_DOUBLE_EQ(LatGrid.sradius, radius);
    EXPECT_DOUBLE_EQ(LatGrid.sradius, 0.5);
    /*
    // minimal value of x, y, z
    EXPECT_DOUBLE_EQ(LatGrid.true_cell_x, 1);
    EXPECT_DOUBLE_EQ(LatGrid.true_cell_y, 1);
    EXPECT_DOUBLE_EQ(LatGrid.true_cell_z, 1);
    // number of cells in x, y, z
    EXPECT_EQ(LatGrid.cell_nx, 3);
    EXPECT_EQ(LatGrid.cell_ny, 3);
    EXPECT_EQ(LatGrid.cell_nz, 3);
    */
    ofs.close();
    remove("test.out");
}

TEST_F(SltkGridTest, OpenMPThreadCountKeepsNeighborSet)
{
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

    std::ofstream ofs_one("test_one.out");
    Grid grid_one(PARAM.input.test_grid);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    grid_one.init(ofs_one,
                  *ucell,
                  radius,
                  pbc,
                  Grid::NeighborBuildMode::AtomPackAndLegacy,
                  Grid::NeighborSearchMode::Half14);
    ofs_one.close();

    std::ofstream ofs_many("test_many.out");
    Grid grid_many(PARAM.input.test_grid);
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    grid_many.init(ofs_many,
                   *ucell,
                   radius,
                   pbc,
                   Grid::NeighborBuildMode::AtomPackAndLegacy,
                   Grid::NeighborSearchMode::Half14);
    ofs_many.close();

    EXPECT_EQ(collect_neighbor_pair_keys(grid_many), collect_neighbor_pair_keys(grid_one));

    remove("test_one.out");
    remove("test_many.out");
}

TEST_F(SltkGridTest, AtomPackNeighborPairsMatchLegacyGridAfterInit)
{
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

    std::ofstream ofs_pbc("test_pair_pbc.out");
    Grid grid_pbc(PARAM.input.test_grid);
    grid_pbc.init(ofs_pbc,
                  *ucell,
                  radius,
                  true,
                  Grid::NeighborBuildMode::AtomPackAndLegacy,
                  Grid::NeighborSearchMode::Half14,
                  Grid::NeighborReferenceMode::Full27);
    ofs_pbc.close();
    EXPECT_EQ(grid_pbc.neighbor_pairs, grid_pbc.neighbor_pairs_27);
    EXPECT_EQ(grid_pbc.neighbor_pairs_27, collect_legacy_neighbor_pairs(grid_pbc));

    std::ofstream ofs_non_pbc("test_pair_non_pbc.out");
    Grid grid_non_pbc(PARAM.input.test_grid);
    grid_non_pbc.init(ofs_non_pbc,
                      *ucell,
                      0.5,
                      false,
                      Grid::NeighborBuildMode::AtomPackAndLegacy,
                      Grid::NeighborSearchMode::Half14,
                      Grid::NeighborReferenceMode::Full27);
    ofs_non_pbc.close();
    EXPECT_EQ(grid_non_pbc.neighbor_pairs, grid_non_pbc.neighbor_pairs_27);
    EXPECT_EQ(grid_non_pbc.neighbor_pairs_27, collect_legacy_neighbor_pairs(grid_non_pbc));

    remove("test_pair_pbc.out");
    remove("test_pair_non_pbc.out");
}

TEST_F(SltkGridTest, Half14DefaultMatchesFull27AndLegacyGrid)
{
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

    std::ofstream ofs_half("test_pair_half14.out");
    Grid grid_half(PARAM.input.test_grid);
    grid_half.init(ofs_half,
                   *ucell,
                   radius,
                   true,
                   Grid::NeighborBuildMode::AtomPackAndLegacy,
                   Grid::NeighborSearchMode::Half14,
                   Grid::NeighborReferenceMode::Full27);
    ofs_half.close();

    std::ofstream ofs_full("test_pair_full27.out");
    Grid grid_full(PARAM.input.test_grid);
    grid_full.init(ofs_full,
                   *ucell,
                   radius,
                   true,
                   Grid::NeighborBuildMode::AtomPackAndLegacy,
                   Grid::NeighborSearchMode::Full27,
                   Grid::NeighborReferenceMode::Full27);
    ofs_full.close();

    EXPECT_EQ(grid_half.neighbor_pairs, grid_half.neighbor_pairs_27);
    EXPECT_EQ(grid_full.neighbor_pairs, grid_full.neighbor_pairs_27);
    EXPECT_EQ(grid_half.neighbor_pairs, grid_full.neighbor_pairs);
    EXPECT_EQ(grid_half.neighbor_pairs, collect_legacy_neighbor_pairs(grid_half));

    std::ofstream ofs_atom_pack_only("test_pair_half14_atom_pack_only.out");
    Grid grid_atom_pack_only(PARAM.input.test_grid);
    grid_atom_pack_only.init(ofs_atom_pack_only,
                             *ucell,
                             radius,
                             true,
                             Grid::NeighborBuildMode::AtomPackOnly,
                             Grid::NeighborSearchMode::Half14);
    ofs_atom_pack_only.close();

    EXPECT_TRUE(grid_atom_pack_only.atoms_in_box.empty());
    EXPECT_TRUE(grid_atom_pack_only.all_adj_info.empty());
    EXPECT_TRUE(grid_atom_pack_only.neighbor_pairs_27.empty());
    EXPECT_FALSE(grid_atom_pack_only.neighbor_pairs.empty());
    EXPECT_FALSE(grid_atom_pack_only.neighbor_pair_indices.empty());

    remove("test_pair_half14.out");
    remove("test_pair_full27.out");
    remove("test_pair_half14_atom_pack_only.out");
}

TEST(SltkGridSyntheticTest, SyntheticCellsKeepHalf14Full27AndLegacyConsistent)
{
    SetGlobalV();
    for (SyntheticNeighborCase test_case: make_synthetic_neighbor_cases())
    {
        UnitCell* ucell = test_case.prepare.SetUcellInfo();
        unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

        for (const double radius: test_case.radii)
        {
            expect_half14_full27_legacy_consistent(*ucell, radius, true, test_case.name);
            expect_half14_full27_legacy_consistent(*ucell, radius, false, test_case.name);
        }

        delete ucell;
    }
}

TEST(SltkGridSyntheticTest, Half14Full27LightweightMetrics)
{
    SetGlobalV();
    std::cout << "[neighbor-metrics] case,pbc,radius,half_pairs,full_pairs,pack_atoms,boxes,half14_us,"
                 "half14_with_ref_us,full27_us,half14_with_legacy_us\n";

    for (SyntheticNeighborCase test_case: make_synthetic_neighbor_cases())
    {
        UnitCell* ucell = test_case.prepare.SetUcellInfo();
        unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

        for (const double radius: test_case.radii)
        {
            for (const bool pbc_case: std::vector<bool>{true, false})
            {
                Grid grid_half_only(PARAM.input.test_grid);
                Grid grid_half_with_reference(PARAM.input.test_grid);
                Grid grid_full_only(PARAM.input.test_grid);
                Grid grid_half_with_legacy(PARAM.input.test_grid);

                std::ofstream ofs_half("metric_half14.out");
                const long long half_us = elapsed_microseconds([&]() {
                    grid_half_only.init(ofs_half,
                                        *ucell,
                                        radius,
                                        pbc_case,
                                        Grid::NeighborBuildMode::AtomPackOnly,
                                        Grid::NeighborSearchMode::Half14);
                });
                ofs_half.close();

                std::ofstream ofs_half_ref("metric_half14_with_reference.out");
                const long long half_ref_us = elapsed_microseconds([&]() {
                    grid_half_with_reference.init(ofs_half_ref,
                                                  *ucell,
                                                  radius,
                                                  pbc_case,
                                                  Grid::NeighborBuildMode::AtomPackOnly,
                                                  Grid::NeighborSearchMode::Half14,
                                                  Grid::NeighborReferenceMode::Full27);
                });
                ofs_half_ref.close();

                std::ofstream ofs_full("metric_full27.out");
                const long long full_us = elapsed_microseconds([&]() {
                    grid_full_only.init(ofs_full,
                                        *ucell,
                                        radius,
                                        pbc_case,
                                        Grid::NeighborBuildMode::AtomPackOnly,
                                        Grid::NeighborSearchMode::Full27);
                });
                ofs_full.close();

                std::ofstream ofs_legacy("metric_half14_with_legacy.out");
                const long long with_legacy_us = elapsed_microseconds([&]() {
                    grid_half_with_legacy.init(ofs_legacy,
                                               *ucell,
                                               radius,
                                               pbc_case,
                                               Grid::NeighborBuildMode::AtomPackAndLegacy,
                                               Grid::NeighborSearchMode::Half14,
                                               Grid::NeighborReferenceMode::Full27);
                });
                ofs_legacy.close();

                SCOPED_TRACE(test_case.name + ", pbc=" + std::to_string(pbc_case)
                             + ", radius=" + std::to_string(radius));
                EXPECT_TRUE(grid_half_only.neighbor_pairs_27.empty());
                EXPECT_TRUE(grid_full_only.neighbor_pairs_27.empty());
                EXPECT_EQ(grid_half_with_reference.neighbor_pairs, grid_half_with_reference.neighbor_pairs_27);
                EXPECT_EQ(grid_half_only.neighbor_pairs, grid_full_only.neighbor_pairs);
                EXPECT_EQ(grid_half_only.neighbor_pairs, grid_half_with_reference.neighbor_pairs);
                EXPECT_EQ(grid_half_only.neighbor_pairs, grid_half_with_legacy.neighbor_pairs);

                std::cout << "[neighbor-metrics] " << test_case.name << "," << pbc_case << "," << radius << ","
                          << grid_half_only.neighbor_pairs.size() << "," << grid_full_only.neighbor_pairs.size()
                          << "," << grid_half_only.atom_pack.size() << "," << grid_half_only.grid_storage.box_size()
                          << "," << half_us << "," << half_ref_us << "," << full_us << "," << with_legacy_us
                          << "\n";

                remove("metric_half14.out");
                remove("metric_half14_with_reference.out");
                remove("metric_full27.out");
                remove("metric_half14_with_legacy.out");
            }
        }

        delete ucell;
    }
}

/*
// This test cannot pass because setAtomLinkArray() is unsuccessful
// if expand_flag is false
TEST_F(SltkGridTest, InitNoExpand)
{
    ofs.open("test.out");
    unitcell::check_dtau(ucell->atoms,ucell->ntype, ucell->lat0, ucell->latvec);
    test_atom_in = 2;
    PARAM.input.test_grid = 1;
    double radius = 1e-1000;
    Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in);
    Grid LatGrid(PARAM.input.test_grid);
    LatGrid.init(ofs, *ucell, Atom_inp);
    ofs.close();
}
*/
