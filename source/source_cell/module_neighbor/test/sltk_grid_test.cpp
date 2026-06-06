#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <tuple>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#define private public
#include "../sltk_grid.h"
#include "prepare_unitcell.h"
#include "source_io/module_parameter/parameter.h"
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
} // namespace

TEST_F(SltkGridTest, Init)
{
    ofs.open("test.out");
    unitcell::check_dtau(ucell->atoms,ucell->ntype, ucell->lat0, ucell->latvec);
    test_atom_in = 2;
    PARAM.input.test_grid = 1;
    Grid LatGrid(PARAM.input.test_grid);
    LatGrid.init(ofs, *ucell, radius, pbc);
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
    grid_one.init(ofs_one, *ucell, radius, pbc);
    ofs_one.close();

    std::ofstream ofs_many("test_many.out");
    Grid grid_many(PARAM.input.test_grid);
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    grid_many.init(ofs_many, *ucell, radius, pbc);
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
    grid_pbc.init(ofs_pbc, *ucell, radius, true);
    ofs_pbc.close();
    EXPECT_EQ(grid_pbc.neighbor_pairs_27, collect_legacy_neighbor_pairs(grid_pbc));

    std::ofstream ofs_non_pbc("test_pair_non_pbc.out");
    Grid grid_non_pbc(PARAM.input.test_grid);
    grid_non_pbc.init(ofs_non_pbc, *ucell, 0.5, false);
    ofs_non_pbc.close();
    EXPECT_EQ(grid_non_pbc.neighbor_pairs_27, collect_legacy_neighbor_pairs(grid_non_pbc));

    remove("test_pair_pbc.out");
    remove("test_pair_non_pbc.out");
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
