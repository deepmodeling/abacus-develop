#include "gmock/gmock.h"
#include "gtest/gtest.h"

#define private public
#include "../sltk_grid.h"
#include "prepare_unitcell.h"
#include "source_io/module_parameter/parameter.h"
#undef private
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>
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
    PARAM.input.test_deconstructor = 0;
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
void SetLattice(UnitCell* ucell,
                const double lat0,
                const ModuleBase::Matrix3& latvec)
{
    ucell->lat0 = lat0;
    ucell->lat0_angstrom = lat0 * ModuleBase::BOHR_TO_A;
    ucell->tpiba = ModuleBase::TWO_PI / lat0;
    ucell->tpiba2 = ucell->tpiba * ucell->tpiba;
    ucell->latvec = latvec;
    ucell->a1.set(latvec.e11, latvec.e12, latvec.e13);
    ucell->a2.set(latvec.e21, latvec.e22, latvec.e23);
    ucell->a3.set(latvec.e31, latvec.e32, latvec.e33);
    ucell->GT = latvec.Inverse();
    ucell->G = ucell->GT.Transpose();
    ucell->GGT = ucell->G * ucell->GT;
    ucell->invGGT = ucell->GGT.Inverse();
    ucell->omega = std::abs(latvec.Det()) * lat0 * lat0 * lat0;
}

ModuleBase::Matrix3 IdentityLatvec()
{
    ModuleBase::Matrix3 latvec;
    latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
    latvec.e21 = 0.0; latvec.e22 = 1.0; latvec.e23 = 0.0;
    latvec.e31 = 0.0; latvec.e32 = 0.0; latvec.e33 = 1.0;
    return latvec;
}

ModuleBase::Matrix3 SkewedLatvec()
{
    ModuleBase::Matrix3 latvec;
    latvec.e11 = 1.0; latvec.e12 = 0.3; latvec.e13 = 0.0;
    latvec.e21 = 0.1; latvec.e22 = 1.0; latvec.e23 = 0.0;
    latvec.e31 = 0.0; latvec.e32 = 0.0; latvec.e33 = 1.0;
    return latvec;
}

void SetTwoAtomPositions(UnitCell* ucell,
                         const ModuleBase::Vector3<double>& tau0,
                         const ModuleBase::Vector3<double>& tau1)
{
    ASSERT_EQ(ucell->ntype, 1);
    ASSERT_GE(ucell->atoms[0].na, 2);
    ucell->atoms[0].tau[0] = tau0;
    ucell->atoms[0].tau[1] = tau1;
}


std::size_t CountExpandedAtomsInBoxes(const Grid& grid)
{
    std::size_t count = 0;
    for (const auto& boxes_yz: grid.atoms_in_box)
    {
        for (const auto& boxes_z: boxes_yz)
        {
            for (const auto& atoms: boxes_z)
            {
                count += atoms.size();
            }
        }
    }
    return count;
}

std::size_t ExpectedExpandedAtomCount(const Grid& grid, const UnitCell& ucell)
{
    const std::size_t layer_x = static_cast<std::size_t>(grid.getGlayerX() + grid.getGlayerX_minus());
    const std::size_t layer_y = static_cast<std::size_t>(grid.getGlayerY() + grid.getGlayerY_minus());
    const std::size_t layer_z = static_cast<std::size_t>(grid.getGlayerZ() + grid.getGlayerZ_minus());
    return static_cast<std::size_t>(ucell.nat) * layer_x * layer_y * layer_z;
}

std::size_t SaturatedAllocatedBoxCount(const Grid& grid)
{
    if (grid.box_nx <= 0 || grid.box_ny <= 0 || grid.box_nz <= 0)
    {
        return 0;
    }
    const std::size_t sx = static_cast<std::size_t>(grid.box_nx);
    const std::size_t sy = static_cast<std::size_t>(grid.box_ny);
    const std::size_t sz = static_cast<std::size_t>(grid.box_nz);
    const std::size_t max_size = std::numeric_limits<std::size_t>::max();
    if (sx > max_size / sy)
    {
        return max_size;
    }
    const std::size_t xy = sx * sy;
    if (xy > max_size / sz)
    {
        return max_size;
    }
    return xy * sz;
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
    EXPECT_GE(LatGrid.box_edge_length, LatGrid.sradius);
    EXPECT_GT(LatGrid.box_nx, 0);
    EXPECT_GT(LatGrid.box_ny, 0);
    EXPECT_GT(LatGrid.box_nz, 0);
    EXPECT_EQ(LatGrid.atoms_in_box.size(), static_cast<size_t>(LatGrid.box_nx));
    EXPECT_EQ(CountExpandedAtomsInBoxes(LatGrid), ExpectedExpandedAtomCount(LatGrid, *ucell));
    EXPECT_EQ(LatGrid.all_adj_info.size(), static_cast<size_t>(ucell->ntype));
    ASSERT_EQ(LatGrid.all_adj_info[0].size(), static_cast<size_t>(ucell->atoms[0].na));
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
    EXPECT_GT(LatGrid.box_nx, 0);
    EXPECT_GT(LatGrid.box_ny, 0);
    EXPECT_GT(LatGrid.box_nz, 0);
    EXPECT_EQ(LatGrid.atoms_in_box.size(), static_cast<size_t>(LatGrid.box_nx));
    ASSERT_FALSE(LatGrid.atoms_in_box.empty());
    EXPECT_EQ(LatGrid.atoms_in_box[0].size(), static_cast<size_t>(LatGrid.box_ny));
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

TEST_F(SltkGridTest, InitOrthogonalCellSetsExpectedLayersAndBoxCounts)
{
    SetLattice(ucell, 1.0, IdentityLatvec());
    SetTwoAtomPositions(ucell,
                        ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                        ModuleBase::Vector3<double>(0.5, 0.0, 0.0));

    ofs.open("orthogonal_grid.out");
    Grid lat_grid(PARAM.input.test_grid);
    lat_grid.init(ofs, *ucell, 1.0, true);

    EXPECT_EQ(lat_grid.getGlayerX(), 2);
    EXPECT_EQ(lat_grid.getGlayerY(), 2);
    EXPECT_EQ(lat_grid.getGlayerZ(), 2);
    EXPECT_EQ(lat_grid.getGlayerX_minus(), 1);
    EXPECT_EQ(lat_grid.getGlayerY_minus(), 1);
    EXPECT_EQ(lat_grid.getGlayerZ_minus(), 1);
    EXPECT_EQ(lat_grid.box_nx, 3);
    EXPECT_EQ(lat_grid.box_ny, 3);
    EXPECT_EQ(lat_grid.box_nz, 3);
    EXPECT_DOUBLE_EQ(lat_grid.sradius, 1.0);
    EXPECT_DOUBLE_EQ(lat_grid.sradius2, 1.0);
    EXPECT_EQ(lat_grid.atoms_in_box.size(), 3u);

    ofs.close();
    remove("orthogonal_grid.out");
}

TEST_F(SltkGridTest, InitNonOrthogonalCellBuildsPositiveExpansion)
{
    SetLattice(ucell, 1.0, SkewedLatvec());
    SetTwoAtomPositions(ucell,
                        ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                        ModuleBase::Vector3<double>(0.45, 0.20, 0.0));

    ofs.open("skewed_grid.out");
    Grid lat_grid(PARAM.input.test_grid);
    lat_grid.init(ofs, *ucell, 0.75, true);

    EXPECT_GE(lat_grid.getGlayerX(), 1);
    EXPECT_GE(lat_grid.getGlayerY(), 1);
    EXPECT_GE(lat_grid.getGlayerZ(), 1);
    EXPECT_GE(lat_grid.box_edge_length, lat_grid.sradius);
    EXPECT_GT(lat_grid.box_nx, 0);
    EXPECT_GT(lat_grid.box_ny, 0);
    EXPECT_GT(lat_grid.box_nz, 0);
    EXPECT_GT(lat_grid.atoms_in_box.size(), 0u);
    EXPECT_EQ(CountExpandedAtomsInBoxes(lat_grid), ExpectedExpandedAtomCount(lat_grid, *ucell));
    EXPECT_EQ(lat_grid.all_adj_info.size(), static_cast<size_t>(ucell->ntype));

    ofs.close();
    remove("skewed_grid.out");
}

TEST_F(SltkGridTest, InitWithPbcDisabledKeepsNoPeriodicExpansion)
{
    SetLattice(ucell, 1.0, IdentityLatvec());
    SetTwoAtomPositions(ucell,
                        ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                        ModuleBase::Vector3<double>(0.5, 0.0, 0.0));

    ofs.open("no_pbc_grid.out");
    Grid lat_grid(PARAM.input.test_grid);
    lat_grid.init(ofs, *ucell, 1.0, false);

    EXPECT_FALSE(lat_grid.pbc);
    EXPECT_EQ(lat_grid.getGlayerX(), 0);
    EXPECT_EQ(lat_grid.getGlayerY(), 0);
    EXPECT_EQ(lat_grid.getGlayerZ(), 0);
    EXPECT_EQ(lat_grid.getGlayerX_minus(), 0);
    EXPECT_EQ(lat_grid.getGlayerY_minus(), 0);
    EXPECT_EQ(lat_grid.getGlayerZ_minus(), 0);
    EXPECT_EQ(lat_grid.box_nx, 0);
    EXPECT_EQ(lat_grid.box_ny, 0);
    EXPECT_EQ(lat_grid.box_nz, 0);
    EXPECT_EQ(lat_grid.all_adj_info.size(), static_cast<size_t>(ucell->ntype));

    ofs.close();
    remove("no_pbc_grid.out");
}

TEST_F(SltkGridTest, CellListCoarsensBinsForSparseHugeAabb)
{
    ModuleBase::Matrix3 large_latvec;
    large_latvec.e11 = 100000.0; large_latvec.e12 = 0.0; large_latvec.e13 = 0.0;
    large_latvec.e21 = 0.0; large_latvec.e22 = 100000.0; large_latvec.e23 = 0.0;
    large_latvec.e31 = 0.0; large_latvec.e32 = 0.0; large_latvec.e33 = 100000.0;
    SetLattice(ucell, 1.0, large_latvec);
    SetTwoAtomPositions(ucell,
                        ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                        ModuleBase::Vector3<double>(0.5, 0.5, 0.5));

    ofs.open("sparse_huge_aabb_grid.out");
    Grid lat_grid(PARAM.input.test_grid);
    lat_grid.init(ofs, *ucell, 0.01, true);

    EXPECT_GT(lat_grid.box_edge_length, 0.1);
    EXPECT_GT(lat_grid.box_nx, 0);
    EXPECT_GT(lat_grid.box_ny, 0);
    EXPECT_GT(lat_grid.box_nz, 0);
    EXPECT_LE(SaturatedAllocatedBoxCount(lat_grid), 1000000u);
    EXPECT_EQ(CountExpandedAtomsInBoxes(lat_grid), ExpectedExpandedAtomCount(lat_grid, *ucell));

    ofs.close();
    remove("sparse_huge_aabb_grid.out");
}

TEST_F(SltkGridTest, ReinitWithPbcDisabledResetsExpansionState)
{
    SetLattice(ucell, 1.0, IdentityLatvec());
    SetTwoAtomPositions(ucell,
                        ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                        ModuleBase::Vector3<double>(0.5, 0.0, 0.0));

    Grid lat_grid(PARAM.input.test_grid);
    ofs.open("reinit_pbc_grid.out");
    lat_grid.init(ofs, *ucell, 1.0, true);
    ofs.close();
    EXPECT_GT(lat_grid.getGlayerX(), 0);
    EXPECT_GT(lat_grid.box_nx, 0);

    ofs.open("reinit_no_pbc_grid.out");
    lat_grid.init(ofs, *ucell, 1.0, false);
    ofs.close();

    EXPECT_FALSE(lat_grid.pbc);
    EXPECT_EQ(lat_grid.getGlayerX(), 0);
    EXPECT_EQ(lat_grid.getGlayerY(), 0);
    EXPECT_EQ(lat_grid.getGlayerZ(), 0);
    EXPECT_EQ(lat_grid.getGlayerX_minus(), 0);
    EXPECT_EQ(lat_grid.getGlayerY_minus(), 0);
    EXPECT_EQ(lat_grid.getGlayerZ_minus(), 0);
    EXPECT_EQ(lat_grid.box_nx, 0);
    EXPECT_EQ(lat_grid.box_ny, 0);
    EXPECT_EQ(lat_grid.box_nz, 0);
    EXPECT_EQ(CountExpandedAtomsInBoxes(lat_grid), 0u);

    remove("reinit_pbc_grid.out");
    remove("reinit_no_pbc_grid.out");
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
