#include "../sltk_atom_arrange.h"

#define private public
#include "source_io/module_parameter/parameter.h"
#undef private
#include <algorithm>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "prepare_unitcell.h"
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

namespace
{
using AdjacentKey = std::tuple<int, int, int, int, int, double, double, double>;

std::vector<AdjacentKey> collect_adjacent_keys(const AdjacentAtomInfo& adjs)
{
    std::vector<AdjacentKey> keys;
    for (int i = 0; i < static_cast<int>(adjs.ntype.size()); ++i)
    {
        keys.push_back(AdjacentKey(adjs.ntype[i],
                                   adjs.natom[i],
                                   adjs.box[i].x,
                                   adjs.box[i].y,
                                   adjs.box[i].z,
                                   adjs.adjacent_tau[i].x,
                                   adjs.adjacent_tau[i].y,
                                   adjs.adjacent_tau[i].z));
    }
    std::sort(keys.begin(), keys.end());
    return keys;
}

void expect_self_is_last(const AdjacentAtomInfo& adjs,
                         const UnitCell& ucell,
                         const int ntype,
                         const int nnumber)
{
    ASSERT_FALSE(adjs.ntype.empty());
    const int last = static_cast<int>(adjs.ntype.size()) - 1;
    EXPECT_EQ(adjs.ntype[last], ntype);
    EXPECT_EQ(adjs.natom[last], nnumber);
    EXPECT_EQ(adjs.box[last].x, 0);
    EXPECT_EQ(adjs.box[last].y, 0);
    EXPECT_EQ(adjs.box[last].z, 0);
    EXPECT_DOUBLE_EQ(adjs.adjacent_tau[last].x, ucell.atoms[ntype].tau[nnumber].x);
    EXPECT_DOUBLE_EQ(adjs.adjacent_tau[last].y, ucell.atoms[ntype].tau[nnumber].y);
    EXPECT_DOUBLE_EQ(adjs.adjacent_tau[last].z, ucell.atoms[ntype].tau[nnumber].z);
}

void expect_grid_driver_default_path_matches_legacy(const UnitCell& ucell, const double radius, const bool pbc)
{
    std::ofstream ofs("grid_driver_atom_pack_compare.out");
    Grid_Driver grid_d(PARAM.input.test_deconstructor, PARAM.input.test_grid);
    grid_d.init(ofs, ucell, radius, pbc);
    ofs.close();

    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            SCOPED_TRACE("pbc=" + std::to_string(pbc) + ", radius=" + std::to_string(radius) + ", type="
                         + std::to_string(it) + ", atom=" + std::to_string(ia));

            AdjacentAtomInfo default_adjs;
            AdjacentAtomInfo legacy_adjs;
            AdjacentAtomInfo atom_pack_adjs;
            AdjacentAtomInfo deprecated_adjs;
            grid_d.Find_atom(ucell, it, ia, &default_adjs);
            grid_d.Find_atom_from_legacy(ucell, it, ia, &legacy_adjs);
            grid_d.Find_atom_from_atom_pack(ucell, it, ia, &atom_pack_adjs);
            grid_d.Find_atom(ucell, ucell.atoms[it].tau[ia], it, ia, &deprecated_adjs);

            EXPECT_EQ(default_adjs.adj_num, legacy_adjs.adj_num);
            EXPECT_EQ(default_adjs.ntype.size(), legacy_adjs.ntype.size());
            EXPECT_EQ(default_adjs.natom.size(), legacy_adjs.natom.size());
            EXPECT_EQ(default_adjs.box.size(), legacy_adjs.box.size());
            EXPECT_EQ(default_adjs.adjacent_tau.size(), legacy_adjs.adjacent_tau.size());
            EXPECT_EQ(collect_adjacent_keys(default_adjs), collect_adjacent_keys(legacy_adjs));
            EXPECT_EQ(atom_pack_adjs.adj_num, legacy_adjs.adj_num);
            EXPECT_EQ(atom_pack_adjs.ntype.size(), legacy_adjs.ntype.size());
            EXPECT_EQ(atom_pack_adjs.natom.size(), legacy_adjs.natom.size());
            EXPECT_EQ(atom_pack_adjs.box.size(), legacy_adjs.box.size());
            EXPECT_EQ(atom_pack_adjs.adjacent_tau.size(), legacy_adjs.adjacent_tau.size());
            EXPECT_EQ(collect_adjacent_keys(atom_pack_adjs), collect_adjacent_keys(legacy_adjs));
            EXPECT_EQ(deprecated_adjs.adj_num, legacy_adjs.adj_num);
            EXPECT_EQ(deprecated_adjs.ntype.size(), legacy_adjs.ntype.size());
            EXPECT_EQ(deprecated_adjs.natom.size(), legacy_adjs.natom.size());
            EXPECT_EQ(deprecated_adjs.box.size(), legacy_adjs.box.size());
            EXPECT_EQ(deprecated_adjs.adjacent_tau.size(), legacy_adjs.adjacent_tau.size());
            EXPECT_EQ(collect_adjacent_keys(deprecated_adjs), collect_adjacent_keys(legacy_adjs));
            expect_self_is_last(default_adjs, ucell, it, ia);
            expect_self_is_last(legacy_adjs, ucell, it, ia);
            expect_self_is_last(atom_pack_adjs, ucell, it, ia);
            expect_self_is_last(deprecated_adjs, ucell, it, ia);
        }
    }

    remove("grid_driver_atom_pack_compare.out");
}

void expect_grid_driver_atom_pack_only_path(const UnitCell& ucell, const double radius, const bool pbc)
{
    std::ofstream ofs("grid_driver_atom_pack_only.out");
    Grid_Driver grid_d(PARAM.input.test_deconstructor, PARAM.input.test_grid);
    grid_d.init(ofs, ucell, radius, pbc, Grid::NeighborBuildMode::AtomPackOnly);
    ofs.close();

    EXPECT_TRUE(grid_d.atoms_in_box.empty());
    EXPECT_TRUE(grid_d.all_adj_info.empty());
    EXPECT_FALSE(grid_d.atom_pack.empty());
    EXPECT_FALSE(grid_d.neighbor_pair_indices.empty());

    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            SCOPED_TRACE("AtomPackOnly pbc=" + std::to_string(pbc) + ", radius=" + std::to_string(radius)
                         + ", type=" + std::to_string(it) + ", atom=" + std::to_string(ia));

            AdjacentAtomInfo default_adjs;
            AdjacentAtomInfo atom_pack_adjs;
            AdjacentAtomInfo deprecated_adjs;
            grid_d.Find_atom(ucell, it, ia, &default_adjs);
            grid_d.Find_atom_from_atom_pack(ucell, it, ia, &atom_pack_adjs);
            grid_d.Find_atom(ucell, ucell.atoms[it].tau[ia], it, ia, &deprecated_adjs);

            EXPECT_EQ(default_adjs.adj_num, atom_pack_adjs.adj_num);
            EXPECT_EQ(default_adjs.ntype.size(), atom_pack_adjs.ntype.size());
            EXPECT_EQ(default_adjs.natom.size(), atom_pack_adjs.natom.size());
            EXPECT_EQ(default_adjs.box.size(), atom_pack_adjs.box.size());
            EXPECT_EQ(default_adjs.adjacent_tau.size(), atom_pack_adjs.adjacent_tau.size());
            EXPECT_EQ(collect_adjacent_keys(default_adjs), collect_adjacent_keys(atom_pack_adjs));
            EXPECT_EQ(collect_adjacent_keys(deprecated_adjs), collect_adjacent_keys(atom_pack_adjs));
            expect_self_is_last(default_adjs, ucell, it, ia);
            expect_self_is_last(atom_pack_adjs, ucell, it, ia);
            expect_self_is_last(deprecated_adjs, ucell, it, ia);
            EXPECT_THROW(grid_d.Find_atom_from_legacy(ucell, it, ia, nullptr), std::runtime_error);
        }
    }

    remove("grid_driver_atom_pack_only.out");
}
} // namespace

/************************************************
 *  unit test of atom_arrange
 ***********************************************/

/**
 * - Tested Functions:
 *   - atom_arrange::delete_vector(void)
 *     - delete vector
 *   - atom_arrange::set_sr_NL
 * 	   - set the sr: search radius including nonlocal beta
 *   - filter_adjs function
 *     - filter AdjacentAtomInfo to the minimized adjacent atoms
 */

void SetGlobalV()
{
    PARAM.input.test_grid = false;
}

class SltkAtomArrangeTest : public testing::Test
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

TEST_F(SltkAtomArrangeTest, setsrNL)
{
    atom_arrange test;
    const std::string teststring = "m";
    double rcutmax_Phi = 1;
    double rcutmax_Beta = 2;
    bool gamma_only_local = true;
    double test_sr = 0;
    std::ofstream ofs;
    ofs.open("./to_test_arrange.txt");
    test_sr = test.set_sr_NL(ofs, teststring, rcutmax_Phi, rcutmax_Beta, gamma_only_local);
    EXPECT_DOUBLE_EQ(test_sr, 2.001);

    gamma_only_local = false;
    test_sr = test.set_sr_NL(ofs, teststring, rcutmax_Phi, rcutmax_Beta, gamma_only_local);
    EXPECT_DOUBLE_EQ(test_sr, 6.001);

    const std::string teststring2 = "no";
    test_sr = test.set_sr_NL(ofs, teststring2, rcutmax_Phi, rcutmax_Beta, gamma_only_local);
    ofs.close();
    std::ifstream ifs;
    std::string test2, s;
    ifs.open("./to_test_arrange.txt");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("Orbital max radius cutoff (Bohr) = 1"));
    EXPECT_THAT(str, testing::HasSubstr("Nonlocal proj. max radius cutoff (Bohr) = 2"));
    ifs.close();
    //remove("./to_test_arrange");
}

TEST_F(SltkAtomArrangeTest, Search)
{
    unitcell::check_dtau(ucell->atoms,ucell->ntype, ucell->lat0, ucell->latvec);
    Grid_Driver grid_d(PARAM.input.test_deconstructor, PARAM.input.test_grid);
    ofs.open("test.out");
    bool test_only = true;
    atom_arrange::search(pbc, ofs, grid_d, *ucell, radius, test_atom_in, test_only);
    EXPECT_EQ(grid_d.getType(0),0);
    EXPECT_EQ(grid_d.getNatom(0), 1); // adjacent atom is 1
    ofs.close();
    ifs.open("test.out");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("search neighboring atoms done."));
    remove("test.out");
}

TEST_F(SltkAtomArrangeTest, Filteradjs)
{
    unitcell::check_dtau(ucell->atoms,ucell->ntype, ucell->lat0, ucell->latvec);
    Grid_Driver grid_d(PARAM.input.test_deconstructor, PARAM.input.test_grid);
    ofs.open("test.out");
    bool test_only = true;
    atom_arrange::search(pbc, ofs, grid_d, *ucell, radius, test_atom_in, test_only);
    EXPECT_EQ(grid_d.getType(0),0);
    EXPECT_EQ(grid_d.getNatom(0), 1); // adjacent atom is 1
    ofs.close();
    ifs.open("test.out");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("search neighboring atoms done."));
    remove("test.out");

    AdjacentAtomInfo adjs;
    grid_d.Find_atom(*ucell, ucell->atoms[0].tau[0], 0, 0, &adjs);
    EXPECT_EQ(adjs.adj_num, 0);
    // add one adjacent atom
    adjs.adj_num++;
    adjs.adjacent_tau.push_back(ModuleBase::Vector3<double>(0,0,0));
    adjs.box.push_back(ModuleBase::Vector3<int>(0,0,0));
    adjs.natom.push_back(1);
    adjs.ntype.push_back(0);
    EXPECT_EQ(adjs.adj_num, 1);
    // filter adjs to no adjacent status
    std::vector<bool> is_adjs(adjs.adj_num + 1, false);
    is_adjs[0] = true;
    filter_adjs(is_adjs, adjs);
    EXPECT_EQ(adjs.adj_num, 0);
}

TEST_F(SltkAtomArrangeTest, GridDriverAtomPackPathMatchesLegacyPath)
{
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

    expect_grid_driver_default_path_matches_legacy(*ucell, radius, true);
    expect_grid_driver_default_path_matches_legacy(*ucell, 0.5, false);
}

TEST_F(SltkAtomArrangeTest, GridDriverAtomPackOnlyPathWorksWithoutLegacyBuild)
{
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

    expect_grid_driver_atom_pack_only_path(*ucell, radius, true);
    expect_grid_driver_atom_pack_only_path(*ucell, 0.5, false);
}
