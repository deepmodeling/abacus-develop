#include "../sltk_atom_arrange.h"
#include "../neighbor_rebuild_state.h"

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
#include "synthetic_neighbor_unitcell.h"
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
    grid_d.init(ofs,
                ucell,
                radius,
                pbc,
                Grid::NeighborBuildMode::AtomPackAndLegacy,
                Grid::NeighborSearchMode::Half14,
                Grid::NeighborReferenceMode::Full27);
    ofs.close();
    EXPECT_EQ(grid_d.neighbor_pairs, grid_d.neighbor_pairs_27);

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
    EXPECT_TRUE(grid_d.neighbor_pairs_27.empty());
    EXPECT_FALSE(grid_d.paged_neighbor_list.type_offset.empty());
    EXPECT_EQ(grid_d.paged_neighbor_list.total_neighbors(),
              static_cast<int>(grid_d.neighbor_pairs.size()));

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

void expect_atom_pack_only_half14_matches_full27_reference(const UnitCell& ucell,
                                                           const double radius,
                                                           const bool pbc,
                                                           const std::string& case_name)
{
    std::ofstream ofs_candidate("grid_driver_synthetic_half14_atom_pack_only.out");
    Grid_Driver candidate(PARAM.input.test_deconstructor, PARAM.input.test_grid);
    candidate.init(ofs_candidate,
                   ucell,
                   radius,
                   pbc,
                   Grid::NeighborBuildMode::AtomPackOnly,
                   Grid::NeighborSearchMode::Half14);
    ofs_candidate.close();

    std::ofstream ofs_reference("grid_driver_synthetic_full27_reference.out");
    Grid_Driver reference(PARAM.input.test_deconstructor, PARAM.input.test_grid);
    reference.init(ofs_reference,
                   ucell,
                   radius,
                   pbc,
                   Grid::NeighborBuildMode::AtomPackAndLegacy,
                   Grid::NeighborSearchMode::Full27,
                   Grid::NeighborReferenceMode::Full27);
    ofs_reference.close();

    EXPECT_TRUE(candidate.atoms_in_box.empty());
    EXPECT_TRUE(candidate.all_adj_info.empty());
    EXPECT_TRUE(candidate.neighbor_pairs_27.empty());
    EXPECT_FALSE(candidate.neighbor_pairs.empty());
    EXPECT_EQ(reference.neighbor_pairs, reference.neighbor_pairs_27);
    EXPECT_EQ(candidate.neighbor_pairs, reference.neighbor_pairs);

    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            SCOPED_TRACE(case_name + ", pbc=" + std::to_string(pbc) + ", radius=" + std::to_string(radius)
                         + ", type=" + std::to_string(it) + ", atom=" + std::to_string(ia));

            AdjacentAtomInfo candidate_adjs;
            AdjacentAtomInfo reference_adjs;
            AdjacentAtomInfo legacy_adjs;
            candidate.Find_atom(ucell, it, ia, &candidate_adjs);
            reference.Find_atom(ucell, it, ia, &reference_adjs);
            reference.Find_atom_from_legacy(ucell, it, ia, &legacy_adjs);

            EXPECT_EQ(candidate_adjs.adj_num, reference_adjs.adj_num);
            EXPECT_EQ(candidate_adjs.ntype.size(), reference_adjs.ntype.size());
            EXPECT_EQ(candidate_adjs.natom.size(), reference_adjs.natom.size());
            EXPECT_EQ(candidate_adjs.box.size(), reference_adjs.box.size());
            EXPECT_EQ(candidate_adjs.adjacent_tau.size(), reference_adjs.adjacent_tau.size());
            EXPECT_EQ(collect_adjacent_keys(candidate_adjs), collect_adjacent_keys(reference_adjs));
            EXPECT_EQ(collect_adjacent_keys(candidate_adjs), collect_adjacent_keys(legacy_adjs));
            expect_self_is_last(candidate_adjs, ucell, it, ia);
            expect_self_is_last(reference_adjs, ucell, it, ia);
            EXPECT_THROW(candidate.Find_atom_from_legacy(ucell, it, ia, nullptr), std::runtime_error);
        }
    }

    remove("grid_driver_synthetic_half14_atom_pack_only.out");
    remove("grid_driver_synthetic_full27_reference.out");
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

TEST(SltkAtomArrangeSyntheticTest, ExplicitHalf14MatchesFull27ReferenceForSyntheticCells)
{
    SetGlobalV();
    for (SyntheticNeighborCase test_case: make_synthetic_neighbor_cases())
    {
        UnitCell* ucell = test_case.prepare.SetUcellInfo();
        unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

        for (const double radius: test_case.radii)
        {
            expect_atom_pack_only_half14_matches_full27_reference(*ucell, radius, true, test_case.name);
            expect_atom_pack_only_half14_matches_full27_reference(*ucell, radius, false, test_case.name);
        }

        delete ucell;
    }
}

TEST(NeighborRebuildStateTest, AppliesSkinThresholdAndInvalidatesMetadata)
{
    SyntheticNeighborCase test_case = make_synthetic_neighbor_cases()[0];
    UnitCell* ucell = test_case.prepare.SetUcellInfo();
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

    ModuleNeighbor::NeighborRebuildState state(3.0);
    const double physical_radius_bohr = 9.5;
    EXPECT_TRUE(state.needs_rebuild(*ucell, physical_radius_bohr, true));
    EXPECT_EQ(state.last_reason(), ModuleNeighbor::NeighborRebuildReason::FirstBuild);
    state.mark_rebuilt(*ucell, physical_radius_bohr, true);
    EXPECT_EQ(state.stats().rebuild_count, 1);

    ucell->atoms[0].tau[0].x += 0.10;
    EXPECT_FALSE(state.needs_rebuild(*ucell, physical_radius_bohr, true));
    EXPECT_EQ(state.last_reason(), ModuleNeighbor::NeighborRebuildReason::Reuse);
    state.mark_reused();
    EXPECT_EQ(state.stats().reuse_count, 1);
    EXPECT_NEAR(state.stats().last_max_displacement_bohr, 1.0, 1.0e-12);

    ucell->atoms[0].tau[0].x += 0.05;
    EXPECT_FALSE(state.needs_rebuild(*ucell, physical_radius_bohr, true));
    ucell->atoms[0].tau[0].x += 0.001;
    EXPECT_TRUE(state.needs_rebuild(*ucell, physical_radius_bohr, true));
    EXPECT_EQ(state.last_reason(), ModuleNeighbor::NeighborRebuildReason::DisplacementExceeded);

    state.mark_rebuilt(*ucell, physical_radius_bohr, true);
    EXPECT_TRUE(state.needs_rebuild(*ucell, physical_radius_bohr + 0.1, true));
    EXPECT_EQ(state.last_reason(), ModuleNeighbor::NeighborRebuildReason::SearchRadiusChanged);
    EXPECT_TRUE(state.needs_rebuild(*ucell, physical_radius_bohr, false));
    EXPECT_EQ(state.last_reason(), ModuleNeighbor::NeighborRebuildReason::BoundaryChanged);
    ucell->latvec.e11 += 0.01;
    EXPECT_TRUE(state.needs_rebuild(*ucell, physical_radius_bohr, true));
    EXPECT_EQ(state.last_reason(), ModuleNeighbor::NeighborRebuildReason::LatticeChanged);

    delete ucell;
}

TEST(NeighborRebuildStateTest, SkinChangeAndZeroSkinForceRebuild)
{
    SyntheticNeighborCase test_case = make_synthetic_neighbor_cases()[0];
    UnitCell* ucell = test_case.prepare.SetUcellInfo();
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

    ModuleNeighbor::NeighborRebuildState state(3.0);
    state.mark_rebuilt(*ucell, 9.5, true);
    state.set_skin_bohr(2.0);
    EXPECT_TRUE(state.needs_rebuild(*ucell, 9.5, true));
    EXPECT_EQ(state.last_reason(), ModuleNeighbor::NeighborRebuildReason::SkinChanged);

    state.mark_rebuilt(*ucell, 9.5, true);
    state.set_skin_bohr(0.0);
    EXPECT_TRUE(state.needs_rebuild(*ucell, 9.5, true));
    EXPECT_EQ(state.last_reason(), ModuleNeighbor::NeighborRebuildReason::SkinDisabled);
    state.mark_rebuilt(*ucell, 9.5, true);
    EXPECT_TRUE(state.needs_rebuild(*ucell, 9.5, true));
    EXPECT_EQ(state.last_reason(), ModuleNeighbor::NeighborRebuildReason::SkinDisabled);

    EXPECT_THROW(state.set_skin_bohr(-0.1), std::invalid_argument);
    delete ucell;
}

TEST(NeighborRebuildStateTest, UsesMinimumImageForPeriodicDisplacement)
{
    SyntheticNeighborCase test_case = make_synthetic_neighbor_cases()[0];
    UnitCell* ucell = test_case.prepare.SetUcellInfo();
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);
    ucell->atoms[0].tau[0] = ModuleBase::Vector3<double>(0.99, 0.0, 0.0);

    ModuleNeighbor::NeighborRebuildState state(3.0);
    state.mark_rebuilt(*ucell, 9.5, true);
    ucell->atoms[0].tau[0].x = 0.01;
    EXPECT_FALSE(state.needs_rebuild(*ucell, 9.5, true));
    EXPECT_NEAR(state.stats().last_max_displacement_bohr, 0.2, 1.0e-12);

    delete ucell;
}

TEST(NeighborRebuildStateTest, ReusedGridUsesCurrentCoordinatesAndPhysicalCutoff)
{
    SetGlobalV();
    std::vector<SyntheticNeighborCase> cases = make_synthetic_neighbor_cases();
    for (int case_index = 0; case_index < 2; ++case_index)
    {
        for (const bool pbc: {false, true})
        {
            UnitCell* ucell = cases[case_index].prepare.SetUcellInfo();
            unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

            const double physical_radius = cases[case_index].radii[0];
            const double skin = 0.30;
            std::ofstream cached_output("grid_driver_verlet_cached.out");
            Grid_Driver cached(PARAM.input.test_deconstructor, PARAM.input.test_grid);
            cached.init(cached_output, *ucell, physical_radius + skin, pbc);
            cached.set_query_radius(physical_radius);
            cached_output.close();
            const std::vector<ModuleNeighbor::NeighborPair> cached_pairs = cached.neighbor_pairs;

            ucell->atoms[0].tau[1].x += 0.05;
            std::ofstream rebuilt_output("grid_driver_verlet_rebuilt.out");
            Grid_Driver rebuilt(PARAM.input.test_deconstructor, PARAM.input.test_grid);
            rebuilt.init(rebuilt_output, *ucell, physical_radius, pbc);
            rebuilt_output.close();

            for (int atom = 0; atom < ucell->atoms[0].na; ++atom)
            {
                AdjacentAtomInfo cached_adjs;
                AdjacentAtomInfo rebuilt_adjs;
                cached.Find_atom(*ucell, 0, atom, &cached_adjs);
                rebuilt.Find_atom(*ucell, 0, atom, &rebuilt_adjs);
                EXPECT_EQ(collect_adjacent_keys(cached_adjs), collect_adjacent_keys(rebuilt_adjs));
                expect_self_is_last(cached_adjs, *ucell, 0, atom);
            }
            EXPECT_EQ(cached.neighbor_pairs, cached_pairs);

            remove("grid_driver_verlet_cached.out");
            remove("grid_driver_verlet_rebuilt.out");
            delete ucell;
        }
    }
}

TEST(NeighborRebuildStateTest, ReusedGridTracksPeriodicCoordinateWrapping)
{
    SetGlobalV();
    SyntheticNeighborCase test_case = make_synthetic_neighbor_cases()[0];
    UnitCell* ucell = test_case.prepare.SetUcellInfo();
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);
    ucell->atoms[0].tau[0].x = 0.01;

    std::ofstream cached_output("grid_driver_verlet_wrapped_cached.out");
    Grid_Driver cached(PARAM.input.test_deconstructor, PARAM.input.test_grid);
    cached.init(cached_output, *ucell, 1.25, true);
    cached.set_query_radius(0.95);
    cached_output.close();

    ucell->atoms[0].tau[0].x = 0.99;
    std::ofstream rebuilt_output("grid_driver_verlet_wrapped_rebuilt.out");
    Grid_Driver rebuilt(PARAM.input.test_deconstructor, PARAM.input.test_grid);
    rebuilt.init(rebuilt_output, *ucell, 0.95, true);
    rebuilt_output.close();

    for (int atom = 0; atom < ucell->atoms[0].na; ++atom)
    {
        AdjacentAtomInfo cached_adjs;
        AdjacentAtomInfo rebuilt_adjs;
        cached.Find_atom(*ucell, 0, atom, &cached_adjs);
        rebuilt.Find_atom(*ucell, 0, atom, &rebuilt_adjs);
        EXPECT_EQ(collect_adjacent_keys(cached_adjs), collect_adjacent_keys(rebuilt_adjs));
    }

    remove("grid_driver_verlet_wrapped_cached.out");
    remove("grid_driver_verlet_wrapped_rebuilt.out");
    delete ucell;
}
