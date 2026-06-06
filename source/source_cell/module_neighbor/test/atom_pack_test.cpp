#include "../atom_pack.h"
#include "../sltk_grid.h"

#include "gtest/gtest.h"
#include "prepare_unitcell.h"

#include "source_cell/read_stru.h"

#include <algorithm>
#include <fstream>
#include <vector>

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
int count_origin_cell_atoms(const ModuleNeighbor::AtomPack& pack)
{
    int count = 0;
    for (int i = 0; i < pack.size(); ++i)
    {
        if (pack.cell_x[i] == 0 && pack.cell_y[i] == 0 && pack.cell_z[i] == 0)
        {
            ++count;
        }
    }
    return count;
}

std::vector<int> collect_grid_box_counts(const Grid& grid)
{
    std::vector<int> counts;
    counts.reserve(grid.box_nx * grid.box_ny * grid.box_nz);
    for (int bx = 0; bx < grid.box_nx; ++bx)
    {
        for (int by = 0; by < grid.box_ny; ++by)
        {
            for (int bz = 0; bz < grid.box_nz; ++bz)
            {
                counts.push_back(static_cast<int>(grid.atoms_in_box[bx][by][bz].size()));
            }
        }
    }
    return counts;
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
                pairs.push_back(ModuleNeighbor::NeighborPair{center_type,
                                                             center_natom,
                                                             atom->type,
                                                             atom->natom,
                                                             atom->cell_x,
                                                             atom->cell_y,
                                                             atom->cell_z});
            }
        }
    }
    std::sort(pairs.begin(), pairs.end());
    return pairs;
}

void expect_neighbor_pair_indices_match_pack(const std::vector<ModuleNeighbor::NeighborPair>& pairs,
                                             const ModuleNeighbor::AtomPack& pack)
{
    for (const ModuleNeighbor::NeighborPair& pair: pairs)
    {
        ASSERT_GE(pair.center_index, 0);
        ASSERT_LT(pair.center_index, pack.size());
        ASSERT_GE(pair.neighbor_index, 0);
        ASSERT_LT(pair.neighbor_index, pack.size());
        EXPECT_EQ(pack.type[pair.center_index], pair.center_type);
        EXPECT_EQ(pack.natom[pair.center_index], pair.center_natom);
        EXPECT_EQ(pack.cell_x[pair.center_index], 0);
        EXPECT_EQ(pack.cell_y[pair.center_index], 0);
        EXPECT_EQ(pack.cell_z[pair.center_index], 0);
        EXPECT_EQ(pack.type[pair.neighbor_index], pair.neighbor_type);
        EXPECT_EQ(pack.natom[pair.neighbor_index], pair.neighbor_natom);
        EXPECT_EQ(pack.cell_x[pair.neighbor_index], pair.cell_x);
        EXPECT_EQ(pack.cell_y[pair.neighbor_index], pair.cell_y);
        EXPECT_EQ(pack.cell_z[pair.neighbor_index], pair.cell_z);
    }
}
} // namespace

TEST(AtomPackTest, BuildsNonPeriodicPackWithoutImages)
{
    UcellTestPrepare utp = UcellTestLib["Si"];
    UnitCell* ucell = utp.SetUcellInfo();
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

    const ModuleNeighbor::AtomPack pack = ModuleNeighbor::build_atom_pack_from_unitcell(*ucell, 0.5, false);

    EXPECT_EQ(pack.size(), ucell->nat);
    EXPECT_EQ(count_origin_cell_atoms(pack), ucell->nat);
    EXPECT_FALSE(pack.empty());
    for (int i = 0; i < pack.size(); ++i)
    {
        EXPECT_FALSE(pack.is_ghost[i]);
        EXPECT_EQ(pack.cell_x[i], 0);
        EXPECT_EQ(pack.cell_y[i], 0);
        EXPECT_EQ(pack.cell_z[i], 0);
    }

    delete ucell;
}

TEST(AtomPackTest, BuildsPeriodicPackWithImages)
{
    UcellTestPrepare utp = UcellTestLib["Si"];
    UnitCell* ucell = utp.SetUcellInfo();
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

    const ModuleNeighbor::AtomPack pack = ModuleNeighbor::build_atom_pack_from_unitcell(*ucell, 0.5, true);

    EXPECT_GT(pack.size(), ucell->nat);
    EXPECT_EQ(count_origin_cell_atoms(pack), ucell->nat);
    EXPECT_EQ(pack.type.size(), pack.x.size());
    EXPECT_EQ(pack.global_index.size(), pack.x.size());

    bool saw_shifted_image = false;
    for (int i = 0; i < pack.size(); ++i)
    {
        saw_shifted_image = saw_shifted_image || pack.cell_x[i] != 0 || pack.cell_y[i] != 0 || pack.cell_z[i] != 0;
        EXPECT_FALSE(pack.is_ghost[i]);
        EXPECT_GE(pack.global_index[i], 0);
        EXPECT_LT(pack.global_index[i], ucell->nat);
    }
    EXPECT_TRUE(saw_shifted_image);

    delete ucell;
}

TEST(AtomPackTest, AppendsMpiGhostRecord)
{
    ModuleNeighbor::AtomPack pack;
    pack.append_atom(0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, false);

    const ModuleNeighbor::MpiAtomRecord ghost(1.0, 2.0, 3.0, 7, std::array<int, 3>{{1, 0, -1}}, true);
    pack.append_mpi_record(ghost);

    ASSERT_EQ(pack.size(), 2);
    EXPECT_TRUE(pack.is_ghost[1]);
    EXPECT_EQ(pack.global_index[1], 7);
    EXPECT_EQ(pack.cell_x[1], 1);
    EXPECT_EQ(pack.cell_y[1], 0);
    EXPECT_EQ(pack.cell_z[1], -1);
    EXPECT_DOUBLE_EQ(pack.x[1], 1.0);
    EXPECT_DOUBLE_EQ(pack.y[1], 2.0);
    EXPECT_DOUBLE_EQ(pack.z[1], 3.0);
}

TEST(AtomPackTest, BuildsFlatGridStorageWithoutLossOrDuplication)
{
    UcellTestPrepare utp = UcellTestLib["Si"];
    UnitCell* ucell = utp.SetUcellInfo();
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

    const double radius = 0.5;
    const ModuleNeighbor::AtomPack pack = ModuleNeighbor::build_atom_pack_from_unitcell(*ucell, radius, true);
    const ModuleNeighbor::GridStorage storage = ModuleNeighbor::build_grid_storage_from_atom_pack(pack, radius + 0.1);

    ASSERT_EQ(storage.atoms_in_box.size(), static_cast<std::size_t>(pack.size()));
    ASSERT_EQ(storage.box_offset.size(), static_cast<std::size_t>(storage.box_size()));
    ASSERT_EQ(storage.box_count.size(), static_cast<std::size_t>(storage.box_size()));
    EXPECT_GT(storage.box_nx, 0);
    EXPECT_GT(storage.box_ny, 0);
    EXPECT_GT(storage.box_nz, 0);

    std::vector<bool> visited(pack.size(), false);
    for (int box_id = 0; box_id < storage.box_size(); ++box_id)
    {
        const int begin = storage.box_offset[box_id];
        const int end = begin + storage.box_count[box_id];
        ASSERT_GE(begin, 0);
        ASSERT_GE(end, begin);
        ASSERT_LE(end, static_cast<int>(storage.atoms_in_box.size()));

        for (int offset = begin; offset < end; ++offset)
        {
            const int atom_index = storage.atoms_in_box[offset];
            ASSERT_GE(atom_index, 0);
            ASSERT_LT(atom_index, pack.size());
            EXPECT_FALSE(visited[atom_index]);
            visited[atom_index] = true;
            EXPECT_EQ(storage.get_box_id_from_coord(pack.x[atom_index], pack.y[atom_index], pack.z[atom_index]), box_id);
        }
    }

    EXPECT_TRUE(std::all_of(visited.begin(), visited.end(), [](const bool value) { return value; }));

    delete ucell;
}

TEST(AtomPackTest, GridStorageClampsOutOfRangeCoordinates)
{
    ModuleNeighbor::AtomPack pack;
    pack.append_atom(0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, false);
    pack.append_atom(1.0, 1.0, 1.0, 0, 1, 0, 0, 0, 1, false);

    const ModuleNeighbor::GridStorage storage = ModuleNeighbor::build_grid_storage_from_atom_pack(pack, 0.6);

    EXPECT_EQ(storage.get_box_id(-100, -100, -100), 0);
    EXPECT_EQ(storage.get_box_id_from_coord(storage.x_min - 100.0, storage.y_min - 100.0, storage.z_min - 100.0), 0);
    EXPECT_EQ(storage.get_box_id(100, 100, 100), storage.box_size() - 1);
    EXPECT_EQ(storage.get_box_id_from_coord(storage.x_max + 100.0, storage.y_max + 100.0, storage.z_max + 100.0),
              storage.box_size() - 1);
}

TEST(AtomPackTest, FlatGridStorageMatchesLegacyGridBoxCountsWithoutPbc)
{
    UcellTestPrepare utp = UcellTestLib["Si"];
    UnitCell* ucell = utp.SetUcellInfo();
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

    const double radius = 0.5;
    const ModuleNeighbor::AtomPack pack = ModuleNeighbor::build_atom_pack_from_unitcell(*ucell, radius, false);
    const ModuleNeighbor::GridStorage storage = ModuleNeighbor::build_grid_storage_from_atom_pack(pack, radius + 0.1);

    std::ofstream ofs("atom_pack_grid_compare.out");
    Grid grid(0);
    grid.init(ofs, *ucell, radius, false);
    ofs.close();

    EXPECT_EQ(storage.box_nx, grid.box_nx);
    EXPECT_EQ(storage.box_ny, grid.box_ny);
    EXPECT_EQ(storage.box_nz, grid.box_nz);
    EXPECT_EQ(storage.box_count, collect_grid_box_counts(grid));

    remove("atom_pack_grid_compare.out");
    delete ucell;
}

TEST(AtomPackTest, FlatGridSearchMatchesLegacyGridPairsWithoutPbc)
{
    UcellTestPrepare utp = UcellTestLib["Si"];
    UnitCell* ucell = utp.SetUcellInfo();
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

    const double radius = 0.5;
    const ModuleNeighbor::AtomPack pack = ModuleNeighbor::build_atom_pack_from_unitcell(*ucell, radius, false);
    const ModuleNeighbor::GridStorage storage = ModuleNeighbor::build_grid_storage_from_atom_pack(pack, radius + 0.1);
    const std::vector<ModuleNeighbor::NeighborPair> flat_pairs
        = ModuleNeighbor::build_neighbor_pairs_27(pack, storage, radius);

    std::ofstream ofs("atom_pack_neighbor_compare_non_pbc.out");
    Grid grid(0);
    grid.init(ofs, *ucell, radius, false);
    ofs.close();
    const std::vector<ModuleNeighbor::NeighborPair> legacy_pairs = collect_legacy_neighbor_pairs(grid);

    EXPECT_EQ(flat_pairs, legacy_pairs);
    EXPECT_EQ(std::adjacent_find(flat_pairs.begin(), flat_pairs.end()), flat_pairs.end());

    remove("atom_pack_neighbor_compare_non_pbc.out");
    delete ucell;
}

TEST(AtomPackTest, FlatGridSearchMatchesLegacyGridPairsWithPbc)
{
    UcellTestPrepare utp = UcellTestLib["Si"];
    UnitCell* ucell = utp.SetUcellInfo();
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

    const double radius = 0.5;
    const ModuleNeighbor::AtomPack pack = ModuleNeighbor::build_atom_pack_from_unitcell(*ucell, radius, true);
    const ModuleNeighbor::GridStorage storage = ModuleNeighbor::build_grid_storage_from_atom_pack(pack, radius + 0.1);
    const std::vector<ModuleNeighbor::NeighborPair> flat_pairs
        = ModuleNeighbor::build_neighbor_pairs_27(pack, storage, radius);

    std::ofstream ofs("atom_pack_neighbor_compare_pbc.out");
    Grid grid(0);
    grid.init(ofs, *ucell, radius, true);
    ofs.close();
    const std::vector<ModuleNeighbor::NeighborPair> legacy_pairs = collect_legacy_neighbor_pairs(grid);

    EXPECT_EQ(flat_pairs, legacy_pairs);
    EXPECT_EQ(std::adjacent_find(flat_pairs.begin(), flat_pairs.end()), flat_pairs.end());

    remove("atom_pack_neighbor_compare_pbc.out");
    delete ucell;
}

TEST(AtomPackTest, HalfDomainSearchRestoresFullPairsWithoutPbc)
{
    UcellTestPrepare utp = UcellTestLib["Si"];
    UnitCell* ucell = utp.SetUcellInfo();
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

    const double radius = 0.5;
    const ModuleNeighbor::AtomPack pack = ModuleNeighbor::build_atom_pack_from_unitcell(*ucell, radius, false);
    const ModuleNeighbor::GridStorage storage = ModuleNeighbor::build_grid_storage_from_atom_pack(pack, radius + 0.1);
    const std::vector<ModuleNeighbor::NeighborPair> full_pairs
        = ModuleNeighbor::build_neighbor_pairs_27(pack, storage, radius);
    const std::vector<ModuleNeighbor::NeighborPair> restored_pairs
        = ModuleNeighbor::build_neighbor_pairs_14(pack, storage, radius);

    EXPECT_EQ(restored_pairs, full_pairs);
    EXPECT_EQ(std::adjacent_find(restored_pairs.begin(), restored_pairs.end()), restored_pairs.end());
    expect_neighbor_pair_indices_match_pack(restored_pairs, pack);

    delete ucell;
}

TEST(AtomPackTest, HalfDomainSearchRestoresFullPairsWithPbc)
{
    UcellTestPrepare utp = UcellTestLib["Si"];
    UnitCell* ucell = utp.SetUcellInfo();
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

    const double radius = 0.5;
    const ModuleNeighbor::AtomPack pack = ModuleNeighbor::build_atom_pack_from_unitcell(*ucell, radius, true);
    const ModuleNeighbor::GridStorage storage = ModuleNeighbor::build_grid_storage_from_atom_pack(pack, radius + 0.1);
    const std::vector<ModuleNeighbor::NeighborPair> full_pairs
        = ModuleNeighbor::build_neighbor_pairs_27(pack, storage, radius);
    const std::vector<ModuleNeighbor::NeighborPair> restored_pairs
        = ModuleNeighbor::build_neighbor_pairs_14(pack, storage, radius);

    EXPECT_EQ(restored_pairs, full_pairs);
    EXPECT_EQ(std::adjacent_find(restored_pairs.begin(), restored_pairs.end()), restored_pairs.end());
    expect_neighbor_pair_indices_match_pack(restored_pairs, pack);

    delete ucell;
}

TEST(AtomPackTest, FlatGridSearchHandlesSingleAtomWithoutSelfNeighbor)
{
    ModuleNeighbor::AtomPack pack;
    pack.append_atom(0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, false);

    const ModuleNeighbor::GridStorage storage = ModuleNeighbor::build_grid_storage_from_atom_pack(pack, 0.6);
    const std::vector<ModuleNeighbor::NeighborPair> pairs
        = ModuleNeighbor::build_neighbor_pairs_27(pack, storage, 0.5);

    EXPECT_TRUE(pairs.empty());
}

TEST(AtomPackTest, HalfDomainSearchHandlesSingleAtomWithoutSelfNeighbor)
{
    ModuleNeighbor::AtomPack pack;
    pack.append_atom(0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, false);

    const ModuleNeighbor::GridStorage storage = ModuleNeighbor::build_grid_storage_from_atom_pack(pack, 0.6);
    const std::vector<ModuleNeighbor::NeighborPair> pairs
        = ModuleNeighbor::build_neighbor_pairs_14(pack, storage, 0.5);

    EXPECT_TRUE(pairs.empty());
}
