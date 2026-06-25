#include "../sltk_grid_driver.h"
#include "../neighbor_rebuild_state.h"

#include "gtest/gtest.h"

#ifdef __MPI
#include <mpi.h>
#endif

#include "source_cell/read_stru.h"
#include "source_io/module_parameter/parameter.h"
#include "synthetic_neighbor_unitcell.h"

#include <algorithm>
#include <array>
#include <fstream>
#include <string>
#include <tuple>
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
using PairKey = std::tuple<int, int, int, int, int, int, int>;

ModuleNeighbor::NeighborMpiComm world_comm()
{
#ifdef __MPI
    return MPI_COMM_WORLD;
#else
    return ModuleNeighbor::kSerialMpiCommWorld;
#endif
}

std::vector<PairKey> collect_pair_keys(const std::vector<ModuleNeighbor::NeighborPair>& pairs)
{
    std::vector<PairKey> keys;
    keys.reserve(pairs.size());
    for (const ModuleNeighbor::NeighborPair& pair: pairs)
    {
        keys.push_back(PairKey(pair.center_type,
                               pair.center_natom,
                               pair.neighbor_type,
                               pair.neighbor_natom,
                               pair.cell_x,
                               pair.cell_y,
                               pair.cell_z));
    }
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    return keys;
}

std::vector<PairKey> gather_pair_keys(const std::vector<PairKey>& local_keys,
                                      ModuleNeighbor::NeighborMpiComm communicator)
{
    std::vector<int> send_buffer;
    send_buffer.reserve(local_keys.size() * 7);
    for (const PairKey& key: local_keys)
    {
        send_buffer.push_back(std::get<0>(key));
        send_buffer.push_back(std::get<1>(key));
        send_buffer.push_back(std::get<2>(key));
        send_buffer.push_back(std::get<3>(key));
        send_buffer.push_back(std::get<4>(key));
        send_buffer.push_back(std::get<5>(key));
        send_buffer.push_back(std::get<6>(key));
    }

#ifdef __MPI
    int size = 1;
    MPI_Comm_size(communicator, &size);
    const int send_count = static_cast<int>(send_buffer.size());
    std::vector<int> counts(size, 0);
    MPI_Allgather(&send_count, 1, MPI_INT, counts.data(), 1, MPI_INT, communicator);
    std::vector<int> displacements(size, 0);
    int total_count = 0;
    for (int rank = 0; rank < size; ++rank)
    {
        displacements[rank] = total_count;
        total_count += counts[rank];
    }
    std::vector<int> recv_buffer(total_count, 0);
    MPI_Allgatherv(send_buffer.empty() ? nullptr : send_buffer.data(),
                   send_count,
                   MPI_INT,
                   recv_buffer.empty() ? nullptr : recv_buffer.data(),
                   counts.data(),
                   displacements.data(),
                   MPI_INT,
                   communicator);
#else
    (void)communicator;
    const std::vector<int>& recv_buffer = send_buffer;
#endif

    std::vector<PairKey> keys;
    for (int offset = 0; offset + 6 < static_cast<int>(recv_buffer.size()); offset += 7)
    {
        keys.push_back(PairKey(recv_buffer[offset + 0],
                               recv_buffer[offset + 1],
                               recv_buffer[offset + 2],
                               recv_buffer[offset + 3],
                               recv_buffer[offset + 4],
                               recv_buffer[offset + 5],
                               recv_buffer[offset + 6]));
    }
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    return keys;
}

std::vector<PairKey> collect_query_keys(const Grid_Driver& grid, const UnitCell& ucell)
{
    std::vector<PairKey> keys;
    for (int type = 0; type < ucell.ntype; ++type)
    {
        for (int natom = 0; natom < ucell.atoms[type].na; ++natom)
        {
            if (!grid.is_local_center(type, natom))
            {
                continue;
            }
            AdjacentAtomInfo adjs;
            grid.Find_atom(ucell, type, natom, &adjs);
            for (int index = 0; index < adjs.adj_num; ++index)
            {
                keys.push_back(PairKey(type,
                                       natom,
                                       adjs.ntype[index],
                                       adjs.natom[index],
                                       adjs.box[index].x,
                                       adjs.box[index].y,
                                       adjs.box[index].z));
            }
        }
    }
    return keys;
}

int mpi_rank()
{
#ifdef __MPI
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
#else
    return 0;
#endif
}

int mpi_size()
{
#ifdef __MPI
    int size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
#else
    return 1;
#endif
}

void expect_mpi_grid_matches_reference(SyntheticNeighborCase test_case,
                                       const double radius,
                                       const bool pbc)
{
    UnitCell* ucell = test_case.prepare.SetUcellInfo();
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

    const std::string suffix = std::to_string(mpi_rank());
    std::ofstream mpi_output("mpi_grid_" + test_case.name + "_" + suffix + ".out");
    Grid_Driver mpi_grid(0, 0);
    ModuleNeighbor::MpiGhostExchangeStats stats;
    mpi_grid.init_mpi(mpi_output, *ucell, radius, pbc, world_comm(), &stats);
    mpi_output.close();

    std::ofstream reference_output("mpi_grid_reference_" + test_case.name + "_" + suffix + ".out");
    Grid_Driver reference_grid(0, 0);
    reference_grid.init(reference_output, *ucell, radius, pbc);
    reference_output.close();

    const std::vector<PairKey> distributed_keys
        = gather_pair_keys(collect_pair_keys(mpi_grid.neighbor_pairs), mpi_grid.mpi_domain().cart_comm());
    EXPECT_EQ(distributed_keys, collect_pair_keys(reference_grid.neighbor_pairs));
    EXPECT_TRUE(mpi_grid.mpi_mode());
    EXPECT_TRUE(mpi_grid.neighbor_pairs_27.empty());
    const int candidate_count = static_cast<int>(
        std::count(mpi_grid.atom_pack.is_ghost.begin(), mpi_grid.atom_pack.is_ghost.end(), true));
    EXPECT_LE(stats.received_ghost_count, candidate_count);

    int global_index = 0;
    for (int type = 0; type < ucell->ntype; ++type)
    {
        for (int natom = 0; natom < ucell->atoms[type].na; ++natom)
        {
            const int local_owner = mpi_grid.is_local_center(type, natom) ? 1 : 0;
            int owner_count = local_owner;
#ifdef __MPI
            MPI_Allreduce(&local_owner,
                          &owner_count,
                          1,
                          MPI_INT,
                          MPI_SUM,
                          mpi_grid.mpi_domain().cart_comm());
#endif
            EXPECT_EQ(owner_count, 1) << "global atom " << global_index;

            AdjacentAtomInfo adjs;
            if (local_owner)
            {
                EXPECT_NO_THROW(mpi_grid.Find_atom(*ucell, type, natom, &adjs));
                EXPECT_FALSE(adjs.ntype.empty());
                if (!adjs.ntype.empty())
                {
                    const int last = static_cast<int>(adjs.ntype.size()) - 1;
                    EXPECT_EQ(adjs.ntype[last], type);
                    EXPECT_EQ(adjs.natom[last], natom);
                    EXPECT_EQ(adjs.box[last].x, 0);
                    EXPECT_EQ(adjs.box[last].y, 0);
                    EXPECT_EQ(adjs.box[last].z, 0);
                }
            }
            else
            {
                EXPECT_THROW(mpi_grid.Find_atom(*ucell, type, natom, &adjs), std::runtime_error);
            }
            ++global_index;
        }
    }

    remove(("mpi_grid_" + test_case.name + "_" + suffix + ".out").c_str());
    remove(("mpi_grid_reference_" + test_case.name + "_" + suffix + ".out").c_str());
    delete ucell;
}
} // namespace

TEST(MpiGridTest, OrthogonalNonPeriodicMatchesReference)
{
    const SyntheticNeighborCase test_case = make_synthetic_neighbor_cases()[0];
    expect_mpi_grid_matches_reference(test_case, test_case.radii[0], false);
}

TEST(MpiGridTest, VerletRebuildDecisionIsCollective)
{
    SyntheticNeighborCase test_case = make_synthetic_neighbor_cases()[0];
    UnitCell* ucell = test_case.prepare.SetUcellInfo();
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

    ModuleNeighbor::NeighborRebuildState state(3.0);
    state.mark_rebuilt(*ucell, 9.5, true);
    if (mpi_rank() == 0)
    {
        ucell->atoms[0].tau[0].x += 0.2;
    }
#ifdef __MPI
    EXPECT_TRUE(state.needs_rebuild_mpi(*ucell, 9.5, true, MPI_COMM_WORLD));
#else
    EXPECT_TRUE(state.needs_rebuild(*ucell, 9.5, true));
#endif

    delete ucell;
}

TEST(MpiGridTest, VerletRebuildReasonUsesExplicitCollectivePriority)
{
    if (mpi_size() < 2)
    {
        GTEST_SKIP() << "Collective reason priority requires at least two MPI ranks.";
    }

    SyntheticNeighborCase test_case = make_synthetic_neighbor_cases()[0];
    UnitCell* ucell = test_case.prepare.SetUcellInfo();
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

    ModuleNeighbor::NeighborRebuildState state(3.0);
    state.mark_rebuilt(*ucell, 9.5, true);
    const bool pbc = mpi_rank() == 0 ? false : true;
    const double radius = mpi_rank() == 1 ? 9.6 : 9.5;
#ifdef __MPI
    EXPECT_TRUE(state.needs_rebuild_mpi(*ucell, radius, pbc, MPI_COMM_WORLD));
    EXPECT_EQ(state.last_reason(), ModuleNeighbor::NeighborRebuildReason::BoundaryChanged);
#endif

    delete ucell;
}

TEST(MpiGridTest, ReusedVerletGridRemainsValidAcrossDomainBoundary)
{
    SyntheticNeighborCase test_case = make_synthetic_neighbor_cases()[0];
    UnitCell* ucell = test_case.prepare.SetUcellInfo();
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);
    ucell->atoms[0].tau[0].x = 0.49;

    std::ofstream cached_output("mpi_grid_verlet_cached_" + std::to_string(mpi_rank()) + ".out");
    Grid_Driver cached(0, 0);
    cached.init_mpi(cached_output, *ucell, 1.25, true, world_comm());
    cached.set_query_radius(0.95);
    cached_output.close();

    ucell->atoms[0].tau[0].x = 0.51;
    std::ofstream rebuilt_output("mpi_grid_verlet_rebuilt_" + std::to_string(mpi_rank()) + ".out");
    Grid_Driver rebuilt(0, 0);
    rebuilt.init_mpi(rebuilt_output, *ucell, 0.95, true, world_comm());
    rebuilt_output.close();

    const std::vector<PairKey> cached_keys
        = gather_pair_keys(collect_query_keys(cached, *ucell), cached.mpi_domain().cart_comm());
    const std::vector<PairKey> rebuilt_keys
        = gather_pair_keys(collect_query_keys(rebuilt, *ucell), rebuilt.mpi_domain().cart_comm());
    EXPECT_EQ(cached_keys, rebuilt_keys);

    remove(("mpi_grid_verlet_cached_" + std::to_string(mpi_rank()) + ".out").c_str());
    remove(("mpi_grid_verlet_rebuilt_" + std::to_string(mpi_rank()) + ".out").c_str());
    delete ucell;
}

TEST(MpiGridTest, OrthogonalPeriodicMatchesReference)
{
    const SyntheticNeighborCase test_case = make_synthetic_neighbor_cases()[0];
    expect_mpi_grid_matches_reference(test_case, test_case.radii[0], true);
}

TEST(MpiGridTest, TriclinicPeriodicMatchesReference)
{
    const SyntheticNeighborCase test_case = make_synthetic_neighbor_cases()[1];
    expect_mpi_grid_matches_reference(test_case, test_case.radii[0], true);
}

TEST(MpiGridTest, SerialReinitializationClearsMpiState)
{
    SyntheticNeighborCase test_case = make_synthetic_neighbor_cases()[0];
    UnitCell* ucell = test_case.prepare.SetUcellInfo();
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);

    const std::string output_name
        = "mpi_grid_reinitialize_" + std::to_string(mpi_rank()) + ".out";
    std::ofstream output(output_name);
    Grid_Driver grid(0, 0);
    grid.init_mpi(output, *ucell, test_case.radii[0], true, world_comm());
    ASSERT_TRUE(grid.mpi_mode());

    grid.init(output, *ucell, test_case.radii[0], true);
    EXPECT_FALSE(grid.mpi_mode());
    EXPECT_FALSE(grid.mpi_domain().initialized());
    EXPECT_TRUE(grid.atoms_in_box.empty());
    EXPECT_TRUE(grid.all_adj_info.empty());

    for (int type = 0; type < ucell->ntype; ++type)
    {
        for (int natom = 0; natom < ucell->atoms[type].na; ++natom)
        {
            AdjacentAtomInfo adjs;
            EXPECT_NO_THROW(grid.Find_atom(*ucell, type, natom, &adjs));
        }
    }

    output.close();
    remove(output_name.c_str());
    delete ucell;
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
#ifdef __MPI
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized)
    {
        MPI_Init(&argc, &argv);
    }
#endif
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (!finalized)
    {
        MPI_Finalize();
    }
#endif
    return result;
}
