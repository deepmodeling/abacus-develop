#include "../atom_pack.h"
#include "../mpi_domain.h"

#include "gtest/gtest.h"

#ifdef __MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <array>
#include <tuple>
#include <vector>

namespace
{
ModuleNeighbor::NeighborMpiComm world_comm()
{
#ifdef __MPI
    return MPI_COMM_WORLD;
#else
    return ModuleNeighbor::kSerialMpiCommWorld;
#endif
}

void ensure_mpi_initialized()
{
#ifdef __MPI
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized)
    {
        int argc = 0;
        char** argv = nullptr;
        MPI_Init(&argc, &argv);
    }
#endif
}

int count_ghost_atoms(const ModuleNeighbor::AtomPack& pack)
{
    int count = 0;
    for (int i = 0; i < pack.size(); ++i)
    {
        if (pack.is_ghost[i])
        {
            ++count;
        }
    }
    return count;
}

using PairKey = std::tuple<int, int, int, int, int, int, int>;

PairKey make_pair_key(const ModuleNeighbor::NeighborPair& pair)
{
    return PairKey(pair.center_type,
                   pair.center_natom,
                   pair.neighbor_type,
                   pair.neighbor_natom,
                   pair.cell_x,
                   pair.cell_y,
                   pair.cell_z);
}

std::vector<PairKey> collect_pair_keys(const std::vector<ModuleNeighbor::NeighborPair>& pairs)
{
    std::vector<PairKey> keys;
    keys.reserve(pairs.size());
    for (const ModuleNeighbor::NeighborPair& pair: pairs)
    {
        keys.push_back(make_pair_key(pair));
    }
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    return keys;
}

ModuleNeighbor::AtomPack build_pack_from_records(const std::vector<ModuleNeighbor::MpiAtomRecord>& local_records,
                                                 const std::vector<ModuleNeighbor::MpiAtomRecord>& ghost_records)
{
    ModuleNeighbor::AtomPack pack;
    pack.reserve(static_cast<int>(local_records.size() + ghost_records.size()));
    for (const ModuleNeighbor::MpiAtomRecord& record: local_records)
    {
        pack.append_mpi_record(record);
    }
    for (const ModuleNeighbor::MpiAtomRecord& record: ghost_records)
    {
        pack.append_mpi_record(record);
    }
    return pack;
}

std::vector<PairKey> gather_pair_keys(const std::vector<PairKey>& local_keys, ModuleNeighbor::NeighborMpiComm comm)
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
    int comm_size = 1;
    MPI_Comm_size(comm, &comm_size);
    const int send_count = static_cast<int>(send_buffer.size());
    std::vector<int> counts(comm_size, 0);
    MPI_Allgather(&send_count, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);

    std::vector<int> displacements(comm_size, 0);
    int total_count = 0;
    for (int i = 0; i < comm_size; ++i)
    {
        displacements[i] = total_count;
        total_count += counts[i];
    }

    std::vector<int> recv_buffer(total_count, 0);
    MPI_Allgatherv(send_buffer.empty() ? 0 : send_buffer.data(),
                   send_count,
                   MPI_INT,
                   recv_buffer.empty() ? 0 : recv_buffer.data(),
                   counts.data(),
                   displacements.data(),
                   MPI_INT,
                   comm);
#else
    (void)comm;
    std::vector<int> recv_buffer = send_buffer;
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

std::vector<ModuleNeighbor::NeighborPair> build_pairs_from_pack(const ModuleNeighbor::AtomPack& pack,
                                                                const double cutoff)
{
    const ModuleNeighbor::GridStorage storage = ModuleNeighbor::build_grid_storage_from_atom_pack(pack, cutoff + 0.1);
    return ModuleNeighbor::build_neighbor_pairs_14(pack, storage, cutoff);
}

std::array<double, 3> local_cartesian_point(const ModuleNeighbor::MpiDomain& domain,
                                            const double upper_x_offset)
{
    const double fx = std::max(domain.local_bounds().lower[0],
                               domain.local_bounds().upper[0] - upper_x_offset);
    const double fy = 0.5 * (domain.local_bounds().lower[1] + domain.local_bounds().upper[1]);
    const double fz = 0.5 * (domain.local_bounds().lower[2] + domain.local_bounds().upper[2]);
    return domain.fractional_to_cartesian(fx, fy, fz);
}
} // namespace

TEST(MpiNeighborPrototypeTest, AppendsGhostsToAtomPack)
{
    ensure_mpi_initialized();

    const double cutoff = 0.75;
    ModuleNeighbor::MpiDomain domain;
    domain.initialize(world_comm(),
                      std::array<double, 3>{{0.0, 0.0, 0.0}},
                      std::array<double, 3>{{8.0, 4.0, 2.0}},
                      cutoff,
                      true);

    const std::array<double, 3> cartesian
        = local_cartesian_point(domain, 0.25 * domain.fractional_ghost_padding()[0]);

    std::vector<ModuleNeighbor::MpiAtomRecord> local_records;
    local_records.push_back(ModuleNeighbor::MpiAtomRecord(cartesian[0],
                                                          cartesian[1],
                                                          cartesian[2],
                                                          domain.rank()));

    ModuleNeighbor::AtomPack pack;
    pack.append_mpi_record(local_records.front(), 0, domain.rank());

    ModuleNeighbor::MpiGhostExchangeStats stats;
    const std::vector<ModuleNeighbor::MpiAtomRecord> ghosts = domain.exchange_ghost_atoms(local_records, &stats);
    for (const ModuleNeighbor::MpiAtomRecord& ghost: ghosts)
    {
        pack.append_mpi_record(ghost);
    }

    EXPECT_EQ(pack.size(), static_cast<int>(ghosts.size()) + 1);
    EXPECT_EQ(count_ghost_atoms(pack), static_cast<int>(ghosts.size()));
    EXPECT_EQ(stats.received_ghost_count, static_cast<int>(ghosts.size()));
    EXPECT_EQ(stats.sent_payload_count, stats.sent_atom_count * 9);
    EXPECT_EQ(stats.received_payload_count, stats.received_ghost_count * 9);

#ifdef __MPI
    int local_ghost_count = static_cast<int>(ghosts.size());
    int global_ghost_count = 0;
    MPI_Allreduce(&local_ghost_count, &global_ghost_count, 1, MPI_INT, MPI_SUM, domain.cart_comm());
    if (domain.size() == 1)
    {
        EXPECT_EQ(global_ghost_count, 0);
    }
    else
    {
        EXPECT_GT(global_ghost_count, 0);
    }
#else
    EXPECT_TRUE(ghosts.empty());
#endif
}

TEST(MpiNeighborPrototypeTest, GhostsParticipateInNonPeriodicAtomPackSearch)
{
    ensure_mpi_initialized();

    const double cutoff = 0.75;
    ModuleNeighbor::MpiDomain domain;
    domain.initialize(world_comm(),
                      std::array<double, 3>{{0.0, 0.0, 0.0}},
                      std::array<double, 3>{{8.0, 2.0, 2.0}},
                      cutoff,
                      false);

    if (domain.size() != 2)
    {
        GTEST_SKIP() << "This MPI neighbor prototype scenario requires exactly 2 ranks.";
    }

    std::vector<ModuleNeighbor::MpiAtomRecord> local_records;
    if (domain.rank() == 0)
    {
        local_records.push_back(ModuleNeighbor::MpiAtomRecord(3.8,
                                                              1.0,
                                                              1.0,
                                                              0,
                                                              std::array<int, 3>{{0, 0, 0}},
                                                              false,
                                                              0,
                                                              0));
    }
    else
    {
        local_records.push_back(ModuleNeighbor::MpiAtomRecord(4.2,
                                                              1.0,
                                                              1.0,
                                                              1,
                                                              std::array<int, 3>{{0, 0, 0}},
                                                              false,
                                                              0,
                                                              1));
    }

    ModuleNeighbor::MpiGhostExchangeStats stats;
    const std::vector<ModuleNeighbor::MpiAtomRecord> ghosts = domain.exchange_ghost_atoms(local_records, &stats);
    const ModuleNeighbor::AtomPack pack = build_pack_from_records(local_records, ghosts);
    const std::vector<PairKey> local_keys = collect_pair_keys(build_pairs_from_pack(pack, cutoff));
    const std::vector<PairKey> global_keys = gather_pair_keys(local_keys, domain.cart_comm());

    std::vector<PairKey> expected;
    expected.push_back(PairKey(0, 0, 0, 1, 0, 0, 0));
    expected.push_back(PairKey(0, 1, 0, 0, 0, 0, 0));
    std::sort(expected.begin(), expected.end());

    EXPECT_EQ(global_keys, expected);
    EXPECT_EQ(count_ghost_atoms(pack), 1);
    EXPECT_EQ(stats.received_ghost_count, 1);
}

TEST(MpiNeighborPrototypeTest, GhostsPreservePeriodicShiftInAtomPackSearch)
{
    ensure_mpi_initialized();

    const double cutoff = 0.75;
    ModuleNeighbor::MpiDomain domain;
    domain.initialize(world_comm(),
                      std::array<double, 3>{{0.0, 0.0, 0.0}},
                      std::array<double, 3>{{8.0, 2.0, 2.0}},
                      cutoff,
                      true);

    if (domain.size() != 2)
    {
        GTEST_SKIP() << "This MPI neighbor prototype scenario requires exactly 2 ranks.";
    }

    std::vector<ModuleNeighbor::MpiAtomRecord> local_records;
    if (domain.rank() == 0)
    {
        local_records.push_back(ModuleNeighbor::MpiAtomRecord(0.2,
                                                              1.0,
                                                              1.0,
                                                              0,
                                                              std::array<int, 3>{{0, 0, 0}},
                                                              false,
                                                              0,
                                                              0));
    }
    else
    {
        local_records.push_back(ModuleNeighbor::MpiAtomRecord(7.8,
                                                              1.0,
                                                              1.0,
                                                              1,
                                                              std::array<int, 3>{{0, 0, 0}},
                                                              false,
                                                              0,
                                                              1));
    }

    ModuleNeighbor::MpiGhostExchangeStats stats;
    const std::vector<ModuleNeighbor::MpiAtomRecord> ghosts = domain.exchange_ghost_atoms(local_records, &stats);
    const ModuleNeighbor::AtomPack pack = build_pack_from_records(local_records, ghosts);
    const std::vector<PairKey> local_keys = collect_pair_keys(build_pairs_from_pack(pack, cutoff));
    const std::vector<PairKey> global_keys = gather_pair_keys(local_keys, domain.cart_comm());

    std::vector<PairKey> expected;
    expected.push_back(PairKey(0, 0, 0, 1, -1, 0, 0));
    expected.push_back(PairKey(0, 1, 0, 0, 1, 0, 0));
    std::sort(expected.begin(), expected.end());

    EXPECT_EQ(global_keys, expected);
    EXPECT_EQ(count_ghost_atoms(pack), 1);
    EXPECT_EQ(stats.received_ghost_count, 1);
}

TEST(MpiNeighborPrototypeTest, FourRanksKeepLocalCentersInAtomPackSearch)
{
    ensure_mpi_initialized();

    const double cutoff = 0.75;
    ModuleNeighbor::MpiDomain domain;
    domain.initialize(world_comm(),
                      std::array<double, 3>{{0.0, 0.0, 0.0}},
                      std::array<double, 3>{{8.0, 4.0, 2.0}},
                      cutoff,
                      false);

    if (domain.size() != 4 || domain.dims()[0] != 2)
    {
        GTEST_SKIP() << "This MPI neighbor prototype scenario requires 4 ranks split along x.";
    }

    const double fx = domain.coords()[0] == 0 ? domain.local_bounds().upper[0] - 0.025
                                              : domain.local_bounds().lower[0] + 0.025;
    const double fy = 0.5 * (domain.local_bounds().lower[1] + domain.local_bounds().upper[1]);
    const double fz = 0.5 * (domain.local_bounds().lower[2] + domain.local_bounds().upper[2]);
    const std::array<double, 3> cartesian = domain.fractional_to_cartesian(fx, fy, fz);

    std::vector<ModuleNeighbor::MpiAtomRecord> local_records;
    local_records.push_back(ModuleNeighbor::MpiAtomRecord(cartesian[0],
                                                          cartesian[1],
                                                          cartesian[2],
                                                          domain.rank(),
                                                          std::array<int, 3>{{0, 0, 0}},
                                                          false,
                                                          0,
                                                          domain.rank()));

    ModuleNeighbor::MpiGhostExchangeStats stats;
    const std::vector<ModuleNeighbor::MpiAtomRecord> ghosts = domain.exchange_ghost_atoms(local_records, &stats);
    const ModuleNeighbor::AtomPack pack = build_pack_from_records(local_records, ghosts);
    const std::vector<ModuleNeighbor::NeighborPair> local_pairs = build_pairs_from_pack(pack, cutoff);
    const std::vector<PairKey> local_keys = collect_pair_keys(local_pairs);
    const std::vector<PairKey> global_keys = gather_pair_keys(local_keys, domain.cart_comm());

    for (const ModuleNeighbor::NeighborPair& pair: local_pairs)
    {
        EXPECT_EQ(pair.center_type, 0);
        EXPECT_EQ(pair.center_natom, domain.rank());
        ASSERT_GE(pair.center_index, 0);
        EXPECT_FALSE(pack.is_ghost[pair.center_index]);
    }
    EXPECT_EQ(stats.received_ghost_count, static_cast<int>(ghosts.size()));
    EXPECT_EQ(stats.sent_payload_count, stats.sent_atom_count * 9);
    EXPECT_EQ(stats.received_payload_count, stats.received_ghost_count * 9);
    EXPECT_FALSE(global_keys.empty());
}

TEST(MpiNeighborPrototypeTest, TriclinicPeriodicSearchMatchesFullCellReference)
{
    ensure_mpi_initialized();

    const double cutoff = 0.75;
    const std::array<std::array<double, 3>, 3> lattice{{
        std::array<double, 3>{{8.0, 0.0, 0.0}},
        std::array<double, 3>{{2.0, 4.0, 0.0}},
        std::array<double, 3>{{1.0, 0.5, 3.0}},
    }};
    ModuleNeighbor::MpiDomain domain;
    domain.initialize_lattice(world_comm(), std::array<double, 3>{{0.0, 0.0, 0.0}}, lattice, cutoff, true);

    if (domain.size() != 2)
    {
        GTEST_SKIP() << "This triclinic neighbor scenario requires exactly 2 ranks.";
    }

    const double fx = domain.rank() == 0 ? 0.025 : 0.975;
    const std::array<double, 3> cartesian = domain.fractional_to_cartesian(fx, 0.5, 0.5);
    std::vector<ModuleNeighbor::MpiAtomRecord> local_records;
    local_records.push_back(ModuleNeighbor::MpiAtomRecord(cartesian[0],
                                                          cartesian[1],
                                                          cartesian[2],
                                                          domain.rank(),
                                                          std::array<int, 3>{{0, 0, 0}},
                                                          false,
                                                          0,
                                                          domain.rank()));

    const std::vector<ModuleNeighbor::MpiAtomRecord> ghosts = domain.exchange_ghost_atoms(local_records);
    const ModuleNeighbor::AtomPack pack = build_pack_from_records(local_records, ghosts);
    const std::vector<ModuleNeighbor::NeighborPair> local_pairs = build_pairs_from_pack(pack, cutoff);
    const std::vector<PairKey> global_keys = gather_pair_keys(collect_pair_keys(local_pairs), domain.cart_comm());

    // Build the same two-atom periodic cell without domain decomposition. The
    // explicit images make this a direct search reference for the distributed
    // local-center plus ghost result.
    ModuleNeighbor::AtomPack reference_pack;
    const std::array<double, 3> atom0 = domain.fractional_to_cartesian(0.025, 0.5, 0.5);
    const std::array<double, 3> atom1 = domain.fractional_to_cartesian(0.975, 0.5, 0.5);
    const std::array<double, 3> atom0_plus = domain.fractional_to_cartesian(1.025, 0.5, 0.5);
    const std::array<double, 3> atom1_minus = domain.fractional_to_cartesian(-0.025, 0.5, 0.5);
    reference_pack.append_atom(atom0[0], atom0[1], atom0[2], 0, 0, 0, 0, 0, 0, false);
    reference_pack.append_atom(atom1[0], atom1[1], atom1[2], 0, 1, 0, 0, 0, 1, false);
    reference_pack.append_atom(atom0_plus[0], atom0_plus[1], atom0_plus[2], 0, 0, 1, 0, 0, 0, false);
    reference_pack.append_atom(atom1_minus[0], atom1_minus[1], atom1_minus[2], 0, 1, -1, 0, 0, 1, false);
    const std::vector<PairKey> reference_keys = collect_pair_keys(build_pairs_from_pack(reference_pack, cutoff));
    EXPECT_EQ(global_keys, reference_keys);
    ASSERT_EQ(ghosts.size(), 1);
    for (const ModuleNeighbor::NeighborPair& pair: local_pairs)
    {
        ASSERT_GE(pair.center_index, 0);
        EXPECT_FALSE(pack.is_ghost[pair.center_index]);
        EXPECT_TRUE(pack.is_ghost[pair.neighbor_index]);
    }
}

TEST(MpiNeighborPrototypeTest, FourRankTriclinicSearchKeepsTopologyPairs)
{
    ensure_mpi_initialized();

    const double cutoff = 0.75;
    const std::array<std::array<double, 3>, 3> lattice{{
        std::array<double, 3>{{8.0, 0.0, 0.0}},
        std::array<double, 3>{{2.0, 4.0, 0.0}},
        std::array<double, 3>{{1.0, 0.5, 3.0}},
    }};
    ModuleNeighbor::MpiDomain domain;
    domain.initialize_lattice(world_comm(), std::array<double, 3>{{0.0, 0.0, 0.0}}, lattice, cutoff, true);

    if (domain.size() != 4 || domain.dims()[0] != 2)
    {
        GTEST_SKIP() << "This triclinic topology scenario requires 4 ranks split along x.";
    }

    const double fx = domain.coords()[0] == 0 ? domain.local_bounds().upper[0] - 0.025
                                              : domain.local_bounds().lower[0] + 0.025;
    const double fy = 0.5 * (domain.local_bounds().lower[1] + domain.local_bounds().upper[1]);
    const double fz = 0.5 * (domain.local_bounds().lower[2] + domain.local_bounds().upper[2]);
    const std::array<double, 3> cartesian = domain.fractional_to_cartesian(fx, fy, fz);
    std::vector<ModuleNeighbor::MpiAtomRecord> local_records;
    local_records.push_back(ModuleNeighbor::MpiAtomRecord(cartesian[0],
                                                          cartesian[1],
                                                          cartesian[2],
                                                          domain.rank(),
                                                          std::array<int, 3>{{0, 0, 0}},
                                                          false,
                                                          0,
                                                          domain.rank()));

    const std::vector<ModuleNeighbor::MpiAtomRecord> ghosts = domain.exchange_ghost_atoms(local_records);
    const ModuleNeighbor::AtomPack pack = build_pack_from_records(local_records, ghosts);
    const std::vector<ModuleNeighbor::NeighborPair> local_pairs = build_pairs_from_pack(pack, cutoff);
    const std::vector<PairKey> global_keys = gather_pair_keys(collect_pair_keys(local_pairs), domain.cart_comm());

    std::vector<PairKey> expected;
#ifdef __MPI
    for (int cy = 0; cy < domain.dims()[1]; ++cy)
    {
        for (int cz = 0; cz < domain.dims()[2]; ++cz)
        {
            std::array<int, 3> left_coords{{0, cy, cz}};
            std::array<int, 3> right_coords{{1, cy, cz}};
            int left_rank = -1;
            int right_rank = -1;
            MPI_Cart_rank(domain.cart_comm(), left_coords.data(), &left_rank);
            MPI_Cart_rank(domain.cart_comm(), right_coords.data(), &right_rank);
            expected.push_back(PairKey(0, left_rank, 0, right_rank, 0, 0, 0));
            expected.push_back(PairKey(0, right_rank, 0, left_rank, 0, 0, 0));
        }
    }
#endif
    std::sort(expected.begin(), expected.end());
    EXPECT_EQ(global_keys, expected);
    for (const ModuleNeighbor::NeighborPair& pair: local_pairs)
    {
        ASSERT_GE(pair.center_index, 0);
        ASSERT_GE(pair.neighbor_index, 0);
        EXPECT_FALSE(pack.is_ghost[pair.center_index]);
        EXPECT_TRUE(pack.is_ghost[pair.neighbor_index]);
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
#ifdef __MPI
    // These tests call MPI collectives inside the test body. Manage MPI here
    // instead of relying on gtest_main, otherwise mpirun can report an abnormal
    // exit even when all gtest assertions pass.
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
