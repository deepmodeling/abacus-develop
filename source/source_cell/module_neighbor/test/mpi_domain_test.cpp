#include "../mpi_domain.h"

#include "gtest/gtest.h"

#ifdef __MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <array>
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

std::array<double, 3> local_fractional_point(const ModuleNeighbor::MpiDomain& domain,
                                             const double upper_x_offset = 0.0)
{
    std::array<double, 3> fractional;
    fractional[0] = upper_x_offset > 0.0
                        ? std::max(domain.local_bounds().lower[0],
                                   domain.local_bounds().upper[0] - upper_x_offset)
                        : 0.5 * (domain.local_bounds().lower[0] + domain.local_bounds().upper[0]);
    fractional[1] = 0.5 * (domain.local_bounds().lower[1] + domain.local_bounds().upper[1]);
    fractional[2] = 0.5 * (domain.local_bounds().lower[2] + domain.local_bounds().upper[2]);
    return fractional;
}
} // namespace

TEST(MpiDomainTest, CreatesCartesianDomainAndOwnsLocalPoint)
{
    ensure_mpi_initialized();

    ModuleNeighbor::MpiDomain domain;
    domain.initialize(world_comm(),
                      std::array<double, 3>{{0.0, 0.0, 0.0}},
                      std::array<double, 3>{{8.0, 4.0, 2.0}},
                      0.5,
                      true);

    EXPECT_TRUE(domain.initialized());
    EXPECT_GE(domain.size(), 1);
    EXPECT_EQ(domain.periods()[0], 1);
    EXPECT_EQ(domain.periods()[1], 1);
    EXPECT_EQ(domain.periods()[2], 1);
    EXPECT_GT(domain.local_bounds().upper[0], domain.local_bounds().lower[0]);
    EXPECT_GT(domain.local_bounds().upper[1], domain.local_bounds().lower[1]);
    EXPECT_GT(domain.local_bounds().upper[2], domain.local_bounds().lower[2]);

    const std::array<double, 3> fractional = local_fractional_point(domain);
    const std::array<double, 3> cartesian
        = domain.fractional_to_cartesian(fractional[0], fractional[1], fractional[2]);
    EXPECT_TRUE(domain.owns(cartesian[0], cartesian[1], cartesian[2]));
}

TEST(MpiDomainTest, SelectsLocalAtomsAndCountsGhostCandidates)
{
    ensure_mpi_initialized();

    const double cutoff = 0.5;
    ModuleNeighbor::MpiDomain domain;
    domain.initialize(world_comm(),
                      std::array<double, 3>{{0.0, 0.0, 0.0}},
                      std::array<double, 3>{{8.0, 4.0, 2.0}},
                      cutoff,
                      true);

    const std::array<double, 3> fractional
        = local_fractional_point(domain, 0.25 * domain.fractional_ghost_padding()[0]);
    const std::array<double, 3> cartesian
        = domain.fractional_to_cartesian(fractional[0], fractional[1], fractional[2]);

    std::vector<ModuleNeighbor::MpiAtomRecord> atoms;
    atoms.push_back(ModuleNeighbor::MpiAtomRecord(cartesian[0], cartesian[1], cartesian[2], domain.rank()));

    const std::vector<int> local_indices = domain.select_local_atoms(atoms);
    EXPECT_EQ(local_indices.size(), 1);

    ModuleNeighbor::MpiGhostExchangeStats stats;
    const std::vector<ModuleNeighbor::MpiAtomRecord> ghosts = domain.exchange_ghost_atoms(atoms, &stats);
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

TEST(MpiDomainTest, ReportsNeighborExchangeStats)
{
    ensure_mpi_initialized();

    const double cutoff = 0.5;
    ModuleNeighbor::MpiDomain domain;
    domain.initialize(world_comm(),
                      std::array<double, 3>{{0.0, 0.0, 0.0}},
                      std::array<double, 3>{{8.0, 4.0, 2.0}},
                      cutoff,
                      false);

    const std::array<double, 3> fractional
        = local_fractional_point(domain, 0.25 * domain.fractional_ghost_padding()[0]);
    const std::array<double, 3> cartesian
        = domain.fractional_to_cartesian(fractional[0], fractional[1], fractional[2]);
    std::vector<ModuleNeighbor::MpiAtomRecord> atoms;
    atoms.push_back(ModuleNeighbor::MpiAtomRecord(cartesian[0],
                                                  cartesian[1],
                                                  cartesian[2],
                                                  domain.rank(),
                                                  std::array<int, 3>{{0, 0, 0}},
                                                  false,
                                                  0,
                                                  domain.rank()));

    ModuleNeighbor::MpiGhostExchangeStats stats;
    const std::vector<ModuleNeighbor::MpiAtomRecord> ghosts = domain.exchange_ghost_atoms(atoms, &stats);

    EXPECT_EQ(stats.received_ghost_count, static_cast<int>(ghosts.size()));
    EXPECT_EQ(stats.sent_payload_count, stats.sent_atom_count * 9);
    EXPECT_EQ(stats.received_payload_count, stats.received_ghost_count * 9);
    EXPECT_LE(stats.neighbor_rank_count, std::min(26, domain.size() - 1));
    EXPECT_LE(stats.nonempty_send_rank_count, stats.neighbor_rank_count);

#ifdef __MPI
    int local_sent_atoms = stats.sent_atom_count;
    int local_received_ghosts = stats.received_ghost_count;
    int global_sent_atoms = 0;
    int global_received_ghosts = 0;
    MPI_Allreduce(&local_sent_atoms, &global_sent_atoms, 1, MPI_INT, MPI_SUM, domain.cart_comm());
    MPI_Allreduce(&local_received_ghosts, &global_received_ghosts, 1, MPI_INT, MPI_SUM, domain.cart_comm());
    EXPECT_EQ(global_sent_atoms, global_received_ghosts);
    if (domain.size() > 1)
    {
        EXPECT_GT(global_sent_atoms, 0);
    }
#else
    EXPECT_EQ(stats.neighbor_rank_count, 0);
    EXPECT_TRUE(ghosts.empty());
#endif
}

TEST(MpiDomainTest, ConvertsAndOwnsTriclinicCoordinates)
{
    ensure_mpi_initialized();

    const std::array<double, 3> origin{{1.0, -2.0, 0.5}};
    const std::array<std::array<double, 3>, 3> lattice{{
        std::array<double, 3>{{8.0, 0.0, 0.0}},
        std::array<double, 3>{{2.0, 4.0, 0.0}},
        std::array<double, 3>{{1.0, 0.5, 3.0}},
    }};
    ModuleNeighbor::MpiDomain domain;
    domain.initialize_lattice(world_comm(), origin, lattice, 0.75, true);

    const std::array<double, 3> fractional = local_fractional_point(domain);
    const std::array<double, 3> cartesian
        = domain.fractional_to_cartesian(fractional[0], fractional[1], fractional[2]);
    const std::array<double, 3> recovered
        = domain.cartesian_to_fractional(cartesian[0], cartesian[1], cartesian[2]);

    EXPECT_NEAR(recovered[0], fractional[0], 1.0e-12);
    EXPECT_NEAR(recovered[1], fractional[1], 1.0e-12);
    EXPECT_NEAR(recovered[2], fractional[2], 1.0e-12);
    EXPECT_TRUE(domain.owns(cartesian[0], cartesian[1], cartesian[2]));
    EXPECT_GT(domain.fractional_ghost_padding()[0], 0.0);
    EXPECT_GT(domain.fractional_ghost_padding()[1], 0.0);
    EXPECT_GT(domain.fractional_ghost_padding()[2], 0.0);
}

TEST(MpiDomainTest, ExchangesTriclinicPeriodicImages)
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
        GTEST_SKIP() << "This triclinic image scenario requires exactly 2 ranks.";
    }

    const double fx = domain.rank() == 0 ? 0.025 : 0.975;
    const std::array<double, 3> cartesian = domain.fractional_to_cartesian(fx, 0.5, 0.5);
    std::vector<ModuleNeighbor::MpiAtomRecord> atoms;
    atoms.push_back(ModuleNeighbor::MpiAtomRecord(cartesian[0],
                                                  cartesian[1],
                                                  cartesian[2],
                                                  domain.rank(),
                                                  std::array<int, 3>{{0, 0, 0}},
                                                  false,
                                                  0,
                                                  domain.rank()));

    ModuleNeighbor::MpiGhostExchangeStats stats;
    const std::vector<ModuleNeighbor::MpiAtomRecord> ghosts = domain.exchange_ghost_atoms(atoms, &stats);
    ASSERT_EQ(ghosts.size(), 1);
    const int expected_shift = domain.rank() == 0 ? -1 : 1;
    EXPECT_EQ(ghosts[0].pbc_shift[0], expected_shift);
    EXPECT_EQ(ghosts[0].pbc_shift[1], 0);
    EXPECT_EQ(ghosts[0].pbc_shift[2], 0);
    EXPECT_EQ(stats.sent_atom_count, 1);
    EXPECT_EQ(stats.received_ghost_count, 1);

    const std::array<double, 3> ghost_fractional
        = domain.cartesian_to_fractional(ghosts[0].x, ghosts[0].y, ghosts[0].z);
    EXPECT_NEAR(ghost_fractional[0], domain.rank() == 0 ? -0.025 : 1.025, 1.0e-12);
    EXPECT_NEAR(ghost_fractional[1], 0.5, 1.0e-12);
    EXPECT_NEAR(ghost_fractional[2], 0.5, 1.0e-12);
}

TEST(MpiDomainTest, RejectsInvalidLattice)
{
    ensure_mpi_initialized();

    ModuleNeighbor::MpiDomain domain;
    const std::array<std::array<double, 3>, 3> singular{{
        std::array<double, 3>{{1.0, 0.0, 0.0}},
        std::array<double, 3>{{2.0, 0.0, 0.0}},
        std::array<double, 3>{{0.0, 0.0, 1.0}},
    }};
    EXPECT_THROW(domain.initialize_lattice(world_comm(),
                                           std::array<double, 3>{{0.0, 0.0, 0.0}},
                                           singular,
                                           0.5,
                                           true),
                 std::invalid_argument);
    EXPECT_THROW(domain.initialize_lattice(world_comm(),
                                           std::array<double, 3>{{0.0, 0.0, 0.0}},
                                           std::array<std::array<double, 3>, 3>{{
                                               std::array<double, 3>{{1.0, 0.0, 0.0}},
                                               std::array<double, 3>{{0.0, 1.0, 0.0}},
                                               std::array<double, 3>{{0.0, 0.0, 1.0}},
                                           }},
                                           -0.5,
                                           true),
                 std::invalid_argument);
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
