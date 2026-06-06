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

    const double x_mid = 0.5 * (domain.local_bounds().lower[0] + domain.local_bounds().upper[0]);
    const double y_mid = 0.5 * (domain.local_bounds().lower[1] + domain.local_bounds().upper[1]);
    const double z_mid = 0.5 * (domain.local_bounds().lower[2] + domain.local_bounds().upper[2]);
    EXPECT_TRUE(domain.owns(x_mid, y_mid, z_mid));
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

    const double x_local = std::max(domain.local_bounds().lower[0], domain.local_bounds().upper[0] - 0.25 * cutoff);
    const double y_local = 0.5 * (domain.local_bounds().lower[1] + domain.local_bounds().upper[1]);
    const double z_local = 0.5 * (domain.local_bounds().lower[2] + domain.local_bounds().upper[2]);

    std::vector<ModuleNeighbor::MpiAtomRecord> atoms;
    atoms.push_back(ModuleNeighbor::MpiAtomRecord(x_local, y_local, z_local, domain.rank()));

    const std::vector<int> local_indices = domain.select_local_atoms(atoms);
    EXPECT_EQ(local_indices.size(), 1);

    const std::vector<ModuleNeighbor::MpiAtomRecord> ghosts = domain.exchange_ghost_atoms(atoms);

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
