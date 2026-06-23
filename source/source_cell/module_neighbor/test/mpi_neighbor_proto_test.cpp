#include "../atom_pack.h"
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

    const double x_local = std::max(domain.local_bounds().lower[0], domain.local_bounds().upper[0] - 0.25 * cutoff);
    const double y_local = 0.5 * (domain.local_bounds().lower[1] + domain.local_bounds().upper[1]);
    const double z_local = 0.5 * (domain.local_bounds().lower[2] + domain.local_bounds().upper[2]);

    std::vector<ModuleNeighbor::MpiAtomRecord> local_records;
    local_records.push_back(ModuleNeighbor::MpiAtomRecord(x_local, y_local, z_local, domain.rank()));

    ModuleNeighbor::AtomPack pack;
    pack.append_mpi_record(local_records.front(), 0, domain.rank());

    const std::vector<ModuleNeighbor::MpiAtomRecord> ghosts = domain.exchange_ghost_atoms(local_records);
    for (const ModuleNeighbor::MpiAtomRecord& ghost: ghosts)
    {
        pack.append_mpi_record(ghost);
    }

    EXPECT_EQ(pack.size(), static_cast<int>(ghosts.size()) + 1);
    EXPECT_EQ(count_ghost_atoms(pack), static_cast<int>(ghosts.size()));

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
