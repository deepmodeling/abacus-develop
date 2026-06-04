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
class ScopedMpiInit
{
  public:
    ScopedMpiInit() : initialized_before_(false), should_finalize_(false)
    {
#ifdef __MPI
        int initialized = 0;
        MPI_Initialized(&initialized);
        initialized_before_ = initialized != 0;
        if (!initialized_before_)
        {
            int argc = 0;
            char** argv = nullptr;
            MPI_Init(&argc, &argv);
            should_finalize_ = true;
        }
#endif
    }

    ~ScopedMpiInit()
    {
#ifdef __MPI
        if (should_finalize_)
        {
            MPI_Finalize();
        }
#endif
    }

  private:
    bool initialized_before_;
    bool should_finalize_;
};
} // namespace

TEST(MpiDomainTest, CreatesCartesianDomainAndOwnsLocalPoint)
{
    ScopedMpiInit mpi;

    ModuleNeighbor::MpiDomain domain;
    domain.initialize(MPI_COMM_WORLD,
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
    ScopedMpiInit mpi;

    const double cutoff = 0.5;
    ModuleNeighbor::MpiDomain domain;
    domain.initialize(MPI_COMM_WORLD,
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
