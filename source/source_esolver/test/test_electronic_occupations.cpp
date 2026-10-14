#include "../electronic_occupations.h"
#include "source_base/matrix.h"
#include "source_base/parallel_global.h"

#include "gtest/gtest.h"

#include <limits>
#include <stdexcept>

namespace
{
int process_count = 1;
int process_rank = 0;

ModuleBase::matrix fractional_weights()
{
    ModuleBase::matrix weights(2, 3);
    weights(0, 0) = 2.0;
    weights(0, 1) = 2.0;
    weights(0, 2) = 2.0;
    weights(1, 0) = 2.0;
    weights(1, 1) = 0.2;
    return weights;
}

TEST(ElectronicOccupations, FractionalPopulationIncludesEveryLocalRowAndBand)
{
    const ModuleBase::matrix weights = fractional_weights();
    EXPECT_NEAR(8.2, ModuleESolver::checked_occupation_electron_count(weights, 8.2, 1, process_count), 1e-14);
}

TEST(ElectronicOccupations, ReturnsOccupationSumInsteadOfEchoingNearbyRequest)
{
    const ModuleBase::matrix weights = fractional_weights();
    const double requested = 8.2 + 5e-9;
    const double actual = ModuleESolver::checked_occupation_electron_count(weights, requested, 1, process_count);
    EXPECT_NEAR(8.2, actual, 1e-14);
    EXPECT_GT(requested - actual, 4e-9);
}

TEST(ElectronicOccupations, RejectsPopulationMismatchOnEitherSide)
{
    const ModuleBase::matrix weights = fractional_weights();
    EXPECT_THROW(ModuleESolver::checked_occupation_electron_count(weights, 8.2 + 2e-8, 1, process_count),
                 std::runtime_error);
    EXPECT_THROW(ModuleESolver::checked_occupation_electron_count(weights, 8.2 - 2e-8, 1, process_count),
                 std::runtime_error);
}

TEST(ElectronicOccupations, ToleranceIsAbsoluteEvenForLargeElectronCounts)
{
    ModuleBase::matrix weights(1, 1);
    weights(0, 0) = 1e6;
    EXPECT_THROW(ModuleESolver::checked_occupation_electron_count(weights, 1e6 + 1e-6, 1, process_count),
                 std::runtime_error);
}

TEST(ElectronicOccupations, GeneralizedSmearingCanContainNegativeWeights)
{
    ModuleBase::matrix weights(1, 3);
    weights(0, 0) = 2.125;
    weights(0, 1) = -0.25;
    weights(0, 2) = 0.125;
    EXPECT_DOUBLE_EQ(2.0, ModuleESolver::checked_occupation_electron_count(weights, 2.0, 1, process_count));
}

TEST(ElectronicOccupations, RejectsNaNOnLastRankTogether)
{
    ModuleBase::matrix weights = fractional_weights();
    if (process_rank == process_count - 1)
    {
        weights(1, 2) = std::numeric_limits<double>::quiet_NaN();
    }
    EXPECT_THROW(ModuleESolver::checked_occupation_electron_count(weights, 8.2, 1, process_count),
                 std::runtime_error);
}

TEST(ElectronicOccupations, RejectsInfinityOnLastRankTogether)
{
    ModuleBase::matrix weights = fractional_weights();
    if (process_rank == process_count - 1)
    {
        weights(1, 2) = std::numeric_limits<double>::infinity();
    }
    EXPECT_THROW(ModuleESolver::checked_occupation_electron_count(weights, 8.2, 1, process_count),
                 std::runtime_error);
}

TEST(ElectronicOccupations, RejectsLocalOverflowBeforePoolReduction)
{
    ModuleBase::matrix weights(1, 2);
    weights(0, 0) = std::numeric_limits<double>::max();
    weights(0, 1) = std::numeric_limits<double>::max();
    EXPECT_THROW(ModuleESolver::checked_occupation_electron_count(weights, 8.2, process_count, 1),
                 std::runtime_error);
}

TEST(ElectronicOccupations, RejectsInvalidRequestOnLastRankTogether)
{
    const ModuleBase::matrix weights = fractional_weights();
    const double requested = process_rank == process_count - 1 ? std::numeric_limits<double>::quiet_NaN() : 8.2;
    EXPECT_THROW(ModuleESolver::checked_occupation_electron_count(weights, requested, 1, process_count),
                 std::runtime_error);
}

TEST(ElectronicOccupations, RejectsInvalidPartitionOnLastRankTogether)
{
    const ModuleBase::matrix weights = fractional_weights();
    const int local_size = process_rank == process_count - 1 ? 0 : process_count;
    EXPECT_THROW(ModuleESolver::checked_occupation_electron_count(weights, 8.2, 1, local_size), std::runtime_error);
    EXPECT_THROW(ModuleESolver::checked_occupation_electron_count(weights, 8.2, 0, process_count), std::runtime_error);
}

TEST(ElectronicOccupations, SinglePoolDoesNotMultiplyReplicatedPopulationByRankCount)
{
    const ModuleBase::matrix weights = fractional_weights();
    EXPECT_NEAR(8.2, ModuleESolver::checked_occupation_electron_count(weights, 8.2, 1, process_count), 1e-14);
}

TEST(ElectronicOccupations, SinglePoolMismatchOnOneReplicaFailsEveryRank)
{
    ModuleBase::matrix weights = fractional_weights();
    if (process_rank == process_count - 1)
    {
        weights(1, 1) += 0.125;
    }
    EXPECT_THROW(ModuleESolver::checked_occupation_electron_count(weights, 8.2, 1, process_count),
                 std::runtime_error);
}

TEST(ElectronicOccupations, KPointPoolsCombineDistinctRowsAndRemovePoolReplicas)
{
    const int pools = process_count % 2 == 0 ? 2 : 1;
    const int ranks_per_pool = process_count / pools;
    const int pool = process_rank / ranks_per_pool;
    ModuleBase::matrix weights(1, 2);
    const double population = pools == 1 ? 8.2 : (pool == 0 ? 2.2 : 6.0);
    weights(0, 0) = 0.25 * population;
    weights(0, 1) = 0.75 * population;
    EXPECT_NEAR(8.2, ModuleESolver::checked_occupation_electron_count(weights, 8.2, pools, ranks_per_pool), 1e-14);
}

TEST(ElectronicOccupations, BandShardsUseLocalColumnCount)
{
    // Independent band shards are reduction groups, just like native BNDPAR.
    const int local_bands = process_rank + 1;
    ModuleBase::matrix weights(2, local_bands);
    const double total_bands = process_count * (process_count + 1) / 2.0;
    for (int band = 0; band < local_bands; ++band)
    {
        weights(0, band) = 3.0 / total_bands;
        weights(1, band) = 5.2 / total_bands;
    }
    EXPECT_NEAR(8.2, ModuleESolver::checked_occupation_electron_count(weights, 8.2, process_count, 1), 1e-14);
}

TEST(ElectronicOccupations, EmptyLocalBlocksStillParticipateInCollectives)
{
    ModuleBase::matrix weights;
    if (process_rank == 0)
    {
        weights = fractional_weights();
    }
    EXPECT_NEAR(8.2, ModuleESolver::checked_occupation_electron_count(weights, 8.2, process_count, 1), 1e-14);
}

TEST(ElectronicOccupations, EmptyPopulationCannotSatisfyPositiveRequest)
{
    const ModuleBase::matrix weights;
    EXPECT_THROW(ModuleESolver::checked_occupation_electron_count(weights, 8.2, process_count, 1),
                 std::runtime_error);
}

TEST(ElectronicOccupations, UnequalPoolSizesRemoveTheirOwnReplicaCounts)
{
    // One rank in the first pool and all remaining ranks in the second.
    const int pools = process_count > 1 ? 2 : 1;
    const int ranks_per_pool = process_rank == 0 ? 1 : process_count - 1;
    ModuleBase::matrix weights(1, 1);
    weights(0, 0) = pools == 1 ? 8.2 : (process_rank == 0 ? 2.2 : 6.0);
    EXPECT_NEAR(8.2, ModuleESolver::checked_occupation_electron_count(weights, 8.2, pools, ranks_per_pool), 1e-14);
}

TEST(ElectronicOccupations, RemoteInvalidEntryCannotBeHiddenByAnEmptyRootBlock)
{
    ModuleBase::matrix weights;
    if (process_rank == process_count - 1)
    {
        weights.create(1, 1);
        weights(0, 0) = -std::numeric_limits<double>::infinity();
    }
    EXPECT_THROW(ModuleESolver::checked_occupation_electron_count(weights, 8.2, process_count, 1),
                 std::runtime_error);
}
} // namespace

int main(int argc, char** argv)
{
    int thread_count = 1;
    Parallel_Global::read_pal_param(argc, argv, process_count, thread_count, process_rank);
#ifdef __MPI
    // Follow source_base/test/mpi_test_main.cpp: no pool communicators are
    // created by these tests, so native finalization must not free them.
    POOL_WORLD = MPI_COMM_NULL;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif
    return result;
}
