/**
 * @file diago_mpi_test.cpp
 * @brief Unit tests for MPI communication helpers (nbcast, nreduce_pool).
 *
 * Tests:
 *   1. MPI communication correctness (broadcast, reduce, edge cases)
 *   2. CommStrategy configuration
 */

#include "source_hsolver/mpi_comm_helper.h"
#include "gtest/gtest.h"
#include "mpi.h"

#include <vector>
#include <complex>
#include <cstdlib>

using namespace hsolver;

// =========================================================================
// Helper: Get MPI rank/size
// =========================================================================

static void getMpiInfo(int& rank, int& nproc) {
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#else
    rank = 0;
    nproc = 1;
#endif
}

// =========================================================================
// Test Fixture
// =========================================================================

class DiagoMPICorrectnessTest : public ::testing::Test {
protected:
    void SetUp() override {
        getMpiInfo(rank_, nproc_);
    }

    int rank_ = 0;
    int nproc_ = 1;
};

// =========================================================================
// Test 1: MPI communication correctness
// =========================================================================

TEST_F(DiagoMPICorrectnessTest, CommunicationCorrectness) {
#ifdef __MPI
    // 1. Broadcast
    {
        double val = (rank_ == 0) ? 42.0 : 0.0;
        MPIRequestTracker tracker;
        MPICommHelper::nbcast(&val, 1, 0, MPI_COMM_WORLD, tracker);
        tracker.wait_all();
        EXPECT_EQ(val, 42.0);
    }

    // 2. Reduce (sum)
    {
        const int N = 100;
        std::vector<double> data(N, static_cast<double>(rank_));
        MPIRequestTracker tracker;
        MPICommHelper::nreduce_pool(data.data(), N, MPI_COMM_WORLD, tracker);
        tracker.wait_all();

        double expected = nproc_ * (nproc_ - 1.0) / 2.0;
        for (int i = 0; i < N; i++) {
            EXPECT_NEAR(data[i], expected, 1e-10)
                << "Reduce result mismatch at index " << i;
        }
    }

    // 3. Edge cases: empty operations
    {
        MPIRequestTracker tracker;
        MPICommHelper::nbcast(static_cast<double*>(nullptr), 0, 0,
                              MPI_COMM_WORLD, tracker);
        tracker.wait_all();
        EXPECT_FALSE(tracker.has_pending());
    }
    {
        MPIRequestTracker tracker;
        std::complex<double> dummy;
        MPICommHelper::nreduce_pool(&dummy, 0, MPI_COMM_WORLD, tracker);
        tracker.wait_all();
        EXPECT_FALSE(tracker.has_pending());
    }
#endif
}

// =========================================================================
// Test 2: CommStrategy Configuration
// =========================================================================

TEST_F(DiagoMPICorrectnessTest, CommStrategyConfiguration) {
    // Adaptive: small problem -> kNonBlocking
    EXPECT_EQ(hsolver::resolve_comm_strategy(hsolver::CommStrategy::kAdaptive,
                                              100, 10),
              hsolver::CommStrategy::kNonBlocking);

    // Adaptive: large problem -> kPipelined
    EXPECT_EQ(hsolver::resolve_comm_strategy(hsolver::CommStrategy::kAdaptive,
                                              1000, 500),
              hsolver::CommStrategy::kPipelined);

    // Explicit override
    EXPECT_EQ(hsolver::resolve_comm_strategy(hsolver::CommStrategy::kBlocking,
                                              1000, 500),
              hsolver::CommStrategy::kBlocking);

    // Default non-blocking
    EXPECT_EQ(hsolver::resolve_comm_strategy(
                  hsolver::CommStrategy::kNonBlocking, 100, 10),
              hsolver::CommStrategy::kNonBlocking);
}

// =========================================================================
// Main
// =========================================================================

int main(int argc, char** argv) {
#ifdef __MPI
    const char* ompi_size = getenv("OMPI_COMM_WORLD_SIZE");
    const char* pmi_size  = getenv("PMI_SIZE");
    if (!ompi_size && !pmi_size) {
        std::cout << "MPI test skipped: not running under mpirun" << std::endl;
        return 0;
    }
    MPI_Init(&argc, &argv);
#endif

    ::testing::InitGoogleTest(&argc, argv);

    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    return result;
}
