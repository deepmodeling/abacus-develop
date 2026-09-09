#ifndef __MPI
#define MODULE_BASE_TEST_DEFINED_MPI
#define __MPI
#endif
#include "source_base/global_variable.h"
#include "source_base/parallel_global.h"
#ifdef MODULE_BASE_TEST_DEFINED_MPI
#undef __MPI
#endif

#include "gtest/gtest.h"

int main(int argc, char** argv)
{
    int process_count = 1;
    int process_rank = 0;
    int thread_count = 1;
    Parallel_Global::read_pal_param(argc, argv, process_count, thread_count, process_rank);
    POOL_WORLD = MPI_COMM_NULL;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    Parallel_Global::finalize_mpi();
    return result;
}
