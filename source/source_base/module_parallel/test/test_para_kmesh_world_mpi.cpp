#include "gtest/gtest.h"

#include "../para_bdiff_ksame_world.h"
#include "../para_kmesh_world.h"

// Run with: mpirun -np 4 ./MODULE_BASE_para_kmesh_world_mpi
//
// The sum reduction protocol has two layers:
//   1. ParaBdiffKsameWorld::reduce_across_bdiff_ksame (band dimension, BPCG shards)
//   2. ParaKmeshWorld::reduce_across_pools      (k dimension, one
//      contribution per k-pool: the first rank of each k-pool injects the
//      partial sum, everyone else injects zero)
// The band layer must run first so that the k layer receives one complete
// per-k-pool partial sum.

TEST(ParaKmeshWorldMpiTest, ReduceAcrossBandGroupsBndpar2)
{
    int nprocs = 0;
    int myrank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    ASSERT_EQ(nprocs, 4);

    // kpar=1, bndpar=2, 4 ranks: band group = myrank/2, rank position
    // inside the band group = myrank%2. Reproduce the BP_WORLD layout,
    // which links the same rank position of every band group.
    MPI_Comm bp_world = MPI_COMM_NULL;
    MPI_Comm_split(MPI_COMM_WORLD, myrank % 2, myrank / 2, &bp_world);
    Parallel::ParaBdiffKsameWorld bdiff(MPI_COMM_WORLD, bp_world, 2);

    // Each band group holds a partial occupation sum of 14 (28 electrons
    // split into two band windows).
    double sumk = 14.0;
    bdiff.reduce_across_bdiff_ksame(sumk);
    EXPECT_DOUBLE_EQ(sumk, 28.0);

    // max/min stay world-wide (idempotent) and must span the band groups
    // as well: the two shards see different eigenvalue windows.
    Parallel::ParaKmeshWorld kmesh(MPI_COMM_WORLD, 1, 0, 0, 1);
    double eup = (myrank < 2) ? 40.0 : 45.0;
    kmesh.reduce_max_across_pools(eup);
    EXPECT_DOUBLE_EQ(eup, 45.0);

    double elw = (myrank < 2) ? -1.0 : -5.0;
    kmesh.reduce_min_across_pools(elw);
    EXPECT_DOUBLE_EQ(elw, -5.0);

    MPI_Comm_free(&bp_world);
}

TEST(ParaKmeshWorldMpiTest, ReduceAcrossKpoolsKpar2)
{
    int nprocs = 0;
    int myrank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    ASSERT_EQ(nprocs, 4);

    // kpar=2, bndpar=1: k-pool 0 = ranks {0,1}, k-pool 1 = ranks {2,3}
    // (consecutive rank blocks, divide_mpi_groups layout).
    const int my_pool = myrank / 2;
    Parallel::ParaKmeshWorld kmesh(MPI_COMM_WORLD, 2, my_pool, 4, 1);
    EXPECT_EQ(kmesh.startpro_pool(0), 0);
    EXPECT_EQ(kmesh.startpro_pool(1), 2);
    EXPECT_EQ(kmesh.kpool_root(), (myrank % 2 == 0));

    // Every process holds its k-pool's partial sum; the reduction must
    // count each pool exactly once.
    double sumk = 3.5;
    kmesh.reduce_across_pools(sumk);
    EXPECT_DOUBLE_EQ(sumk, 7.0);
}

TEST(ParaKmeshWorldMpiTest, UnevenKpoolsKpar3)
{
    int nprocs = 0;
    int myrank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    ASSERT_EQ(nprocs, 4);

    // nproc=4, kpar=3 (the 007_PW_UPF201_USPP_Fe layout): k-pool sizes
    // are [2,1,1]. divide_mpi_groups puts ranks {0,1} in pool 0, rank 2
    // in pool 1 and rank 3 in pool 2.
    int my_pool = 0;
    if (myrank >= 3)
    {
        my_pool = 2;
    }
    else if (myrank >= 2)
    {
        my_pool = 1;
    }
    Parallel::ParaKmeshWorld kmesh(MPI_COMM_WORLD, 3, my_pool, 0, 1);
    EXPECT_EQ(kmesh.startpro_pool(0), 0);
    EXPECT_EQ(kmesh.startpro_pool(1), 2);
    EXPECT_EQ(kmesh.startpro_pool(2), 3);
    EXPECT_EQ(kmesh.kpool_root(), (myrank != 1));

    // Every process holds its k-pool's partial sum. The reduction must
    // count each pool exactly once even though the pools are uneven:
    // the legacy average-pool-size division (4/3 = 1) double-counted
    // pool 0 here and corrupted the electron count / Fermi level.
    const double pool_sum = (my_pool == 0) ? 10.0 : ((my_pool == 1) ? 20.0 : 30.0);
    double sumk = pool_sum;
    kmesh.reduce_across_pools(sumk);
    EXPECT_DOUBLE_EQ(sumk, 60.0);
}

TEST(ParaKmeshWorldMpiTest, SinglePoolIsNoOp)
{
    int nprocs = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ASSERT_EQ(nprocs, 4);

    // kpar == 1: the sum reduction must be a no-op regardless of the
    // world size (the band dimension is handled by ParaBdiffKsameWorld).
    Parallel::ParaKmeshWorld kmesh(MPI_COMM_WORLD, 1, 0, 0, 1);
    double sumk = 42.0;
    kmesh.reduce_across_pools(sumk);
    EXPECT_DOUBLE_EQ(sumk, 42.0);
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
