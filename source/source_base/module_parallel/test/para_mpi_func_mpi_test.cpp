#include "gtest/gtest.h"

#include "../para_mpi_func.h"
#include "../para_world.h"

TEST(ParaMpiFuncMpiTest, BcastIntFromRoot)
{
    auto world = Parallel::ParaWorld::serial("test");
    int v = (world.rank() == 0) ? 99 : 0;
    Parallel::bcast_int(v, world);
    EXPECT_EQ(v, 99);
}

TEST(ParaMpiFuncMpiTest, ReduceAllSum)
{
    auto world = Parallel::ParaWorld::serial("test");
    int v = 1; // each rank contributes 1
    Parallel::reduce_all(v, world);
    EXPECT_EQ(v, 1); // serial: size=1
}

TEST(ParaMpiFuncMpiTest, GatherIntAll)
{
    auto world = Parallel::ParaWorld::serial("test");
    int v = world.rank();
    int all[1] = {0};
    Parallel::gather_int(v, all, world);
    EXPECT_EQ(all[0], 0);
}

TEST(ParaMpiFuncMpiTest, BarrierNoHang)
{
    auto world = Parallel::ParaWorld::serial("test");
    Parallel::barrier(world);
    SUCCEED();
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
