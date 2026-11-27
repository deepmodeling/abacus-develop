#include "gtest/gtest.h"

#include "../para_world.h"

TEST(ParaWorldMpiTest, SerialFactory)
{
    const Parallel::ParaWorld world = Parallel::ParaWorld::serial("pw");
    EXPECT_EQ(world.tag(), "pw");
    EXPECT_EQ(world.rank(), 0);
    EXPECT_EQ(world.size(), 1);
    EXPECT_TRUE(world.valid());
}

TEST(ParaWorldMpiTest, WrapCommunicator)
{
    const Parallel::ParaWorld world = Parallel::ParaWorld::serial("pw");
    EXPECT_EQ(world.tag(), "pw");
    EXPECT_TRUE(world.valid());
    EXPECT_EQ(world.rank(), 0);
    EXPECT_EQ(world.size(), 1);
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
