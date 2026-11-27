#include "gtest/gtest.h"

#include "../para_world.h"

TEST(ParaWorldTest, SerialFactory)
{
    const Parallel::ParaWorld world = Parallel::ParaWorld::serial("pw");
    EXPECT_EQ(world.tag(), "pw");
    EXPECT_EQ(world.rank(), 0);
    EXPECT_EQ(world.size(), 1);
    EXPECT_TRUE(world.valid());
}

TEST(ParaWorldTest, EmptyTagIsInvalid)
{
    const Parallel::ParaWorld world = Parallel::ParaWorld::serial("");
    EXPECT_TRUE(world.tag().empty());
    EXPECT_EQ(world.rank(), 0);
    EXPECT_EQ(world.size(), 1);
    EXPECT_FALSE(world.valid());
}
