#include "gtest/gtest.h"

#include "../para_pw_world.h"

TEST(ParaPwWorldTest, SerialMode)
{
    const Parallel::ParaPwWorld world(128);
    EXPECT_EQ(world.tag(), "pw");
    EXPECT_EQ(world.npw(), 128);
    EXPECT_EQ(world.npwtot(), 128);
    EXPECT_EQ(world.poolnproc(), 1);
    EXPECT_EQ(world.poolrank(), 0);
    EXPECT_EQ(world.npw_per(0), 128);
}

TEST(ParaPwWorldTest, SerialZero)
{
    const Parallel::ParaPwWorld world(0);
    EXPECT_EQ(world.npw(), 0);
    EXPECT_EQ(world.npwtot(), 0);
}

TEST(ParaPwWorldTest, Validity)
{
    const Parallel::ParaPwWorld world(64);
    EXPECT_TRUE(world.valid());
    EXPECT_EQ(world.size(), 1);
}
