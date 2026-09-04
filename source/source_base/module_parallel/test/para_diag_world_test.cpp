#include "gtest/gtest.h"

#include "../para_diag_world.h"

TEST(ParaDiagWorldTest, SerialMode)
{
    const Parallel::ParaDiagWorld world;
    EXPECT_EQ(world.tag(), "diag");
    EXPECT_EQ(world.drank(), 0);
    EXPECT_EQ(world.dsize(), 1);
    EXPECT_EQ(world.dcolor(), 0);
    EXPECT_TRUE(world.valid());
}

TEST(ParaDiagWorldTest, AliasesMatchBase)
{
    const Parallel::ParaDiagWorld world;
    EXPECT_EQ(world.drank(), world.rank());
    EXPECT_EQ(world.dsize(), world.size());
}
