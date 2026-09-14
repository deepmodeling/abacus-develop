#include "gtest/gtest.h"

#include "../para_bdiff_ksame_world.h"

TEST(ParaBdiffKsameWorldTest, SerialMode)
{
    const Parallel::ParaBdiffKsameWorld world;
    EXPECT_EQ(world.tag(), "bdiff_ksame");
    EXPECT_EQ(world.my_bndgroup(), 0);
    EXPECT_EQ(world.nbndgroup(), 1);
    EXPECT_EQ(world.rank_in_bpgroup(), 0);
    EXPECT_EQ(world.nproc_in_bndgroup(), 1);
    EXPECT_TRUE(world.valid());
}

TEST(ParaBdiffKsameWorldTest, AliasesMatchBase)
{
    const Parallel::ParaBdiffKsameWorld world;
    EXPECT_EQ(world.rank_in_bpgroup(), world.rank());
    EXPECT_EQ(world.nproc_in_bndgroup(), world.size());
}
