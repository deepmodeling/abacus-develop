#include "gtest/gtest.h"

#include "../para_matrix_world.h"

TEST(ParaMatrixWorldTest, SerialMode)
{
    const Parallel::ParaMatrixWorld world;
    EXPECT_EQ(world.tag(), "matrix");
    EXPECT_EQ(world.dim0(), 1);
    EXPECT_EQ(world.dim1(), 1);
    EXPECT_EQ(world.coord_row(), 0);
    EXPECT_EQ(world.coord_col(), 0);
    EXPECT_TRUE(world.valid());
}

TEST(ParaMatrixWorldTest, GridProductMatchesSize)
{
    const Parallel::ParaMatrixWorld world;
    EXPECT_EQ(world.dim0() * world.dim1(), world.size());
}
