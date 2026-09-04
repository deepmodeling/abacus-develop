#include "gtest/gtest.h"

#include "../para_setup.h"
#include "../para_world.h"
#include "../para_collection.h"

TEST(ParaSetupTest, DivideMpiGroupsSerial)
{
    int procs_in_group, my_group, rank_in_group;
    Parallel::divide_mpi_groups(1, 1, 0, false, procs_in_group, my_group, rank_in_group);
    EXPECT_EQ(procs_in_group, 1);
    EXPECT_EQ(my_group, 0);
    EXPECT_EQ(rank_in_group, 0);
}

TEST(ParaSetupTest, DivideMpiGroupsEven)
{
    // 8 procs, 4 groups -> 2 procs per group
    int procs_in_group, my_group, rank_in_group;
    Parallel::divide_mpi_groups(8, 4, 3, true, procs_in_group, my_group, rank_in_group);
    EXPECT_EQ(procs_in_group, 2);
    EXPECT_EQ(my_group, 1);
    EXPECT_EQ(rank_in_group, 1);
}

TEST(ParaSetupTest, DivideMpiGroupsUneven)
{
    // 7 procs, 3 groups -> group 0,1 have 3 procs, group 2 has 1
    // rank 0-5 -> group 0,1 (3 procs each)
    // rank 6 -> group 2 (1 proc)
    int procs_in_group, my_group, rank_in_group;

    // rank 2: first group (procs_in_group+1=3), 2/3=0, 2%3=2
    Parallel::divide_mpi_groups(7, 3, 2, false, procs_in_group, my_group, rank_in_group);
    EXPECT_EQ(procs_in_group, 3);
    EXPECT_EQ(my_group, 0);
    EXPECT_EQ(rank_in_group, 2);

    // rank 6: (6-1)/2=2, (6-1)%2=1 -> wait, extra_procs=1, procs_in_group=2
    // rank 6 >= 1*3 = 3, so: (6-1)/2=2, (6-1)%2=1
    Parallel::divide_mpi_groups(7, 3, 6, false, procs_in_group, my_group, rank_in_group);
    EXPECT_EQ(procs_in_group, 2);
    EXPECT_EQ(my_group, 2);
    EXPECT_EQ(rank_in_group, 1);
}

TEST(ParaSetupTest, DivideMpiGroupsRank0)
{
    int procs_in_group, my_group, rank_in_group;
    Parallel::divide_mpi_groups(12, 4, 0, true, procs_in_group, my_group, rank_in_group);
    EXPECT_EQ(procs_in_group, 3);
    EXPECT_EQ(my_group, 0);
    EXPECT_EQ(rank_in_group, 0);
}
