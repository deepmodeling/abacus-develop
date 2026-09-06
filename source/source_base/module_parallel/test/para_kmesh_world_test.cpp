#include "gtest/gtest.h"

#include "../para_kmesh_world.h"

TEST(ParaKmeshWorldTest, SerialSinglePool)
{
    const Parallel::ParaKmeshWorld world(4, 1);
    EXPECT_EQ(world.tag(), "kmesh");
    EXPECT_EQ(world.kpar(), 1);
    EXPECT_EQ(world.my_pool(), 0);
    EXPECT_EQ(world.rank_in_pool(), 0);
    EXPECT_EQ(world.nproc(), 1);
    EXPECT_EQ(world.nspin(), 1);
    EXPECT_EQ(world.nkstot(), 4);
    EXPECT_EQ(world.nks_local(), 4);
    EXPECT_EQ(world.startk_global(), 0);
}

TEST(ParaKmeshWorldTest, EvenDistribution)
{
    const Parallel::ParaKmeshWorld world(6, 1);
    // serial: kpar=1, so all 6 k-points in pool 0
    EXPECT_EQ(world.nks_pool(0), 6);
    EXPECT_EQ(world.startk_pool(0), 0);
    EXPECT_EQ(world.max_nks_pool(), 6);
}

TEST(ParaKmeshWorldTest, WhichPool)
{
    const Parallel::ParaKmeshWorld world(5, 1);
    for (int ik = 0; ik < 5; ++ik)
    {
        EXPECT_EQ(world.which_pool(ik), 0);
    }
}

TEST(ParaKmeshWorldTest, PoolCollectionSerial)
{
    const Parallel::ParaKmeshWorld world(3, 1);
    const double wk[] = {0.5, 0.3, 0.2};
    double value = 0.0;
    world.pool_collection(value, wk, 1);
    EXPECT_DOUBLE_EQ(value, 0.3);
}

TEST(ParaKmeshWorldTest, PoolCollectionArraySerial)
{
    const Parallel::ParaKmeshWorld world(2, 1);
    // 2 k-points, 3 elements each: k0={1,2,3}, k1={4,5,6}
    const double w[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double value[3] = {0.0, 0.0, 0.0};
    world.pool_collection(value, w, 3, 1);
    EXPECT_DOUBLE_EQ(value[0], 4.0);
    EXPECT_DOUBLE_EQ(value[1], 5.0);
    EXPECT_DOUBLE_EQ(value[2], 6.0);
}

TEST(ParaKmeshWorldTest, GatherKvecSerial)
{
    const Parallel::ParaKmeshWorld world(2, 1);
    const std::vector<double> local = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0};
    std::vector<double> global;
    world.gather_kvec(local, global);
    ASSERT_EQ(global.size(), 6u);
    EXPECT_DOUBLE_EQ(global[0], 1.0);
    EXPECT_DOUBLE_EQ(global[4], 1.0);
}

TEST(ParaKmeshWorldTest, Nspin2Restriction)
{
    // nspin=2 forces pool=0 for pool_collection (legacy behavior)
    const Parallel::ParaKmeshWorld world(4, 2);
    EXPECT_EQ(world.nspin(), 2);
}
