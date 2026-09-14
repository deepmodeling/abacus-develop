#include "gtest/gtest.h"

#include "../para_rgrid_world.h"
#include "../para_world.h"

#include <vector>

TEST(ParaRgridWorldTest, SerialMode)
{
    const Parallel::ParaRgridWorld world(4, 5, 6);
    EXPECT_EQ(world.tag(), "rgrid");
    EXPECT_EQ(world.ncx(), 4);
    EXPECT_EQ(world.ncy(), 5);
    EXPECT_EQ(world.ncz(), 6);
    EXPECT_EQ(world.nczp(), 6);
    EXPECT_EQ(world.nrxx(), 120);
    EXPECT_EQ(world.numz(0), 6);
    EXPECT_EQ(world.startz(0), 0);
}

TEST(ParaRgridWorldTest, WhichproSerial)
{
    const Parallel::ParaRgridWorld world(2, 2, 8);
    for (int iz = 0; iz < 8; ++iz)
    {
        EXPECT_EQ(world.whichpro(iz), 0);
    }
}

TEST(ParaRgridWorldTest, Validity)
{
    const Parallel::ParaRgridWorld world(1, 1, 1);
    EXPECT_TRUE(world.valid());
}

// ===== Cross-domain operation tests (serial mode) =====

TEST(ParaRgridWorldTest, ReduceAcrossPoolsSerial)
{
    Parallel::ParaRgridWorld rgrid(2, 2, 4);
    auto fake_kmesh = Parallel::ParaWorld::serial("kmesh");

    std::vector<double> data(rgrid.nrxx(), 1.0);
    rgrid.reduce_across_pools(data.data(), fake_kmesh);

    // Serial mode: no-op, data unchanged
    for (int i = 0; i < rgrid.nrxx(); ++i)
    {
        EXPECT_DOUBLE_EQ(data[i], 1.0);
    }
}

TEST(ParaRgridWorldTest, BcastDataSerial)
{
    const int ncx = 2, ncy = 2, ncz = 4;
    Parallel::ParaRgridWorld rgrid(ncx, ncy, ncz);
    auto fake_comm = Parallel::ParaWorld::serial("comm");

    // Build global grid: value = ixy * ncz + iz
    std::vector<double> global(ncx * ncy * ncz);
    for (int ixy = 0; ixy < ncx * ncy; ++ixy)
    {
        for (int iz = 0; iz < ncz; ++iz)
        {
            global[ixy * ncz + iz] = ixy * ncz + iz;
        }
    }

    std::vector<double> local(rgrid.nrxx(), -1.0);
    rgrid.bcast_data(global.data(), local.data(), fake_comm);

    // Serial: local should have all z-planes
    for (int ixy = 0; ixy < ncx * ncy; ++ixy)
    {
        for (int iz = 0; iz < ncz; ++iz)
        {
            EXPECT_DOUBLE_EQ(local[ixy * ncz + iz], ixy * ncz + iz);
        }
    }
}

TEST(ParaRgridWorldTest, ReduceDataSerial)
{
    const int ncx = 2, ncy = 2, ncz = 4;
    Parallel::ParaRgridWorld rgrid(ncx, ncy, ncz);
    auto fake_comm = Parallel::ParaWorld::serial("comm");

    // Local grid: value = ixy * nczp + iz
    std::vector<double> local(rgrid.nrxx());
    for (int i = 0; i < rgrid.nrxx(); ++i)
    {
        local[i] = static_cast<double>(i);
    }

    std::vector<double> global(ncx * ncy * ncz, -1.0);
    rgrid.reduce_data(global.data(), local.data(), fake_comm);

    // Serial: global should match local (single process owns all z)
    for (int ixy = 0; ixy < ncx * ncy; ++ixy)
    {
        for (int iz = 0; iz < ncz; ++iz)
        {
            EXPECT_DOUBLE_EQ(global[ixy * ncz + iz], local[ixy * ncz + iz]);
        }
    }
}

TEST(ParaRgridWorldTest, ReduceAcrossPoolsInvalidWorld)
{
    Parallel::ParaRgridWorld rgrid(2, 2, 4);
    auto invalid = Parallel::ParaWorld::serial("");

    std::vector<double> data(rgrid.nrxx(), 5.0);
    rgrid.reduce_across_pools(data.data(), invalid);

    // Invalid world: no-op
    for (int i = 0; i < rgrid.nrxx(); ++i)
    {
        EXPECT_DOUBLE_EQ(data[i], 5.0);
    }
}
