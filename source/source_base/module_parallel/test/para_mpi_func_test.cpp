#include "gtest/gtest.h"

#include "../para_mpi_func.h"
#include "../para_world.h"

TEST(ParaMpiFuncTest, BcastIntSerial)
{
    auto world = Parallel::ParaWorld::serial("test");
    int v = 42;
    Parallel::bcast_int(v, world);
    EXPECT_EQ(v, 42);
}

TEST(ParaMpiFuncTest, BcastDoubleSerial)
{
    auto world = Parallel::ParaWorld::serial("test");
    double v = 3.14;
    Parallel::bcast_double(v, world);
    EXPECT_DOUBLE_EQ(v, 3.14);
}

TEST(ParaMpiFuncTest, BcastBoolSerial)
{
    auto world = Parallel::ParaWorld::serial("test");
    bool v = true;
    Parallel::bcast_bool(v, world);
    EXPECT_TRUE(v);
}

TEST(ParaMpiFuncTest, BcastStringSerial)
{
    auto world = Parallel::ParaWorld::serial("test");
    std::string s = "hello";
    Parallel::bcast_string(s, world);
    EXPECT_EQ(s, "hello");
}

TEST(ParaMpiFuncTest, BcastIntArraySerial)
{
    auto world = Parallel::ParaWorld::serial("test");
    int v[] = {1, 2, 3};
    Parallel::bcast_int(v, 3, world);
    EXPECT_EQ(v[0], 1);
    EXPECT_EQ(v[1], 2);
    EXPECT_EQ(v[2], 3);
}

TEST(ParaMpiFuncTest, BcastComplexSerial)
{
    auto world = Parallel::ParaWorld::serial("test");
    std::complex<double> v(1.0, 2.0);
    Parallel::bcast_complex(v, world);
    EXPECT_DOUBLE_EQ(v.real(), 1.0);
    EXPECT_DOUBLE_EQ(v.imag(), 2.0);
}

TEST(ParaMpiFuncTest, BcastCharArraySerial)
{
    auto world = Parallel::ParaWorld::serial("test");
    char buf[] = "abc";
    Parallel::bcast_char(buf, 4, world);
    EXPECT_EQ(buf[0], 'a');
}

TEST(ParaMpiFuncTest, ReduceAllSerial)
{
    auto world = Parallel::ParaWorld::serial("test");
    double v = 5.0;
    Parallel::reduce_all(v, world);
    EXPECT_DOUBLE_EQ(v, 5.0);
}

TEST(ParaMpiFuncTest, ReduceMinMaxSerial)
{
    auto world = Parallel::ParaWorld::serial("test");
    double v = 7.5;
    Parallel::reduce_min(v, world);
    EXPECT_DOUBLE_EQ(v, 7.5);
    Parallel::reduce_max(v, world);
    EXPECT_DOUBLE_EQ(v, 7.5);
}

TEST(ParaMpiFuncTest, GatherIntSerial)
{
    auto world = Parallel::ParaWorld::serial("test");
    int v = 99;
    int all[1] = {0};
    Parallel::gather_int(v, all, world);
    EXPECT_EQ(all[0], 99);
}

TEST(ParaMpiFuncTest, InvalidWorldIsNoop)
{
    auto world = Parallel::ParaWorld::serial("");
    EXPECT_FALSE(world.valid());
    int v = 42;
    Parallel::bcast_int(v, world);
    EXPECT_EQ(v, 42); // unchanged
    Parallel::reduce_all(v, world);
    EXPECT_EQ(v, 42); // unchanged
    Parallel::barrier(world); // should not hang or crash
}
