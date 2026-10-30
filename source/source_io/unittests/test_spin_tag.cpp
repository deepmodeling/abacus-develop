#include "source_io/module_output/spin_tag.h"

#include "gtest/gtest.h"

#include <string>

/**
 * Tests for ModuleIO::make_spin_tag, the running-log suffix shared by the
 * H(R), DM, Vxc(R), band and cube writers.
 */

TEST(SpinTagTest, CollinearSpinChannels)
{
    EXPECT_EQ(ModuleIO::make_spin_tag(0, 2), " (spin up  )");
    EXPECT_EQ(ModuleIO::make_spin_tag(1, 2), " (spin down)");
}

TEST(SpinTagTest, NoTagForSingleSpinOrNoncollinear)
{
    EXPECT_EQ(ModuleIO::make_spin_tag(0, 1), "");
    EXPECT_EQ(ModuleIO::make_spin_tag(0, 4), "");
    EXPECT_EQ(ModuleIO::make_spin_tag(3, 4), "");
}

TEST(SpinTagTest, NoTagForSpinSummedCall)
{
    EXPECT_EQ(ModuleIO::make_spin_tag(-1, 2), "");
}
