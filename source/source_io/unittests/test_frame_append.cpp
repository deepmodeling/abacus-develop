#include "gtest/gtest.h"
#include "source_io/module_output/frame_append.h"

#include <stdexcept>

TEST(FrameAppend, ModeOneAppendsAfterTheFirstStep)
{
    // Step 0 defines the fresh run: the file is truncated, never appended.
    EXPECT_FALSE(ModuleIO::frame_appends(1, 0));
    EXPECT_TRUE(ModuleIO::frame_appends(1, 1));
    EXPECT_TRUE(ModuleIO::frame_appends(1, 5));
}

TEST(FrameAppend, ModesZeroAndTwoNeverAppend)
{
    EXPECT_FALSE(ModuleIO::frame_appends(0, 0));
    EXPECT_FALSE(ModuleIO::frame_appends(0, 3));
    EXPECT_FALSE(ModuleIO::frame_appends(2, 0));
    EXPECT_FALSE(ModuleIO::frame_appends(2, 3));
}

TEST(FrameAppend, RejectsInvalidArguments)
{
    EXPECT_THROW(ModuleIO::frame_appends(-1, 0), std::invalid_argument);
    EXPECT_THROW(ModuleIO::frame_appends(3, 1), std::invalid_argument);
    EXPECT_THROW(ModuleIO::frame_appends(1, -1), std::invalid_argument);
}
