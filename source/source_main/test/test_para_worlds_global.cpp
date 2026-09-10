#include "gtest/gtest.h"

#include "../para_worlds_global.h"

#include "source_base/module_parallel/para_collection.h"
#include "source_base/module_parallel/para_tag.h"

using namespace Parallel;

// The global holder is a process-wide singleton; reset it before and after
// each case so the state machine can be exercised repeatedly in one binary.
class ParaWorldsGlobalTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        reset_global_para_worlds_for_test();
    }
    void TearDown() override
    {
        reset_global_para_worlds_for_test();
    }
};

// nimage = 1 on a single process: both domains exist and are trivially sized.
TEST_F(ParaWorldsGlobalTest, InitRegistersEsolverAndImagesDomains)
{
    ParaCollection& collection = init_global_para_worlds(1, 0, 1);
    EXPECT_EQ(collection.size(), 2u);

    const ParaCollection& fetched = global_para_worlds();
    EXPECT_EQ(&fetched, &collection);

    const ParaWorld& esolver = fetched.find(ParaTag::esolver);
    const ParaWorld& images = fetched.find(ParaTag::images);
    EXPECT_EQ(esolver.size(), 1);
    EXPECT_EQ(esolver.rank(), 0);
    EXPECT_EQ(images.size(), 1);
    EXPECT_EQ(images.rank(), 0);
}

// A second initialization after reset must succeed and yield a fresh object.
TEST_F(ParaWorldsGlobalTest, ResetAllowsReinitialization)
{
    ParaCollection& first = init_global_para_worlds(1, 0, 1);
    EXPECT_EQ(first.size(), 2u);

    reset_global_para_worlds_for_test();
    // Re-initialization must succeed (the latch was cleared by reset) and the
    // fresh collection must be usable. The container address is not asserted:
    // the allocator is free to reuse the just-freed storage.
    const ParaCollection& second = init_global_para_worlds(1, 0, 1);
    EXPECT_EQ(second.size(), 2u);
    EXPECT_EQ(second.find(ParaTag::esolver).size(), 1);
}

// Invalid image counts are rejected instead of producing a bad split.
TEST_F(ParaWorldsGlobalTest, RejectsNimageOutOfRange)
{
    EXPECT_DEATH(init_global_para_worlds(1, 0, 0), ".*");
    reset_global_para_worlds_for_test();
    EXPECT_DEATH(init_global_para_worlds(2, 0, 4), ".*");
    reset_global_para_worlds_for_test();
}
