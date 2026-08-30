#include "source_lcao/module_ri/sternheimer_runtime_options.h"

#include <cstdlib>
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>

namespace
{

constexpr const char* kTestFlag = "ABACUS_STERNHEIMER_TEST_FLAG";
constexpr const char* kBatchWidth = "ABACUS_STERNHEIMER_CHANNEL_BATCH_WIDTH";

class SternheimerRuntimeOptionsTest : public ::testing::Test
{
  protected:
    void TearDown() override
    {
        unsetenv(kTestFlag);
        unsetenv(kBatchWidth);
    }
};

} // namespace

TEST_F(SternheimerRuntimeOptionsTest, UsesDefaultWhenEnvironmentIsMissing)
{
    unsetenv(kTestFlag);

    EXPECT_TRUE(ModuleRI::sternheimer_environment_flag(kTestFlag, true));
    EXPECT_FALSE(ModuleRI::sternheimer_environment_flag(kTestFlag, false));
}

TEST_F(SternheimerRuntimeOptionsTest, UsesTwoChannelBatchesByDefaultAndOneForRollback)
{
    unsetenv(kBatchWidth);
    EXPECT_EQ(ModuleRI::sternheimer_channel_batch_width(), 2);

    setenv(kBatchWidth, "1", 1);
    EXPECT_EQ(ModuleRI::sternheimer_channel_batch_width(), 1);
    setenv(kBatchWidth, "8", 1);
    EXPECT_EQ(ModuleRI::sternheimer_channel_batch_width(), 8);
}

TEST_F(SternheimerRuntimeOptionsTest, RejectsInvalidChannelBatchWidths)
{
    for (const char* value: {"", "0", "-1", "four", "65"})
    {
        setenv(kBatchWidth, value, 1);
        EXPECT_THROW(ModuleRI::sternheimer_channel_batch_width(), std::invalid_argument) << value;
    }
}

TEST_F(SternheimerRuntimeOptionsTest, ParsesRecognizedBooleanTextCaseInsensitively)
{
    for (const char* value: {"1", "true", "TRUE", "on", "On", "yes", "YES"})
    {
        setenv(kTestFlag, value, 1);
        EXPECT_TRUE(ModuleRI::sternheimer_environment_flag(kTestFlag, false)) << value;
    }
    for (const char* value: {"0", "false", "FALSE", "off", "Off", "no", "NO"})
    {
        setenv(kTestFlag, value, 1);
        EXPECT_FALSE(ModuleRI::sternheimer_environment_flag(kTestFlag, true)) << value;
    }
}

TEST_F(SternheimerRuntimeOptionsTest, RejectsEmptyAndUnknownBooleanText)
{
    for (const char* value: {"", "maybe", "2"})
    {
        setenv(kTestFlag, value, 1);
        EXPECT_THROW(ModuleRI::sternheimer_environment_flag(kTestFlag, true), std::invalid_argument) << value;
    }
}
