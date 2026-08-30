#include "source_lcao/module_ri/sternheimer_runtime_options.h"

#include <cstdlib>
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>

namespace
{

constexpr const char* kTestFlag = "ABACUS_STERNHEIMER_TEST_FLAG";

class SternheimerRuntimeOptionsTest : public ::testing::Test
{
  protected:
    void TearDown() override
    {
        unsetenv(kTestFlag);
    }
};

} // namespace

TEST_F(SternheimerRuntimeOptionsTest, UsesDefaultWhenEnvironmentIsMissing)
{
    unsetenv(kTestFlag);

    EXPECT_TRUE(ModuleRI::sternheimer_environment_flag(kTestFlag, true));
    EXPECT_FALSE(ModuleRI::sternheimer_environment_flag(kTestFlag, false));
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
