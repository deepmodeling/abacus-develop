#include "source_lcao/module_ri/sternheimer_runtime_options.h"

#include "source_io/module_parameter/input_parameter.h"

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

TEST(SternheimerMolecularCoulomb, RequiresExplicitIsolatedContractAndMatchingBasis)
{
    Input_para input;
    EXPECT_FALSE(ModuleRI::validate_sternheimer_molecular_coulomb(input));
    input.sternheimer_ao_potential_file = "ao.dat";
    EXPECT_THROW(ModuleRI::validate_sternheimer_molecular_coulomb(input), std::invalid_argument);
    input.sternheimer_molecular_coulomb = "isolated_ri";
    input.basis_type = "lcao";
    input.sternheimer_delta = true;
    input.out_sternheimer_librpa = true;
    input.symmetry = "-1";
    input.exx_singularity_correction = "limits";
    EXPECT_TRUE(ModuleRI::validate_sternheimer_molecular_coulomb(input));
    for (const auto member: {&Input_para::exx_rotate_abfs,
                             &Input_para::exx_coul_moment,
                             &Input_para::out_sternheimer_siab,
                             &Input_para::cal_force})
    {
        input.*member = true;
        EXPECT_THROW(ModuleRI::validate_sternheimer_molecular_coulomb(input), std::invalid_argument);
        input.*member = false;
    }
    input.sternheimer_q_index = 1;
    EXPECT_THROW(ModuleRI::validate_sternheimer_molecular_coulomb(input), std::invalid_argument);
    input.sternheimer_q_index = 0;
    input.exx_singularity_correction = "massidda";
    EXPECT_THROW(ModuleRI::validate_sternheimer_molecular_coulomb(input), std::invalid_argument);
    input.exx_singularity_correction = "limits";
    input.out_librpa_3d_coulomb_method = "direct_reciprocal";
    EXPECT_THROW(ModuleRI::validate_sternheimer_molecular_coulomb(input), std::invalid_argument);
    input.out_librpa_3d_coulomb_method = "ewald";
    input.shrink_abfs_pca_thr = 1.e-6;
    EXPECT_THROW(ModuleRI::validate_sternheimer_molecular_coulomb(input), std::invalid_argument);
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
