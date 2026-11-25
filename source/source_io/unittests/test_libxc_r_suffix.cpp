#include "gtest/gtest.h"
#include "source_io/module_chgpot/libxc_r_suffix.h"

#include <stdexcept>
#include <string>
#include <vector>

namespace
{

std::vector<std::string> suffixes(const std::string& data_name, const int nspin)
{
    return ModuleIO::libxc_r_channel_suffixes(data_name, nspin);
}

} // namespace

TEST(LibxcRSuffix, SpinQuantities)
{
    EXPECT_EQ(suffixes("rho", 1), (std::vector<std::string>{"s1"}));
    EXPECT_EQ(suffixes("rho", 2), (std::vector<std::string>{"s1", "s2"}));
    EXPECT_EQ(suffixes("vrho", 2), (std::vector<std::string>{"s1", "s2"}));
}

TEST(LibxcRSuffix, ScalarQuantities)
{
    EXPECT_EQ(suffixes("exc", 1), (std::vector<std::string>{""}));
    EXPECT_EQ(suffixes("amag", 2), (std::vector<std::string>{""}));
}

TEST(LibxcRSuffix, Tensor3Quantities)
{
    const std::vector<std::string> polarized = {"00", "01", "11"};
    EXPECT_EQ(suffixes("sigma", 1), (std::vector<std::string>{"00"}));
    EXPECT_EQ(suffixes("sigma", 2), polarized);
    EXPECT_EQ(suffixes("vsigma", 2), polarized);
    EXPECT_EQ(suffixes("v2rho2", 2), polarized);
}

TEST(LibxcRSuffix, SecondOrderQuantities)
{
    EXPECT_EQ(suffixes("v2rhosigma", 1), (std::vector<std::string>{"0_00"}));
    EXPECT_EQ(suffixes("v2rhosigma", 2),
              (std::vector<std::string>{"0_00", "0_01", "0_11", "1_00", "1_01", "1_11"}));
    EXPECT_EQ(suffixes("v2sigma2", 1), (std::vector<std::string>{"00_00"}));
    EXPECT_EQ(suffixes("v2sigma2", 2),
              (std::vector<std::string>{"00_00", "00_01", "00_11", "01_01", "01_11", "11_11"}));
}

TEST(LibxcRSuffix, ThirdOrderQuantities)
{
    EXPECT_EQ(suffixes("v3rho3", 1), (std::vector<std::string>{"0_0_0"}));
    EXPECT_EQ(suffixes("v3rho3", 2),
              (std::vector<std::string>{"0_0_0", "0_0_1", "0_1_1", "1_1_1"}));

    const std::vector<std::string> rho2sigma = suffixes("v3rho2sigma", 2);
    ASSERT_EQ(rho2sigma.size(), 9);
    EXPECT_EQ(rho2sigma.front(), "00_00");
    EXPECT_EQ(rho2sigma.back(), "11_11");

    const std::vector<std::string> rhosigma2 = suffixes("v3rhosigma2", 2);
    ASSERT_EQ(rhosigma2.size(), 12);
    EXPECT_EQ(rhosigma2.front(), "0_00_00");
    EXPECT_EQ(rhosigma2.back(), "1_11_11");

    const std::vector<std::string> sigma3 = suffixes("v3sigma3", 2);
    ASSERT_EQ(sigma3.size(), 10);
    EXPECT_EQ(sigma3.front(), "00_00_00");
    EXPECT_EQ(sigma3.back(), "11_11_11");

    EXPECT_EQ(suffixes("v3rho2sigma", 1), (std::vector<std::string>{"00_00"}));
    EXPECT_EQ(suffixes("v3rhosigma2", 1), (std::vector<std::string>{"0_00_00"}));
    EXPECT_EQ(suffixes("v3sigma3", 1), (std::vector<std::string>{"00_00_00"}));
}

TEST(LibxcRSuffix, FourthOrderQuantities)
{
    EXPECT_EQ(suffixes("v4rho4", 1), (std::vector<std::string>{"0_0_0_0"}));
    EXPECT_EQ(suffixes("v4rho4", 2).size(), 5);

    const std::vector<std::string> rho3sigma = suffixes("v4rho3sigma", 2);
    ASSERT_EQ(rho3sigma.size(), 12);
    EXPECT_EQ(rho3sigma.front(), "000_00");
    EXPECT_EQ(rho3sigma.back(), "111_11");

    const std::vector<std::string> rho2sigma2 = suffixes("v4rho2sigma2", 2);
    ASSERT_EQ(rho2sigma2.size(), 18);
    EXPECT_EQ(rho2sigma2.front(), "00_00_00");
    EXPECT_EQ(rho2sigma2.back(), "11_11_11");

    const std::vector<std::string> rhosigma3 = suffixes("v4rhosigma3", 2);
    ASSERT_EQ(rhosigma3.size(), 20);
    EXPECT_EQ(rhosigma3.front(), "0_00_00_00");
    EXPECT_EQ(rhosigma3.back(), "1_11_11_11");

    const std::vector<std::string> sigma4 = suffixes("v4sigma4", 2);
    ASSERT_EQ(sigma4.size(), 15);
    EXPECT_EQ(sigma4.front(), "00_00_00_00");
    EXPECT_EQ(sigma4.back(), "11_11_11_11");

    EXPECT_EQ(suffixes("v4rho3sigma", 1), (std::vector<std::string>{"000_00"}));
    EXPECT_EQ(suffixes("v4rho2sigma2", 1), (std::vector<std::string>{"00_00_00"}));
    EXPECT_EQ(suffixes("v4rhosigma3", 1), (std::vector<std::string>{"0_00_00_00"}));
    EXPECT_EQ(suffixes("v4sigma4", 1), (std::vector<std::string>{"00_00_00_00"}));
}

TEST(LibxcRSuffix, RejectsInvalidArguments)
{
    EXPECT_THROW(ModuleIO::libxc_r_channel_suffixes("rho", 0), std::invalid_argument);
    EXPECT_THROW(ModuleIO::libxc_r_channel_suffixes("rho", 3), std::invalid_argument);
    EXPECT_THROW(ModuleIO::libxc_r_channel_suffixes("unknown", 1), std::invalid_argument);
}
