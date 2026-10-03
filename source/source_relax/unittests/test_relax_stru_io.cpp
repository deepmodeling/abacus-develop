#include "../relax_stru_io.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>

using testing::HasSubstr;
using testing::Not;

namespace
{
ModuleBase::matrix identity_stress()
{
    ModuleBase::matrix stress(3, 3);
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            stress(i, j) = (i == j) ? 1.0 : 0.0;
        }
    }
    return stress;
}
} // namespace

// build_stru_header -------------------------------------------------------

TEST(RelaxStruIO, HeaderContainsVersionAndStep)
{
    const std::string header = relax_stru_io::build_stru_header(0, 1.0, identity_stress(), false);
    EXPECT_THAT(header, HasSubstr("# ABACUS version:"));
    EXPECT_THAT(header, HasSubstr("# Written at"));
    EXPECT_THAT(header, HasSubstr("# RELAX STEP 1,"));
    EXPECT_THAT(header, Not(HasSubstr("(FINAL)")));
}

TEST(RelaxStruIO, HeaderFinalStepLabel)
{
    const std::string header = relax_stru_io::build_stru_header(4, 1.0, identity_stress(), true);
    EXPECT_THAT(header, HasSubstr("# RELAX STEP 5 (FINAL),"));
}

TEST(RelaxStruIO, HeaderConvertsEnergyToEv)
{
    // etot in Ry, printed in eV via ModuleBase::Ry_to_eV
    const std::string header = relax_stru_io::build_stru_header(0, 2.0, identity_stress(), false);
    char expected[64];
    std::snprintf(expected, sizeof(expected), "%.8f eV", 2.0 * ModuleBase::Ry_to_eV);
    EXPECT_THAT(header, HasSubstr(expected));
}

TEST(RelaxStruIO, HeaderWritesThreeStressRowsInKbar)
{
    ModuleBase::matrix stress(3, 3);
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            stress(i, j) = 0.0;
        }
    }
    stress(0, 0) = 1.0;
    stress(1, 1) = 2.0;
    stress(2, 2) = 3.0;

    const std::string header = relax_stru_io::build_stru_header(0, 1.0, stress, false);
    // three "# Stress (kbar):" lines
    int count = 0;
    size_t pos = 0;
    while ((pos = header.find("# Stress (kbar):", pos)) != std::string::npos)
    {
        count++;
        pos += 1;
    }
    EXPECT_EQ(count, 3);
}

// need_orbital ------------------------------------------------------------

TEST(RelaxStruIO, NeedOrbitalForLcao)
{
    Input_para inp;
    inp.basis_type = "lcao";
    EXPECT_TRUE(relax_stru_io::need_orbital(inp));
}

TEST(RelaxStruIO, NeedOrbitalForLcaoInPw)
{
    Input_para inp;
    inp.basis_type = "lcao_in_pw";
    EXPECT_TRUE(relax_stru_io::need_orbital(inp));
}

TEST(RelaxStruIO, NeedOrbitalForPwOnlyWithNaoWfc)
{
    Input_para inp;
    inp.basis_type = "pw";
    inp.init_wfc = "atomic";
    EXPECT_FALSE(relax_stru_io::need_orbital(inp));

    inp.init_wfc = "nao";
    EXPECT_TRUE(relax_stru_io::need_orbital(inp));
}
