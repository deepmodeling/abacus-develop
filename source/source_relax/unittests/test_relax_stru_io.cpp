#include "../relax_stru_io.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <sstream>
#include <vector>

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

Input_para fully_evaluated_inp()
{
    Input_para inp;
    inp.cal_force = true;
    inp.cal_stress = true;
    return inp;
}

} // namespace

TEST(RelaxStruIO, HeaderContainsVersionAndStep)
{
    const Input_para inp = fully_evaluated_inp();
    const std::string header = relax_stru_io::build_stru_header(0, 1.0, identity_stress(), inp, false, true);
    EXPECT_THAT(header, HasSubstr("# ABACUS version:"));
    EXPECT_THAT(header, HasSubstr("# RELAX STEP 1,"));
    EXPECT_THAT(header, Not(HasSubstr("(FINAL)")));
}

TEST(RelaxStruIO, HeaderFinalStepLabel)
{
    const Input_para inp = fully_evaluated_inp();
    const std::string header = relax_stru_io::build_stru_header(4, 1.0, identity_stress(), inp, true, true);
    EXPECT_THAT(header, HasSubstr("# RELAX STEP 5 (FINAL),"));
}

TEST(RelaxStruIO, HeaderConvertsEnergyToEv)
{
    // etot in Ry, printed in eV via ModuleBase::Ry_to_eV
    const Input_para inp = fully_evaluated_inp();
    const std::string header = relax_stru_io::build_stru_header(0, 2.0, identity_stress(), inp, false, true);
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

    const Input_para inp = fully_evaluated_inp();
    const std::string header = relax_stru_io::build_stru_header(0, 1.0, stress, inp, false, true);
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

// cal_force / cal_stress header combinations -------------------------------

// Fixture providing a populated stress matrix and the helpers used by the
// header behavior tests. These tests check only build_stru_header output;
// driver-loop control flow is covered by the relax sync/nsync suites.
class HeaderBehaviorTest : public testing::Test
{
protected:
    Input_para inp;
    ModuleBase::matrix stress;
    const double etot = -10.5; // Ry

    void SetUp() override
    {
        stress.create(3, 3);
        stress(0, 0) = 1.0;
        stress(0, 1) = 0.1;
        stress(0, 2) = 0.2;
        stress(1, 0) = 0.1;
        stress(1, 1) = 1.0;
        stress(1, 2) = 0.3;
        stress(2, 0) = 0.2;
        stress(2, 1) = 0.3;
        stress(2, 2) = 1.0;
    }

    static std::vector<std::string> split_lines(const std::string& s)
    {
        std::vector<std::string> lines;
        std::istringstream iss(s);
        std::string line;
        while (std::getline(iss, line))
        {
            lines.push_back(line);
        }
        return lines;
    }

    static int count_substr(const std::string& s, const std::string& sub)
    {
        int count = 0;
        size_t pos = 0;
        while ((pos = s.find(sub, pos)) != std::string::npos)
        {
            ++count;
            pos += 1;
        }
        return count;
    }
};

TEST_F(HeaderBehaviorTest, CalForceStressTrue)
{
    inp.cal_force = true;
    inp.cal_stress = true;

    const std::string header = relax_stru_io::build_stru_header(0, etot, stress, inp, false, true);
    const auto lines = split_lines(header);

    ASSERT_EQ(lines.size(), 7u);
    EXPECT_THAT(lines[0], testing::HasSubstr("# ABACUS version:"));
    EXPECT_THAT(lines[1], testing::HasSubstr("# Written at"));
    EXPECT_THAT(lines[2], testing::HasSubstr("# RELAX STEP 1, Energy:"));
    EXPECT_THAT(lines[3], testing::HasSubstr("# Stress (kbar): "));
    EXPECT_THAT(lines[4], testing::HasSubstr("# Stress (kbar): "));
    EXPECT_THAT(lines[5], testing::HasSubstr("# Stress (kbar): "));
    EXPECT_THAT(lines[6], testing::HasSubstr("# NOTE: stress and forces computed for this geometry"));
    EXPECT_THAT(header, testing::Not(testing::HasSubstr("N/A")));
}

TEST_F(HeaderBehaviorTest, CalStressFalse)
{
    inp.cal_force = true;
    inp.cal_stress = false;

    const std::string header = relax_stru_io::build_stru_header(0, etot, stress, inp, false, true);
    const auto lines = split_lines(header);

    ASSERT_EQ(lines.size(), 7u);
    // Stress lines are pure N/A, no inline reason comment
    EXPECT_EQ(count_substr(header, "# Stress (kbar): N/A N/A N/A"), 3);
    EXPECT_THAT(lines[3], testing::Not(testing::HasSubstr("cal_stress")));
    EXPECT_THAT(lines[6], testing::HasSubstr("# NOTE: stress not computed (cal_stress=0)"));
}

TEST_F(HeaderBehaviorTest, CalForceFalse)
{
    inp.cal_force = false;
    inp.cal_stress = true;

    const std::string header = relax_stru_io::build_stru_header(0, etot, stress, inp, false, true);
    const auto lines = split_lines(header);

    ASSERT_EQ(lines.size(), 7u);
    EXPECT_THAT(header, testing::Not(testing::HasSubstr("N/A")));
    EXPECT_THAT(lines[6], testing::HasSubstr("# NOTE: forces not computed (cal_force=0); per-atom f fields omitted intentionally"));
}

TEST_F(HeaderBehaviorTest, CalForceStressBothFalse)
{
    inp.cal_force = false;
    inp.cal_stress = false;

    const std::string header = relax_stru_io::build_stru_header(0, etot, stress, inp, false, true);
    const auto lines = split_lines(header);

    ASSERT_EQ(lines.size(), 7u);
    EXPECT_EQ(count_substr(header, "# Stress (kbar): N/A N/A N/A"), 3);
    EXPECT_THAT(lines[6], testing::HasSubstr("# NOTE: stress not computed (cal_stress=0); forces not computed (cal_force=0)"));
}

TEST_F(HeaderBehaviorTest, FinalStep)
{
    inp.cal_force = true;
    inp.cal_stress = true;

    const std::string header = relax_stru_io::build_stru_header(4, etot, stress, inp, true, true);
    const auto lines = split_lines(header);

    ASSERT_EQ(lines.size(), 7u);
    EXPECT_THAT(lines[2], testing::HasSubstr("# RELAX STEP 5 (FINAL), Energy:"));
}

TEST_F(HeaderBehaviorTest, GeometryNotEvaluated)
{
    inp.cal_force = true;
    inp.cal_stress = true;

    const std::string header = relax_stru_io::build_stru_header(4, etot, stress, inp, true, false);
    const auto lines = split_lines(header);

    ASSERT_EQ(lines.size(), 7u);
    EXPECT_THAT(lines[2], testing::HasSubstr("Energy: N/A"));
    EXPECT_EQ(count_substr(header, "# Stress (kbar): N/A N/A N/A"), 3);
    EXPECT_THAT(lines[6], testing::HasSubstr("# NOTE: geometry proposed by optimizer but not evaluated"));
    EXPECT_THAT(lines[6], testing::HasSubstr("energy N/A"));
    EXPECT_THAT(lines[6], testing::HasSubstr("stress N/A"));
    EXPECT_THAT(lines[6], testing::HasSubstr("forces omitted"));
    EXPECT_THAT(lines[6], testing::Not(testing::HasSubstr("energy above belongs to the last evaluated geometry")));
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
