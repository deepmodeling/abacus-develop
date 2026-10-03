#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "../relax_driver.h"
#include "source_base/matrix.h"
#include "source_io/module_parameter/input_parameter.h"

#include <sstream>
#include <vector>

// Test build_stru_header for stru_out / final_out header generation.
// The header logic is the only part of Relax_Driver that does not depend on
// PARAM.globalv or MPI rank, so it can be unit-tested in isolation.
//
// The header is always exactly 7 lines:
//   1: version
//   2: timestamp
//   3: relax step + energy
//   4-6: stress (3 rows, values or N/A)
//   7: NOTE describing the state of stress/forces/geometry

class RelaxDriverHeaderTest : public testing::Test
{
protected:
    Input_para inp;
    ModuleBase::matrix stress;
    const double etot = -10.5; // Ry

    void SetUp() override
    {
        stress.create(3, 3);
        stress(0, 0) = 1.0; stress(0, 1) = 0.1; stress(0, 2) = 0.2;
        stress(1, 0) = 0.1; stress(1, 1) = 1.0; stress(1, 2) = 0.3;
        stress(2, 0) = 0.2; stress(2, 1) = 0.3; stress(2, 2) = 1.0;
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

TEST_F(RelaxDriverHeaderTest, CalForceStressTrue)
{
    inp.cal_force = true;
    inp.cal_stress = true;

    const std::string header = Relax_Driver::build_stru_header(0, etot, stress, inp, false, true);
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

TEST_F(RelaxDriverHeaderTest, CalStressFalse)
{
    inp.cal_force = true;
    inp.cal_stress = false;

    const std::string header = Relax_Driver::build_stru_header(0, etot, stress, inp, false, true);
    const auto lines = split_lines(header);

    ASSERT_EQ(lines.size(), 7u);
    // Stress lines are pure N/A, no inline reason comment
    EXPECT_EQ(count_substr(header, "# Stress (kbar): N/A N/A N/A"), 3);
    EXPECT_THAT(lines[3], testing::Not(testing::HasSubstr("cal_stress")));
    EXPECT_THAT(lines[6], testing::HasSubstr("# NOTE: stress not computed (cal_stress=0)"));
}

TEST_F(RelaxDriverHeaderTest, CalForceFalse)
{
    inp.cal_force = false;
    inp.cal_stress = true;

    const std::string header = Relax_Driver::build_stru_header(0, etot, stress, inp, false, true);
    const auto lines = split_lines(header);

    ASSERT_EQ(lines.size(), 7u);
    EXPECT_THAT(header, testing::Not(testing::HasSubstr("N/A")));
    EXPECT_THAT(lines[6], testing::HasSubstr("# NOTE: forces not computed (cal_force=0); per-atom f fields omitted intentionally"));
}

TEST_F(RelaxDriverHeaderTest, CalForceStressBothFalse)
{
    inp.cal_force = false;
    inp.cal_stress = false;

    const std::string header = Relax_Driver::build_stru_header(0, etot, stress, inp, false, true);
    const auto lines = split_lines(header);

    ASSERT_EQ(lines.size(), 7u);
    EXPECT_EQ(count_substr(header, "# Stress (kbar): N/A N/A N/A"), 3);
    EXPECT_THAT(lines[6], testing::HasSubstr("# NOTE: stress not computed (cal_stress=0); forces not computed (cal_force=0)"));
}

TEST_F(RelaxDriverHeaderTest, FinalStep)
{
    inp.cal_force = true;
    inp.cal_stress = true;

    const std::string header = Relax_Driver::build_stru_header(4, etot, stress, inp, true, true);
    const auto lines = split_lines(header);

    ASSERT_EQ(lines.size(), 7u);
    EXPECT_THAT(lines[2], testing::HasSubstr("# RELAX STEP 5 (FINAL), Energy:"));
}

TEST_F(RelaxDriverHeaderTest, GeometryNotEvaluated)
{
    inp.cal_force = true;
    inp.cal_stress = true;

    const std::string header = Relax_Driver::build_stru_header(4, etot, stress, inp, true, false);
    const auto lines = split_lines(header);

    ASSERT_EQ(lines.size(), 7u);
    EXPECT_EQ(count_substr(header, "# Stress (kbar): N/A N/A N/A"), 3);
    EXPECT_THAT(lines[6], testing::HasSubstr("# NOTE: geometry proposed by optimizer but not evaluated"));
    EXPECT_THAT(lines[6], testing::HasSubstr("stress N/A"));
    EXPECT_THAT(lines[6], testing::HasSubstr("forces omitted"));
    EXPECT_THAT(lines[6], testing::HasSubstr("energy above belongs to the last evaluated geometry"));
}

TEST_F(RelaxDriverHeaderTest, PreviewFinalHeaders)
{
    inp.cal_force = true;
    inp.cal_stress = true;

    const std::string header_ok = Relax_Driver::build_stru_header(4, etot, stress, inp, true, true);
    const std::string header_bad = Relax_Driver::build_stru_header(4, etot, stress, inp, true, false);

    std::cout << "\n========== STRU_FINAL header (converged, geometry evaluated) ==========\n"
              << header_ok
              << "\n========== STRU_FINAL header (early exit, geometry NOT evaluated) ==========\n"
              << header_bad
              << std::endl;

    EXPECT_EQ(split_lines(header_ok).size(), 7u);
    EXPECT_EQ(split_lines(header_bad).size(), 7u);
}
