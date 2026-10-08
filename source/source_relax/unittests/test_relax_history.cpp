#include "../relax_history.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using testing::HasSubstr;
using testing::Not;

TEST(RelaxHistory, EmptyHistory)
{
    const std::vector<double> hist;
    EXPECT_EQ(format_relax_history(hist), "\n");
}

TEST(RelaxHistory, FewStepsAllShown)
{
    const std::vector<double> hist = {1.0, 2.0, 3.0};
    const std::string out = format_relax_history(hist);
    EXPECT_THAT(out, HasSubstr("1.000e+00"));
    EXPECT_THAT(out, HasSubstr("2.000e+00"));
    EXPECT_THAT(out, HasSubstr("3.000e+00"));
    EXPECT_THAT(out, Not(HasSubstr("omitted")));
}

TEST(RelaxHistory, LineBreakEveryFive)
{
    const std::vector<double> hist(7, 1.0);
    const std::string out = format_relax_history(hist);
    // 7 values with per_line = 5 -> two leading "\n  " line starts
    EXPECT_THAT(out, HasSubstr("\n  1.000e+00 1.000e+00 1.000e+00 1.000e+00 1.000e+00\n  1.000e+00"));
}

TEST(RelaxHistory, LongHistoryTruncated)
{
    std::vector<double> hist(150, 0.0);
    hist.front() = 1.0;   // first value kept
    hist.back() = 2.0;    // last value kept
    hist[75] = 9.0;       // middle value must be omitted
    const std::string out = format_relax_history(hist);
    EXPECT_THAT(out, HasSubstr("... (omitted 130 step(s)) ..."));
    EXPECT_THAT(out, HasSubstr("1.000e+00"));
    EXPECT_THAT(out, HasSubstr("2.000e+00"));
    EXPECT_THAT(out, Not(HasSubstr("9.000e+00")));
}
