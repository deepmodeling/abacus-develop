#include "../snap_clebsch_gordan.h"

#include <gtest/gtest.h>

#include <cmath>

namespace
{

TEST(SnapClebschGordan, HighestWeightCoefficientIsOne)
{
    const ModuleMD::SNAP::ClebschGordanTable table(2);
    EXPECT_NEAR(table.coefficient(1, 1, 1, 1, 2, 2), 1.0, 1.0e-14);
}

TEST(SnapClebschGordan, SingletHasExpectedSigns)
{
    const ModuleMD::SNAP::ClebschGordanTable table(2);
    const double value = 1.0 / std::sqrt(2.0);
    EXPECT_NEAR(table.coefficient(1, 1, 1, -1, 0, 0), value, 1.0e-14);
    EXPECT_NEAR(table.coefficient(1, -1, 1, 1, 0, 0), -value, 1.0e-14);
}

TEST(SnapClebschGordan, SelectionRulesReturnZero)
{
    const ModuleMD::SNAP::ClebschGordanTable table(2);
    EXPECT_DOUBLE_EQ(table.coefficient(1, 1, 1, 1, 0, 0), 0.0);
    EXPECT_DOUBLE_EQ(table.coefficient(1, 1, 1, -1, 2, 1), 0.0);
    EXPECT_DOUBLE_EQ(table.h_coefficient(1, 1, 1, 1, -1, -1, 0, 0, 0), 0.5);
}

} // namespace
