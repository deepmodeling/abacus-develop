#include <gtest/gtest.h>
#include "../state_group.h"
#include <cmath>

TEST(RootGroup, RotatedDegenerateSpaceBeatsUnrelatedRoot)
{
    const double c = std::sqrt(0.5);
    const std::vector<std::complex<double>> overlap{c, c, 0.0, -c, c, 0.0};
    const auto result = LR::match_state_group(overlap, 2, 3, {{0, 1}, {2}}, {c, c, 0.0}, false);
    EXPECT_EQ(result.group, 0);
    EXPECT_NEAR(result.old_coverage, 1.0, 1e-14);
    ASSERT_EQ(result.singular_values.size(), 2);
    EXPECT_NEAR(result.singular_values[0], 1.0, 1e-14);
    EXPECT_NEAR(result.singular_values[1], 1.0, 1e-14);
}

TEST(RootGroup, SplitUsesReferenceToBreakSubspaceTie)
{
    const auto result = LR::match_state_group({1.0, 0.0, 0.0, 1.0}, 2, 2,
        {{0}, {1}}, {0.1, 0.9}, false);
    EXPECT_EQ(result.group, 1);
    EXPECT_DOUBLE_EQ(result.old_coverage, 0.5);
    EXPECT_DOUBLE_EQ(result.new_coverage, 1.0);
    EXPECT_DOUBLE_EQ(result.singular_values[0], 1.0);
}

TEST(RootGroup, MergeReportsUnrelatedAddedDirection)
{
    const auto result = LR::match_state_group({1.0, 0.0}, 1, 2, {{0, 1}}, {1.0, 0.0}, false);
    EXPECT_DOUBLE_EQ(result.old_coverage, 1.0);
    EXPECT_DOUBLE_EQ(result.new_coverage, 0.5);
}

TEST(RootGroup, HandlesMoreThanThreeRootsAndComplexGauge)
{
    std::vector<std::complex<double>> overlap(16, 0.0);
    for (int i = 0; i < 4; ++i) { overlap[i * 4 + i] = std::complex<double>(0.0, 1.0); }
    const auto result = LR::match_state_group(overlap, 4, 4, {{0, 1, 2, 3}}, {1.0, 0.0, 0.0, 0.0}, false);
    ASSERT_EQ(result.singular_values.size(), 4);
    for (const double value : result.singular_values) { EXPECT_DOUBLE_EQ(value, 1.0); }
}

TEST(RootGroup, PenalizesUnrelatedExtraMembers)
{
    const auto result = LR::match_state_group({0.8, 0.0, 0.7}, 1, 3,
        {{0, 1}, {2}}, {0.8, 0.0, 0.7}, false);
    EXPECT_EQ(result.group, 1);
}

TEST(RootGroup, RejectsInvalidDimensions)
{
    EXPECT_THROW(LR::match_state_group({}, 0, 0, {}, {}, false), std::invalid_argument);
}

TEST(RootGroup, ComplexRotationPreservesBothDirections)
{
    const double c = std::sqrt(0.5);
    const std::complex<double> phase(0.0, c);
    const auto result = LR::match_state_group({c, phase, phase, c}, 2, 2, {{0, 1}}, {c, c}, false);
    EXPECT_NEAR(result.singular_values[0], 1.0, 1e-14);
    EXPECT_NEAR(result.singular_values[1], 1.0, 1e-14);
}

TEST(RootGroup, ProjectionLossIsNotNormalizedAway)
{
    const auto result = LR::match_state_group({0.4, 0.0, 0.0, 0.4}, 2, 2, {{0, 1}}, {0.4, 0.0}, false);
    EXPECT_NEAR(result.old_coverage, 0.16, 1e-14);
    EXPECT_NEAR(result.new_coverage, 0.16, 1e-14);
    EXPECT_NEAR(result.singular_values[0], 0.4, 1e-14);
    EXPECT_NEAR(result.singular_values[1], 0.4, 1e-14);
}

TEST(RootGroup, JTTracksSingleBranchAfterThreeToOnePlusTwoSplit)
{
    const std::vector<std::complex<double>> overlap{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
    const std::vector<std::vector<int>> groups{{0}, {1, 2}};
    const std::vector<double> reference{1.0, 0.0, 0.0};
    const auto average = LR::match_state_group(overlap, 3, 3, groups, reference, false);
    const auto jt = LR::match_state_group(overlap, 3, 3, groups, reference, true);
    EXPECT_EQ(average.group, 1);
    EXPECT_EQ(jt.group, 0);
    EXPECT_DOUBLE_EQ(jt.reference_weights[0], 1.0);
    EXPECT_DOUBLE_EQ(jt.reference_weights[1], 0.0);
    EXPECT_NEAR(jt.old_coverage, 1.0 / 3.0, 1e-14);
    EXPECT_DOUBLE_EQ(jt.new_coverage, 1.0);
}

TEST(RootGroup, JTGroupWeightIsInvariantUnderRotationWithinNewGroup)
{
    const double c = std::sqrt(0.5);
    const std::vector<std::complex<double>> overlap{c, c, 0.0, -c, c, 0.0, 0.0, 0.0, 1.0};
    const double outside = std::sqrt(0.395);
    const auto result = LR::match_state_group(overlap, 3, 3, {{0, 1}, {2}}, {0.55, 0.55, outside}, true);
    const std::vector<std::complex<double>> identity{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
    const double inside = std::sqrt(0.605);
    const auto unrotated = LR::match_state_group(identity, 3, 3, {{0, 1}, {2}}, {inside, 0.0, outside}, true);
    EXPECT_EQ(result.group, 0);
    EXPECT_EQ(result.group, unrotated.group);
    EXPECT_NEAR(result.reference_weights[0], 0.605, 1e-14);
    EXPECT_NEAR(result.reference_weights[1], 0.395, 1e-14);
    EXPECT_NEAR(result.reference_weights[0], unrotated.reference_weights[0], 1e-14);
}
