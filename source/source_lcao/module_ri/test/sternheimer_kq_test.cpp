#include "source_lcao/module_ri/sternheimer_kq.h"

#include <array>
#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

namespace
{

using KPoint = ModuleRI::SternheimerReducedKPoint;

void expect_kpoint_near(const KPoint& actual, const KPoint& expected, const double tolerance = 1.0e-14)
{
    for (int direction = 0; direction != 3; ++direction)
    {
        EXPECT_NEAR(actual[static_cast<std::size_t>(direction)],
                    expected[static_cast<std::size_t>(direction)],
                    tolerance);
    }
}

} // namespace

TEST(SternheimerKQ, FoldsReducedCoordinatesAndRecordsReciprocalShift)
{
    const auto result = ModuleRI::fold_sternheimer_kpoint({0.6, -0.6, 1.5});

    expect_kpoint_near(result.kpoint, {-0.4, 0.4, -0.5});
    EXPECT_EQ(result.reciprocal_shift, (std::array<int, 3>{1, -1, 2}));
}

TEST(SternheimerKQ, MeasuresDistanceAcrossBrillouinZoneBoundary)
{
    const double distance = ModuleRI::periodic_sternheimer_kpoint_distance({0.49, 0.1, -0.2}, {-0.49, 0.1, -0.2});
    EXPECT_NEAR(distance, 0.02, 1.0e-14);
}

TEST(SternheimerKQ, MapsEverySourcePointToUniqueFoldedKPlusQ)
{
    const std::vector<KPoint> kpoints = {
        {-0.25, -0.25, 0.0},
        {0.25, -0.25, 0.0},
        {-0.25, 0.25, 0.0},
        {0.25, 0.25, 0.0},
    };

    const auto mapping = ModuleRI::build_sternheimer_kq_map(kpoints, {0.5, 0.0, 0.0}, 1.0e-12);

    ASSERT_EQ(mapping.size(), kpoints.size());
    EXPECT_EQ(mapping[0].source_index, 0);
    EXPECT_EQ(mapping[0].target_index, 1);
    EXPECT_EQ(mapping[1].source_index, 1);
    EXPECT_EQ(mapping[1].target_index, 0);
    EXPECT_EQ(mapping[2].target_index, 3);
    EXPECT_EQ(mapping[3].target_index, 2);
    expect_kpoint_near(mapping[0].folded_k_plus_q, kpoints[1]);
    expect_kpoint_near(mapping[1].folded_k_plus_q, kpoints[0]);
    EXPECT_EQ(mapping[0].reciprocal_shift, (std::array<int, 3>{0, 0, 0}));
    EXPECT_EQ(mapping[1].reciprocal_shift, (std::array<int, 3>{1, 0, 0}));
}

TEST(SternheimerKQ, RecordsShiftRelativeToTheSelectedABACUSBoundaryRepresentative)
{
    const std::vector<KPoint> kpoints = {
        {0.0, 0.0, 0.0},
        {0.5, 0.0, 0.0},
    };

    const auto mapping = ModuleRI::build_sternheimer_kq_map(kpoints, {0.5, 0.0, 0.0}, 1.0e-12);

    ASSERT_EQ(mapping.size(), 2);
    EXPECT_EQ(mapping[0].target_index, 1);
    expect_kpoint_near(mapping[0].folded_k_plus_q, kpoints[1]);
    EXPECT_EQ(mapping[0].reciprocal_shift, (std::array<int, 3>{0, 0, 0}));
    EXPECT_EQ(mapping[1].target_index, 0);
    expect_kpoint_near(mapping[1].folded_k_plus_q, kpoints[0]);
    EXPECT_EQ(mapping[1].reciprocal_shift, (std::array<int, 3>{1, 0, 0}));
}

TEST(SternheimerKQ, RejectsQThatIsNotCommensurateWithKMesh)
{
    const std::vector<KPoint> kpoints = {
        {-0.25, 0.0, 0.0},
        {0.25, 0.0, 0.0},
    };

    EXPECT_THROW(ModuleRI::build_sternheimer_kq_map(kpoints, {0.125, 0.0, 0.0}, 1.0e-12), std::invalid_argument);
}

TEST(SternheimerKQ, RejectsDuplicatePeriodicTargetPoints)
{
    const std::vector<KPoint> kpoints = {
        {-0.5, 0.0, 0.0},
        {0.5, 0.0, 0.0},
    };

    EXPECT_THROW(ModuleRI::build_sternheimer_kq_map(kpoints, {0.0, 0.0, 0.0}, 1.0e-12), std::invalid_argument);
}

TEST(SternheimerKQ, UsesPositiveABACUSBlochPhaseConvention)
{
    const std::complex<double> phase = ModuleRI::sternheimer_bloch_phase({0.25, 0.125, 0.0}, {1, 1, 0});
    const double inverse_sqrt_two = 1.0 / std::sqrt(2.0);

    EXPECT_NEAR(phase.real(), -inverse_sqrt_two, 1.0e-14);
    EXPECT_NEAR(phase.imag(), inverse_sqrt_two, 1.0e-14);
}

TEST(SternheimerKQ, BlochPhasePreservesGammaAndComposesTranslations)
{
    const KPoint gamma{0.0, 0.0, 0.0};
    EXPECT_EQ(ModuleRI::sternheimer_bloch_phase(gamma, {7, -4, 3}), std::complex<double>(1.0, 0.0));

    const KPoint kpoint{0.125, -0.25, 0.375};
    const std::array<int, 3> first{1, -2, 0};
    const std::array<int, 3> second{-3, 1, 2};
    const std::array<int, 3> sum{
        first[0] + second[0],
        first[1] + second[1],
        first[2] + second[2],
    };
    const std::complex<double> composed
        = ModuleRI::sternheimer_bloch_phase(kpoint, first) * ModuleRI::sternheimer_bloch_phase(kpoint, second);
    const std::complex<double> direct = ModuleRI::sternheimer_bloch_phase(kpoint, sum);

    EXPECT_NEAR(composed.real(), direct.real(), 1.0e-14);
    EXPECT_NEAR(composed.imag(), direct.imag(), 1.0e-14);
}
