#include <gtest/gtest.h>

#if __has_include("source_lcao/module_ri/sternheimer_weak_q_unit.h")
#include "source_lcao/module_ri/sternheimer_weak_q_unit.h"
#include "source_lcao/module_ri/sternheimer_weak_augmented.h"

#include <limits>

using ModuleRI::SternheimerWeakQUnit;

TEST(SternheimerWeakQUnit, DefaultsCoverSelectedFullKAllBandsAndFrequencies)
{
    SternheimerWeakQUnit::Raw raw;
    raw.source_k = "2";
    const auto unit = SternheimerWeakQUnit::parse(raw, {20, 18, 20}, 6, 1113);
    EXPECT_EQ(unit.source_k, 2);
    EXPECT_EQ(unit.band_begin, 1);
    EXPECT_EQ(unit.band_end, 18);
    EXPECT_EQ(unit.frequency_begin, 1);
    EXPECT_EQ(unit.frequency_end, 6);
    EXPECT_TRUE(unit.all_source_bands(18));
    EXPECT_TRUE(unit.all_frequencies(6));
    EXPECT_EQ(unit.columns().size(), 1113U);
    EXPECT_EQ(unit.expected_equations(), 2ULL * 18 * 6 * 1113);
}

TEST(SternheimerWeakQUnit, OneBandPilotAndUnevenColumnShards)
{
    SternheimerWeakQUnit::Raw raw;
    raw.source_k = "16";
    raw.band_begin = "20";
    raw.band_end = "20";
    raw.frequency_begin = "2";
    raw.frequency_end = "3";
    raw.shard_index = "15";
    raw.shard_count = "16";
    const auto unit = SternheimerWeakQUnit::parse(raw, std::vector<int>(16, 20), 6, 1113);
    EXPECT_FALSE(unit.all_source_bands(20));
    EXPECT_FALSE(unit.all_frequencies(6));
    EXPECT_EQ(unit.columns().front(), 15);
    EXPECT_EQ(unit.columns().back(), 1103);
    EXPECT_EQ(unit.columns().size(), 69U);
    EXPECT_EQ(unit.expected_equations(), 276U);
}

TEST(SternheimerWeakQUnit, RejectsMalformedAndOutOfBoundsSelections)
{
    for (const char* value : {"", "0", "-1", "+1", " 1", "1x", "2147483648", "3"})
    {
        SternheimerWeakQUnit::Raw raw;
        raw.source_k = value;
        EXPECT_THROW(SternheimerWeakQUnit::parse(raw, {2, 3}, 6, 5), std::invalid_argument);
    }
    SternheimerWeakQUnit::Raw raw;
    raw.band_begin = "3";
    EXPECT_THROW(SternheimerWeakQUnit::parse(raw, {2}, 6, 5), std::invalid_argument);
    raw.band_begin = "2";
    raw.band_end = "1";
    EXPECT_THROW(SternheimerWeakQUnit::parse(raw, {2}, 6, 5), std::invalid_argument);
    raw = {};
    raw.frequency_end = "7";
    EXPECT_THROW(SternheimerWeakQUnit::parse(raw, {2}, 6, 5), std::invalid_argument);
    raw.frequency_end = "0";
    EXPECT_THROW(SternheimerWeakQUnit::parse(raw, {2}, 6, 5), std::invalid_argument);
    raw = {};
    raw.shard_index = "2";
    raw.shard_count = "2";
    EXPECT_THROW(SternheimerWeakQUnit::parse(raw, {2}, 6, 5), std::invalid_argument);
    EXPECT_THROW(SternheimerWeakQUnit::parse({}, {}, 6, 5), std::invalid_argument);
    EXPECT_THROW(SternheimerWeakQUnit::parse({}, {0}, 6, 5), std::invalid_argument);
    EXPECT_THROW(SternheimerWeakQUnit::parse({}, {2}, 0, 5), std::invalid_argument);
}

TEST(SternheimerWeakQUnit, CanonicalTargetAdmitsIntegerFoldWithoutChangingFields)
{
    using Point = std::array<double, 3>;
    using Shift = std::array<int, 3>;
    EXPECT_NO_THROW(ModuleRI::validate_sternheimer_weak_q_fold(
        Point{0.25, -0.5, 0.0}, Point{0.25, -0.25, 0.0}, Point{-0.5, 0.25, 0.0}, Shift{1, -1, 0}));
    EXPECT_NO_THROW(ModuleRI::validate_sternheimer_weak_q_fold(
        Point{0, 0, 0}, Point{0.25, 0, 0}, Point{0.25, 0, 0}, Shift{0, 0, 0}));
    EXPECT_THROW(ModuleRI::validate_sternheimer_weak_q_fold(
        Point{0.25, 0, 0}, Point{0.25, 0, 0}, Point{-0.5, 0, 0}, Shift{0, 0, 0}), std::invalid_argument);
    EXPECT_THROW(ModuleRI::validate_sternheimer_weak_q_fold(
        Point{0.25, 0, 0}, Point{0.25, 0, 0}, Point{-0.49, 0, 0}, Shift{1, 0, 0}), std::invalid_argument);
    EXPECT_THROW(ModuleRI::validate_sternheimer_weak_q_fold(
        Point{std::numeric_limits<double>::quiet_NaN(), 0, 0}, Point{}, Point{}, Shift{}), std::invalid_argument);
}

TEST(SternheimerWeakQUnit, HaRyAndSpinWeightedColumnMajorSignedSum)
{
    using Complex = std::complex<double>;
    EXPECT_DOUBLE_EQ(ModuleRI::sternheimer_weak_q_omega_ry(0.25), 0.5);
    EXPECT_DOUBLE_EQ(ModuleRI::sternheimer_weak_q_left_vertex_factor(), 0.5);
    const double weight = ModuleRI::sternheimer_weak_q_weight(2.0 / 16, 1.0);
    EXPECT_DOUBLE_EQ(weight, 0.125);
    std::vector<Complex> result(6, Complex(0.0));
    const std::vector<Complex> plus{{1, 2}, {3, -4}, {5, 6}};
    const std::vector<Complex> minus{{7, 8}, {9, 10}, {11, -12}};
    ModuleRI::accumulate_sternheimer_weak_q_column(result, plus, 3, 1, weight);
    ModuleRI::accumulate_sternheimer_weak_q_column(result, minus, 3, 1, weight);
    for (int row = 0; row < 3; ++row)
    {
        EXPECT_EQ(result[row], Complex(0));
        EXPECT_EQ(result[row + 3], weight * (plus[row] + minus[row]));
    }
    EXPECT_THROW(ModuleRI::sternheimer_weak_q_weight(-0.1, 1), std::invalid_argument);
    EXPECT_THROW(ModuleRI::sternheimer_weak_q_weight(0.1, 0.5), std::invalid_argument);
    EXPECT_THROW(ModuleRI::sternheimer_weak_q_omega_ry(0), std::invalid_argument);
    EXPECT_THROW(ModuleRI::accumulate_sternheimer_weak_q_column(result, plus, 3, 2, weight), std::invalid_argument);
    auto bad = plus;
    bad[1] = {std::numeric_limits<double>::infinity(), 0};
    EXPECT_THROW(ModuleRI::accumulate_sternheimer_weak_q_column(result, bad, 3, 1, weight), std::invalid_argument);
}

TEST(SternheimerWeakQUnit, OriginalResidualMustBeFiniteAndWithinTolerance)
{
    EXPECT_NO_THROW(ModuleRI::validate_sternheimer_weak_q_residual(true, 1e-8, 2e-8));
    EXPECT_THROW(ModuleRI::validate_sternheimer_weak_q_residual(false, 0, 0), std::runtime_error);
    EXPECT_THROW(ModuleRI::validate_sternheimer_weak_q_residual(true, 1.01e-8, 1e-9), std::runtime_error);
    EXPECT_THROW(ModuleRI::validate_sternheimer_weak_q_residual(
        true, std::numeric_limits<double>::quiet_NaN(), 0), std::runtime_error);
    EXPECT_THROW(ModuleRI::validate_sternheimer_weak_q_residual(true, 0, -1), std::runtime_error);
}

TEST(SternheimerWeakQUnit, ActualCompactDensitySupportIncludesLastInterpolationInterval)
{
    EXPECT_DOUBLE_EQ(ModuleRI::sternheimer_weak_q_radial_support({0, 1, 2, 3, 4}, {1, 0.2, 0, 0, 0}), 2.0);
    EXPECT_DOUBLE_EQ(ModuleRI::sternheimer_weak_q_radial_support({0, 1, 2}, {1, 0.2, 0.1}), 2.0);
    EXPECT_TRUE(ModuleRI::sternheimer_weak_q_z_support_contained(10, 2, 20));
    EXPECT_FALSE(ModuleRI::sternheimer_weak_q_z_support_contained(1, 2, 20));
    EXPECT_FALSE(ModuleRI::sternheimer_weak_q_z_support_contained(19, 2, 20));
    EXPECT_FALSE(ModuleRI::sternheimer_weak_q_z_support_contained(-1, 0.5, 20));
    EXPECT_THROW(ModuleRI::sternheimer_weak_q_radial_support({0, 2, 1}, {1, 0, 0}), std::invalid_argument);
}

TEST(SternheimerWeakQUnit, BothSignedSolvesMatchDirectResolventWithUnrotatedSourceEnergies)
{
    using Blocks = ModuleRI::SternheimerWeakAugmented;
    using Complex = Blocks::Complex;
    Blocks::Data data;
    data.nocc = 1;
    data.nvirtual = 2;
    data.ncoarse = 0;
    data.hu = {-2.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 3.0};
    auto blocks = std::make_shared<const Blocks>(std::move(data));
    const std::array<double, 2> source_energies{-0.7, -1.2}; // Deliberately not energy-sorted.
    const std::array<double, 2> target_energies{1.5, 3.0};
    const double omega_ha = 0.3;
    const double omega_ry = ModuleRI::sternheimer_weak_q_omega_ry(omega_ha);
    const double weight = ModuleRI::sternheimer_weak_q_weight(2.0 / 16, 1.0);
    Blocks::Vector actual(4, 0.0), expected(4, 0.0);
    for (int band = 0; band != 2; ++band)
    {
        const double amplitude = band + 1.0;
        const std::array<Blocks::Vector, 2> vertices{{
            {amplitude * Complex(0.3, 0.7), Complex(-0.8, 0.2)},
            {Complex(0.6, -0.5), amplitude * Complex(0.4, 0.9)}}};
        for (int sign : {1, -1})
        {
            Blocks::Worker worker(blocks, [](const auto&, auto& y) { y.clear(); },
                                  source_energies[band], sign * omega_ry);
            for (int j = 0; j != 2; ++j)
            {
                ModuleRI::SternheimerRPA::SolverOptions options;
                options.residual_tol = 1e-8;
                const auto solved = worker.solve({vertices[j], {}}, options);
                ASSERT_TRUE(solved.converged);
                ModuleRI::validate_sternheimer_weak_q_residual(
                    solved.converged, solved.relative_residual, solved.absolute_residual);
                Blocks::Vector column(2, 0.0);
                for (int i = 0; i != 2; ++i)
                    for (int a = 0; a != 2; ++a)
                        column[i] += ModuleRI::sternheimer_weak_q_left_vertex_factor()
                            * std::conj(vertices[i][a]) * solved.coefficients[a];
                ModuleRI::accumulate_sternheimer_weak_q_column(actual, column, 2, j, weight);
            }
        }
        for (int j = 0; j != 2; ++j)
            for (int i = 0; i != 2; ++i)
                for (int a = 0; a != 2; ++a)
                {
                    const double gap = target_energies[a] - source_energies[band];
                    expected[i + 2 * j] -= weight * gap / (gap * gap + omega_ry * omega_ry)
                        * std::conj(vertices[i][a]) * vertices[j][a];
                }
    }
    for (int i = 0; i != 4; ++i) EXPECT_NEAR(std::abs(actual[i] - expected[i]), 0.0, 1e-12);
    EXPECT_GT(std::abs(expected[1].imag()), 1e-3);
    EXPECT_LT(actual[0].real(), 0.0);
    // A second B+B^dagger would fail this direct normalization test by a factor of two.
    EXPECT_GT(std::abs(2.0 * actual[0] - expected[0]), 1e-3);
}
#else
TEST(SternheimerWeakQUnit, MissingImplementationIsRed)
{
    FAIL() << "Missing sternheimer_weak_q_unit.h: test-first remote RED snapshot";
}
#endif
