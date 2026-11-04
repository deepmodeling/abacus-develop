#include "../smooth_cut_coulomb.h"
#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>

TEST(SmoothCutCoulomb, GammaMomentAndSharpLimit)
{
    const double cut_width = 1.091868692359;
    const double cut_sigma = std::log(cut_width);
    const double cut_z = 1.e-7;
    EXPECT_NEAR(Conv_Coulomb_Pot_K::log_erfc_cut_factor(cut_z, 1.0, cut_width) / (cut_z * cut_z),
                0.5 * std::exp(cut_sigma * cut_sigma), 1.e-13);
    EXPECT_DOUBLE_EQ(Conv_Coulomb_Pot_K::log_erfc_cut_factor(0.0, 32.0, cut_width), 0.0);
    EXPECT_NEAR(Conv_Coulomb_Pot_K::log_erfc_cut_factor(0.2, 32.0, 1.00000001),
                1.0 - std::cos(6.4), 1.e-13);
    for (int cut_index = 1; cut_index <= 100; ++cut_index)
    {
        const double cut_factor = Conv_Coulomb_Pot_K::log_erfc_cut_factor(cut_index, 20.0, cut_width);
        EXPECT_GE(cut_factor, 0.0);
        EXPECT_LE(cut_factor, 2.0 + 1.e-14);
    }
}

TEST(SmoothCutCoulomb, InscribedRadiusUsesShortestCombination)
{
    const std::array<std::array<double, 3>, 3> cut_cubic = {{{{8., 0., 0.}}, {{0., 10., 0.}}, {{0., 0., 12.}}}};
    const std::array<std::array<double, 3>, 3> cut_fcc = {{{{0., 4., 4.}}, {{4., 0., 4.}}, {{4., 4., 0.}}}};
    const std::array<std::array<double, 3>, 3> cut_skew = {{{{10., 0., 0.}}, {{9.9, 0.1, 0.}}, {{0., 0., 5.}}}};
    EXPECT_DOUBLE_EQ(Conv_Coulomb_Pot_K::bvk_inscribed_radius(cut_cubic), 4.0);
    EXPECT_NEAR(Conv_Coulomb_Pot_K::bvk_inscribed_radius(cut_fcc), std::sqrt(8.0), 1.e-14);
    EXPECT_NEAR(Conv_Coulomb_Pot_K::bvk_inscribed_radius(cut_skew), std::sqrt(0.02) / 2., 1.e-14);
    EXPECT_THROW(Conv_Coulomb_Pot_K::log_erfc_cut_factor(1., 1., 1.), std::invalid_argument);
    EXPECT_THROW(Conv_Coulomb_Pot_K::log_erfc_cut_factor(-1., 1., 1.1), std::invalid_argument);
}
