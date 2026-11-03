#include "../utils/spectrum_weights.h"
#include <gtest/gtest.h>
#include <cmath>
#include <complex>

TEST(SpectrumWeights, RemovesOnlyUnpolarizedSpinDegeneracy)
{
    EXPECT_DOUBLE_EQ(LR_Util::spatial_k_weight(0.002, 1), 0.001);
    EXPECT_DOUBLE_EQ(LR_Util::spatial_k_weight(0.001, 2), 0.001);
    EXPECT_DOUBLE_EQ(LR_Util::spatial_k_weight(0.001, 4), 0.001);
}

TEST(SpectrumWeights, PreservesNonuniformSpatialWeights)
{
    const double ks_weights[] = {0.2, 0.6, 1.2};
    const double expected[] = {0.1, 0.3, 0.6};
    double sum = 0.0;
    for (int k = 0; k < 3; ++k)
    {
        const double w = LR_Util::spatial_k_weight(ks_weights[k], 1);
        EXPECT_DOUBLE_EQ(w, expected[k]);
        sum += w;
    }
    EXPECT_DOUBLE_EQ(sum, 1.0);
}

TEST(SpectrumWeights, SingletIntensityMatchesExplicitSpinChannels)
{
    // A unit-normalized singlet has X_up=X_down=X/sqrt(2).
    // Compare that explicit contraction with the spatial-orbital sqrt(2) path.
    const int nk = 2;
    const std::complex<double> dipoles[] = {{1.0, 0.5}, {-0.2, 0.7}};
    const double amplitudes[] = {0.6, 0.8};
    std::complex<double> singlet = 0.0;
    std::complex<double> explicit_spins = 0.0;
    for (int k = 0; k < nk; ++k)
    {
        const double unpolarized_w = LR_Util::spatial_k_weight(1.0, 1);
        const double polarized_w = LR_Util::spatial_k_weight(0.5, 2);
        singlet += std::sqrt(2.0) * dipoles[k] * amplitudes[k] * unpolarized_w * double(nk);
        for (int spin = 0; spin < 2; ++spin)
        {
            explicit_spins += dipoles[k] * amplitudes[k] / std::sqrt(2.0) * polarized_w * double(nk);
        }
    }
    EXPECT_NEAR(std::abs(singlet - explicit_spins), 0.0, 1e-14);
    EXPECT_NEAR(std::norm(singlet), 1.8664, 1e-14);
}

TEST(SpectrumWeights, LRConversionPreservesSharedKSWeights)
{
    const std::vector<double> ks_weights = {0.2, 0.6, 1.2};
    const std::vector<double> spatial = LR_Util::spatial_k_weights(ks_weights, 1);
    EXPECT_EQ(spatial, (std::vector<double>{0.1, 0.3, 0.6}));
    EXPECT_EQ(ks_weights, (std::vector<double>{0.2, 0.6, 1.2}));
    EXPECT_EQ(LR_Util::spatial_k_weights(spatial, 2), spatial);
}

TEST(SpectrumWeights, LRAndBSEIndependentTransitionHaveSameSingletStrength)
{
    // BSE RI_kRlist supplies 1/Nk directly; KS supplies 2/Nk for nspin=1.
    // An IPA eigenvector has one unit entry, with no compensating spin scale.
    for (const int nk : {1, 2, 8})
    {
        const std::vector<double> ks_weights(nk, 2.0 / nk);
        const std::vector<double> lr_weights = LR_Util::spatial_k_weights(ks_weights, 1);
        const std::vector<double> bse_weights(nk, 1.0 / nk);
        const std::complex<double> dipole(1.0, 0.5);
        const std::complex<double> lr_amplitude = std::sqrt(2.0) * LR_Util::weighted_velocity_transition(
            dipole, 1.0, 0.0, 1.0, true, lr_weights, 0, 0, nk);
        const std::complex<double> bse_amplitude = std::sqrt(2.0) * LR_Util::weighted_velocity_transition(
            dipole, 1.0, 0.0, 1.0, true, bse_weights, 0, 0, nk);
        EXPECT_NEAR(std::norm(lr_amplitude), 2.5, 1e-14);
        EXPECT_NEAR(std::abs(lr_amplitude - bse_amplitude), 0.0, 1e-14);
    }
}

TEST(SpectrumWeights, ProductionContractionSelectsEachSpinBlock)
{
    const std::vector<double> ks_weights = {0.1, 0.9, 0.3, 0.7};
    const std::vector<double> weights = LR_Util::spatial_k_weights(ks_weights, 2);
    const std::complex<double> velocity(2.0, 4.0);
    const double expected[] = {0.2, 1.8, 0.6, 1.4};
    for (int spin = 0; spin < 2; ++spin)
    {
        for (int ik = 0; ik < 2; ++ik)
        {
            const std::complex<double> actual = LR_Util::weighted_velocity_transition(
                velocity, 1.0, 0.0, 2.0, true, weights, spin, ik, 2);
            EXPECT_NEAR(std::abs(actual - std::complex<double>(1.0, 2.0) * expected[spin * 2 + ik]), 0.0, 1e-14);
        }
    }
}

TEST(SpectrumWeights, ProductionFullBSEContractionKeepsDeexcitationSign)
{
    const std::vector<double> bse_weights = {1.0};
    const std::complex<double> velocity(2.0, 4.0);
    const std::complex<double> ks_gap = LR_Util::weighted_velocity_transition(
        velocity, 3.0, 1.0, 2.0, true, bse_weights, 0, 0, 1);
    const std::complex<double> lr_energy = LR_Util::weighted_velocity_transition(
        velocity, 3.0, 1.0, 2.0, false, bse_weights, 0, 0, 1);
    EXPECT_EQ(ks_gap, std::complex<double>(2.0, 8.0));
    EXPECT_EQ(lr_energy, std::complex<double>(8.0, 8.0));
}
