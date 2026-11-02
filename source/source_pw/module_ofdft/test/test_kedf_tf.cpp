#include "../kedf_tf.h"

#include <gtest/gtest.h>

class TestKEDFTF : public ::testing::Test
{
protected:
    KEDF_TF kedf;
};

TEST_F(TestKEDFTF, IntegratesNonuniformDensity)
{
    // For rho = {0, 1, 8}, integral(rho^(5/3)) = (0 + 1 + 32) * dV.
    const double density[] = {0.0, 1.0, 8.0};
    const double* rho[] = {density};
    kedf.set_para(3, 0.25, 1.0);

    const double energy = kedf.get_energy(rho);
    const double coefficient = 5.742468000376382;
    const double expected = coefficient * 33.0 * 0.25;
    EXPECT_NEAR(energy, expected, 1e-12);
    EXPECT_DOUBLE_EQ(kedf.tf_energy, energy);
}

TEST_F(TestKEDFTF, ZeroDensityHasZeroEnergy)
{
    const double density[] = {0.0, 0.0};
    const double* rho[] = {density};
    kedf.set_para(2, 0.5, 1.0);

    EXPECT_DOUBLE_EQ(kedf.get_energy(rho), 0.0);
    EXPECT_DOUBLE_EQ(kedf.get_energy_density(rho, 0, 1), 0.0);
}

TEST_F(TestKEDFTF, EnergyDensityUsesWeightAndSelectedSpin)
{
    const double spin_up[] = {1.0, 8.0};
    const double spin_down[] = {8.0, 1.0};
    const double* rho[] = {spin_up, spin_down};
    const double coefficient = 5.742468000376382;
    kedf.set_para(2, 0.25, 0.5);

    EXPECT_NEAR(kedf.get_energy_density(rho, 0, 1), coefficient * 16.0, 1e-12);
    EXPECT_NEAR(kedf.get_energy_density(rho, 1, 1), coefficient * 0.5, 1e-12);

    kedf.set_para(2, 0.25, 0.0);
    EXPECT_DOUBLE_EQ(kedf.get_energy_density(rho, 0, 1), 0.0);
}

TEST_F(TestKEDFTF, StressIsIsotropicAndScalesWithVolume)
{
    const double density[] = {1.0};
    const double* rho[] = {density};
    kedf.set_para(1, 0.5, 1.0);
    const double energy = kedf.get_energy(rho);
    kedf.get_stress(2.0);
    const double expected = energy / 3.0;

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            if (i == j)
            {
                EXPECT_NEAR(kedf.stress(i, j), expected, 1e-12);
            }
            else
            {
                EXPECT_DOUBLE_EQ(kedf.stress(i, j), 0.0);
            }
        }
    }

    kedf.get_stress(4.0);
    for (int i = 0; i < 3; ++i)
    {
        EXPECT_NEAR(kedf.stress(i, i), expected * 0.5, 1e-12);
    }
}

TEST_F(TestKEDFTF, WeightedEnergyMatchesPointwiseIntegralAndStress)
{
    const double density[] = {1.0, 8.0};
    const double* rho[] = {density};
    const double coefficient = 5.742468000376382;
    const double dV = 0.25;
    const double volume = 2.0 * dV;
    const double weights[] = {0.0, 0.4, 1.0, 1.1};
    for (const double weight : weights)
    {
        SCOPED_TRACE(weight);
        kedf.set_para(2, dV, weight);
        const double expected = weight * coefficient * 33.0 * dV;
        const double energy = kedf.get_energy(rho);
        EXPECT_NEAR(energy, expected, 1e-12);

        double integral = 0.0;
        for (int ir = 0; ir < 2; ++ir)
        {
            integral += kedf.get_energy_density(rho, 0, ir) * dV;
        }
        EXPECT_NEAR(energy, integral, 1e-12);

        ModuleBase::matrix potential(1, 2);
        kedf.tf_potential(rho, potential);
        EXPECT_NEAR(kedf.tf_energy, expected, 1e-12);
        kedf.get_stress(volume);
        const double expected_stress = 2.0 * expected / (3.0 * volume);
        for (int i = 0; i < 3; ++i)
        {
            EXPECT_NEAR(kedf.stress(i, i), expected_stress, 1e-12);
        }
    }
}

TEST_F(TestKEDFTF, WeightedPotentialMatchesConservedDensityDerivative)
{
    const double density[] = {1.0, 8.0};
    const double direction[] = {1.0, -1.0};
    const double* rho[] = {density};
    const double dV = 0.25;
    const double weights[] = {0.0, 0.4, 1.0, 1.1};
    const double steps[] = {1e-3, 1e-4};
    for (const double weight : weights)
    {
        SCOPED_TRACE(weight);
        kedf.set_para(2, dV, weight);
        ModuleBase::matrix potential(1, 2);
        kedf.tf_potential(rho, potential);
        const double product = dV * (potential(0, 0) - potential(0, 1));
        for (const double step : steps)
        {
            SCOPED_TRACE(step);
            double plus_density[2];
            double minus_density[2];
            for (int ir = 0; ir < 2; ++ir)
            {
                plus_density[ir] = density[ir] + step * direction[ir];
                minus_density[ir] = density[ir] - step * direction[ir];
            }
            const double* plus_rho[] = {plus_density};
            const double* minus_rho[] = {minus_density};
            const double plus_energy = kedf.get_energy(plus_rho);
            const double minus_energy = kedf.get_energy(minus_rho);
            const double derivative = (plus_energy - minus_energy) / (2.0 * step);
            EXPECT_NEAR(derivative, product, 1e-6);
        }
    }
}
