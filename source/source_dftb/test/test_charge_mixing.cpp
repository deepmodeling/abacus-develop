#include "source_dftb/charge_mixing.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace ModuleDFTB
{

TEST(DftbNativeChargeMixerTest, ModifiedBroydenConvergesAffineFixedPointAndConservesCharge)
{
    DftbChargeMixerParameters parameters;
    parameters.method = "broyden";
    parameters.mixing_parameter = 0.2;
    parameters.history = 6;
    DftbChargeMixer mixer(parameters);

    const std::vector<double> fixed_point = {0.2, 0.1};
    std::vector<double> charges = {0.8, -0.5};
    constexpr double total_charge = 0.3;
    bool used_broyden = false;
    for (int iteration = 0; iteration < 40; ++iteration)
    {
        std::vector<double> residual(charges.size(), 0.0);
        for (std::size_t atom = 0; atom < charges.size(); ++atom)
        {
            const double output_charge = fixed_point[atom] + 0.7 * (charges[atom] - fixed_point[atom]);
            residual[atom] = output_charge - charges[atom];
        }
        const double max_residual = std::max(std::abs(residual[0]), std::abs(residual[1]));
        if (max_residual < 1.0e-9) break;

        charges = mixer.mix(charges, residual, total_charge);
        used_broyden = used_broyden || mixer.last_step() == "Broyden";
        EXPECT_NEAR(charges[0] + charges[1], total_charge, 1.0e-14);
    }

    EXPECT_TRUE(used_broyden);
    EXPECT_NEAR(charges[0], fixed_point[0], 1.0e-7);
    EXPECT_NEAR(charges[1], fixed_point[1], 1.0e-7);
}

TEST(DftbNativeChargeMixerTest, PulayStillConvergesAndConservesCharge)
{
    DftbChargeMixerParameters parameters;
    parameters.method = "pulay";
    parameters.mixing_parameter = 0.2;
    parameters.history = 6;
    DftbChargeMixer mixer(parameters);

    const std::vector<double> fixed_point = {0.2, 0.1};
    std::vector<double> charges = {0.8, -0.5};
    constexpr double total_charge = 0.3;
    for (int iteration = 0; iteration < 40; ++iteration)
    {
        std::vector<double> residual(charges.size(), 0.0);
        for (std::size_t atom = 0; atom < charges.size(); ++atom)
        {
            const double output_charge = fixed_point[atom] + 0.7 * (charges[atom] - fixed_point[atom]);
            residual[atom] = output_charge - charges[atom];
        }
        if (std::max(std::abs(residual[0]), std::abs(residual[1])) < 1.0e-9) break;
        charges = mixer.mix(charges, residual, total_charge);
        EXPECT_NEAR(charges[0] + charges[1], total_charge, 1.0e-14);
    }
    EXPECT_NEAR(charges[0], fixed_point[0], 1.0e-7);
    EXPECT_NEAR(charges[1], fixed_point[1], 1.0e-7);
}

TEST(DftbNativeChargeMixerTest, FirstBroydenIterationUsesLinearDamping)
{
    DftbChargeMixerParameters parameters;
    parameters.method = "broyden";
    parameters.mixing_parameter = 0.25;
    DftbChargeMixer mixer(parameters);
    const std::vector<double> next = mixer.mix({0.2, -0.2}, {-0.4, 0.4}, 0.0);
    EXPECT_EQ(mixer.last_step(), "linear startup");
    EXPECT_NEAR(next[0], 0.1, 1.0e-14);
    EXPECT_NEAR(next[1], -0.1, 1.0e-14);
}

} // namespace ModuleDFTB
