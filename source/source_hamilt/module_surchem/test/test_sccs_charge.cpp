#include "../sccs/sccs_charge.h"

#include <gtest/gtest.h>

#include <vector>

namespace
{

TEST(SccsCharge, SumsOnlyChargeBearingSpinChannels)
{
    const std::vector<std::vector<double>> spin_density = {{0.2, 0.3},
                                                           {0.1, 0.4},
                                                           {8.0, 9.0}};
    const std::vector<double> density = ModuleSccs::sum_electron_density(spin_density, 2);
    ASSERT_EQ(density.size(), 2);
    EXPECT_DOUBLE_EQ(density[0], 0.3);
    EXPECT_DOUBLE_EQ(density[1], 0.7);
}

TEST(SccsCharge, PreservesPositiveAndNegativeNetCharge)
{
    const ModuleSurchem::SerialChargeReduction reduction;
    const std::vector<double> ionic_density(4, 0.75);

    const ModuleSccs::ChargeDensity cation
        = ModuleSccs::assemble_charge_density(std::vector<double>(4, 0.5),
                                              ionic_density,
                                              1.0,
                                              2.0,
                                              3.0,
                                              1.0e-12,
                                              reduction);
    EXPECT_DOUBLE_EQ(cation.net_charge, 1.0);
    for (std::size_t index = 0; index < cation.solute.size(); ++index)
    {
        EXPECT_DOUBLE_EQ(cation.solute[index], 0.25);
    }

    const ModuleSccs::ChargeDensity anion
        = ModuleSccs::assemble_charge_density(std::vector<double>(4, 1.0),
                                              ionic_density,
                                              1.0,
                                              4.0,
                                              3.0,
                                              1.0e-12,
                                              reduction);
    EXPECT_DOUBLE_EQ(anion.net_charge, -1.0);
}

TEST(SccsCharge, RejectsNormalizationMismatch)
{
    const ModuleSurchem::SerialChargeReduction reduction;
    EXPECT_THROW(ModuleSccs::assemble_charge_density(std::vector<double>(2, 1.0),
                                                     std::vector<double>(2, 1.0),
                                                     1.0,
                                                     3.0,
                                                     2.0,
                                                     1.0e-12,
                                                     reduction),
                 std::runtime_error);
}

} // namespace
