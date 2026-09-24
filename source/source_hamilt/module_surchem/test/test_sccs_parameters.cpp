#include "../sccs_parameters.h"

#include "gtest/gtest.h"

#include <stdexcept>

namespace
{

TEST(SccsParameters, ParsesNamesCaseInsensitively)
{
    EXPECT_EQ(ModuleSccs::parse_preset("WATER-CATION"), ModuleSccs::Preset::WaterCation);
    EXPECT_EQ(ModuleSccs::parse_boundary("PCC_0D"), ModuleSccs::Boundary::Pcc0d);
    EXPECT_EQ(ModuleSccs::parse_boundary("PCC-2D"), ModuleSccs::Boundary::Pcc2d);
    EXPECT_THROW(ModuleSccs::parse_preset("automatic"), std::invalid_argument);
}

TEST(SccsParameters, ReproducesPinnedEnvironWaterPresets)
{
    const ModuleSccs::SccsConfig neutral = ModuleSccs::water_preset(ModuleSccs::Preset::WaterNeutral);
    const ModuleSccs::SccsConfig cation = ModuleSccs::water_preset(ModuleSccs::Preset::WaterCation);
    const ModuleSccs::SccsConfig anion = ModuleSccs::water_preset(ModuleSccs::Preset::WaterAnion);

    EXPECT_DOUBLE_EQ(neutral.cavity.density_min, 1.0e-4);
    EXPECT_DOUBLE_EQ(neutral.cavity.density_max, 5.0e-3);
    EXPECT_DOUBLE_EQ(cation.cavity.density_min, 2.0e-4);
    EXPECT_DOUBLE_EQ(cation.cavity.density_max, 3.5e-3);
    EXPECT_DOUBLE_EQ(anion.cavity.density_min, 2.4e-3);
    EXPECT_DOUBLE_EQ(anion.cavity.density_max, 1.55e-2);
    EXPECT_DOUBLE_EQ(neutral.cavity.epsilon_bulk, 78.3);
    EXPECT_GT(neutral.surface_tension, cation.surface_tension);
    EXPECT_DOUBLE_EQ(anion.surface_tension, 0.0);
    EXPECT_LT(neutral.pressure, 0.0);
    EXPECT_GT(cation.pressure, 0.0);
    EXPECT_GT(anion.pressure, cation.pressure);
}

TEST(SccsParameters, ConvertsPublishedInputUnitsToAtomicUnits)
{
    EXPECT_NEAR(ModuleSccs::dyn_per_cm_to_hartree_per_bohr2(1.0), 6.423048558616410e-7, 1.0e-18);
    EXPECT_NEAR(ModuleSccs::gpa_to_hartree_per_bohr3(1.0), 3.398930921743166e-5, 1.0e-16);
}

TEST(SccsParameters, RejectsInvalidSolverControls)
{
    ModuleSccs::SccsConfig config = ModuleSccs::water_preset(ModuleSccs::Preset::WaterNeutral);
    config.mixing = 0.0;
    EXPECT_THROW(ModuleSccs::validate_config(config), std::invalid_argument);
    config.mixing = 0.5;
    config.max_iterations = 0;
    EXPECT_THROW(ModuleSccs::validate_config(config), std::invalid_argument);
    config.max_iterations = 200;
    config.adaptive_mixing = true;
    config.mixing_min = 0.6;
    config.mixing_max = 0.8;
    EXPECT_THROW(ModuleSccs::validate_config(config), std::invalid_argument);
    config.mixing = 0.7;
    EXPECT_NO_THROW(ModuleSccs::validate_config(config));
}

} // namespace
