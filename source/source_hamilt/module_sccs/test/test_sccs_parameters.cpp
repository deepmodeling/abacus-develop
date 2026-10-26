#include "../sccs_parameters.h"
#include <gtest/gtest.h>

TEST(SccsParameters, ExactPresetNamesOnly)
{
    EXPECT_EQ(ModuleSccs::parse_preset("water-neutral"), ModuleSccs::Preset::WaterNeutral);
    EXPECT_EQ(ModuleSccs::parse_preset("water-cation"), ModuleSccs::Preset::WaterCation);
    EXPECT_EQ(ModuleSccs::parse_preset("water-anion"), ModuleSccs::Preset::WaterAnion);
    EXPECT_EQ(ModuleSccs::parse_preset("vacuum"), ModuleSccs::Preset::Vacuum);
    EXPECT_EQ(ModuleSccs::parse_preset("custom"), ModuleSccs::Preset::Custom);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ModuleSccs::parse_preset("unknown"), ::testing::ExitedWithCode(1), "");
    EXPECT_EXIT(ModuleSccs::parse_preset("Water-Neutral"), ::testing::ExitedWithCode(1), "");
    testing::internal::GetCapturedStdout();
}

// The water presets in Hartree atomic units, which also fixes the dyn/cm and
// GPa conversions, and the vacuum preset.
TEST(SccsParameters, PresetValues)
{
    const ModuleSccs::Preset presets[] = {ModuleSccs::Preset::WaterNeutral,
                                         ModuleSccs::Preset::WaterCation,
                                         ModuleSccs::Preset::WaterAnion};
    const double minima[] = {1e-4, 2e-4, 2.4e-3};
    const double maxima[] = {5e-3, 3.5e-3, 1.55e-2};
    const double tensions[] = {3.076638359192053e-5, 3.211522295607572e-6, 0.0};
    const double pressures[] = {-1.223613888352364e-5, 4.248659334556821e-6, 1.529517360440455e-5};
    for (int index = 0; index < 3; ++index)
    {
        const ModuleSccs::SccsConfig config = ModuleSccs::preset_config(presets[index]);
        EXPECT_DOUBLE_EQ(config.cavity.density_min, minima[index]);
        EXPECT_DOUBLE_EQ(config.cavity.density_max, maxima[index]);
        EXPECT_DOUBLE_EQ(config.cavity.epsilon_bulk, 78.3);
        EXPECT_DOUBLE_EQ(config.surface_regularization, 1e-8);
        EXPECT_NEAR(config.surface_tension, tensions[index], 1e-18);
        EXPECT_NEAR(config.pressure, pressures[index], 1e-18);
    }
    const ModuleSccs::SccsConfig vacuum = ModuleSccs::preset_config(ModuleSccs::Preset::Vacuum);
    EXPECT_DOUBLE_EQ(vacuum.cavity.epsilon_bulk, 1.0);
    EXPECT_DOUBLE_EQ(vacuum.surface_tension, 0.0);
    EXPECT_DOUBLE_EQ(vacuum.pressure, 0.0);
}
