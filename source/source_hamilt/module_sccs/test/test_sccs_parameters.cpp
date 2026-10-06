#include "../sccs_parameters.h"
#include <gtest/gtest.h>

TEST(SccsParameters, PresetNames)
{
    EXPECT_EQ(ModuleSccs::parse_preset("water-neutral"), ModuleSccs::Preset::WaterNeutral);
    EXPECT_EQ(ModuleSccs::parse_preset("water-cation"), ModuleSccs::Preset::WaterCation);
    EXPECT_EQ(ModuleSccs::parse_preset("water-anion"), ModuleSccs::Preset::WaterAnion);
    EXPECT_EQ(ModuleSccs::parse_preset("vacuum"), ModuleSccs::Preset::Vacuum);
    EXPECT_EQ(ModuleSccs::parse_preset("custom"), ModuleSccs::Preset::Custom);
}

TEST(SccsParameters, UnknownPresetStops)
{
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ModuleSccs::parse_preset("unknown"), ::testing::ExitedWithCode(1), "");
    EXPECT_EXIT(ModuleSccs::parse_preset("Water-Neutral"), ::testing::ExitedWithCode(1), "");
    testing::internal::GetCapturedStdout();
}

TEST(SccsParameters, OriginalWaterPresets)
{
    const ModuleSccs::Preset presets[] = {ModuleSccs::Preset::WaterNeutral,
                                         ModuleSccs::Preset::WaterCation,
                                         ModuleSccs::Preset::WaterAnion};
    const double minima[] = {1e-4, 2e-4, 2.4e-3};
    const double maxima[] = {5e-3, 3.5e-3, 1.55e-2};
    const double tensions[] = {3.076640259577260e-5, 3.211524279308205e-6, 0.0};
    const double pressures[] = {-1.223615131827540e-5, 4.248663652178958e-6, 1.529518914784425e-5};
    for (int index = 0; index < 3; ++index)
    {
        const ModuleSccs::SccsConfig config = ModuleSccs::make_sccs_config(presets[index]);
        EXPECT_DOUBLE_EQ(config.cavity.density_min, minima[index]);
        EXPECT_DOUBLE_EQ(config.cavity.density_max, maxima[index]);
        EXPECT_DOUBLE_EQ(config.cavity.epsilon_bulk, 78.3);
        EXPECT_DOUBLE_EQ(config.surface_regularization, 1e-8);
        EXPECT_NEAR(config.surface_tension, tensions[index], 1e-18);
        EXPECT_NEAR(config.pressure, pressures[index], 1e-18);
    }
}

TEST(SccsParameters, Vacuum)
{
    const ModuleSccs::SccsConfig config = ModuleSccs::make_sccs_config(ModuleSccs::Preset::Vacuum);
    EXPECT_DOUBLE_EQ(config.cavity.epsilon_bulk, 1.0);
    EXPECT_DOUBLE_EQ(config.surface_tension, 0.0);
    EXPECT_DOUBLE_EQ(config.pressure, 0.0);
}

TEST(SccsParameters, PhysicalUnitConversions)
{
    const double tension = ModuleSccs::dyn_per_cm_to_hartree_per_bohr2(1.0);
    const double pressure = ModuleSccs::gpa_to_hartree_per_bohr3(1.0);
    EXPECT_NEAR(tension, 6.423048558616410e-7, 1e-20);
    EXPECT_NEAR(pressure, 3.398930921743166e-5, 1e-18);
}
