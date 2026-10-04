#include "../sccs_parameters.h"
#include <gtest/gtest.h>

#include <limits>

TEST(SccsParameters, PresetNamesAndFailure)
{
    std::string error;
    ModuleSccs::Preset preset = ModuleSccs::Preset::Custom;
    ASSERT_TRUE(ModuleSccs::parse_preset("Water-Neutral", preset, error));
    EXPECT_EQ(preset, ModuleSccs::Preset::WaterNeutral);
    ASSERT_TRUE(ModuleSccs::parse_preset("WATER-CATION", preset, error));
    EXPECT_EQ(preset, ModuleSccs::Preset::WaterCation);
    ASSERT_TRUE(ModuleSccs::parse_preset("water-anion", preset, error));
    EXPECT_EQ(preset, ModuleSccs::Preset::WaterAnion);
    EXPECT_FALSE(ModuleSccs::parse_preset("unknown", preset, error));
    EXPECT_FALSE(error.empty());
    EXPECT_EQ(preset, ModuleSccs::Preset::WaterAnion);
    ASSERT_TRUE(ModuleSccs::parse_preset("vacuum", preset, error));
    EXPECT_EQ(preset, ModuleSccs::Preset::Vacuum);
    EXPECT_TRUE(error.empty());
    ASSERT_TRUE(ModuleSccs::parse_preset("custom", preset, error));
    EXPECT_EQ(preset, ModuleSccs::Preset::Custom);
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
    std::string error;
    for (int index = 0; index < 3; ++index)
    {
        ModuleSccs::SccsConfig config;
        ASSERT_TRUE(ModuleSccs::make_sccs_config(presets[index], config, error)) << error;
        ASSERT_TRUE(ModuleSccs::validate_config(config, error)) << error;
        EXPECT_DOUBLE_EQ(config.cavity.density_min, minima[index]);
        EXPECT_DOUBLE_EQ(config.cavity.density_max, maxima[index]);
        EXPECT_DOUBLE_EQ(config.cavity.epsilon_bulk, 78.3);
        EXPECT_DOUBLE_EQ(config.surface_regularization, 1e-8);
        EXPECT_NEAR(config.surface_tension, tensions[index], 1e-18);
        EXPECT_NEAR(config.pressure, pressures[index], 1e-18);
    }
}

TEST(SccsParameters, VacuumAndExplicitCustom)
{
    ModuleSccs::SccsConfig config;
    std::string error;
    ASSERT_TRUE(ModuleSccs::make_sccs_config(ModuleSccs::Preset::Vacuum, config, error));
    ASSERT_TRUE(ModuleSccs::validate_config(config, error));
    EXPECT_DOUBLE_EQ(config.cavity.epsilon_bulk, 1.0);
    EXPECT_DOUBLE_EQ(config.surface_tension, 0.0);
    EXPECT_DOUBLE_EQ(config.pressure, 0.0);
    config.cavity.epsilon_bulk = 12.0;
    EXPECT_FALSE(ModuleSccs::make_sccs_config(ModuleSccs::Preset::Custom, config, error));
    EXPECT_FALSE(error.empty());
    EXPECT_DOUBLE_EQ(config.cavity.epsilon_bulk, 12.0);
    EXPECT_TRUE(ModuleSccs::validate_config(config, error));
}

TEST(SccsParameters, PhysicalUnitConversions)
{
    const double tension = ModuleSccs::dyn_per_cm_to_hartree_per_bohr2(1.0);
    const double pressure = ModuleSccs::gpa_to_hartree_per_bohr3(1.0);
    EXPECT_NEAR(tension, 6.423048558616410e-7, 1e-20);
    EXPECT_NEAR(pressure, 3.398930921743166e-5, 1e-18);
}

TEST(SccsParameters, InvalidConfiguration)
{
    ModuleSccs::SccsConfig config;
    std::string error;
    ASSERT_TRUE(ModuleSccs::make_sccs_config(ModuleSccs::Preset::WaterNeutral, config, error));
    config.surface_regularization = 0.0;
    EXPECT_FALSE(ModuleSccs::validate_config(config, error));
    EXPECT_FALSE(error.empty());
    config.surface_regularization = 1e-8;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    config.pressure = nan;
    EXPECT_FALSE(ModuleSccs::validate_config(config, error));
    config.pressure = -0.1;
    config.surface_tension = nan;
    EXPECT_FALSE(ModuleSccs::validate_config(config, error));
    config.surface_tension = 0.0;
    config.cavity.density_min = -1.0;
    EXPECT_FALSE(ModuleSccs::validate_config(config, error));
    config.cavity.density_min = 1e-4;
    EXPECT_TRUE(ModuleSccs::validate_config(config, error));
    EXPECT_TRUE(error.empty());
}
