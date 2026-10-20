#include "sccs_parameters.h"

#include "source_base/constants.h"
#include "source_base/tool_quit.h"

namespace ModuleSccs
{

double dyn_cm_to_au(double value)
{
    const double joule_per_square_metre = value * 1e-3;
    return joule_per_square_metre * ModuleBase::BOHR_RADIUS_SI * ModuleBase::BOHR_RADIUS_SI / ModuleBase::HARTREE_SI;
}

double gpa_to_au(double value)
{
    const double pascal = value * 1e9;
    return pascal * ModuleBase::BOHR_RADIUS_SI * ModuleBase::BOHR_RADIUS_SI * ModuleBase::BOHR_RADIUS_SI
           / ModuleBase::HARTREE_SI;
}

Preset parse_preset(const std::string& name)
{
    if (name == "custom") { return Preset::Custom; }
    if (name == "vacuum") { return Preset::Vacuum; }
    if (name == "water-neutral") { return Preset::WaterNeutral; }
    if (name == "water-cation") { return Preset::WaterCation; }
    if (name == "water-anion") { return Preset::WaterAnion; }
    const std::string message = "Unknown SCCS preset: " + name;
    ModuleBase::WARNING_QUIT("ModuleSccs::parse_preset", message);
}

SccsConfig preset_config(Preset preset)
{
    if (preset == Preset::Custom)
    {
        ModuleBase::WARNING_QUIT("ModuleSccs::preset_config",
                                 "SCCS custom configurations must be supplied explicitly");
    }
    SccsConfig config;
    config.cavity.epsilon_bulk = 78.3;
    config.surface_regularization = 1e-8;
    if (preset == Preset::Vacuum || preset == Preset::WaterNeutral)
    {
        config.cavity.density_min = 1e-4;
        config.cavity.density_max = 5e-3;
        if (preset == Preset::Vacuum)
        {
            config.cavity.epsilon_bulk = 1.0;
        }
        else
        {
            config.surface_tension = dyn_cm_to_au(47.9);
            config.pressure = gpa_to_au(-0.36);
        }
    }
    else if (preset == Preset::WaterCation)
    {
        config.cavity.density_min = 2e-4;
        config.cavity.density_max = 3.5e-3;
        config.surface_tension = dyn_cm_to_au(5.0);
        config.pressure = gpa_to_au(0.125);
    }
    else
    {
        config.cavity.density_min = 2.4e-3;
        config.cavity.density_max = 1.55e-2;
        config.pressure = gpa_to_au(0.45);
    }
    return config;
}
} // namespace ModuleSccs
