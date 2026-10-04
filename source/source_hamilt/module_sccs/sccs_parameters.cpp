#include "sccs_parameters.h"

#include <cctype>
#include <cmath>

namespace ModuleSccs
{
namespace
{
// CODATA 2018 SI values retain the original SCCS preset conversion.
// ModuleBase's older SI values differ, so using them would change these presets.
const double hartree_joule = 4.3597447222071e-18;
const double bohr_metre = 5.29177210903e-11;
}

double dyn_per_cm_to_hartree_per_bohr2(double value)
{
    const double joule_per_square_metre = value * 1e-3;
    return joule_per_square_metre * bohr_metre * bohr_metre / hartree_joule;
}

double gpa_to_hartree_per_bohr3(double value)
{
    const double pascal = value * 1e9;
    return pascal * bohr_metre * bohr_metre * bohr_metre / hartree_joule;
}

bool parse_preset(const std::string& name, Preset& preset, std::string& error)
{
    error.clear();
    std::string normalized = name;
    for (char& character : normalized)
    {
        const unsigned char byte = static_cast<unsigned char>(character);
        const int lower = std::tolower(byte);
        character = static_cast<char>(lower);
    }
    if (normalized == "custom")
    {
        preset = Preset::Custom;
    }
    else if (normalized == "vacuum")
    {
        preset = Preset::Vacuum;
    }
    else if (normalized == "water-neutral")
    {
        preset = Preset::WaterNeutral;
    }
    else if (normalized == "water-cation")
    {
        preset = Preset::WaterCation;
    }
    else if (normalized == "water-anion")
    {
        preset = Preset::WaterAnion;
    }
    else
    {
        error = "Unknown SCCS preset: " + name;
        return false;
    }
    return true;
}

bool make_sccs_config(Preset preset, SccsConfig& config, std::string& error)
{
    error.clear();
    SccsConfig candidate;
    candidate.cavity.epsilon_bulk = 78.3;
    candidate.surface_regularization = 1e-8;
    if (preset == Preset::Vacuum || preset == Preset::WaterNeutral)
    {
        candidate.cavity.density_min = 1e-4;
        candidate.cavity.density_max = 5e-3;
        if (preset == Preset::Vacuum)
        {
            candidate.cavity.epsilon_bulk = 1.0;
        }
        else
        {
            candidate.surface_tension = dyn_per_cm_to_hartree_per_bohr2(47.9);
            candidate.pressure = gpa_to_hartree_per_bohr3(-0.36);
        }
    }
    else if (preset == Preset::WaterCation)
    {
        candidate.cavity.density_min = 2e-4;
        candidate.cavity.density_max = 3.5e-3;
        candidate.surface_tension = dyn_per_cm_to_hartree_per_bohr2(5.0);
        candidate.pressure = gpa_to_hartree_per_bohr3(0.125);
    }
    else if (preset == Preset::WaterAnion)
    {
        candidate.cavity.density_min = 2.4e-3;
        candidate.cavity.density_max = 1.55e-2;
        candidate.pressure = gpa_to_hartree_per_bohr3(0.45);
    }
    else
    {
        error = "SCCS custom configurations must be supplied explicitly";
        return false;
    }
    config = candidate;
    return true;
}

bool validate_config(const SccsConfig& config, std::string& error)
{
    if (!validate_cavity_parameters(config.cavity, error))
    {
        return false;
    }
    if (!std::isfinite(config.surface_tension) || !std::isfinite(config.pressure))
    {
        error = "SCCS requires finite surface tension and pressure";
        return false;
    }
    if (!std::isfinite(config.surface_regularization) || config.surface_regularization <= 0.0)
    {
        error = "SCCS requires finite positive surface regularization";
        return false;
    }
    return true;
}
} // namespace ModuleSccs
