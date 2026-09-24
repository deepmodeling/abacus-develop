#include "sccs_parameters.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>

namespace ModuleSccs
{
namespace
{

const double HARTREE_JOULE = 4.3597447222071e-18;
const double BOHR_METRE = 5.29177210903e-11;

std::string lower_case(const std::string& value)
{
    std::string result = value;
    std::transform(result.begin(), result.end(), result.begin(), [](const char character) {
        return static_cast<char>(std::tolower(static_cast<unsigned char>(character)));
    });
    return result;
}

SccsConfig common_water_config()
{
    SccsConfig config;
    config.cavity.epsilon_bulk = 78.3;
    config.surface_regularization = 1.0e-8;
    config.boundary = Boundary::Periodic;
    config.max_iterations = 200;
    config.mixing_method = "linear";
    config.mixing_history = 8;
    config.mixing = 0.5;
    config.adaptive_mixing = false;
    config.mixing_min = 0.1;
    config.mixing_max = 0.8;
    config.tolerance_rms = 1.0e-10;
    config.tolerance_max = 1.0e-8;
    return config;
}

} // namespace

Preset parse_preset(const std::string& value)
{
    const std::string normalized = lower_case(value);
    if (normalized == "custom")
    {
        return Preset::Custom;
    }
    if (normalized == "water-neutral")
    {
        return Preset::WaterNeutral;
    }
    if (normalized == "water-cation")
    {
        return Preset::WaterCation;
    }
    if (normalized == "water-anion")
    {
        return Preset::WaterAnion;
    }
    throw std::invalid_argument("unknown SCCS preset: " + value);
}

Boundary parse_boundary(const std::string& value)
{
    const std::string normalized = lower_case(value);
    if (normalized == "periodic")
    {
        return Boundary::Periodic;
    }
    if (normalized == "pcc_0d" || normalized == "pcc-0d")
    {
        return Boundary::Pcc0d;
    }
    if (normalized == "pcc_2d" || normalized == "pcc-2d")
    {
        return Boundary::Pcc2d;
    }
    throw std::invalid_argument("unknown SCCS boundary: " + value);
}

SccsConfig water_preset(const Preset preset)
{
    SccsConfig config = common_water_config();
    if (preset == Preset::WaterNeutral)
    {
        config.cavity.density_min = 1.0e-4;
        config.cavity.density_max = 5.0e-3;
        config.surface_tension = dyn_per_cm_to_hartree_per_bohr2(47.9);
        config.pressure = gpa_to_hartree_per_bohr3(-0.36);
    }
    else if (preset == Preset::WaterCation)
    {
        config.cavity.density_min = 2.0e-4;
        config.cavity.density_max = 3.5e-3;
        config.surface_tension = dyn_per_cm_to_hartree_per_bohr2(5.0);
        config.pressure = gpa_to_hartree_per_bohr3(0.125);
    }
    else if (preset == Preset::WaterAnion)
    {
        config.cavity.density_min = 2.4e-3;
        config.cavity.density_max = 1.55e-2;
        config.surface_tension = 0.0;
        config.pressure = gpa_to_hartree_per_bohr3(0.45);
    }
    else
    {
        throw std::invalid_argument("custom SCCS parameters cannot be generated from a water preset");
    }
    validate_config(config);
    return config;
}

void validate_config(const SccsConfig& config)
{
    validate_cavity_parameters(config.cavity);
    if (!std::isfinite(config.surface_tension) || !std::isfinite(config.pressure))
    {
        throw std::invalid_argument("SCCS surface tension and pressure must be finite");
    }
    if (!std::isfinite(config.surface_regularization) || config.surface_regularization <= 0.0)
    {
        throw std::invalid_argument("SCCS surface regularization must be positive and finite");
    }
    if (config.max_iterations <= 0)
    {
        throw std::invalid_argument("SCCS maximum iteration count must be positive");
    }
    if (config.mixing_method != "linear" && config.mixing_method != "pulay"
        && config.mixing_method != "anderson")
    {
        throw std::invalid_argument("unknown SCCS polarization mixing method: "
                                    + config.mixing_method);
    }
    if (config.mixing_history < 2)
    {
        throw std::invalid_argument("SCCS accelerated-mixing history must be at least two");
    }
    if (!std::isfinite(config.mixing) || config.mixing <= 0.0 || config.mixing > 1.0)
    {
        throw std::invalid_argument("SCCS mixing must be in the interval (0, 1]");
    }
    if (!std::isfinite(config.mixing_min) || !std::isfinite(config.mixing_max)
        || config.mixing_min <= 0.0 || config.mixing_max > 1.0
        || config.mixing_min > config.mixing_max)
    {
        throw std::invalid_argument("SCCS adaptive-mixing bounds must satisfy 0 < min <= max <= 1");
    }
    if (config.adaptive_mixing
        && (config.mixing < config.mixing_min || config.mixing > config.mixing_max))
    {
        throw std::invalid_argument("SCCS initial mixing must lie within the adaptive-mixing bounds");
    }
    if (!std::isfinite(config.tolerance_rms) || config.tolerance_rms <= 0.0
        || !std::isfinite(config.tolerance_max) || config.tolerance_max <= 0.0)
    {
        throw std::invalid_argument("SCCS residual tolerances must be positive and finite");
    }
}

double dyn_per_cm_to_hartree_per_bohr2(const double value)
{
    const double joule_per_square_metre = value * 1.0e-3;
    return joule_per_square_metre * BOHR_METRE * BOHR_METRE / HARTREE_JOULE;
}

double gpa_to_hartree_per_bohr3(const double value)
{
    const double pascal = value * 1.0e9;
    return pascal * BOHR_METRE * BOHR_METRE * BOHR_METRE / HARTREE_JOULE;
}

} // namespace ModuleSccs
