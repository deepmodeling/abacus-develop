#include "sccs_cavity.h"
#include "source_base/constants.h"

#include <cmath>

namespace ModuleSccs
{
bool validate_cavity_parameters(const CavityParameters& parameters, std::string& error)
{
    error.clear();
    if (!std::isfinite(parameters.density_min) || !std::isfinite(parameters.density_max)
        || parameters.density_min <= 0.0 || parameters.density_max <= parameters.density_min)
    {
        error = "SCCS requires finite density thresholds with 0 < density_min < density_max";
        return false;
    }
    if (!std::isfinite(parameters.epsilon_bulk) || parameters.epsilon_bulk < 1.0)
    {
        error = "SCCS requires a finite bulk dielectric constant >= 1";
        return false;
    }
    return true;
}

bool evaluate_cavity(double density,
                     const CavityParameters& parameters,
                     CavityPoint& result,
                     std::string& error)
{
    if (!validate_cavity_parameters(parameters, error))
    {
        return false;
    }
    if (!std::isfinite(density))
    {
        error = "SCCS cavity requires a finite electronic density";
        return false;
    }

    CavityPoint point;
    if (density <= parameters.density_min)
    {
        point.epsilon = parameters.epsilon_bulk;
    }
    else if (density >= parameters.density_max)
    {
        point.solute = 1.0;
        point.epsilon = 1.0;
    }
    else
    {
        const double density_ratio = parameters.density_max / parameters.density_min;
        if (!std::isfinite(density_ratio))
        {
            error = "SCCS cavity density ratio exceeds the numerical range";
            return false;
        }
        const double log_width = std::log(density_ratio);
        const double relative_density = parameters.density_max / density;
        const double log_relative_density = std::log(relative_density);
        const double x = log_relative_density / log_width;
        const double angle = ModuleBase::TWO_PI * x;
        const double solvent = x - std::sin(angle) / ModuleBase::TWO_PI;
        const double dsolvent_drho = -(1.0 - std::cos(angle)) / (log_width * density);
        const double log_epsilon = std::log(parameters.epsilon_bulk);
        const double log_epsilon_point = log_epsilon * solvent;
        point.solute = 1.0 - solvent;
        point.dsolute_drho = -dsolvent_drho;
        point.epsilon = std::exp(log_epsilon_point);
        point.depsilon_drho = point.epsilon * log_epsilon * dsolvent_drho;
        if (!std::isfinite(point.solute) || !std::isfinite(point.dsolute_drho)
            || !std::isfinite(point.epsilon) || !std::isfinite(point.depsilon_drho))
        {
            error = "SCCS cavity transition exceeds the numerical range";
            return false;
        }
    }
    result = point;
    return true;
}
} // namespace ModuleSccs
