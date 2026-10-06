#include "sccs_cavity.h"
#include "source_base/constants.h"

#include <cmath>

namespace ModuleSccs
{
CavityPoint evaluate_cavity(double density, const CavityParameters& parameters)
{
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
    }
    return point;
}
} // namespace ModuleSccs
