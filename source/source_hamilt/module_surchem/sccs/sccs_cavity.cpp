#include "sccs_cavity.h"

#include "source_base/constants.h"

#include <cmath>
#include <stdexcept>

namespace ModuleSccs
{

void validate_cavity_parameters(const CavityParameters& parameters)
{
    if (!std::isfinite(parameters.density_min) || !std::isfinite(parameters.density_max)
        || !std::isfinite(parameters.epsilon_bulk))
    {
        throw std::invalid_argument("SCCS cavity parameters must be finite");
    }
    if (parameters.density_min <= 0.0 || parameters.density_max <= parameters.density_min)
    {
        throw std::invalid_argument("SCCS cavity requires 0 < density_min < density_max");
    }
    if (parameters.epsilon_bulk < 1.0)
    {
        throw std::invalid_argument("SCCS bulk permittivity must be at least one");
    }
}

CavityPoint evaluate_cavity(const double density, const CavityParameters& parameters)
{
    validate_cavity_parameters(parameters);
    if (!std::isfinite(density))
    {
        throw std::domain_error("SCCS cavity density must be finite");
    }

    CavityPoint result;
    // LCAO-to-grid transforms can produce negative Fourier ringing in the
    // vacuum. Every such value is below density_min and therefore belongs to
    // the constant bulk-solvent branch, whose density derivative is zero.
    if (density <= parameters.density_min)
    {
        result.epsilon = parameters.epsilon_bulk;
        return result;
    }
    if (density >= parameters.density_max)
    {
        result.solute = 1.0;
        result.epsilon = 1.0;
        return result;
    }

    const double log_width = std::log(parameters.density_max / parameters.density_min);
    const double x = std::log(parameters.density_max / density) / log_width;
    const double solvent = x - std::sin(ModuleBase::TWO_PI * x) / ModuleBase::TWO_PI;
    const double dsolvent_drho
        = -(1.0 - std::cos(ModuleBase::TWO_PI * x)) / (log_width * density);
    const double log_epsilon = std::log(parameters.epsilon_bulk);

    result.solute = 1.0 - solvent;
    result.dsolute_drho = -dsolvent_drho;
    result.epsilon = std::exp(log_epsilon * solvent);
    result.depsilon_drho = result.epsilon * log_epsilon * dsolvent_drho;
    return result;
}

} // namespace ModuleSccs
