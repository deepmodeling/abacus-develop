#include "sccs_poisson.h"

#include <cmath>
#include <stdexcept>

namespace ModuleSccs
{

void PolarizationReduction::reduce_sum(double& value) const
{
    if (!std::isfinite(value))
    {
        throw std::domain_error("SCCS reduction input must be finite");
    }
}

void SerialPolarizationReduction::reduce_residual(double& square_sum,
                                                  double& maximum,
                                                  double& point_count) const
{
    if (!std::isfinite(square_sum) || !std::isfinite(maximum) || !std::isfinite(point_count))
    {
        throw std::domain_error("SCCS residual values must be finite before reduction");
    }
}

} // namespace ModuleSccs
