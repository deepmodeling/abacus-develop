#include "pcc_moments.h"

#include "../common/charge_reduction.h"

#include <cmath>
#include <stdexcept>

namespace ModulePcc
{

MultipoleMoments reduce_pcc_moments(MultipoleMoments moments,
                                     const ModuleSurchem::ChargeReduction& reduction)
{
    double values[5] = {moments.charge,
                        moments.dipole.x,
                        moments.dipole.y,
                        moments.dipole.z,
                        moments.quadrupole_trace};
    reduction.reduce_sum(values, 5);
    for (int index = 0; index < 5; ++index)
    {
        if (!std::isfinite(values[index]))
        {
            throw std::domain_error("zero-dimensional PCC reduced moments must be finite");
        }
    }
    moments.charge = values[0];
    moments.dipole.x = values[1];
    moments.dipole.y = values[2];
    moments.dipole.z = values[3];
    moments.quadrupole_trace = values[4];
    return moments;
}

MultipoleMoments reduced_pcc_density_moments(
    const std::vector<double>& density,
    const std::vector<ModuleBase::Vector3<double>>& positions,
    const double volume_element,
    const PccGeometry& geometry,
    const ModuleSurchem::ChargeReduction& reduction)
{
    const MultipoleMoments moments = density_moments(density, positions, volume_element, geometry);
    return reduce_pcc_moments(moments, reduction);
}

Pcc2dMoments reduce_pcc_2d_moments(Pcc2dMoments moments,
                                    const ModuleSurchem::ChargeReduction& reduction)
{
    double values[3] = {moments.charge, moments.dipole, moments.quadrupole};
    reduction.reduce_sum(values, 3);
    if (!std::isfinite(values[0]) || !std::isfinite(values[1])
        || !std::isfinite(values[2]))
    {
        throw std::domain_error("two-dimensional PCC reduced moments must be finite");
    }
    moments.charge = values[0];
    moments.dipole = values[1];
    moments.quadrupole = values[2];
    return moments;
}

Pcc2dMoments reduced_pcc_2d_density_moments(
    const std::vector<double>& density,
    const std::vector<ModuleBase::Vector3<double>>& positions,
    const double volume_element,
    const Pcc2dGeometry& geometry,
    const ModuleSurchem::ChargeReduction& reduction)
{
    const Pcc2dMoments moments = pcc_2d_density_moments(density, positions, volume_element, geometry);
    return reduce_pcc_2d_moments(moments, reduction);
}

} // namespace ModulePcc
