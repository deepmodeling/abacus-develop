#include "pcc_0d.h"

#include "source_base/constants.h"
#include "source_cell/cell_geometry.h"

#include <cmath>

namespace elecstate
{

bool make_pcc_0d_parameters(const unitcell::OrthogonalCell& cell,
                            const double relative_tolerance,
                            Pcc0dParameters& parameters,
                            std::string& error)
{
    error.clear();
    if (!std::isfinite(relative_tolerance) || relative_tolerance <= 0.0)
    {
        error = "PCC 0D requires a positive finite geometry tolerance";
        return false;
    }
    const double length = cell.lengths[0];
    for (int axis = 0; axis < 3; ++axis)
    {
        const double edge = cell.lengths[axis];
        const double difference = edge - length;
        if (!std::isfinite(edge) || edge <= 0.0 || std::abs(difference) > relative_tolerance * length)
        {
            error = "PCC 0D requires an equal-edge cubic cell";
            return false;
        }
    }
    Pcc0dParameters candidate;
    candidate.length = length;
    parameters = candidate;
    return true;
}

double pcc_0d_potential(const ChargeMoments& moments,
                        const ModuleBase::Vector3<double>& position,
                        const Pcc0dParameters& parameters)
{
    const double length = parameters.length;
    const double volume = length * length * length;
    const double dipole_projection = moments.dipole * position;
    const double parabolic = moments.charge * position.norm2()
                             - 2.0 * dipole_projection + moments.second_moment;
    return parameters.madelung * moments.charge / length
           - 2.0 * ModuleBase::PI * parabolic / (3.0 * volume);
}

ModuleBase::Vector3<double> pcc_0d_gradient(const ChargeMoments& moments,
                                          const ModuleBase::Vector3<double>& position,
                                          const Pcc0dParameters& parameters)
{
    const double length = parameters.length;
    const double volume = length * length * length;
    const double factor = -4.0 * ModuleBase::PI / (3.0 * volume);
    const ModuleBase::Vector3<double> field = position * moments.charge - moments.dipole;
    return field * factor;
}

double pcc_0d_bilinear_energy(const ChargeMoments& left,
                             const ChargeMoments& right,
                             const Pcc0dParameters& parameters)
{
    const double length = parameters.length;
    const double volume = length * length * length;
    const double monopole = parameters.madelung * left.charge * right.charge / length;
    const double dipole_product = left.dipole * right.dipole;
    const double multipole = left.second_moment * right.charge
                             + left.charge * right.second_moment - 2.0 * dipole_product;
    return monopole - 2.0 * ModuleBase::PI * multipole / (3.0 * volume);
}

double pcc_0d_energy(const ChargeMoments& moments, const Pcc0dParameters& parameters)
{
    const double bilinear = pcc_0d_bilinear_energy(moments, moments, parameters);
    return 0.5 * bilinear;
}

ModuleBase::Vector3<double> pcc_0d_force(const ChargeMoments& moments,
                                       const double ionic_charge,
                                       const ModuleBase::Vector3<double>& position,
                                       const Pcc0dParameters& parameters)
{
    const ModuleBase::Vector3<double> gradient = pcc_0d_gradient(moments, position, parameters);
    return gradient * (-ionic_charge);
}

} // namespace elecstate
