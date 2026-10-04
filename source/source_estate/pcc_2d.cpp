#include "pcc_2d.h"

#include "source_base/constants.h"
#include "source_cell/cell_geometry.h"

#include <cmath>

namespace elecstate
{
bool make_pcc_2d_parameters(const unitcell::SlabCell& cell,
                            Pcc2dParameters& parameters,
                            std::string& error)
{
    error.clear();
    if (!std::isfinite(cell.area) || cell.area <= 0.0 || !std::isfinite(cell.length) || cell.length <= 0.0)
    {
        error = "PCC 2D requires positive finite periodic area and normal length";
        return false;
    }
    parameters.area = cell.area;
    parameters.length = cell.length;
    return true;
}

double pcc_2d_potential(const ChargeMoments& moments,
                        const double coordinate,
                        const ModuleBase::Vector3<double>& normal,
                        const Pcc2dParameters& parameters)
{
    const double dipole = moments.dipole * normal;
    const double parabolic = moments.charge * coordinate * coordinate
                             - 2.0 * dipole * coordinate + moments.second_moment;
    const double constant = -ModuleBase::PI * parameters.length / (3.0 * parameters.area);
    const double factor = 2.0 * ModuleBase::PI / (parameters.area * parameters.length);
    return constant * moments.charge - factor * parabolic;
}

ModuleBase::Vector3<double> pcc_2d_gradient(const ChargeMoments& moments,
                                           const double coordinate,
                                           const ModuleBase::Vector3<double>& normal,
                                           const Pcc2dParameters& parameters)
{
    const double dipole = moments.dipole * normal;
    const double factor = -4.0 * ModuleBase::PI / (parameters.area * parameters.length);
    const double derivative = factor * (moments.charge * coordinate - dipole);
    return normal * derivative;
}

double pcc_2d_energy(const ChargeMoments& moments, const Pcc2dParameters& parameters)
{
    const double constant = -ModuleBase::PI * parameters.length / (3.0 * parameters.area);
    const double factor = 2.0 * ModuleBase::PI / (parameters.area * parameters.length);
    const double multipole = 2.0 * (moments.charge * moments.second_moment - moments.dipole.norm2());
    return 0.5 * (constant * moments.charge * moments.charge - factor * multipole);
}

ModuleBase::Vector3<double> pcc_2d_force(const ChargeMoments& moments,
                                        const double ionic_charge,
                                        const double coordinate,
                                        const ModuleBase::Vector3<double>& normal,
                                        const Pcc2dParameters& parameters)
{
    const ModuleBase::Vector3<double> gradient = pcc_2d_gradient(moments, coordinate, normal, parameters);
    return gradient * (-ionic_charge);
}
} // namespace elecstate
