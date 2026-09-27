#include "sccs_pcc_coulomb.h"

#include "../sccs/sccs_charge.h"

#include "source_basis/module_pw/pw_basis.h"
#include "source_base/constants.h"

#include <cmath>
#include <stdexcept>

namespace ModuleSccs
{
namespace
{

MultipoleMoments reduce_pcc_moments(MultipoleMoments moments,
                                     const ChargeReduction& reduction)
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

} // namespace

MultipoleMoments reduced_pcc_density_moments(
    const std::vector<double>& density,
    const std::vector<ModuleBase::Vector3<double>>& positions,
    const double volume_element,
    const PccGeometry& geometry,
    const ChargeReduction& reduction)
{
    return reduce_pcc_moments(
        density_moments(density, positions, volume_element, geometry),
        reduction);
}

PccCoulombOperator::PccCoulombOperator(
    const ModulePW::PW_Basis& basis,
    const double tpiba,
    const std::vector<ModuleBase::Vector3<double>>& positions,
    const double volume_element,
    const PccGeometry& geometry,
    const ChargeReduction& reduction)
    : periodic_(basis, tpiba),
      positions_(positions),
      relative_positions_(positions.size()),
      volume_element_(volume_element),
      geometry_(geometry),
      reduction_(reduction)
{
    validate_pcc_geometry(geometry_);
    if (positions_.size() != static_cast<std::size_t>(basis.nrxx))
    {
        throw std::invalid_argument("SCCS PCC positions must match the local PW real-space grid");
    }
    if (!std::isfinite(volume_element_) || volume_element_ <= 0.0)
    {
        throw std::invalid_argument("SCCS PCC integration geometry must be positive and finite");
    }
    for (std::size_t index = 0; index < positions_.size(); ++index)
    {
        relative_positions_[index] = pcc_relative_position(positions_[index], geometry_);
    }
}

void PccCoulombOperator::apply(const std::vector<double>& charge,
                               ElectrostaticField& field) const
{
    periodic_.apply(charge, field);
    const MultipoleMoments moments = reduce_pcc_moments(
        density_moments_from_relative_positions(charge,
                                                relative_positions_,
                                                volume_element_),
        reduction_);
    for (std::size_t index = 0; index < charge.size(); ++index)
    {
        field.potential[index]
            += pcc_potential(moments,
                             relative_positions_[index],
                             geometry_.parameters);
        const ModuleBase::Vector3<double> correction
            = pcc_potential_gradient(moments,
                                     relative_positions_[index],
                                     geometry_.parameters);
        field.gradient[index].x += correction.x;
        field.gradient[index].y += correction.y;
        field.gradient[index].z += correction.z;
    }
}

void PccCoulombOperator::apply_gradient(
    const std::vector<double>& charge,
    std::vector<ModuleBase::Vector3<double>>& gradient) const
{
    periodic_.apply_gradient(charge, gradient);
    const MultipoleMoments local_moments
        = density_moments_from_relative_positions(charge, relative_positions_, volume_element_);
    const MultipoleMoments moments = reduce_pcc_moments(local_moments, reduction_);
    for (std::size_t index = 0; index < charge.size(); ++index)
    {
        const ModuleBase::Vector3<double> correction
            = pcc_potential_gradient(moments,
                                     relative_positions_[index],
                                     geometry_.parameters);
        gradient[index].x += correction.x;
        gradient[index].y += correction.y;
        gradient[index].z += correction.z;
    }
}

void PccCoulombOperator::apply_gradient_adjoint(
    const std::vector<ModuleBase::Vector3<double>>& field,
    std::vector<double>& result) const
{
    periodic_.apply_gradient_adjoint(field, result);
    double moments[4] = {};
    for (std::size_t i = 0; i < field.size(); ++i)
    {
        for (int d = 0; d < 3; ++d)
        {
            moments[d] += field[i][d] * volume_element_;
            moments[3] += field[i][d] * relative_positions_[i][d] * volume_element_;
        }
    }
    reduction_.reduce_sum(moments, 4);
    const double length = geometry_.parameters.cube_length;
    const double factor = -ModuleBase::FOUR_PI / (3.0 * length * length * length);
    for (std::size_t i = 0; i < field.size(); ++i)
    {
        const double projection = relative_positions_[i].x * moments[0]
                                  + relative_positions_[i].y * moments[1]
                                  + relative_positions_[i].z * moments[2];
        result[i] += factor * (moments[3] - projection);
    }
}

} // namespace ModuleSccs
