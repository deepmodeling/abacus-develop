#include "sccs_pcc_coulomb.h"

#include "../common/charge_reduction.h"
#include "../pcc/pcc_moments.h"

#include "source_basis/module_pw/pw_basis.h"
#include "source_base/constants.h"

#include <cmath>
#include <stdexcept>

namespace ModuleSccs
{

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

// Scalar clients retain the same PCC moments and gauge without computing
// three unused gradient transforms.
void PccCoulombOperator::apply_potential(const std::vector<double>& charge,
                                           std::vector<double>& potential) const
{
    periodic_.apply_potential(charge, potential);
    const MultipoleMoments local_moments = density_moments_from_relative_positions(
        charge, relative_positions_, volume_element_);
    const MultipoleMoments moments = reduce_pcc_moments(local_moments, reduction_);
    for (std::size_t index = 0; index < charge.size(); ++index)
    {
        potential[index] += pcc_potential(moments, relative_positions_[index], geometry_.parameters);
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

// Apply the exact discrete transpose of apply_gradient, including PCC.
// Moment reductions account for all three Cartesian components; they preserve
// <u, Gq> = <G^T u, q> with the uniform-grid volume weight.
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
