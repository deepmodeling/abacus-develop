#include "sccs_pcc_2d_coulomb.h"

#include "../common/charge_reduction.h"
#include "../pcc/pcc_moments.h"

#include "source_basis/module_pw/pw_basis.h"
#include "source_base/constants.h"

#include <cmath>
#include <stdexcept>

namespace ModuleSccs
{

Pcc2dCoulombOperator::Pcc2dCoulombOperator(
    const ModulePW::PW_Basis& basis,
    const double tpiba,
    const std::vector<ModuleBase::Vector3<double>>& positions,
    const double volume_element,
    const Pcc2dGeometry& geometry,
    const ChargeReduction& reduction)
    : periodic_(basis, tpiba),
      positions_(positions),
      relative_y_(positions.size()),
      volume_element_(volume_element),
      geometry_(geometry),
      reduction_(reduction)
{
    validate_pcc_2d_parameters(geometry_.parameters);
    if (positions_.size() != static_cast<std::size_t>(basis.nrxx))
    {
        throw std::invalid_argument(
            "two-dimensional SCCS PCC positions must match the local PW grid");
    }
    if (!std::isfinite(volume_element_) || volume_element_ <= 0.0
        || !std::isfinite(geometry_.origin_y))
    {
        throw std::invalid_argument(
            "two-dimensional SCCS PCC integration geometry must be finite and positive");
    }
    for (std::size_t index = 0; index < positions_.size(); ++index)
    {
        relative_y_[index] = pcc_2d_relative_y(positions_[index].y, geometry_);
    }
}

void Pcc2dCoulombOperator::apply(const std::vector<double>& charge,
                                 ElectrostaticField& field) const
{
    periodic_.apply(charge, field);
    const Pcc2dMoments moments
        = reduce_pcc_2d_moments(
            pcc_2d_density_moments_from_relative_y(charge,
                                                   relative_y_,
                                                   volume_element_),
            reduction_);
    for (std::size_t index = 0; index < charge.size(); ++index)
    {
        field.potential[index]
            += pcc_2d_potential(moments, relative_y_[index], geometry_.parameters);
        field.gradient[index].y
            += pcc_2d_potential_gradient(moments,
                                         relative_y_[index],
                                         geometry_.parameters).y;
    }
}

// Scalar clients retain the same PCC moments and gauge without computing
// three unused gradient transforms.
void Pcc2dCoulombOperator::apply_potential(const std::vector<double>& charge,
                                           std::vector<double>& potential) const
{
    periodic_.apply_potential(charge, potential);
    const Pcc2dMoments local_moments = pcc_2d_density_moments_from_relative_y(
        charge, relative_y_, volume_element_);
    const Pcc2dMoments moments = reduce_pcc_2d_moments(local_moments, reduction_);
    for (std::size_t index = 0; index < charge.size(); ++index)
    {
        potential[index] += pcc_2d_potential(moments, relative_y_[index], geometry_.parameters);
    }
}

void Pcc2dCoulombOperator::apply_gradient(
    const std::vector<double>& charge,
    std::vector<ModuleBase::Vector3<double>>& gradient) const
{
    periodic_.apply_gradient(charge, gradient);
    const Pcc2dMoments local_moments
        = pcc_2d_density_moments_from_relative_y(charge, relative_y_, volume_element_);
    const Pcc2dMoments moments = reduce_pcc_2d_moments(local_moments, reduction_);
    for (std::size_t index = 0; index < charge.size(); ++index)
    {
        gradient[index].y
            += pcc_2d_potential_gradient(moments,
                                         relative_y_[index],
                                         geometry_.parameters).y;
    }
}

// Apply the exact discrete transpose of apply_gradient, including PCC.
// Moment reductions account for the open y component only; they preserve
// <u, Gq> = <G^T u, q> with the uniform-grid volume weight.
void Pcc2dCoulombOperator::apply_gradient_adjoint(
    const std::vector<ModuleBase::Vector3<double>>& field,
    std::vector<double>& result) const
{
    periodic_.apply_gradient_adjoint(field, result);
    double moments[2] = {};
    for (std::size_t i = 0; i < field.size(); ++i)
    {
        moments[0] += field[i].y * volume_element_;
        moments[1] += field[i].y * relative_y_[i] * volume_element_;
    }
    reduction_.reduce_sum(moments, 2);
    const double volume = geometry_.parameters.periodic_area * geometry_.parameters.cell_length_y;
    const double factor = -ModuleBase::FOUR_PI / volume;
    for (std::size_t i = 0; i < field.size(); ++i)
    {
        result[i] += factor * (moments[1] - relative_y_[i] * moments[0]);
    }
}

} // namespace ModuleSccs
