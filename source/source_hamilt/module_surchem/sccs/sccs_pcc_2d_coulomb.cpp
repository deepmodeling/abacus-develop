#include "sccs_pcc_2d_coulomb.h"

#include "../common/charge_reduction.h"
#include "../pcc/pcc_moments.h"

#include "source_basis/module_pw/pw_basis.h"

#include <cmath>
#include <stdexcept>

namespace ModuleSccs
{

Pcc2dCoulombOperator::Pcc2dCoulombOperator(
    const ModulePW::PW_Basis& basis,
    const double tpiba,
    const std::vector<ModuleBase::Vector3<double>>& positions,
    const double volume_element,
    const ModulePcc::Pcc2dGeometry& geometry,
    const ModuleSurchem::ChargeReduction& reduction)
    : periodic_(basis, tpiba),
      relative_y_(positions.size()),
      volume_element_(volume_element),
      geometry_(geometry),
      reduction_(reduction)
{
    ModulePcc::validate_pcc_2d_geometry(geometry_);
    if (positions.size() != static_cast<std::size_t>(basis.nrxx))
    {
        throw std::invalid_argument(
            "two-dimensional SCCS PCC positions must match the local PW grid");
    }
    if (!std::isfinite(volume_element_) || volume_element_ <= 0.0)
    {
        throw std::invalid_argument(
            "two-dimensional SCCS PCC integration geometry must be finite and positive");
    }
    for (std::size_t index = 0; index < positions.size(); ++index)
    {
        relative_y_[index] = ModulePcc::pcc_2d_relative_y(positions[index].y, geometry_);
    }
}

// Periodic potential plus the PCC2D term of the current charge's moments.
void Pcc2dCoulombOperator::apply_potential(const std::vector<double>& charge,
                                           std::vector<double>& potential) const
{
    periodic_.apply_potential(charge, potential);
    const ModulePcc::Pcc2dMoments local_moments = ModulePcc::pcc_2d_density_moments_from_relative_y(
        charge, relative_y_, volume_element_);
    const ModulePcc::Pcc2dMoments moments = ModulePcc::reduce_pcc_2d_moments(local_moments, reduction_);
    for (std::size_t index = 0; index < charge.size(); ++index)
    {
        potential[index] += ModulePcc::pcc_2d_potential(moments, relative_y_[index], geometry_.parameters);
    }
}

} // namespace ModuleSccs
