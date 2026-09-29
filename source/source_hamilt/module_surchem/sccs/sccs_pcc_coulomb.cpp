#include "sccs_pcc_coulomb.h"

#include "../common/charge_reduction.h"
#include "../pcc/pcc_moments.h"

#include "source_basis/module_pw/pw_basis.h"

#include <cmath>
#include <stdexcept>

namespace ModuleSccs
{

PccCoulombOperator::PccCoulombOperator(
    const ModulePW::PW_Basis& basis,
    const double tpiba,
    const std::vector<ModuleBase::Vector3<double>>& positions,
    const double volume_element,
    const ModulePcc::PccGeometry& geometry,
    const ModuleSurchem::ChargeReduction& reduction)
    : periodic_(basis, tpiba),
      relative_positions_(positions.size()),
      volume_element_(volume_element),
      geometry_(geometry),
      reduction_(reduction)
{
    ModulePcc::validate_pcc_geometry(geometry_);
    if (positions.size() != static_cast<std::size_t>(basis.nrxx))
    {
        throw std::invalid_argument("SCCS PCC positions must match the local PW real-space grid");
    }
    if (!std::isfinite(volume_element_) || volume_element_ <= 0.0)
    {
        throw std::invalid_argument("SCCS PCC integration geometry must be positive and finite");
    }
    for (std::size_t index = 0; index < positions.size(); ++index)
    {
        relative_positions_[index] = ModulePcc::pcc_relative_position(positions[index], geometry_);
    }
}

// Periodic potential plus the PCC term of the current charge's moments.
void PccCoulombOperator::apply_potential(const std::vector<double>& charge,
                                           std::vector<double>& potential) const
{
    periodic_.apply_potential(charge, potential);
    const ModulePcc::MultipoleMoments local_moments = ModulePcc::density_moments_from_relative_positions(
        charge, relative_positions_, volume_element_);
    const ModulePcc::MultipoleMoments moments = ModulePcc::reduce_pcc_moments(local_moments, reduction_);
    for (std::size_t index = 0; index < charge.size(); ++index)
    {
        potential[index] += ModulePcc::pcc_potential(moments, relative_positions_[index], geometry_.parameters);
    }
}

} // namespace ModuleSccs
