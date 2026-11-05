#include "sccs_coul_0d.h"

#include "sccs_moments.h"

#include "source_basis/module_pw/pw_basis.h"

namespace ModuleSccs
{
Pcc0dCoulombOperator::Pcc0dCoulombOperator(
    const ModulePW::PW_Basis& basis,
    double tpiba,
    const std::vector<ModuleBase::Vector3<double>>& relative_positions,
    const elecstate::Pcc0dParameters& parameters)
    : basis_(basis), periodic_(basis, tpiba), relative_positions_(relative_positions), parameters_(parameters)
{
}

bool Pcc0dCoulombOperator::has_correction() const
{
    return true;
}

void Pcc0dCoulombOperator::apply_potential(const std::vector<double>& charge, std::vector<double>& potential)
{
    const elecstate::ChargeMoments moments = pool_moments(charge, relative_positions_, basis_);
    periodic_.apply_potential(charge, potential);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        potential[ir] += elecstate::pcc_0d_potential(moments, relative_positions_[ir], parameters_);
    }
}
} // namespace ModuleSccs
