#include "sccs_coul_2d.h"

#include "sccs_moments.h"

#include "source_basis/module_pw/pw_basis.h"

namespace ModuleSccs
{
Pcc2dCoulombOperator::Pcc2dCoulombOperator(
    const ModulePW::PW_Basis& basis,
    double tpiba,
    const std::vector<ModuleBase::Vector3<double>>& projected_positions,
    const ModuleBase::Vector3<double>& normal,
    const elecstate::Pcc2dParameters& parameters)
    : basis_(basis), periodic_(basis, tpiba), projected_positions_(projected_positions),
      normal_(normal), parameters_(parameters)
{
}

bool Pcc2dCoulombOperator::has_correction() const
{
    return true;
}

void Pcc2dCoulombOperator::apply_potential(const std::vector<double>& charge, std::vector<double>& potential)
{
    const elecstate::ChargeMoments moments = pool_moments(charge, projected_positions_, basis_);
    periodic_.apply_potential(charge, potential);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        const double coordinate = projected_positions_[ir] * normal_;
        potential[ir] += elecstate::pcc_2d_potential(moments, coordinate, normal_, parameters_);
    }
}
} // namespace ModuleSccs
