#include "sccs_pw_coulomb.h"

#include "source_base/constants.h"
#include "source_basis/module_pw/pw_basis.h"

namespace ModuleSccs
{
PeriodicCoulombOperator::PeriodicCoulombOperator(const ModulePW::PW_Basis& basis, double tpiba)
    : basis_(basis), tpiba_(tpiba)
{
}

bool PeriodicCoulombOperator::has_correction() const
{
    return false;
}

void PeriodicCoulombOperator::apply_potential(const std::vector<double>& charge, std::vector<double>& potential)
{
    reciprocal_work_.resize(basis_.npw);
    basis_.real2recip(charge.data(), reciprocal_work_.data());
    ++counts_.forward_calls;
    const double tpiba2 = tpiba_ * tpiba_;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ig = 0; ig < basis_.npw; ++ig)
    {
        if (ig == basis_.ig_gge0 || basis_.gg[ig] == 0.0)
        {
            reciprocal_work_[ig] = 0.0;
        }
        else
        {
            reciprocal_work_[ig] = ModuleBase::FOUR_PI * reciprocal_work_[ig] / (tpiba2 * basis_.gg[ig]);
        }
    }
    potential.resize(basis_.nrxx);
    basis_.recip2real(reciprocal_work_.data(), potential.data());
    ++counts_.inverse_calls;
}
} // namespace ModuleSccs
