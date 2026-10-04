#include "sccs_pw_coulomb.h"

#include "source_base/constants.h"
#include "source_base/parallel_reduce.h"
#include "source_basis/module_pw/pw_basis.h"

#include <cmath>

namespace ModuleSccs
{
bool validate_pw_grid(const ModulePW::PW_Basis& basis, double tpiba, std::string& error)
{
    error.clear();
    double invalid = 0.0;
    if (!std::isfinite(tpiba) || tpiba <= 0.0 || !std::isfinite(basis.omega) || basis.omega <= 0.0
        || basis.gamma_only || basis.nxyz <= 0 || basis.nrxx < 0 || basis.npw < 0
        || basis.nmaxgr < basis.nrxx || basis.nmaxgr < basis.npw
        || (basis.npw > 0 && (basis.gg == nullptr || basis.gcar == nullptr)))
    {
        invalid = 1.0;
    }
    Parallel_Reduce::reduce_max_pool(basis.poolnproc, invalid);
    if (invalid != 0.0)
    {
        error = "SCCS requires an initialized full-G PW grid and positive finite tpiba and volume on every pool rank";
        return false;
    }
    return true;
}

bool validate_grid_values(const std::vector<double>& values,
                          const ModulePW::PW_Basis& basis,
                          std::string& error)
{
    error.clear();
    const std::size_t local_size = basis.nrxx;
    double invalid = 0.0;
    if (values.size() != local_size)
    {
        invalid = 1.0;
    }
    for (double value : values)
    {
        if (!std::isfinite(value))
        {
            invalid = 1.0;
        }
    }
    Parallel_Reduce::reduce_max_pool(basis.poolnproc, invalid);
    if (invalid != 0.0)
    {
        error = "SCCS requires finite values matching the local PW grid on every pool rank";
        return false;
    }
    return true;
}

PeriodicCoulombOperator::PeriodicCoulombOperator(const ModulePW::PW_Basis& basis, double tpiba)
    : basis_(basis), tpiba_(tpiba)
{
}

bool PeriodicCoulombOperator::apply_potential(const std::vector<double>& charge,
                                             std::vector<double>& potential,
                                             std::string& error)
{
    if (!validate_pw_grid(basis_, tpiba_, error) || !validate_grid_values(charge, basis_, error))
    {
        return false;
    }
    reciprocal_work_.resize(basis_.npw);
    basis_.real2recip(charge.data(), reciprocal_work_.data());
    const double tpiba2 = tpiba_ * tpiba_;
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
    std::vector<double> candidate(basis_.nrxx);
    basis_.recip2real(reciprocal_work_.data(), candidate.data());
    if (!validate_grid_values(candidate, basis_, error))
    {
        return false;
    }
    potential.swap(candidate);
    return true;
}
} // namespace ModuleSccs
