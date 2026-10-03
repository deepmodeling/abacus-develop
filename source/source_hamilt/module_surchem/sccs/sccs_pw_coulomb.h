#ifndef SCCS_PW_COULOMB_H
#define SCCS_PW_COULOMB_H

#include "sccs_poisson.h"

#include <complex>

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSccs
{

std::vector<ModuleBase::Vector3<double>> periodic_gradient(
    const std::vector<double>& values,
    const ModulePW::PW_Basis& basis,
    double tpiba);

class PeriodicCoulombOperator : public CoulombOperator
{
  public:
    PeriodicCoulombOperator(const ModulePW::PW_Basis& basis, double tpiba);

    CoulombTransformCounts transform_counts() const override { return counts_; }

    void apply_potential(const std::vector<double>& charge,
                         std::vector<double>& potential) const override;

  private:
    const ModulePW::PW_Basis& basis_;
    double tpiba_ = 0.0;

    // Scratch only: every application overwrites it before use. Like the
    // underlying PW_Basis FFT workspace, this operator is not reentrant.
    mutable std::vector<std::complex<double>> reciprocal_work_;
    mutable CoulombTransformCounts counts_;
};

} // namespace ModuleSccs

#endif
