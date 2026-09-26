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

std::vector<double> periodic_negative_divergence(
    const std::vector<ModuleBase::Vector3<double>>& field,
    const ModulePW::PW_Basis& basis,
    double tpiba);

class PeriodicCoulombOperator : public CoulombOperator
{
  public:
    PeriodicCoulombOperator(const ModulePW::PW_Basis& basis, double tpiba);

    void apply(const std::vector<double>& charge, ElectrostaticField& field) const override;

    // Adjoint of charge -> electrostatic field gradient under the grid inner product.
    void apply_gradient_adjoint(
        const std::vector<ModuleBase::Vector3<double>>& field,
        std::vector<double>& result) const override;


  private:
    const ModulePW::PW_Basis& basis_;
    double tpiba_ = 0.0;

    // Scratch only: every transform overwrites its inputs before use. Like
    // the underlying PW_Basis FFT workspace, this operator is not reentrant.
    // Sharing buffers between forward and adjoint applications bounds memory.
    mutable std::vector<std::complex<double>> reciprocal_work_;
    mutable std::vector<std::complex<double>> reciprocal_aux_;
    mutable std::vector<double> real_work_;
};

} // namespace ModuleSccs

#endif
