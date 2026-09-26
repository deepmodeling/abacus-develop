#ifndef SCCS_PW_COULOMB_H
#define SCCS_PW_COULOMB_H

#include "sccs_poisson.h"

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

    void apply(const std::vector<double>& charge, ElectrostaticField& field) const override;

  private:
    const ModulePW::PW_Basis& basis_;
    double tpiba_ = 0.0;
};

} // namespace ModuleSccs

#endif
