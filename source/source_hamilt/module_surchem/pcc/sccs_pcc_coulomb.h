#ifndef SCCS_PCC_COULOMB_H
#define SCCS_PCC_COULOMB_H

#include "../sccs/sccs_charge.h"
#include "sccs_pcc.h"
#include "../sccs/sccs_poisson.h"
#include "../sccs/sccs_pw_coulomb.h"

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSccs
{

MultipoleMoments reduced_pcc_density_moments(
    const std::vector<double>& density,
    const std::vector<ModuleBase::Vector3<double>>& positions,
    double volume_element,
    const PccGeometry& geometry,
    const ChargeReduction& reduction);

class PccCoulombOperator : public CoulombOperator
{
  public:
    PccCoulombOperator(const ModulePW::PW_Basis& basis,
                       double tpiba,
                       const std::vector<ModuleBase::Vector3<double>>& positions,
                       double volume_element,
                       const PccGeometry& geometry,
                       const ChargeReduction& reduction);

    void apply(const std::vector<double>& charge, ElectrostaticField& field) const override;

  private:
    PeriodicCoulombOperator periodic_;
    const std::vector<ModuleBase::Vector3<double>>& positions_;
    std::vector<ModuleBase::Vector3<double>> relative_positions_;
    double volume_element_ = 0.0;
    PccGeometry geometry_;
    const ChargeReduction& reduction_;
};

} // namespace ModuleSccs

#endif
