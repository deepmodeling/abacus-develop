#ifndef SCCS_PCC_2D_COULOMB_H
#define SCCS_PCC_2D_COULOMB_H

#include "sccs_pcc_2d.h"
#include "../sccs/sccs_poisson.h"
#include "../sccs/sccs_pw_coulomb.h"

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSccs
{

class ChargeReduction;

Pcc2dMoments reduced_pcc_2d_density_moments(
    const std::vector<double>& density,
    const std::vector<ModuleBase::Vector3<double>>& positions,
    double volume_element,
    const Pcc2dGeometry& geometry,
    const ChargeReduction& reduction);

std::vector<double> pcc_2d_plane_average(const std::vector<double>& values,
                                         const ModulePW::PW_Basis& basis,
                                         const ChargeReduction& reduction);

class Pcc2dCoulombOperator : public CoulombOperator
{
  public:
    Pcc2dCoulombOperator(const ModulePW::PW_Basis& basis,
                         double tpiba,
                         const std::vector<ModuleBase::Vector3<double>>& positions,
                         double volume_element,
                         const Pcc2dGeometry& geometry,
                         const ChargeReduction& reduction);

    void apply(const std::vector<double>& charge, ElectrostaticField& field) const override;

  private:
    PeriodicCoulombOperator periodic_;
    const std::vector<ModuleBase::Vector3<double>>& positions_;
    std::vector<double> relative_y_;
    double volume_element_ = 0.0;
    Pcc2dGeometry geometry_;
    const ChargeReduction& reduction_;
};

} // namespace ModuleSccs

#endif
