#ifndef SCCS_PCC_2D_COULOMB_H
#define SCCS_PCC_2D_COULOMB_H

#include "../pcc/pcc_2d.h"
#include "sccs_poisson.h"
#include "sccs_pw_coulomb.h"

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSurchem
{
class ChargeReduction;
}

namespace ModuleSccs
{

class Pcc2dCoulombOperator : public CoulombOperator
{
  public:
    Pcc2dCoulombOperator(const ModulePW::PW_Basis& basis,
                         double tpiba,
                         const std::vector<ModuleBase::Vector3<double>>& positions,
                         double volume_element,
                         const ModulePcc::Pcc2dGeometry& geometry,
                         const ModuleSurchem::ChargeReduction& reduction);

    bool has_boundary_correction() const override
    {
        return true;
    }

    CoulombTransformCounts transform_counts() const override
    {
        return periodic_.transform_counts();
    }

    void apply_potential(const std::vector<double>& charge,
                         std::vector<double>& potential) const override;

  private:
    PeriodicCoulombOperator periodic_;
    std::vector<double> relative_y_;
    double volume_element_ = 0.0;
    ModulePcc::Pcc2dGeometry geometry_;
    const ModuleSurchem::ChargeReduction& reduction_;
};

} // namespace ModuleSccs

#endif
