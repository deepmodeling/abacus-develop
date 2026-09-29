#ifndef SCCS_PCC_COULOMB_H
#define SCCS_PCC_COULOMB_H

#include "../pcc/pcc_0d.h"
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

class PccCoulombOperator : public CoulombOperator
{
  public:
    PccCoulombOperator(const ModulePW::PW_Basis& basis,
                       double tpiba,
                       const std::vector<ModuleBase::Vector3<double>>& positions,
                       double volume_element,
                       const ModulePcc::PccGeometry& geometry,
                       const ModuleSurchem::ChargeReduction& reduction);

    bool has_boundary_correction() const override
    {
        return true;
    }

    void apply(const std::vector<double>& charge, ElectrostaticField& field) const override;

    void apply_potential(const std::vector<double>& charge,
                         std::vector<double>& potential) const override;

    CoulombTransformCounts transform_counts() const override
    {
        return periodic_.transform_counts();
    }

    void apply_gradient(const std::vector<double>& charge,
                        std::vector<ModuleBase::Vector3<double>>& gradient) const override;

    // Adjoint of charge -> electrostatic field gradient under the grid inner product.
    void apply_gradient_adjoint(
        const std::vector<ModuleBase::Vector3<double>>& field,
        std::vector<double>& result) const override;


  private:
    PeriodicCoulombOperator periodic_;
    const std::vector<ModuleBase::Vector3<double>>& positions_;
    std::vector<ModuleBase::Vector3<double>> relative_positions_;
    double volume_element_ = 0.0;
    ModulePcc::PccGeometry geometry_;
    const ModuleSurchem::ChargeReduction& reduction_;
};

} // namespace ModuleSccs

#endif
