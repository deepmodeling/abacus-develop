#ifndef SCCS_PCC_COULOMB_H
#define SCCS_PCC_COULOMB_H

#include "sccs_pcc.h"
#include "../sccs/sccs_poisson.h"
#include "../sccs/sccs_pw_coulomb.h"

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSccs
{

class ChargeReduction;

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

    bool has_boundary_correction() const override
    {
        return true;
    }

    void apply(const std::vector<double>& charge, ElectrostaticField& field) const override;

    void apply_potential(const std::vector<double>& charge,
                         std::vector<double>& potential) const override;

    CoulombTransformProfile transform_profile() const override
    {
        return periodic_.transform_profile();
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
    PccGeometry geometry_;
    const ChargeReduction& reduction_;
};

} // namespace ModuleSccs

#endif
