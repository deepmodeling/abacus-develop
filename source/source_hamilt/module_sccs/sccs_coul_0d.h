#ifndef SCCS_COUL_0D_H
#define SCCS_COUL_0D_H

#include "sccs_pw_coulomb.h"
#include "source_estate/pcc_0d.h"

namespace ModuleSccs
{
// Grid positions are minimum-image displacements about a common, fixed origin
// supplied by the cell layer. Only grid-charge moments enter this operator;
// the host's point-ion PCC energy must not be added here.
class Pcc0dCoulombOperator : public CoulombOperator
{
public:
    Pcc0dCoulombOperator(const ModulePW::PW_Basis& basis,
                         double tpiba,
                         const std::vector<ModuleBase::Vector3<double>>& relative_positions,
                         const elecstate::Pcc0dParameters& parameters);
    bool has_correction() const override;
    CoulombTransformCounts fft_counts() const override { return periodic_.fft_counts(); }
    void apply_potential(const std::vector<double>& charge, std::vector<double>& potential) override;

private:
    const ModulePW::PW_Basis& basis_;
    PeriodicCoulombOperator periodic_;
    const std::vector<ModuleBase::Vector3<double>> relative_positions_;
    const elecstate::Pcc0dParameters parameters_;
};
} // namespace ModuleSccs

#endif
