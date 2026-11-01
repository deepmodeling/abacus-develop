#ifndef SCCS_COUL_2D_H
#define SCCS_COUL_2D_H

#include "sccs_pw_coulomb.h"
#include "source_estate/pcc_2d.h"

namespace ModuleSccs
{
// Positions are minimum-image normal projections about a fixed origin supplied
// by the cell layer. Only grid-charge moments enter this operator; the host
// supplies the point-ion vacuum PCC energy separately.
class Pcc2dCoulombOperator : public CoulombOperator
{
public:
    Pcc2dCoulombOperator(const ModulePW::PW_Basis& basis,
                         double tpiba,
                         const std::vector<ModuleBase::Vector3<double>>& projected_positions,
                         const ModuleBase::Vector3<double>& normal,
                         const elecstate::Pcc2dParameters& parameters);
    bool has_correction() const override;
    CoulombTransformCounts fft_counts() const override { return periodic_.fft_counts(); }
    void apply_potential(const std::vector<double>& charge, std::vector<double>& potential) override;

private:
    const ModulePW::PW_Basis& basis_;
    PeriodicCoulombOperator periodic_;
    const std::vector<ModuleBase::Vector3<double>> projected_positions_;
    const ModuleBase::Vector3<double> normal_;
    const elecstate::Pcc2dParameters parameters_;
};
} // namespace ModuleSccs

#endif
