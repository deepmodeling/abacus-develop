#ifndef POT_PCC_H
#define POT_PCC_H

#include "pot_base.h"
#include "source_cell/cell_geometry.h"
#include "source_estate/pcc_0d.h"

namespace elecstate
{

/// Self-consistent vacuum PCC. Owns the result of its last density update.
class PotPcc : public PotBase
{
  public:
    explicit PotPcc(const ModulePW::PW_Basis* basis);

    void cal_v_eff(const Charge* charge, const UnitCell* cell, ModuleBase::matrix& potential) override;
    double get_energy() const override;
    const std::vector<double>& electron_potential() const;
    void add_force(const UnitCell& cell, ModuleBase::matrix& force) const;

  private:
    void prepare_ions(const UnitCell& cell, ChargeMoments& ionic_moments);
    ChargeMoments collect_electrons(const Charge& charge,
                                    const UnitCell& cell,
                                    std::vector<ModuleBase::Vector3<double>>& positions) const;

    unitcell::OrthogonalCell geometry_;
    Pcc0dParameters parameters_;
    ChargeMoments moments_;
    std::vector<ModuleBase::Vector3<double>> ionic_positions_;
    std::vector<double> ionic_charges_;
    std::vector<double> electron_potential_;
    double energy_rydberg_ = 0.0;
    bool result_valid_ = false;
};

} // namespace elecstate

#endif
