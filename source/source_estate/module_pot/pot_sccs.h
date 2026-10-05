#ifndef POT_SCCS_H
#define POT_SCCS_H

#include "pot_base.h"
#include "source_hamilt/module_sccs/sccs_parameters.h"
#include "source_hamilt/module_sccs/sccs_response.h"

struct Input_para;
namespace elecstate
{
bool make_sccs_config_from_input(const Input_para& input,
                                 ModuleSccs::SccsConfig& config,
                                 ModuleSccs::PolarizationSolverParameters& solver,
                                 std::string& error);

class PotSccs : public PotBase
{
public:
    PotSccs(const ModulePW::PW_Basis* basis,
            const ModuleSccs::SccsConfig& config,
            const ModuleSccs::PolarizationSolverParameters& solver);
    void cal_v_eff(const Charge* charge, const UnitCell* cell, ModuleBase::matrix& potential) override;
    double get_energy() const override;
    void add_solvation_force(const UnitCell& cell, ModuleBase::matrix& force) const override;
    void get_solvation_energy(double& electrostatic, double& non_electrostatic) const override;
    const std::vector<double>* solvent_electrostatic_potential() const override;

private:
    const ModuleSccs::SccsConfig config_;
    const ModuleSccs::PolarizationSolverParameters solver_;
    double electrostatic_rydberg_ = 0.0;
    double non_electrostatic_rydberg_ = 0.0;
    std::vector<double> electrostatic_potential_;
    std::vector<double> restart_potential_;
};
} // namespace elecstate
#endif
