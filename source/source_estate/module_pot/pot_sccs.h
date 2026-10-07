#ifndef POT_SCCS_H
#define POT_SCCS_H

#include "pot_base.h"
#include "source_hamilt/module_sccs/sccs_parameters.h"

#include <string>

struct Input_para;
namespace elecstate
{
// Map INPUT values, already validated by the INPUT reader, onto the solver.
void make_sccs_config_from_input(const Input_para& input,
                                 ModuleSccs::SccsConfig& config,
                                 ModuleSccs::PolarizationSolverParameters& solver);

// Returns a warning (empty if none) for a charged cell in a dielectric: the
// periodic Poisson solver drops the G = 0 component of the net charge, so the
// energy depends on the cell size. electron_count is INPUT nelec.
std::string check_sccs_charge(const ModuleSccs::SccsConfig& config, const UnitCell& cell, double electron_count);

class PotSccs : public PotBase
{
public:
    // electron_count is the expected number of electrons (INPUT nelec); a grid
    // density that differs from it by more than 1e-6 e is reported in the
    // warning log, and the solute charge uses the grid density.
    PotSccs(const ModulePW::PW_Basis* basis,
            const ModuleSccs::SccsConfig& config,
            const ModuleSccs::PolarizationSolverParameters& solver,
            double electron_count);
    void cal_v_eff(const Charge* charge, const UnitCell* cell, ModuleBase::matrix& potential) override;
    double get_energy() const override;
    void add_solvation_force(const UnitCell& cell, ModuleBase::matrix& force) const override;
    void get_solvation_energy(double& electrostatic, double& non_electrostatic) const override;
    const std::vector<double>* solvent_electrostatic_potential() const override;

private:
    const ModuleSccs::SccsConfig config_;
    const ModuleSccs::PolarizationSolverParameters solver_;
    const double electron_count_;
    double electrostatic_rydberg_ = 0.0;
    double non_electrostatic_rydberg_ = 0.0;
    std::vector<double> electrostatic_potential_;
    std::vector<double> restart_potential_;
};
} // namespace elecstate
#endif
