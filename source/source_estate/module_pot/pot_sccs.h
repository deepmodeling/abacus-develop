#ifndef POT_SCCS_H
#define POT_SCCS_H

#include "pot_base.h"
#include "sccs_output.h"
#include "sccs_fixed_src.h"
#include "source_hamilt/module_sccs/sccs_parameters.h"

#include <iosfwd>
#include <string>

struct Input_para;
namespace elecstate
{
// Map INPUT values, already validated by the INPUT reader, onto the solver.
void sccs_cfg_input(const Input_para& input,
                    ModuleSccs::SccsConfig& config,
                    ModuleSccs::PolarizationSolverParameters& solver);

// Structure-dependent checks that the INPUT reader cannot make: a
// sccs_corespread count other than one or nat stops the run. Returns a warning
// (empty if none) when a single width meets an unknown pseudopotential
// element, whose automatic core exclusion is then not possible.
std::string check_sccs_cell(const ModuleSccs::SccsConfig& config, const UnitCell& cell);

// Returns a warning (empty if none) for a charged cell in a dielectric: the
// periodic Poisson solver drops the G = 0 component of the net charge, so the
// energy depends on the cell size. electron_count is the resolved nelec (the
// INPUT value, or the valence count when it is not set).
std::string check_charge(const ModuleSccs::SccsConfig& config, const UnitCell& cell, double electron_count);

// State carried from the component of the previous ionic step, which is
// rebuilt every ionic step: an activation already reached, the last
// unshifted sqrt-CG solution as a warm start, and the solvent-aware probe
// kernel of the fixed cell and grid. The solver keeps the warm start only when
// it lowers the initial charge residual, so moved atoms are safe.
struct SccsResume
{
    bool active = false;
    std::vector<double> restart_potential;
    bool probe_ready = false; // same on every pool rank
    std::vector<double> probe_kernel; // local G vectors
};

class PotSccs : public PotBase
{
public:
    // electron_count is the expected number of electrons (resolved nelec); a grid
    // density that differs from it by more than 1e-6 e is reported in the
    // warning log, and the solute charge uses the grid density.
    // resume carries the activation and warm start of the previous ionic step.
    PotSccs(const ModulePW::PW_Basis* basis,
            const ModuleSccs::SccsConfig& config,
            const ModuleSccs::PolarizationSolverParameters& solver,
            double electron_count,
            const SccsResume& resume);
    void cal_v_eff(const Charge* charge, const UnitCell* cell, ModuleBase::matrix& potential) override;
    double get_energy() const override;
    void add_sol_force(const UnitCell& cell, ModuleBase::matrix& force) const override;
    void sol_energy(double& electrostatic, double& non_electrostatic) const override;
    const std::vector<double>* sol_elec_pot() const override;
    // eps and the cavity s of the last active evaluation, plus the local s
    // (cavity_local) and the probe solute fraction c (filled_fraction) of a
    // solvent-aware cavity; none before the start.
    void add_sol_fields(std::vector<SolventGridField>& fields) const override;

    // Delayed start: called once per electronic iteration while SCCS waits,
    // with the density residual broadcast to every rank. Returns true only in
    // the iteration that activates SCCS; activation is irreversible.
    bool try_activate(int electronic_iteration, double density_residual);
    bool is_active() const { return sccs_active_; }
    SccsResume resume_state() const;
    bool has_output() const { return output_.valid; }
    void write_iter(std::ostream& output, int level, double residual, double pcc_energy) const;

private:
    bool sccs_active_ = false;
    const ModuleSccs::SccsConfig config_;
    const ModuleSccs::PolarizationSolverParameters solver_;
    const double electron_count_;
    SccsOutput output_;
    SccsFixedSources fixed_sources_;
    double electrostatic_rydberg_ = 0.0;
    double non_electrostatic_rydberg_ = 0.0;
    std::vector<double> electrostatic_potential_;
    std::vector<double> restart_potential_;
    std::vector<double> cavity_potential_;
    std::vector<double> epsilon_;
    std::vector<double> solute_;
    std::vector<double> local_solute_;
    std::vector<double> filled_fraction_;
};
} // namespace elecstate
#endif
