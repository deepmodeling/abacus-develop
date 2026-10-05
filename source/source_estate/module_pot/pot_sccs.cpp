#include "pot_sccs.h"

#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_cell/cell_tools.h"
#include "source_hamilt/module_sccs/sccs_functional.h"
#include "source_hamilt/module_sccs/sccs_ionic_charge.h"
#include "source_io/module_parameter/input_parameter.h"

#include <cmath>
#include <utility>

namespace
{
void require_valid_on_pool(bool valid, const std::string& error)
{
    double invalid = valid ? 0.0 : 1.0;
    Parallel_Reduce::reduce_pool(invalid);
    if (invalid != 0.0)
    {
        const std::string message = error.empty() ? "Invalid SCCS input on another pool rank" : error;
        ModuleBase::WARNING_QUIT("PotSccs", message);
    }
}
}

namespace elecstate
{
bool make_sccs_config_from_input(const Input_para& input,
                                 ModuleSccs::SccsConfig& config,
                                 ModuleSccs::PolarizationSolverParameters& solver,
                                 std::string& error)
{
    ModuleSccs::Preset preset;
    if (!ModuleSccs::parse_preset(input.sccs_preset, preset, error)) { return false; }
    ModuleSccs::SccsConfig candidate;
    if (preset == ModuleSccs::Preset::Custom)
    {
        candidate.cavity.density_min = input.sccs_rho_min;
        candidate.cavity.density_max = input.sccs_rho_max;
        candidate.cavity.epsilon_bulk = input.sccs_epsilon;
        candidate.surface_tension = ModuleSccs::dyn_per_cm_to_hartree_per_bohr2(input.sccs_gamma);
        candidate.pressure = ModuleSccs::gpa_to_hartree_per_bohr3(input.sccs_pressure);
    }
    else if (!ModuleSccs::make_sccs_config(preset, candidate, error)) { return false; }
    candidate.surface_regularization = input.sccs_surface_eta;
    if (!ModuleSccs::validate_config(candidate, error)) { return false; }
    if (input.sccs_maxiter <= 0 || !std::isfinite(input.sccs_tol_rms) || input.sccs_tol_rms <= 0.0
        || !std::isfinite(input.sccs_tol_max) || input.sccs_tol_max <= 0.0)
    {
        error = "SCCS requires a positive iteration limit and finite positive residual tolerances";
        return false;
    }
    config = candidate;
    solver.max_iterations = input.sccs_maxiter;
    solver.tolerance_rms = input.sccs_tol_rms;
    solver.tolerance_max = input.sccs_tol_max;
    return true;
}

PotSccs::PotSccs(const ModulePW::PW_Basis* basis,
                 const ModuleSccs::SccsConfig& config,
                 const ModuleSccs::PolarizationSolverParameters& solver)
    : config_(config), solver_(solver)
{
    this->rho_basis_ = basis;
    this->dynamic_mode = true;
}

void PotSccs::cal_v_eff(const Charge* charge, const UnitCell* cell, ModuleBase::matrix& potential)
{
    ModuleBase::timer::start("PotSccs", "cal_v_eff");
    const bool storage_valid = charge != nullptr && cell != nullptr && this->rho_basis_ != nullptr;
    require_valid_on_pool(storage_valid, "SCCS requires charge, cell and PW basis storage");
    const ModulePW::PW_Basis& basis = *this->rho_basis_;
    const bool grid_valid = (charge->nspin == 1 || charge->nspin == 2)
                            && potential.nr == charge->nspin && potential.nc == basis.nrxx
                            && charge->rho != nullptr && cell->atoms != nullptr && cell->ntype > 0
                            && std::isfinite(cell->lat0) && cell->lat0 > 0.0;
    require_valid_on_pool(grid_valid, "SCCS requires initialized atom/density storage and an nspin=1/2 potential");
    bool density_valid = true;
    for (int spin = 0; spin < charge->nspin; ++spin)
    {
        if (basis.nrxx > 0 && charge->rho[spin] == nullptr) { density_valid = false; }
    }
    require_valid_on_pool(density_valid, "SCCS charge density is not available");
    const std::vector<unitcell::AtomData> atoms = unitcell::get_atom_data(cell->atoms, cell->ntype, cell->lat0);
    const int atom_count = atoms.size();
    const bool count_valid = atom_count == cell->nat;
    require_valid_on_pool(count_valid, "SCCS atom count does not match UnitCell");
    std::vector<double> ions;
    std::string error;
    const bool ions_valid = ModuleSccs::gaussian_ionic_density(atoms, basis, cell->tpiba,
                                                              ModuleSccs::gaussian_ion_spread, ions, error);
    require_valid_on_pool(ions_valid, error);
    std::vector<double> density(basis.nrxx, 0.0);
    std::vector<double> solute_charge(basis.nrxx);
    double ionic_sum = 0.0;
    double net_charge = 0.0;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        for (int spin = 0; spin < charge->nspin; ++spin) { density[ir] += charge->rho[spin][ir]; }
        solute_charge[ir] = ions[ir] - density[ir];
        ionic_sum += ions[ir];
        net_charge += solute_charge[ir];
    }
    Parallel_Reduce::reduce_pool(ionic_sum);
    Parallel_Reduce::reduce_pool(net_charge);
    const double dv = basis.omega / basis.nxyz;
    ionic_sum *= dv;
    net_charge *= dv;
    double expected_ionic_charge = 0.0;
    for (const unitcell::AtomData& atom : atoms) { expected_ionic_charge += atom.valence_charge; }
    const double normalization_error = ionic_sum - expected_ionic_charge;
    const bool normalization_valid = std::isfinite(normalization_error) && std::abs(normalization_error) < 1e-6;
    require_valid_on_pool(normalization_valid, "SCCS Gaussian ionic charge normalization failed");
    const bool neutral = std::isfinite(net_charge) && std::abs(net_charge) < 1e-6;
    require_valid_on_pool(neutral, "Periodic SCCS currently requires a neutral cell");
    ModuleSccs::SccsResponse response;
    const bool response_valid = ModuleSccs::solve_sccs_response(density, solute_charge, config_.cavity,
                                                               solver_, restart_potential_, basis, cell->tpiba,
                                                               response, error);
    require_valid_on_pool(response_valid, error);
    ModuleSccs::FunctionalResult functional;
    const bool functional_valid = ModuleSccs::evaluate_functional(solute_charge, response, config_, basis,
                                                                  cell->tpiba, functional, error);
    require_valid_on_pool(functional_valid, error);
    electrostatic_rydberg_ = 2.0 * functional.reaction_energy;
    non_electrostatic_rydberg_ = 2.0 * (functional.surface_energy + functional.volume_energy);
    electrostatic_potential_.resize(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        electrostatic_potential_[ir] = -2.0 * functional.reaction_potential[ir];
        const double value = 2.0 * functional.electron_potential[ir];
        for (int spin = 0; spin < charge->nspin; ++spin) { potential(spin, ir) += value; }
    }
    restart_potential_ = std::move(response.restart_potential);
    ModuleBase::timer::end("PotSccs", "cal_v_eff");
}

double PotSccs::get_energy() const
{
    return electrostatic_rydberg_ + non_electrostatic_rydberg_;
}

void PotSccs::get_solvation_energy(double& electrostatic, double& non_electrostatic) const
{
    electrostatic = electrostatic_rydberg_;
    non_electrostatic = non_electrostatic_rydberg_;
}

const std::vector<double>* PotSccs::solvent_electrostatic_potential() const
{
    return &electrostatic_potential_;
}
} // namespace elecstate
