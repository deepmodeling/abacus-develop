#include "pot_sccs.h"

#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_cell/cell_tools.h"
#include "source_hamilt/module_sccs/sccs_functional.h"
#include "source_hamilt/module_sccs/sccs_ionic_charge.h"
#include "source_hamilt/module_sccs/sccs_ionic_force.h"
#include "source_hamilt/module_sccs/sccs_response.h"
#include "source_io/module_parameter/input_parameter.h"

#include <cmath>
#include <sstream>
#include <utility>

namespace elecstate
{
void make_sccs_config_from_input(const Input_para& input,
                                 ModuleSccs::SccsConfig& config,
                                 ModuleSccs::PolarizationSolverParameters& solver)
{
    const ModuleSccs::Preset preset = ModuleSccs::parse_preset(input.sccs_preset);
    if (preset == ModuleSccs::Preset::Custom)
    {
        config = ModuleSccs::SccsConfig();
        config.cavity.density_min = input.sccs_rho_min;
        config.cavity.density_max = input.sccs_rho_max;
        config.cavity.epsilon_bulk = input.sccs_epsilon;
        config.surface_tension = ModuleSccs::dyn_per_cm_to_hartree_per_bohr2(input.sccs_gamma);
        config.pressure = ModuleSccs::gpa_to_hartree_per_bohr3(input.sccs_pressure);
    }
    else
    {
        config = ModuleSccs::make_sccs_config(preset);
    }
    config.surface_regularization = input.sccs_surface_eta;
    solver.max_iterations = input.sccs_maxiter;
    solver.tolerance_rms = input.sccs_tol_rms;
    solver.tolerance_max = input.sccs_tol_max;
}

std::string check_sccs_charge(const ModuleSccs::SccsConfig& config, const UnitCell& cell, double electron_count)
{
    double ionic_charge = 0.0;
    for (int it = 0; it < cell.ntype; ++it)
    {
        ionic_charge += cell.atoms[it].ncpp.zv * cell.atoms[it].na;
    }
    const double net_charge = ionic_charge - electron_count;
    if (config.cavity.epsilon_bulk <= 1.0 || std::abs(net_charge) <= 1e-6) { return std::string(); }
    return "charged SCCS with periodic boundaries: the periodic Poisson solver drops the G = 0 component "
           "of the net charge, so the energy depends on the cell size.";
}

PotSccs::PotSccs(const ModulePW::PW_Basis* basis,
                 const ModuleSccs::SccsConfig& config,
                 const ModuleSccs::PolarizationSolverParameters& solver,
                 double electron_count)
    : config_(config), solver_(solver), electron_count_(electron_count)
{
    this->rho_basis_ = basis;
    this->dynamic_mode = true;
}

void PotSccs::cal_v_eff(const Charge* charge, const UnitCell* cell, ModuleBase::matrix& potential)
{
    ModuleBase::timer::start("PotSccs", "cal_v_eff");
    const ModulePW::PW_Basis& basis = *this->rho_basis_;
    const std::vector<unitcell::AtomData> atoms = unitcell::get_atom_data(cell->atoms, cell->ntype, cell->lat0);
    std::vector<double> ions;
    ModuleSccs::gaussian_ionic_density(atoms, basis, cell->tpiba, ModuleSccs::gaussian_ion_spread, ions);
    std::vector<double> density(basis.nrxx, 0.0);
    std::vector<double> solute_charge(basis.nrxx);
    double grid_charges[2] = {0.0, 0.0}; // electrons, ions
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        for (int spin = 0; spin < charge->nspin; ++spin) { density[ir] += charge->rho[spin][ir]; }
        solute_charge[ir] = ions[ir] - density[ir];
        grid_charges[0] += density[ir];
        grid_charges[1] += ions[ir];
    }
    // Pool-reduced, so every rank takes the same decision.
    Parallel_Reduce::reduce_pool(grid_charges, 2);
    const double dv = basis.omega / basis.nxyz;
    const double electron_error = grid_charges[0] * dv - electron_count_;
    double valence_charge = 0.0;
    for (const unitcell::AtomData& atom : atoms) { valence_charge += atom.valence_charge; }
    const double ionic_error = grid_charges[1] * dv - valence_charge;
    if (std::abs(ionic_error) > 1e-6)
    {
        ModuleBase::WARNING_QUIT("PotSccs::cal_v_eff", "SCCS ionic density normalization does not match the valence charge");
    }
    if (std::abs(electron_error) > 1e-6)
    {
        std::ostringstream message;
        message << "SCCS grid electron count differs from the expected value by " << electron_error
                << " e; the solute charge uses the grid density";
        ModuleBase::WARNING("PotSccs::cal_v_eff", message.str());
    }
    ModuleSccs::SccsResponse response;
    ModuleSccs::solve_sccs_response(density, solute_charge, config_.cavity, solver_, restart_potential_, basis,
                                    cell->tpiba, response);
    ModuleSccs::FunctionalResult functional;
    ModuleSccs::evaluate_functional(solute_charge, response, config_, basis, cell->tpiba, functional);
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

void PotSccs::add_solvation_force(const UnitCell& cell, ModuleBase::matrix& force) const
{
    ModuleBase::timer::start("PotSccs", "add_solvation_force");
    const ModulePW::PW_Basis& basis = *this->rho_basis_;
    const std::vector<unitcell::AtomData> atoms = unitcell::get_atom_data(cell.atoms, cell.ntype, cell.lat0);
    std::vector<double> reaction(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        reaction[ir] = -0.5 * electrostatic_potential_[ir];
    }
    std::vector<ModuleBase::Vector3<double>> ionic_force;
    ModuleSccs::gaussian_ionic_force(atoms, reaction, basis, cell.tpiba, ModuleSccs::gaussian_ion_spread,
                                     ionic_force);
    for (int ia = 0; ia < cell.nat; ++ia)
    {
        for (int axis = 0; axis < 3; ++axis)
        {
            force(ia, axis) += 2.0 * ionic_force[ia][axis];
        }
    }
    ModuleBase::timer::end("PotSccs", "add_solvation_force");
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
