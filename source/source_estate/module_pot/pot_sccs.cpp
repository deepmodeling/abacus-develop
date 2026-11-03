#include "pot_sccs.h"
#include "sol_grid_field.h"
#include "sccs_pcc_0d.h"
#include "sccs_pcc_2d.h"

#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_cell/cell_tools.h"
#include "source_hamilt/module_sccs/sccs_functional.h"
#include "source_hamilt/module_sccs/sccs_ion_charge.h"
#include "source_hamilt/module_sccs/sccs_ion_force.h"
#include "source_hamilt/module_sccs/sccs_response.h"
#include "source_hamilt/module_sccs/sccs_pw_coulomb.h"
#include "source_io/module_parameter/input_parameter.h"

#include <cmath>
#include <sstream>
#include <utility>

namespace elecstate
{
namespace
{
// Charge tolerance in e, the same as the PotPcc electron-count check: far
// above the grid round-off (~1e-12 e) and far below any physical net charge.
const double charge_tolerance = 1e-6;
} // namespace

void sccs_cfg_input(const Input_para& input,
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
        config.surface_tension = ModuleSccs::dyn_cm_to_au(input.sccs_gamma);
        config.pressure = ModuleSccs::gpa_to_au(input.sccs_pressure);
    }
    else
    {
        config = ModuleSccs::preset_config(preset);
    }
    config.start_drho = input.sccs_start_drho;
    config.start_nmax = input.sccs_start_nmax;
    config.collect_report = input.sccs_debug >= 2;
    config.core_electrons = input.sccs_solvent_mode == "full";
    config.core_spreads = input.sccs_corespread;
    config.surface_regularization = input.sccs_surface_eta;
    config.cavity.lowpass_p1 = input.sccs_lowpass_p1;
    config.cavity.lowpass_p2 = input.sccs_lowpass_p2;
    if (input.assume_isolated == "pcc_0d")
    {
        config.boundary = ModuleSccs::Boundary::Pcc0d;
    }
    else if (input.assume_isolated == "pcc_2d")
    {
        config.boundary = ModuleSccs::Boundary::Pcc2d;
        config.pcc_2d_axis = input.pcc_2d_axis;
    }
    else
    {
        config.boundary = ModuleSccs::Boundary::Periodic;
    }
    solver.max_iterations = input.sccs_maxiter;
    solver.tolerance_rms = input.sccs_tol_rms;
    solver.tolerance_max = input.sccs_tol_max;
}

std::string check_sccs_cell(const ModuleSccs::SccsConfig& config, const UnitCell& cell)
{
    if (!config.core_electrons)
    {
        return std::string();
    }
    const std::size_t width_count = config.core_spreads.size();
    const std::size_t atom_count = cell.nat;
    if (width_count != 1 && width_count != atom_count)
    {
        ModuleBase::WARNING_QUIT("check_sccs_cell", "sccs_corespread must contain one value or exactly nat values");
    }
    if (width_count != 1)
    {
        return std::string();
    }
    std::ostringstream message;
    for (int it = 0; it < cell.ntype; ++it)
    {
        const int atomic_number = unitcell::pseudo_atom_z(cell.atoms[it]);
        if (atomic_number != 0)
        {
            continue;
        }
        message << "sccs_solvent_mode full: the element of atom type " << cell.atoms[it].label
                << " (pseudopotential element '" << cell.atoms[it].ncpp.psd
                << "') is unknown, so the single sccs_corespread value is applied to it even if its"
                << " pseudopotential has no core electrons; give one value per atom to control it.\n";
    }
    return message.str();
}

std::string check_charge(const ModuleSccs::SccsConfig& config, const UnitCell& cell, double electron_count)
{
    double ionic_charge = 0.0;
    for (int it = 0; it < cell.ntype; ++it)
    {
        ionic_charge += cell.atoms[it].ncpp.zv * cell.atoms[it].na;
    }
    const double net_charge = ionic_charge - electron_count;
    const bool periodic = config.boundary == ModuleSccs::Boundary::Periodic;
    if (!periodic || config.cavity.epsilon_bulk <= 1.0 || std::abs(net_charge) <= charge_tolerance)
    {
        return std::string();
    }
    return "charged SCCS with periodic boundaries: the periodic Poisson solver drops the G = 0 component "
           "of the net charge, so the energy depends on the cell size. Use assume_isolated pcc_0d or pcc_2d "
           "for converged charged energies.";
}

PotSccs::PotSccs(const ModulePW::PW_Basis* basis,
                 const ModuleSccs::SccsConfig& config,
                 const ModuleSccs::PolarizationSolverParameters& solver,
                 double electron_count,
                 const SccsResume& resume)
    : config_(config), solver_(solver), electron_count_(electron_count)
{
    this->rho_basis_ = basis;
    this->dynamic_mode = true;
    sccs_active_ = config_.start_drho <= 0.0 || resume.active;
    restart_potential_ = resume.restart_potential;
}

SccsResume PotSccs::resume_state() const
{
    SccsResume resume;
    resume.active = sccs_active_;
    resume.restart_potential = restart_potential_;
    return resume;
}

bool PotSccs::try_activate(int electronic_iteration, double density_residual)
{
    if (sccs_active_)
    {
        return false;
    }
    if (density_residual > config_.start_drho && electronic_iteration < config_.start_nmax)
    {
        return false;
    }
    sccs_active_ = true;
    restart_potential_.clear();
    return true;
}

void PotSccs::cal_v_eff(const Charge* charge, const UnitCell* cell, ModuleBase::matrix& potential)
{
    ModuleBase::timer::start("PotSccs", "cal_v_eff");
    const ModulePW::PW_Basis& basis = *this->rho_basis_;
    if (!sccs_active_)
    {
        electrostatic_rydberg_ = 0.0;
        non_electrostatic_rydberg_ = 0.0;
        electrostatic_potential_.assign(basis.nrxx, 0.0);
        cavity_potential_.clear();
        ModuleBase::timer::end("PotSccs", "cal_v_eff");
        return;
    }
    const std::vector<unitcell::AtomData> atoms = unitcell::get_atom_data(cell->atoms, cell->ntype, cell->lat0);
    // The ion and core densities are rebuilt only when the atoms changed.
    const bool reused_fixed_sources = fixed_sources_.update(*cell, basis, atoms, config_);
    output_.reused_fixed_sources = reused_fixed_sources;
    const std::vector<double>& ions = fixed_sources_.ionic_density();
    const std::vector<ModuleBase::Vector3<double>>& grid_positions = fixed_sources_.positions();
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
    if (std::abs(ionic_error) > charge_tolerance)
    {
        ModuleBase::WARNING_QUIT("PotSccs::cal_v_eff", "SCCS ionic density normalization does not match the valence charge");
    }
    // With PCC the PCC component reports the same mismatch.
    if (config_.boundary == ModuleSccs::Boundary::Periodic && std::abs(electron_error) > charge_tolerance)
    {
        std::ostringstream message;
        message << "SCCS grid electron count differs from the expected value by " << electron_error
                << " e; the solute charge uses the grid density";
        ModuleBase::WARNING("PotSccs::cal_v_eff", message.str());
    }
    std::unique_ptr<ModuleSccs::CoulombOperator> coulomb;
    if (config_.boundary == ModuleSccs::Boundary::Pcc0d)
    {
        pcc_0d_operator(*cell, basis, atoms, grid_positions, coulomb);
    }
    else if (config_.boundary == ModuleSccs::Boundary::Pcc2d)
    {
        pcc_2d_operator(*cell, basis, atoms, config_.pcc_2d_axis, grid_positions, coulomb);
    }
    else
    {
        coulomb.reset(new ModuleSccs::PeriodicCoulombOperator(basis, cell->tpiba));
    }
    if (config_.core_electrons)
    {
        const std::vector<double>& core_density = fixed_sources_.core_density();
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            density[ir] += core_density[ir];
        }
    }
    ModuleSccs::SccsResponse response;
    ModuleSccs::solve_response(density, solute_charge, config_.cavity, solver_, restart_potential_, basis,
                               cell->tpiba, *coulomb, response);
    // Only the Coulomb-operator FFTs of the sqrt-CG solve, its main cost; the
    // functional below and the cavity derivatives transform outside this count.
    output_.transforms = coulomb->fft_counts();
    ModuleSccs::FunctionalResult functional;
    ModuleSccs::eval_functional(solute_charge, response, config_, basis, cell->tpiba, *coulomb, functional);
    output_.valid = true;
    output_.iterations = response.polarization.iterations;
    output_.warm_started = response.polarization.warm_started;
    output_.residual_rms = response.polarization.residual_rms;
    output_.residual_max = response.polarization.residual_max;
    output_.reaction_energy = functional.reaction_energy;
    output_.volume = functional.volume;
    output_.surface = functional.surface;
    output_.far_field_charge = response.far_field_polarization_charge;
    if (config_.collect_report && config_.boundary != ModuleSccs::Boundary::Periodic)
    {
        collect_output(*cell, basis, atoms, grid_positions, ions, solute_charge, config_, response, output_);
    }
    electrostatic_rydberg_ = 2.0 * functional.reaction_energy;
    non_electrostatic_rydberg_ = 2.0 * (functional.surface_energy + functional.volume_energy);
    electrostatic_potential_.resize(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        electrostatic_potential_[ir] = -2.0 * functional.reaction_potential[ir];
        const double value = 2.0 * functional.electron_potential[ir];
        for (int spin = 0; spin < charge->nspin; ++spin) { potential(spin, ir) += value; }
    }
    cavity_potential_ = std::move(functional.cavity_potential);
    epsilon_ = std::move(response.epsilon);
    solute_ = std::move(response.solute);
    restart_potential_ = std::move(response.restart_potential);
    ModuleBase::timer::end("PotSccs", "cal_v_eff");
}

void PotSccs::add_sol_force(const UnitCell& cell, ModuleBase::matrix& force) const
{
    ModuleBase::timer::start("PotSccs", "add_sol_force");
    if (!sccs_active_)
    {
        ModuleBase::timer::end("PotSccs", "add_sol_force");
        return;
    }
    const ModulePW::PW_Basis& basis = *this->rho_basis_;
    const std::vector<unitcell::AtomData> atoms = unitcell::get_atom_data(cell.atoms, cell.ntype, cell.lat0);
    std::vector<double> reaction(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        reaction[ir] = -0.5 * electrostatic_potential_[ir];
    }
    std::vector<ModuleBase::Vector3<double>> ionic_force;
    ModuleSccs::gauss_ion_force(atoms, reaction, basis, cell.tpiba, ModuleSccs::gaussian_ion_spread,
                                ionic_force);
    if (config_.core_electrons)
    {
        std::vector<ModuleBase::Vector3<double>> core_force;
        ModuleSccs::gauss_core_frc(atoms, cavity_potential_, basis, cell.tpiba, config_.core_spreads,
                                   core_force);
        for (int ia = 0; ia < cell.nat; ++ia)
        {
            ionic_force[ia] += core_force[ia];
        }
    }
    for (int ia = 0; ia < cell.nat; ++ia)
    {
        for (int axis = 0; axis < 3; ++axis)
        {
            force(ia, axis) += 2.0 * ionic_force[ia][axis];
        }
    }
    ModuleBase::timer::end("PotSccs", "add_sol_force");
}

void PotSccs::write_iter(std::ostream& output, int level, double residual, double pcc_energy) const
{
    if (level == 0)
    {
        return;
    }
    if (!sccs_active_)
    {
        write_sccs_wait(output, residual, config_.start_drho, config_.start_nmax);
        return;
    }
    const double energy = get_energy();
    const bool pcc = config_.boundary != ModuleSccs::Boundary::Periodic;
    write_sccs_out(output, output_, level, energy, pcc, pcc_energy);
}

double PotSccs::get_energy() const
{
    return electrostatic_rydberg_ + non_electrostatic_rydberg_;
}

void PotSccs::sol_energy(double& electrostatic, double& non_electrostatic) const
{
    electrostatic = electrostatic_rydberg_;
    non_electrostatic = non_electrostatic_rydberg_;
}

const std::vector<double>* PotSccs::sol_elec_pot() const
{
    return &electrostatic_potential_;
}

void PotSccs::add_sol_fields(std::vector<SolventGridField>& fields) const
{
    // Fields are written collectively: gate on output_.valid, set on every rank,
    // not on epsilon_, which is empty on a rank without grid points.
    if (!sccs_active_ || !output_.valid)
    {
        return;
    }
    SolventGridField dielectric;
    dielectric.name = "eps";
    dielectric.values = epsilon_;
    fields.push_back(dielectric);
    SolventGridField cavity;
    cavity.name = "cavity";
    cavity.values = solute_;
    fields.push_back(cavity);
}
} // namespace elecstate
