#include "relax_nsync.h"
#include "relax_history.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_cell/cell_tools.h"
#include "source_cell/update_cell.h"

void IonCellOptimizer::print_converged_summary(const int istep,
                                               const bool with_stress,
                                               std::ofstream& ofs_running) const
{
    ofs_running << " Relaxation method: " << inp_->relax_method[0] << std::endl;
    ofs_running << " Relaxation converged in " << istep << " step(s)." << std::endl;
    if (!max_force_history_.empty())
    {
        ofs_running << " Largest force per step (eV/Angstrom):" << format_relax_history(max_force_history_);
    }
    if (with_stress && !max_stress_history_.empty())
    {
        ofs_running << " Largest stress per step (kbar):" << format_relax_history(max_stress_history_);
    }
    ofs_running << "\n Relaxation is converged!" << std::endl;
}

/**
 * @brief Initialize relaxation algorithms based on calculation type.
 *
 * Allocates memory and initializes the appropriate relaxation methods:
 * - For "relax" calculation: only initializes Ions_Move_Methods
 * - For "cell-relax" calculation: initializes both Ions_Move_Methods and 
 *   Lattice_Change_Methods
 * 
 * @param natom Number of atoms in the system
 */
void IonCellOptimizer::init_relax(const int& natom, const Input_para& inp)
{
    inp_ = &inp;
    max_force_history_.clear();
    max_stress_history_.clear();

    if (inp_->calculation == "relax")
    {
        IMM.allocate(natom, inp_->relax_method[0], inp_->relax_method[1]);
    }
    if (inp_->calculation == "cell-relax")
    {
        IMM.allocate(natom, inp_->relax_method[0], inp_->relax_method[1]);
        LCM.allocate();
    }
}

/**
 * @brief Perform one step of relaxation (atomic and/or cell).
 * 
 * Main relaxation loop that coordinates atomic and cell relaxation:
 * 1. Check for maximum iteration limit
 * 2. Determine calculation mode (relax vs cell-relax)
 * 3. Perform atomic relaxation if needed and atoms can move
 * 4. If in cell-relax mode and atomic relaxation converged, perform cell relaxation
 * 
 * Convergence behavior:
 * - Returns false if relaxation is still in progress
 * - Returns true if relaxation has converged or maximum iterations reached
 * 
 * @param istep Current total iteration step
 * @param energy Total energy of the system
 * @param ucell Unit cell containing atomic positions and lattice vectors
 * @param force Ionic forces matrix (natoms x 3)
 * @param stress Stress tensor matrix (3 x 3)
 * @param force_step Current step counter for force-based relaxation (output)
 * @param stress_step Current step counter for stress-based relaxation (output)
 * @return true if relaxation is converged, false otherwise
 */
bool IonCellOptimizer::relax_step(const int& istep,
                           const double& energy,
                           UnitCell& ucell,
                           ModuleBase::matrix force,
                           ModuleBase::matrix stress,
                           int& force_step,
                           int& stress_step,
                           std::ofstream& ofs_running)
{
    ModuleBase::TITLE("IonCellOptimizer", "relax_step");

    // Reset update flags at the beginning of each step
    ucell.ionic_position_updated = false;
    ucell.cell_parameter_updated = false;

    // relax_nmax == 0 is a valid dry-run mode: no relaxation step may run, and
    // it must not be reported as a failed relaxation (mirrors the dry-run
    // branch of Relax_Driver::final_out). Terminate immediately and silently.
    if (inp_->relax_nmax == 0)
    {
        return true;
    }

    // Check if we've reached the maximum number of iterations.
    if (istep == inp_->relax_nmax)
    {
        // This step never ran cal_movement, so no force was recorded for it;
        // the actual number of ionic steps taken is the history length.
        const int ionic_steps = static_cast<int>(max_force_history_.size());
        ofs_running << " Relaxation method: " << inp_->relax_method[0] << std::endl;
        ofs_running << " Relaxation stopped after " << ionic_steps << " ionic step(s) (relax_nmax = "
                    << inp_->relax_nmax << " reached)." << std::endl;
        if (!max_force_history_.empty())
        {
            ofs_running << " Largest force per step (eV/Angstrom):" << format_relax_history(max_force_history_);
        }
        if (!max_stress_history_.empty())
        {
            ofs_running << " Largest stress per step (kbar):" << format_relax_history(max_stress_history_);
        }
        ofs_running << " Relaxation is not converged after reaching relax_nmax!" << std::endl;
        return true;
    }

    // Determine calculation mode
    const bool is_cell_relax = (inp_->calculation == "cell-relax");
    const bool is_relax = (inp_->calculation == "relax");

    // In non-cell-relax mode, force_step follows istep
    if (!is_cell_relax)
    {
        force_step = istep;
    }

    // Determine what relaxation steps are needed
    const bool need_atom_relax = (is_relax || is_cell_relax) && unitcell::if_atoms_can_move(ucell.atoms, ucell.ntype);
    const bool need_cell_relax = is_cell_relax && unitcell::if_cell_can_change(ucell.lat_axis_free);

    // Track whether the atomic part has converged this step (used to decide
    // which final marker to print and whether cell relaxation may proceed)
    bool atom_converged = false;
    bool atom_relax_performed = false;

    // Atomic relaxation branch
    if (need_atom_relax)
    {
        assert(inp_->cal_force == 1);
        
        // Calculate and apply atomic movement
        std::vector<std::string> relax_method = inp_->relax_method;

        Relax_Criteria criteria;
        criteria.force_thr = inp_->force_thr;
        criteria.force_thr_ev = inp_->force_thr_ev;
        criteria.stress_thr = inp_->stress_thr;
        criteria.fixed_ibrav = inp_->fixed_ibrav;
        criteria.out_level = inp_->out_level;
        criteria.test_relax_method = inp_->test_relax_method;

        IMM.cal_movement(istep, force_step, force, energy, ucell, ofs_running, relax_method, criteria);
        ++force_step;
        atom_relax_performed = true;
        atom_converged = IMM.get_converged();

        const double max_force_ev = IMM.get_largest_grad() * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A;
        max_force_history_.push_back(max_force_ev);
        ofs_running << " Largest force is " << max_force_ev
                    << " eV/Angstrom while threshold is " << inp_->force_thr_ev << " eV/Angstrom" << std::endl;
        if (!atom_converged)
        {
            ofs_running << "\n Relaxation is not converged yet!" << std::endl;
            ucell.ionic_position_updated = true;
            return false; // not converged
        }
        else if (!is_cell_relax)
        {
            print_converged_summary(istep, false, ofs_running);
            return true; // converged
        }
        // When the lattice is fully fixed (e.g. fixed_axes = abc), there is no
        // cell step to run; the atomic convergence above is the final result.
        else if (!need_cell_relax)
        {
            print_converged_summary(istep, false, ofs_running);
            return true; // converged, ions relaxed with a fixed lattice
        }
        // Otherwise, continue to cell relaxation
    }
    else if (is_relax)
    {
        // Relax mode but no atoms can move - nothing to do. The no-op run is
        // still a valid converged relaxation, so emit the unified summary to
        // give the running log an explicit success marker.
        ModuleBase::WARNING("IonCellOptimizer", "No atoms are allowed to move!");
        print_converged_summary(istep, false, ofs_running);
        return true;
    }

    // Cell relaxation branch (only in cell-relax mode)
    if (need_cell_relax)
    {
        assert(inp_->cal_stress == 1);
        
        // Calculate and apply lattice change
        Relax_Criteria criteria;
        criteria.force_thr = inp_->force_thr;
        criteria.force_thr_ev = inp_->force_thr_ev;
        criteria.stress_thr = inp_->stress_thr;
        criteria.fixed_ibrav = inp_->fixed_ibrav;
        criteria.out_level = inp_->out_level;
        criteria.test_relax_method = inp_->test_relax_method;

        LCM.cal_lattice_change(istep, stress_step, stress, energy, ucell, ofs_running, criteria);
        bool converged = LCM.get_converged();

        const double max_stress_kbar = LCM.get_largest_grad();
        max_stress_history_.push_back(max_stress_kbar);
        ofs_running << " Largest stress is " << max_stress_kbar
                    << " kbar while threshold is " << inp_->stress_thr << " kbar" << std::endl;
        if (converged)
        {
            print_converged_summary(istep, true, ofs_running);
        }
        else
        {
            ofs_running << "\n Relaxation is not converged yet!" << std::endl;
            // Reset force_step counter after cell change for fresh atomic relaxation
            force_step = 1;
            stress_step++;
            IMM.reset_after_cell_change(inp_->relax_method, ofs_running);
            ucell.cell_parameter_updated = true;
            
            // Update cell-related parameters after volume change
            unitcell::setup_cell_after_vc(ucell, ofs_running, inp_->nspin);
            ModuleBase::GlobalFunc::DONE(ofs_running, "SETUP UNITCELL");
        }
        
        return converged;
    }
    else if (is_cell_relax && !unitcell::if_cell_can_change(ucell.lat_axis_free))
    {
        ModuleBase::WARNING("IonCellOptimizer", "Lattice vectors are not allowed to change!");
        return true;
    }

    return true;
}
