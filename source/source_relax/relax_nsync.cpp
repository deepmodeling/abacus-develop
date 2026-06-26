#include "relax_nsync.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_io/module_parameter/parameter.h"
#include "source_cell/update_cell.h"

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
void IonCellOptimizer::init_relax(const int& natom)
{
    if (PARAM.inp.calculation == "relax")
    {
        IMM.allocate(natom, PARAM.inp.relax_method[0], PARAM.inp.relax_method[1]);
    }
    if (PARAM.inp.calculation == "cell-relax")
    {
        IMM.allocate(natom, PARAM.inp.relax_method[0], PARAM.inp.relax_method[1]);
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

    // Check if we've reached the maximum number of iterations
    if (istep == PARAM.inp.relax_nmax)
    {
        return true;
    }

    // Determine calculation mode
    const bool is_cell_relax = (PARAM.inp.calculation == "cell-relax");
    const bool is_relax = (PARAM.inp.calculation == "relax");

    // In non-cell-relax mode, force_step follows istep
    if (!is_cell_relax)
    {
        force_step = istep;
    }

    // Determine what relaxation steps are needed
    const bool need_atom_relax = (is_relax || is_cell_relax) && ucell.if_atoms_can_move();
    const bool need_cell_relax = is_cell_relax && ucell.if_cell_can_change();

    // Atomic relaxation branch
    if (need_atom_relax)
    {
        assert(PARAM.inp.cal_force == 1);
        
        // Calculate and apply atomic movement
        std::vector<std::string> relax_method = PARAM.inp.relax_method;
        IMM.cal_movement(istep, force_step, force, energy, ucell, ofs_running, relax_method);
        ++force_step;
        
        // Check convergence
        bool converged = IMM.get_converged();
        if (!converged)
        {
            ucell.ionic_position_updated = true;
            return false; // not converged
        }
        else if (!is_cell_relax)
        {
            return true; // converged
        }
        // Otherwise, continue to cell relaxation
    }
    else if (is_relax)
    {
        // Relax mode but no atoms can move - nothing to do
        ModuleBase::WARNING("IonCellOptimizer", "No atom is allowed to move!");
        return true;
    }

    // Cell relaxation branch (only in cell-relax mode)
    if (need_cell_relax)
    {
        assert(PARAM.inp.cal_stress == 1);
        
        // Calculate and apply lattice change
        LCM.cal_lattice_change(istep, stress_step, stress, energy, ucell, ofs_running);
        bool converged = LCM.get_converged();
        
        if (!converged)
        {
            // Reset force_step counter after cell change for fresh atomic relaxation
            force_step = 1;
            stress_step++;
            IMM.reset_after_cell_change(PARAM.inp.relax_method, ofs_running);
            ucell.cell_parameter_updated = true;
            
            // Update cell-related parameters after volume change
            unitcell::setup_cell_after_vc(ucell, ofs_running);
            ModuleBase::GlobalFunc::DONE(ofs_running, "SETUP UNITCELL");
        }
        
        return converged;
    }
    else if (is_cell_relax && !ucell.if_cell_can_change())
    {
        ModuleBase::WARNING("IonCellOptimizer", "Lattice vectors are not allowed to change!");
        return true;
    }

    return true;
}
