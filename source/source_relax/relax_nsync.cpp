#include "relax_nsync.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_io/module_parameter/parameter.h"
#include "source_cell/update_cell.h"

void IonCellOptimizer::init_relax(const int& natom)
{
    // Geometry optimization algorithm setup.
    if (PARAM.inp.calculation == "relax")
    {
        // Ions_Move_Methods
        IMM.allocate(natom, PARAM.inp.relax_method[0], PARAM.inp.relax_method[1]);
    }
    if (PARAM.inp.calculation == "cell-relax")
    {
        // Ions_Move_Methods
        IMM.allocate(natom, PARAM.inp.relax_method[0], PARAM.inp.relax_method[1]);
        // allocate arrays related to changes of lattice vectors
        LCM.allocate();
    }
}

// The interface for relaxation
bool IonCellOptimizer::relax_step(const int& istep,
                           const double& energy,
                           UnitCell& ucell,
                           ModuleBase::matrix force,
                           ModuleBase::matrix stress,
                           int& force_step,
                           int& stress_step)
{
    ModuleBase::TITLE("IonCellOptimizer", "relax_step");

    ucell.ionic_position_updated = false;
    ucell.cell_parameter_updated = false;

    if (istep == PARAM.inp.relax_nmax)
    {
        return true;
    }

    const bool is_cell_relax = (PARAM.inp.calculation == "cell-relax");
    const bool is_relax = (PARAM.inp.calculation == "relax");

    if (!is_cell_relax)
    {
        force_step = istep;
    }

    const bool need_atom_relax = (is_relax || is_cell_relax) && ucell.if_atoms_can_move();
    const bool need_cell_relax = is_cell_relax && ucell.if_cell_can_change();

    if (need_atom_relax)
    {
        assert(PARAM.inp.cal_force == 1);
        bool converged = this->do_relax(istep, force, energy, ucell, force_step, GlobalV::ofs_running);
        if (!converged)
        {
            ucell.ionic_position_updated = true;
            return false;
        }
        else if (!is_cell_relax)
        {
            return true;
        }
    }
    else if (is_relax)
    {
        ModuleBase::WARNING("IonCellOptimizer", "No atom is allowed to move!");
        return true;
    }

    if (need_cell_relax)
    {
        if (ucell.if_atoms_can_move() && !IMM.get_converged())
        {
            GlobalV::ofs_running << "Note: Need to wait for atomic relaxation first!" << std::endl;
            return false;
        }

        assert(PARAM.inp.cal_stress == 1);
        bool converged = this->do_cellrelax(istep, stress_step, stress, energy, ucell, GlobalV::ofs_running);
        if (!converged)
        {
            force_step = 1;
            stress_step++;
            ucell.cell_parameter_updated = true;
            unitcell::setup_cell_after_vc(ucell, GlobalV::ofs_running);
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SETUP UNITCELL");
        }
        return converged;
    }
    else if (is_cell_relax)
    {
        ModuleBase::WARNING("IonCellOptimizer", "Lattice vectors are not allowed to change!");
        return true;
    }

    return true;
}
bool IonCellOptimizer::do_relax(const int& istep,
                         const ModuleBase::matrix& ionic_force,
                         const double& total_energy,
                         UnitCell& ucell,
                         int& jstep,
                         std::ofstream& ofs)
{
    ModuleBase::TITLE("IonCellOptimizer", "do_relax");
    IMM.cal_movement(istep, jstep, ionic_force, total_energy, ucell, ofs);
    ++jstep;
    return IMM.get_converged();
}
bool IonCellOptimizer::do_cellrelax(const int& istep,
                             const int& stress_step,
                             const ModuleBase::matrix& stress,
                             const double& total_energy,
                             UnitCell& ucell,
                             std::ofstream& ofs)
{
    ModuleBase::TITLE("IonCellOptimizer", "do_cellrelax");
    LCM.cal_lattice_change(istep, stress_step, stress, total_energy, ucell, ofs);
    return LCM.get_converged();
}
