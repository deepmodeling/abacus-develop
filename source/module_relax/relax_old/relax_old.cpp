#include "relax_old.h"

#include "module_base/global_function.h"
#include "module_base/global_variable.h"

void Relax_old::init_relax(const int& natom)
{
    // Geometry optimization algorithm setup.
    if (GlobalV::CALCULATION == "relax")
    {
        // Ions_Move_Methods
        IMM.allocate(natom);
    }
    if (GlobalV::CALCULATION == "cell-relax")
    {
        // Ions_Move_Methods
        IMM.allocate(natom);
        // allocate arrays related to changes of lattice vectors
        LCM.allocate();
    }
}

// The interface for relaxation
bool Relax_old::relax_step(const int& istep,
                           const double& energy,
                           UnitCell& ucell,
                           ModuleBase::matrix force,
                           ModuleBase::matrix stress,
                           int& force_step,
                           int& stress_step)
{
    ModuleBase::TITLE("Relax_old", "relax_step");

    // should not do it this way, will change after the refactor of ucell class
    ucell.ionic_position_updated = false;
    ucell.cell_parameter_updated = false;

    // stop in last step
    if (istep == GlobalV::RELAX_NMAX)
    {
        return 1;
    }
    // choose what to do next
    if (GlobalV::CALCULATION != "cell-relax")
        force_step = istep;
    if (this->if_do_relax(ucell))
    {
        // do relax calculation and generate next structure
        bool converged = 0;
        converged = this->do_relax(istep, force, energy, ucell, force_step);
        if (!converged)
        {
            ucell.ionic_position_updated = true;
            return converged;
        }
        else if (GlobalV::CALCULATION != "cell-relax")
        {
            return converged;
        }
    }
    if (this->if_do_cellrelax(ucell))
    {
        // do cell relax calculation and generate next structure
        bool converged = 0;
        converged = this->do_cellrelax(istep, stress_step, stress, energy, ucell);
        if (!converged)
        {
            force_step = 1;
            stress_step++;
            ucell.cell_parameter_updated = true;
            ucell.setup_cell_after_vc(GlobalV::ofs_running);
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SETUP UNITCELL");
        }
        return converged;
    }

    return 1;
}

bool Relax_old::if_do_relax(const UnitCell& ucell)
{
    ModuleBase::TITLE("Relax_old", "if_do_relax");
    if (GlobalV::CALCULATION == "relax" || GlobalV::CALCULATION == "cell-relax")
    {
        if (!ucell.if_atoms_can_move())
        {
            ModuleBase::WARNING("Ions", "No atom is allowed to move!");
            return 0;
        }
        //		if(!IMM.get_converged()) return 1;
        else
        {
            assert(GlobalV::CAL_FORCE == 1);
            return 1;
        }
    }
    else
        return 0;
}
bool Relax_old::if_do_cellrelax(const UnitCell& ucell)
{
    ModuleBase::TITLE("Relax_old", "if_do_cellrelax");
    if (GlobalV::CALCULATION == "cell-relax")
    {
        if (!ucell.if_cell_can_change())
        {
            ModuleBase::WARNING("Ions", "Lattice vectors are not allowed to change!");
            return 0;
        }
        else if (ucell.if_atoms_can_move() && !IMM.get_converged())
        {
            GlobalV::ofs_running << "Note: Need to wait for atomic relaxation first!";
            return 0;
        }
        else
        {
            assert(GlobalV::CAL_STRESS == 1);
            return 1;
        }
    }
    else
        return 0;
}
bool Relax_old::do_relax(const int& istep,
                         const ModuleBase::matrix& ionic_force,
                         const double& total_energy,
                         UnitCell& ucell,
                         int& jstep)
{
    ModuleBase::TITLE("Relax_old", "do_relax");
    IMM.cal_movement(istep, jstep, ionic_force, total_energy, ucell);
    ++jstep;
    return IMM.get_converged();
}
bool Relax_old::do_cellrelax(const int& istep,
                             const int& stress_step,
                             const ModuleBase::matrix& stress,
                             const double& total_energy,
                             UnitCell& ucell)
{
    ModuleBase::TITLE("Relax_old", "do_cellrelax");
    LCM.cal_lattice_change(istep, stress_step, stress, total_energy, ucell);
    return LCM.get_converged();
}
