#include "relax_old.h"

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_parameter/parameter.h"
#include "module_cell/update_cell.h"
void Relax_old::init_relax(const int& natom)
{
    // Geometry optimization algorithm setup.
    if (PARAM.inp.calculation == "relax")
    {
        // Ions_Move_Methods
        IMM.allocate(natom);
    }
    if (PARAM.inp.calculation == "cell-relax")
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
    if (istep == PARAM.inp.relax_nmax)
    {
        return true;
    }
    // choose what to do next
    if (PARAM.inp.calculation != "cell-relax") {
        force_step = istep;
}
    if (this->if_do_relax(ucell))
    {
        // do relax calculation and generate next structure
        bool converged = false;
        converged = this->do_relax(istep, force, energy, ucell, force_step);
        if (!converged)
        {
            ucell.ionic_position_updated = true;
            return converged;
        }
        else if (PARAM.inp.calculation != "cell-relax")
        {
            return converged;
        }
    }
    if (this->if_do_cellrelax(ucell))
    {
        // do cell relax calculation and generate next structure
        bool converged = false;
        converged = this->do_cellrelax(istep, stress_step, stress, energy, ucell);
        if (!converged)
        {
            force_step = 1;
            stress_step++;
            ucell.cell_parameter_updated = true;
            unitcell::setup_cell_after_vc(ucell,GlobalV::ofs_running);
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SETUP UNITCELL");
        }
        return converged;
    }

    return true;
}

bool Relax_old::if_do_relax(const UnitCell& ucell)
{
    ModuleBase::TITLE("Relax_old", "if_do_relax");
    if (PARAM.inp.calculation == "relax" || PARAM.inp.calculation == "cell-relax")
    {
        if (!ucell.if_atoms_can_move())
        {
            ModuleBase::WARNING("Ions", "No atom is allowed to move!");
            return false;
        }
        //		if(!IMM.get_converged()) return 1;
        else
        {
            assert(PARAM.inp.cal_force == 1);
            return true;
        }
    }
    else {
        return false;
}
}
bool Relax_old::if_do_cellrelax(const UnitCell& ucell)
{
    ModuleBase::TITLE("Relax_old", "if_do_cellrelax");
    if (PARAM.inp.calculation == "cell-relax")
    {
        if (!ucell.if_cell_can_change())
        {
            ModuleBase::WARNING("Ions", "Lattice vectors are not allowed to change!");
            return false;
        }
        else if (ucell.if_atoms_can_move() && !IMM.get_converged())
        {
            GlobalV::ofs_running << "Note: Need to wait for atomic relaxation first!";
            return false;
        }
        else
        {
            assert(PARAM.inp.cal_stress == 1);
            return true;
        }
    }
    else {
        return false;
}
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
