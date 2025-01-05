#include "run_md.h"

#include "module_parameter/parameter.h"
#include "fire.h"
#include "langevin.h"
#include "md_func.h"
#include "module_base/global_file.h"
#include "module_base/timer.h"
#include "module_io/print_info.h"
#include "msst.h"
#include "nhchain.h"
#include "verlet.h"
#include "module_cell/update_cell.h"
namespace Run_MD
{

void md_line(UnitCell& unit_in, ModuleESolver::ESolver* p_esolver, const Parameter& param_in)
{
    ModuleBase::TITLE("Run_MD", "md_line");
    ModuleBase::timer::tick("Run_MD", "md_line");

    /// determine the md_type
    MD_base* mdrun;
    if (param_in.mdp.md_type == "fire")
    {
        mdrun = new FIRE(param_in, unit_in);
    }
    else if ((param_in.mdp.md_type == "nvt" && param_in.mdp.md_thermostat == "nhc") || param_in.mdp.md_type == "npt")
    {
        mdrun = new Nose_Hoover(param_in, unit_in);
    }
    else if (param_in.mdp.md_type == "nve" || param_in.mdp.md_type == "nvt")
    {
        mdrun = new Verlet(param_in, unit_in);
    }
    else if (param_in.mdp.md_type == "langevin")
    {
        mdrun = new Langevin(param_in, unit_in);
    }
    else if (param_in.mdp.md_type == "msst")
    {
        mdrun = new MSST(param_in, unit_in);
    }
    else
    {
        ModuleBase::WARNING_QUIT("md_line", "no such md_type!");
    }

    /// md cycle
    while ((mdrun->step_ + mdrun->step_rst_) <= param_in.mdp.md_nstep && !mdrun->stop)
    {
        if (mdrun->step_ == 0)
        {
            mdrun->setup(p_esolver, PARAM.globalv.global_readin_dir);
        }
        else
        {
            ModuleIO::print_screen(0, 0, mdrun->step_ + mdrun->step_rst_);
            mdrun->first_half(GlobalV::ofs_running);

            /// update force and virial due to the update of atom positions
            MD_func::force_virial(p_esolver,
                                  mdrun->step_,
                                  unit_in,
                                  mdrun->potential,
                                  mdrun->force,
                                  param_in.inp.cal_stress,
                                  mdrun->virial);

            mdrun->second_half();

            MD_func::compute_stress(unit_in,
                                    mdrun->vel,
                                    mdrun->allmass,
                                    param_in.inp.cal_stress,
                                    mdrun->virial,
                                    mdrun->stress);
            mdrun->t_current = MD_func::current_temp(mdrun->kinetic,
                                                     unit_in.nat,
                                                     mdrun->frozen_freedom_,
                                                     mdrun->allmass,
                                                     mdrun->vel);
        }

        if ((mdrun->step_ + mdrun->step_rst_) % param_in.mdp.md_dumpfreq == 0)
        {
            mdrun->print_md(GlobalV::ofs_running, PARAM.inp.cal_stress);

            MD_func::dump_info(mdrun->step_ + mdrun->step_rst_,
                               PARAM.globalv.global_out_dir,
                               unit_in,
                               param_in,
                               mdrun->virial,
                               mdrun->force,
                               mdrun->vel);
        }

        if ((mdrun->step_ + mdrun->step_rst_) % param_in.mdp.md_restartfreq == 0)
        {
            unitcell::update_vel(mdrun->vel,unit_in.ntype,unit_in.nat,unit_in.atoms);
            std::stringstream file;
            file << PARAM.globalv.global_stru_dir << "STRU_MD_" << mdrun->step_ + mdrun->step_rst_;
            // changelog 20240509
            // because I move out the dependence on GlobalV from UnitCell::print_stru_file
            // so its parameter is calculated here
            bool need_orb = PARAM.inp.basis_type=="pw";
            need_orb = need_orb && PARAM.inp.init_wfc.substr(0, 3)=="nao";
            need_orb = need_orb || PARAM.inp.basis_type=="lcao";
            need_orb = need_orb || PARAM.inp.basis_type=="lcao_in_pw";
            unit_in.print_stru_file(file.str(), 
                                    PARAM.inp.nspin, 
                                    false, // Cartesian coordinates
                                    PARAM.inp.calculation == "md", 
                                    PARAM.inp.out_mul,
                                    need_orb,
                                    PARAM.globalv.deepks_setorb,
                                    GlobalV::MY_RANK);
            mdrun->write_restart(PARAM.globalv.global_out_dir);
        }

        mdrun->step_++;
    }

    delete mdrun;
    ModuleBase::timer::tick("Run_MD", "md_line");
    return;
}

} // namespace Run_MD