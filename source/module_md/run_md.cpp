#include "run_md.h"

#include "fire.h"
#include "langevin.h"
#include "md_func.h"
#include "module_base/timer.h"
#include "module_io/print_info.h"
#include "msst.h"
#include "nhchain.h"
#include "verlet.h"

Run_MD::Run_MD()
{
}

Run_MD::~Run_MD()
{
}

void Run_MD::md_line(UnitCell &unit_in, ModuleESolver::ESolver *p_esolver, MD_parameters &md_para)
{
    ModuleBase::TITLE("Run_MD", "md_line");
    ModuleBase::timer::tick("Run_MD", "md_line");

    // determine the md_type
    MDrun *mdrun;
    if (md_para.md_type == "fire")
    {
        mdrun = new FIRE(md_para, unit_in);
    }
    else if ((md_para.md_type == "nvt" && md_para.md_thermostat == "nhc") || md_para.md_type == "npt")
    {
        mdrun = new Nose_Hoover(md_para, unit_in);
    }
    else if (md_para.md_type == "nve" || md_para.md_type == "nvt")
    {
        mdrun = new Verlet(md_para, unit_in);
    }
    else if (md_para.md_type == "langevin")
    {
        mdrun = new Langevin(md_para, unit_in);
    }
    else if (md_para.md_type == "msst")
    {
        mdrun = new MSST(md_para, unit_in);
    }
    else
    {
        ModuleBase::WARNING_QUIT("md_line", "no such md_type!");
    }

    // md cycle
    while ((mdrun->step_ + mdrun->step_rst_) <= md_para.md_nstep && !mdrun->stop)
    {
        if (mdrun->step_ == 0)
        {
            mdrun->setup(p_esolver, GlobalV::MY_RANK, GlobalV::global_readin_dir);
        }
        else
        {
            Print_Info::print_screen(0, 0, mdrun->step_ + mdrun->step_rst_);
            mdrun->first_half(GlobalV::MY_RANK, GlobalV::ofs_running);

            // update force and virial due to the update of atom positions
            MD_func::force_virial(p_esolver, mdrun->step_, mdrun->ucell, mdrun->potential, mdrun->force, mdrun->virial);

            mdrun->second_half(GlobalV::MY_RANK);

            MD_func::compute_stress(mdrun->ucell, mdrun->vel, mdrun->allmass, mdrun->virial, mdrun->stress);
            mdrun->t_current = MD_func::current_temp(mdrun->kinetic,
                                                     mdrun->ucell.nat,
                                                     mdrun->frozen_freedom_,
                                                     mdrun->allmass,
                                                     mdrun->vel);
        }

        if ((mdrun->step_ + mdrun->step_rst_) % mdrun->mdp.md_dumpfreq == 0)
        {
            mdrun->outputMD(GlobalV::ofs_running, GlobalV::CAL_STRESS, GlobalV::MY_RANK);

            MD_func::MDdump(mdrun->step_ + mdrun->step_rst_,
                            mdrun->ucell,
                            md_para,
                            mdrun->virial,
                            mdrun->force,
                            mdrun->vel);
        }

        if ((mdrun->step_ + mdrun->step_rst_) % mdrun->mdp.md_restartfreq == 0)
        {
            mdrun->ucell.update_vel(mdrun->vel);
            std::stringstream file;
            file << GlobalV::global_stru_dir << "STRU_MD_" << mdrun->step_ + mdrun->step_rst_;
            mdrun->ucell.print_stru_file(file.str(), 1, 1);
            mdrun->write_restart(GlobalV::MY_RANK, GlobalV::global_out_dir);
        }

        mdrun->step_++;
    }

    delete mdrun;
    ModuleBase::timer::tick("Run_MD", "md_line");
    return;
}