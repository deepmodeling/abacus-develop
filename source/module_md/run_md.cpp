#include "run_md.h"

#include "../input.h"
#include "../module_base/timer.h"
#include "../src_io/print_info.h"
#include "FIRE.h"
#include "Langevin.h"
#include "MD_func.h"
#include "MSST.h"
#include "Nose_Hoover.h"
#include "verlet.h"

Run_MD::Run_MD()
{
}

Run_MD::~Run_MD()
{
}

void Run_MD::md_line(UnitCell &unit_in, ModuleESolver::ESolver *p_esolver)
{
    ModuleBase::TITLE("Run_MD", "md_line");
    ModuleBase::timer::tick("Run_MD", "md_line");

    // determine the md_type
    MDrun *mdrun;
    if (INPUT.mdp.md_type == -1)
    {
        mdrun = new FIRE(INPUT.mdp, unit_in);
    }
    else if (INPUT.mdp.md_type == 0)
    {
        mdrun = new Verlet(INPUT.mdp, unit_in);
    }
    else if (INPUT.mdp.md_type == 1)
    {
        mdrun = new Nose_Hoover(INPUT.mdp, unit_in);
    }
    else if (INPUT.mdp.md_type == 2)
    {
        mdrun = new Langevin(INPUT.mdp, unit_in);
    }
    else if (INPUT.mdp.md_type == 4)
    {
        mdrun = new MSST(INPUT.mdp, unit_in);
        unit_in.cell_parameter_updated = true;
    }

    // md cycle
    while ((mdrun->step_ + mdrun->step_rst_) <= GlobalV::MD_NSTEP && !mdrun->stop)
    {
        if (mdrun->step_ == 0)
        {
            mdrun->setup(p_esolver);
        }
        else
        {
            Print_Info::print_screen(0, 0, mdrun->step_ + mdrun->step_rst_);

            ModuleBase::Vector3<double> *vel_tmp;
            if (GlobalV::ESOLVER_TYPE == "tddft" )
            {
                vel_tmp = new ModuleBase::Vector3<double>[mdrun->ucell.nat];
                for (int i = 0; i < mdrun>ucell.nat; ++i)
                {
                    for (int k = 0; k < 3; ++k)
                    {
                        if (mdrun->ionmbl[i][k])
                        {
                            vel_tmp[i][k] = mdrun->vel[i][k];
                        }
                    }
                }
            }

            mdrun->first_half();

            // update force and virial due to the update of atom positions
            if ( GlobalV::ESOLVER_TYPE == "tddft" )
            {
                MD_func::force_virial(p_esolver,
                                      mdrun->step_,
                                      mdrun->ucell,
                                      mdrun->potential,
                                      mdrun->force,
                                      mdrun->virial,
                                      vel_tmp);
            }
            else
            {
                MD_func::force_virial(p_esolver,
                                      mdrun->step_,
                                      mdrun->ucell,
                                      mdrun->potential,
                                      mdrun->force,
                                      mdrun->virial);
            }

            mdrun->second_half();

            MD_func::compute_stress(mdrun->ucell, mdrun->vel, mdrun->allmass, mdrun->virial, mdrun->stress);
            mdrun->t_current = MD_func::current_temp(mdrun->kinetic,
                                                     mdrun->ucell.nat,
                                                     mdrun->frozen_freedom_,
                                                     mdrun->allmass,
                                                     mdrun->vel);
        }

        if ((mdrun->step_ + mdrun->step_rst_) % mdrun->mdp.md_dumpfreq == 0)
        {
            // Print_Info::print_screen(0, 0, mdrun->step_ + mdrun->step_rst_);
            mdrun->outputMD(GlobalV::ofs_running, GlobalV::CAL_STRESS);

            MD_func::MDdump(mdrun->step_ + mdrun->step_rst_, mdrun->ucell, mdrun->virial, mdrun->force);
        }

        if ((mdrun->step_ + mdrun->step_rst_) % mdrun->mdp.md_restartfreq == 0)
        {
            mdrun->ucell.update_vel(mdrun->vel);
            std::stringstream file;
            file << GlobalV::global_stru_dir << "STRU_MD_" << mdrun->step_ + mdrun->step_rst_;
            mdrun->ucell.print_stru_file(file.str(), 1, 1);
            mdrun->write_restart();
        }

        mdrun->step_++;
    }

    delete mdrun;
    ModuleBase::timer::tick("Run_MD", "md_line");
    return;
}