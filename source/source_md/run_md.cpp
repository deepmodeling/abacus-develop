#include "run_md.h"

#include "source_cell/distributed_mdcell_reader.h"
#include "source_cell/md_cell.h"
#include "source_io/module_parameter/parameter.h"
#include "fire.h"
#include "langevin.h"
#include "md_func.h"
#include "source_base/global_file.h"
#include "source_base/timer.h"
#include "source_io/module_output/print_info.h"
#include "msst.h"
#include "nhchain.h"
#include "verlet.h"
#include "source_cell/update_cell.h"
#include "source_cell/print_cell.h"
namespace Run_MD
{
namespace
{
}

void md_line(MdCell& mdcell, ModuleESolver::ESolver* p_esolver, const Parameter& param_in)
{
    ModuleBase::TITLE("Run_MD", "md_line");
    ModuleBase::timer::start("Run_MD", "md_line");

    /// determine the md_type
    MD_base* mdrun = nullptr;
    if (param_in.mdp.md_type == "fire")
    {
        mdrun = new FIRE(param_in, mdcell);
    }
    else if ((param_in.mdp.md_type == "nvt" && param_in.mdp.md_thermostat == "nhc") || param_in.mdp.md_type == "npt")
    {
        mdrun = new Nose_Hoover(param_in, mdcell);
    }
    else if (param_in.mdp.md_type == "nve" || param_in.mdp.md_type == "nvt")
    {
        mdrun = new Verlet(param_in, mdcell);
    }
    else if (param_in.mdp.md_type == "langevin")
    {
        mdrun = new Langevin(param_in, mdcell);
    }
    else if (param_in.mdp.md_type == "msst")
    {
        mdrun = new MSST(param_in, mdcell);
    }
    else
    {
        ModuleBase::WARNING_QUIT("md_line", "no such md_type!");
    }

    /// md cycle, mohan update 2026-01-04, change '<=' to '<'
    while ((mdrun->step_ + mdrun->step_rst_) < param_in.mdp.md_nstep && !mdrun->stop)
    {
        if (mdrun->step_ == 0)
        {
            mdrun->setup(p_esolver, PARAM.globalv.global_readin_dir);
        }
        else
        {
            // mohan add 2026-01-04
            const int stress_step = 0;
            const int force_step = 0;
            const int istep_print = mdrun->step_ + mdrun->step_rst_ + 1;
            ModuleIO::print_screen(stress_step, force_step, istep_print);
            mdrun->first_half(GlobalV::ofs_running);

            /// update force and virial due to the update of atom positions
            MD_func::force_virial(p_esolver,
                                  mdrun->step_,
                                  mdcell,
                                  mdrun->potential,
                                  mdrun->force,
                                  param_in.inp.cal_stress,
                                  mdrun->virial);

            mdrun->second_half();

            MD_func::compute_stress(mdcell,
                                    mdrun->vel,
                                    mdrun->allmass,
                                    param_in.inp.cal_stress,
                                    mdrun->virial,
                                    mdrun->stress);
            mdrun->t_current = MD_func::current_temp(mdrun->kinetic,
                                                     mdcell,
                                                     mdrun->frozen_freedom_,
                                                     mdrun->allmass,
                                                     mdrun->vel);
        }

        if ((mdrun->step_ + mdrun->step_rst_) % param_in.mdp.md_dumpfreq == 0)
        {
            mdrun->print_md(GlobalV::ofs_running, PARAM.inp.cal_stress);

            MD_func::dump_info(mdrun->step_ + mdrun->step_rst_,
                               PARAM.globalv.global_out_dir,
                               mdcell,
                               param_in,
                               mdrun->virial,
                               mdrun->force,
                               mdrun->vel);
        }

        if ((mdrun->step_ + mdrun->step_rst_) % param_in.mdp.md_restartfreq == 0)
        {
            if (mdcell.has_backing_unitcell())
            {
                UnitCell& unit_in = mdcell.backing_unitcell();
                int iat = 0;
                for (int it = 0; it < unit_in.ntype; ++it)
                {
                    for (int ia = 0; ia < unit_in.atoms[it].na; ++ia)
                    {
                        unit_in.atoms[it].tau[ia] = mdcell.owned_atoms()[static_cast<std::size_t>(iat)].cart;
                        unit_in.atoms[it].taud[ia] = mdcell.owned_atoms()[static_cast<std::size_t>(iat)].frac;
                        unit_in.atoms[it].vel[ia] = mdrun->vel[iat];
                        ++iat;
                    }
                }
                unitcell::update_vel(mdrun->vel, unit_in.ntype, unit_in.nat, unit_in.atoms);
                std::stringstream file;
                file << PARAM.globalv.global_stru_dir << "STRU_MD_" << mdrun->step_ + mdrun->step_rst_;
                bool need_orb = PARAM.inp.basis_type=="pw";
                need_orb = need_orb && PARAM.inp.init_wfc.substr(0, 3)=="nao";
                need_orb = need_orb || PARAM.inp.basis_type=="lcao";
                need_orb = need_orb || PARAM.inp.basis_type=="lcao_in_pw";
                unitcell::print_stru_file(unit_in,
                                        unit_in.atoms,
                                        unit_in.latvec,
                                        file.str(),
                                        PARAM.inp.nspin,
                                        false,
                                        PARAM.inp.calculation == "md",
                                        PARAM.inp.out_mul,
                                        need_orb,
                                        PARAM.globalv.deepks_setorb,
                                        GlobalV::MY_RANK);
            }
            mdrun->write_restart(PARAM.globalv.global_out_dir);
        }

        mdrun->step_++;
    }

    delete mdrun;
    ModuleBase::timer::end("Run_MD", "md_line");
    return;
}

} // namespace Run_MD
