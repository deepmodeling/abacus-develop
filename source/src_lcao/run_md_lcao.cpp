#include "run_md_lcao.h"
#include "FORCE_STRESS.h"
#include "../src_pw/global.h"
#include "../module_orbital/parallel_orbitals.h"
#include "../src_pdiag/pdiag_double.h"
#include "../src_io/write_HS.h"
#include "../src_io/cal_r_overlap_R.h"
#include "../src_io/print_info.h"
#include "../module_relaxation/variable_cell.h" // mohan add 2021-02-01
#include "../src_ri/exx_abfs.h"
#include "../src_ri/exx_opt_orb.h"
#include "../module_neighbor/sltk_atom_arrange.h"
#include "../module_md/MD_func.h"
#include "../module_md/FIRE.h"
#include "../module_md/verlet.h"
#include "../module_md/MSST.h"
#include "../module_md/Nose_Hoover.h"
#include "../module_md/Langevin.h"

Run_MD_LCAO::Run_MD_LCAO()
{
    cellchange = false;
    CE.Init_CE();
}

Run_MD_LCAO::~Run_MD_LCAO()
{
}

void Run_MD_LCAO::opt_ions(ModuleESolver::ESolver* p_esolver)
{
    ModuleBase::TITLE("Run_MD_LCAO", "opt_ions");
    ModuleBase::timer::tick("Run_MD_LCAO", "opt_ions");

    if (GlobalV::OUT_LEVEL == "i")
    {
        std::cout << std::setprecision(12);
        std::cout << " " << std::setw(7) << "ISTEP" << std::setw(5) << "NE" << std::setw(18) << "ETOT(eV)"
                  << std::setw(10) << "dE(meV)" << std::setw(10) << "F(eV/A)" << std::setw(10) << "T(MIN)" << std::endl;
    }

    //Charge_Extrapolation
    // CE.allocate_ions();

    // determine the md_type
    MDrun* mdrun;
    if (INPUT.mdp.md_type == -1)
    {
        mdrun = new FIRE(INPUT.mdp, GlobalC::ucell);
    }
    else if (INPUT.mdp.md_type == 0)
    {
        mdrun = new Verlet(INPUT.mdp, GlobalC::ucell);
    }
    else if (INPUT.mdp.md_type == 1)
    {
        mdrun = new Nose_Hoover(INPUT.mdp, GlobalC::ucell);
    }
    else if (INPUT.mdp.md_type == 2)
    {
        mdrun = new Langevin(INPUT.mdp, GlobalC::ucell);
    }
    else if (INPUT.mdp.md_type == 4)
    {
        mdrun = new MSST(INPUT.mdp, GlobalC::ucell);
        cellchange = true;
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
            CE.update_all_pos(GlobalC::ucell);

            mdrun->first_half();

            if (cellchange)
            {
                CE.update_istep(1);
            }
            else
            {
                CE.update_istep(mdrun->step_);
            }

            CE.save_pos_next(GlobalC::ucell);
            CE.extrapolate_charge();

            if (cellchange)
            {
                Variable_Cell::init_after_vc();
            }

            // reset local potential
            GlobalC::pot.init_pot(mdrun->step_, GlobalC::sf.strucFac);

            // update force and virial due to the update of atom positions

            MD_func::force_virial(p_esolver, mdrun->step_, mdrun->mdp, mdrun->ucell, mdrun->potential, mdrun->force, mdrun->virial);

            mdrun->second_half();

            MD_func::compute_stress(mdrun->ucell, mdrun->vel, mdrun->allmass, mdrun->virial, mdrun->stress);
            mdrun->t_current = MD_func::current_temp(mdrun->kinetic, mdrun->ucell.nat, mdrun->frozen_freedom_, mdrun->allmass, mdrun->vel);
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

    if (GlobalC::pot.out_pot == 2)
    {
        std::stringstream ssp;
        std::stringstream ssp_ave;
        ssp << GlobalV::global_out_dir << "ElecStaticPot";
        ssp_ave << GlobalV::global_out_dir << "ElecStaticPot_AVE";
        GlobalC::pot.write_elecstat_pot(ssp.str(), ssp_ave.str(), GlobalC::rhopw); // output 'Hartree + local pseudopot'
    }

#ifdef __MPI
    atom_arrange::delete_vector(GlobalV::ofs_running,
                                GlobalV::SEARCH_PBC,
                                GlobalC::GridD,
                                GlobalC::ucell,
                                GlobalV::SEARCH_RADIUS,
                                GlobalV::test_atom_input);
#endif

    delete mdrun;
    ModuleBase::timer::tick("Run_MD_LCAO","opt_ions"); 
    return;
}
