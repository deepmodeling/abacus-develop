#include "run_md_pw.h"
#include "global.h" // use chr.
#include "../module_relaxation/variable_cell.h" // mohan add 2021-02-01
#include "../module_md/MD_func.h"
#include "../module_md/FIRE.h"
#include "../module_md/NVE.h"
#include "../module_md/MSST.h"
#include "../module_md/NVT_ADS.h"
#include "../module_md/NVT_NHC.h"
#include "../module_md/Langevin.h"
#include "../src_io/print_info.h"

Run_MD_PW::Run_MD_PW()
{
    cellchange = false;
    CE.Init_CE();
}

Run_MD_PW::~Run_MD_PW(){}

void Run_MD_PW::md_ions_pw(ModuleESolver::ESolver *p_esolver)
{
    ModuleBase::TITLE("Run_MD_PW", "md_ions_pw");
    ModuleBase::timer::tick("Run_MD_PW", "md_ions_pw");

    if (GlobalV::OUT_LEVEL == "i")
    {
        std::cout << std::setprecision(12);
        std::cout << " " << std::setw(7) << "ISTEP"
             << std::setw(5) << "NE"
             << std::setw(15) << "ETOT(eV)"
             << std::setw(15) << "EDIFF(eV)"
             << std::setw(15) << "MAX_F(eV/A)"
             << std::setw(15) << "TRADIUS(Bohr)"
             << std::setw(8) << "UPDATE"
             << std::setw(11) << "ETIME(MIN)"
             << std::setw(11) << "FTIME(MIN)"
             << std::endl;
    }

    // allocation for ion movement.
    // CE.allocate_ions();

    // determine the md_type
    Verlet *verlet;
    if(INPUT.mdp.md_type == -1)
    {
        verlet = new FIRE(INPUT.mdp, GlobalC::ucell); 
    }
    else if(INPUT.mdp.md_type == 0)
    {
        verlet = new NVE(INPUT.mdp, GlobalC::ucell); 
    }
    else if(INPUT.mdp.md_type==1)
    {
        verlet = new NVT_NHC(INPUT.mdp, GlobalC::ucell);
    }
    else if(INPUT.mdp.md_type==2)
    {
        verlet = new Langevin(INPUT.mdp, GlobalC::ucell);
    }
    else if(INPUT.mdp.md_type == 3)
    {
        verlet = new NVT_ADS(INPUT.mdp, GlobalC::ucell);
    }
    else if(INPUT.mdp.md_type==4)
    {
        verlet = new MSST(INPUT.mdp, GlobalC::ucell); 
        cellchange = true;
    }

    // md cycle
    while ( (verlet->step_ + verlet->step_rst_) <= GlobalV::MD_NSTEP && !verlet->stop)
    {
        if(verlet->step_ == 0)
        {
            verlet->setup(p_esolver);
        }
        else
        {
            Print_Info::print_screen(0, 0, verlet->step_ + verlet->step_rst_);
            CE.update_all_pos(GlobalC::ucell);

            verlet->first_half();

            if(cellchange)
            {
                CE.update_istep(1);
            }
            else
            {
                CE.update_istep(verlet->step_);
            }

            CE.save_pos_next(GlobalC::ucell);
            CE.extrapolate_charge();

            if(cellchange)
            {
                GlobalC::ucell.cell_parameter_updated = true;
                Variable_Cell::init_after_vc();
            }

            // reset local potential and initial wave function
            GlobalC::pot.init_pot(verlet->step_, GlobalC::sf.strucFac);
            
            // new wave functions
            //GlobalC::wf.wfcinit();

            // update force and virial due to the update of atom positions
            MD_func::force_virial(p_esolver, verlet->step_, verlet->mdp, verlet->ucell, verlet->potential, verlet->force, verlet->virial);

            verlet->second_half();

            MD_func::kinetic_stress(verlet->ucell, verlet->vel, verlet->allmass, verlet->kinetic, verlet->stress);

            verlet->stress += verlet->virial;
        }

        if((verlet->step_ + verlet->step_rst_) % verlet->mdp.md_dumpfreq == 0)
        {
            // Print_Info::print_screen(0, 0, verlet->step_ + verlet->step_rst_);
            verlet->outputMD(GlobalV::ofs_running, GlobalV::CAL_STRESS);

            MD_func::MDdump(verlet->step_ + verlet->step_rst_, verlet->ucell, verlet->virial, verlet->force);
        }

        if((verlet->step_ + verlet->step_rst_) % verlet->mdp.md_restartfreq == 0)
        {
            verlet->ucell.update_vel(verlet->vel);
            std::stringstream file;
            file << GlobalV::global_stru_dir << "STRU_MD_" << verlet->step_ + verlet->step_rst_;
#ifdef __LCAO
            verlet->ucell.print_stru_file(GlobalC::ORB, file.str(), 1, 1);
#else
            verlet->ucell.print_stru_file(file.str(), 1, 1);
#endif
            verlet->write_restart();
        }

        verlet->step_++;
    }

    if (GlobalC::pot.out_pot == 2)
    {
        std::stringstream ssp;
        std::stringstream ssp_ave;
        ssp << GlobalV::global_out_dir << "ElecStaticPot";
        ssp_ave << GlobalV::global_out_dir << "ElecStaticPot_AVE";
        GlobalC::pot.write_elecstat_pot(ssp.str(), ssp_ave.str(),GlobalC::rhopw); //output 'Hartree + local pseudopot'
    }

    delete verlet;
    ModuleBase::timer::tick("Run_MD_PW", "md_ions_pw");
    return;
}
