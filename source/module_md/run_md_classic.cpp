#include "run_md_classic.h"
#include "MD_func.h"
#include "NVE.h"
#include "MSST.h"
#include "FIRE.h"
#include "NVT_ADS.h"
#include "NVT_NHC.h"
#include "Langevin.h"
#include "../input.h"
#include "../src_io/print_info.h"
#include "../module_base/timer.h"

Run_MD_CLASSIC::Run_MD_CLASSIC(){}

Run_MD_CLASSIC::~Run_MD_CLASSIC(){}

void Run_MD_CLASSIC::classic_md_line(UnitCell_pseudo &unit_in, ModuleESolver::ESolver *p_esolver)
{
	ModuleBase::TITLE("Run_MD_CLASSIC", "classic_md_line");
    ModuleBase::timer::tick("Run_MD_CLASSIC", "classic_md_line");

    // determine the md_type
    Verlet *verlet;
    if(INPUT.mdp.md_type == -1)
    {
        verlet = new FIRE(INPUT.mdp, unit_in); 
    }
    else if(INPUT.mdp.md_type == 0)
    {
        verlet = new NVE(INPUT.mdp, unit_in); 
    }
    else if(INPUT.mdp.md_type == 1)
    {
        verlet = new NVT_NHC(INPUT.mdp, unit_in);
    }
    else if(INPUT.mdp.md_type == 2)
    {
        verlet = new Langevin(INPUT.mdp, unit_in);
    }
    else if(INPUT.mdp.md_type == 3)
    {
        verlet = new NVT_ADS(INPUT.mdp, unit_in);
    }
    else if(INPUT.mdp.md_type == 4)
    {
        verlet = new MSST(INPUT.mdp, unit_in); 
    }

    // md cycle
    while ( (verlet->step_ + verlet->step_rst_) <= GlobalV::MD_NSTEP && !verlet->stop )
    {
        if(verlet->step_ == 0)
        {
            verlet->setup(p_esolver);
        }
        else
        {
            Print_Info::print_screen(0, 0, verlet->step_ + verlet->step_rst_);
            verlet->first_half();

            // update force and virial due to the update of atom positions
            MD_func::force_virial(p_esolver, verlet->step_, verlet->mdp, verlet->ucell, verlet->potential, verlet->force, verlet->virial);

            verlet->second_half();

            MD_func::kinetic_stress(verlet->ucell, verlet->vel, verlet->allmass, verlet->kinetic, verlet->stress);

            verlet->stress +=  verlet->virial;
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

	GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << verlet->potential*ModuleBase::Hartree_to_eV << " eV" << std::endl; 
    GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;

    delete verlet;
    ModuleBase::timer::tick("Run_MD_CLASSIC", "classic_md_line");
    return;
}