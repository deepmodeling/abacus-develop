#include "run_md_lcao.h"
#include "LOOP_elec.h"
#include "FORCE_STRESS.h"
#include "../src_pw/global.h"
#include "../src_pw/vdwd2.h"
#include "../src_pw/vdwd2_parameters.h"
#include "../src_pw/vdwd3_parameters.h"
#include "../module_orbital/parallel_orbitals.h"
#include "../src_pdiag/pdiag_double.h"
#include "../src_io/write_HS.h"
#include "../src_io/cal_r_overlap_R.h"
#include "../src_io/print_info.h"
#include "../src_ions/variable_cell.h" // mohan add 2021-02-01
#include "../src_ri/exx_abfs.h"
#include "../src_ri/exx_opt_orb.h"
#include "../module_neighbor/sltk_atom_arrange.h"
#include "../module_md/MD_func.h"
#include "../module_md/FIRE.h"
#include "../module_md/NVE.h"
#include "../module_md/MSST.h"
#include "../module_md/NVT_ADS.h"
#include "../module_md/NVT_NHC.h"
#include "../module_md/Langevin.h"

Run_MD_LCAO::Run_MD_LCAO(Parallel_Orbitals &pv)
{
    cellchange = false;
    this->LM_md.ParaV = &pv;
    // * allocate H and S matrices according to computational resources
	// * set the 'trace' between local H/S and global H/S
	this->LM_md.divide_HS_in_frag(GlobalV::GAMMA_ONLY_LOCAL, pv);
}

Run_MD_LCAO::~Run_MD_LCAO(){}


void Run_MD_LCAO::opt_cell(ORB_control &orb_con, ModuleESolver::ESolver *p_esolver)
{
	ModuleBase::TITLE("Run_MD_LCAO","opt_cell");

    // Initialize the local wave functions.
    // npwx, eigenvalues, and weights
    // npwx may change according to cell change
    // this function belongs to cell LOOP
    GlobalC::wf.allocate_ekb_wg(GlobalC::kv.nks);

    // Initialize the FFT.
    // this function belongs to cell LOOP
    GlobalC::UFFT.allocate();

    // output is GlobalC::ppcell.vloc 3D local pseudopotentials
    // this function belongs to cell LOOP
    GlobalC::ppcell.init_vloc(GlobalC::pw.nggm, GlobalC::ppcell.vloc);

    // Initialize the sum of all local potentials.
    // if ion_step==0, read in/initialize the potentials
    // this function belongs to ions LOOP
    int ion_step=0;
    GlobalC::pot.init_pot(ion_step, GlobalC::pw.strucFac);

    opt_ions(p_esolver);
    orb_con.clear_after_ions(GlobalC::UOT, GlobalC::ORB, GlobalV::deepks_setorb, GlobalC::ucell.infoNL.nproj);
    
    return;
}


void Run_MD_LCAO::opt_ions(ModuleESolver::ESolver *p_esolver)
{
    ModuleBase::TITLE("Run_MD_LCAO","opt_ions"); 
    ModuleBase::timer::tick("Run_MD_LCAO","opt_ions"); 
		
    if(GlobalV::OUT_LEVEL=="i")
    {
        std::cout << std::setprecision(12);
        std::cout<< " " << std::setw(7)<< "ISTEP"
        <<std::setw(5)<< "NE"
        <<std::setw(18)<< "ETOT(eV)"
        <<std::setw(10)<< "dE(meV)"
        <<std::setw(10)<< "F(eV/A)"
        <<std::setw(10)<< "T(MIN)"
        <<std::endl;
    }

    //Charge_Extrapolation
    CE.allocate_ions();

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
    else if(INPUT.mdp.md_type == 1)
    {
        verlet = new NVT_NHC(INPUT.mdp, GlobalC::ucell);
    }
    else if(INPUT.mdp.md_type == 2)
    {
        verlet = new Langevin(INPUT.mdp, GlobalC::ucell);
    }
    else if(INPUT.mdp.md_type == 3)
    {
        verlet = new NVT_ADS(INPUT.mdp, GlobalC::ucell);
    }
    else if(INPUT.mdp.md_type == 4)
    {
        verlet = new MSST(INPUT.mdp, GlobalC::ucell); 
        cellchange = true;
    }

    // md cycle
    while ((verlet->step_ + verlet->step_rst_) <= GlobalV::MD_NSTEP && !verlet->stop)
    {
        if(verlet->step_ == 0)
        {
            MD_func::ParaV = this->LM_md.ParaV;
            verlet->setup(p_esolver);
        }
        else
        {
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

            CE.extrapolate_charge();
            CE.save_pos_next(GlobalC::ucell);

            if(cellchange)
            {
                Variable_Cell::init_after_vc();
            }

            // reset local potential
            GlobalC::pot.init_pot(verlet->step_, GlobalC::pw.strucFac);

            // update force and virial due to the update of atom positions
            MD_func::force_virial(p_esolver, verlet->step_, verlet->mdp, verlet->ucell, verlet->potential, verlet->force, verlet->virial);

            verlet->second_half();

            MD_func::kinetic_stress(verlet->ucell, verlet->vel, verlet->allmass, verlet->kinetic, verlet->stress);

            verlet->stress += verlet->virial;
        }

        if((verlet->step_ + verlet->step_rst_) % verlet->mdp.md_dumpfreq == 0)
        {
            Print_Info::print_screen(0, 0, verlet->step_ + verlet->step_rst_);
            verlet->outputMD();

            MD_func::MDdump(verlet->step_ + verlet->step_rst_, verlet->ucell, verlet->virial, verlet->force);
        }

        if((verlet->step_ + verlet->step_rst_) % verlet->mdp.md_restartfreq == 0)
        {
            verlet->ucell.update_vel(verlet->vel);
            std::stringstream file;
            file << GlobalV::global_out_dir << "STRU_MD_" << verlet->step_ + verlet->step_rst_;
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
        GlobalC::pot.write_elecstat_pot(ssp.str(), ssp_ave.str()); //output 'Hartree + local pseudopot'
    }

    GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << GlobalC::en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl; 
    GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;

	// mohan update 2021-02-10

    ModuleBase::timer::tick("Run_MD_LCAO","opt_ions"); 
    return;
}

void Run_MD_LCAO::md_force_virial(
    ModuleESolver::ESolver *p_esolver,
    const int &istep,
    const int& numIon, 
    double &potential, 
    ModuleBase::Vector3<double>* force, 
    ModuleBase::matrix& virial)
{
    //----------------------------------------------------------
	// about vdw, jiyy add vdwd3 and linpz add vdwd2
	//----------------------------------------------------------	
	if(INPUT.vdw_method=="d2")
	{
		// setup vdwd2 parameters
		GlobalC::vdwd2_para.initial_parameters(INPUT);
	    GlobalC::vdwd2_para.initset(GlobalC::ucell);
    }
    if(INPUT.vdw_method=="d3_0" || INPUT.vdw_method=="d3_bj")
    {
        GlobalC::vdwd3_para.initial_parameters(INPUT);
    }
    // Peize Lin add 2014.04.04, update 2021.03.09
    if(GlobalC::vdwd2_para.flag_vdwd2)
    {
        Vdwd2 vdwd2(GlobalC::ucell,GlobalC::vdwd2_para);
        vdwd2.cal_energy();
        GlobalC::en.evdw = vdwd2.get_energy();
    }
    // jiyy add 2019-05-18, update 2021.05.02
    else if(GlobalC::vdwd3_para.flag_vdwd3)
    {
        Vdwd3 vdwd3(GlobalC::ucell,GlobalC::vdwd3_para);
        vdwd3.cal_energy();
        GlobalC::en.evdw = vdwd3.get_energy();
    }

    Local_Orbital_wfc LOWF_md;
    Local_Orbital_Charge LOC_md;
    LOC_md.ParaV = LOWF_md.ParaV = this->LM_md.ParaV;
    if (GlobalV::GAMMA_ONLY_LOCAL)
    {
        LOWF_md.wfc_gamma.resize(GlobalV::NSPIN);
	}
	else
	{
        LOWF_md.wfc_k.resize(GlobalC::kv.nks);
    }

    LOC_md.init_dm_2d();
    // solve electronic structures in terms of LCAO
    // mohan add 2021-02-09
    LCAO_Hamilt UHM_md;
    UHM_md.genH.LM = UHM_md.LM = &this->LM_md;
    
    Record_adj RA_md;

    LOOP_elec LOE;
    LOE.solve_elec_stru(istep + 1, RA_md, LOC_md, LOWF_md, UHM_md);

    //to call the force of each atom
	ModuleBase::matrix fcs;//temp force matrix
	Force_Stress_LCAO FSL(RA_md);
    FSL.getForceStress(GlobalV::CAL_FORCE, GlobalV::CAL_STRESS,
        GlobalV::TEST_FORCE, GlobalV::TEST_STRESS,
        LOC_md, LOWF_md, UHM_md, fcs, virial);
    RA_md.delete_grid();
	for(int ion=0; ion<numIon; ++ion)
    {
		force[ion].x = fcs(ion, 0)/2.0;
		force[ion].y = fcs(ion, 1)/2.0;
		force[ion].z = fcs(ion, 2)/2.0;
	}

    virial = 0.5 * virial;

    potential = GlobalC::en.etot/2;

#ifdef __MPI //2015-10-01, xiaohui
	atom_arrange::delete_vector(
		GlobalV::ofs_running, 
		GlobalV::SEARCH_PBC, 
		GlobalC::GridD, 
		GlobalC::ucell, 
		GlobalV::SEARCH_RADIUS, 
		GlobalV::test_atom_input);
#endif //2015-10-01, xiaohui
}