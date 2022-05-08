#include "LOOP_ions.h"
#include "../src_pw/global.h"
#include "../module_orbital/parallel_orbitals.h"
#include "../src_pdiag/pdiag_double.h"
#include "../module_base/global_function.h"
#include "../src_io/write_HS.h"
#include "../src_io/print_info.h"
#include "../src_io/cal_r_overlap_R.h"
#include "../src_ions/variable_cell.h" // mohan add 2021-02-01
#include "../src_ri/exx_abfs.h"
#include "../src_ri/exx_opt_orb.h"
#include "ELEC_scf.h"
#include "../module_neighbor/sltk_atom_arrange.h"
#include "../src_pw/vdwd2.h"
#include "../src_pw/vdwd3.h"
#include "../src_pw/vdwd2_parameters.h"
#include "../src_pw/vdwd3_parameters.h"
#include "dmft.h"
#include "src_lcao/LCAO_matrix.h"
#ifdef __DEEPKS
#include "../module_deepks/LCAO_deepks.h"    //caoyu add 2021-07-26
#endif

LOOP_ions::LOOP_ions()
{}

LOOP_ions::~LOOP_ions()
{}

void LOOP_ions::opt_ions(ModuleESolver::ESolver *p_esolver)
{
    ModuleBase::TITLE("LOOP_ions","opt_ions");
    ModuleBase::timer::tick("LOOP_ions","opt_ions");

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

    // Geometry optimization algorithm setup.
    if(GlobalV::CAL_FORCE)
    {
        //Ions_Move_Methods
        IMM.allocate();
        //Charge_Extrapolation
        CE.allocate_ions();
    }

    // pengfei Li 2018-05-14
    if(GlobalV::CAL_STRESS)
    {
        // allocate arrays related to changes of lattice vectors
        LCM.allocate();
    }

    this->istep = 1;
    int force_step = 1;
    int stress_step = 1;
    bool stop = false;
    while(istep <= GlobalV::RELAX_NMAX && !stop)
    {
        time_t estart = time(NULL);

		// xiaohui add "m" option, 2015-09-16
        if(GlobalV::OUT_LEVEL=="ie" || GlobalV::OUT_LEVEL=="m")
        {
            Print_Info::print_screen(stress_step, force_step, istep);
/*          std::cout << " ---------------------------------------------------------" << std::endl;
            if(GlobalV::CALCULATION=="relax")
            {
                std::cout << " RELAX IONS : " << istep << std::endl;
            }
            else if(GlobalV::CALCULATION=="cell-relax")
            {
                std::cout << " RELAX CELL : " << stress_step << std::endl;
                std::cout << " RELAX IONS : " << force_step << " (in total: " << istep << ")" << std::endl;
            }
            else if(GlobalV::CALCULATION=="scf")
            {
                std::cout << " SELF-CONSISTENT : " << std::endl;
            }
            else if(GlobalV::CALCULATION=="nscf")
            {
                std::cout << " NONSELF-CONSISTENT : " << std::endl;
            }
            std::cout << " ---------------------------------------------------------" << std::endl;

            GlobalV::ofs_running << " ---------------------------------------------------------" << std::endl;
            if(GlobalV::CALCULATION=="relax")
            {
                GlobalV::ofs_running << " RELAX IONS : " << istep << std::endl;
            }
            else if(GlobalV::CALCULATION=="cell-relax")
            {
                GlobalV::ofs_running << " RELAX CELL : " << stress_step << std::endl;
                GlobalV::ofs_running << " RELAX IONS : " << force_step << " (in total: " << istep << ")" << std::endl;
            }
            else if(GlobalV::CALCULATION=="scf")
            {
                GlobalV::ofs_running << " SELF-CONSISTENT" << std::endl;
            }
            else if(GlobalV::CALCULATION=="nscf")
            {
                GlobalV::ofs_running << " NONSELF-CONSISTENT" << std::endl;
            }
            GlobalV::ofs_running << " ---------------------------------------------------------" << std::endl;*/
        }

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

		// solve electronic structures in terms of LCAO
        // mohan add 2021-02-09
        p_esolver->Run(this->istep-1, GlobalC::ucell);

		time_t eend = time(NULL);

		//for second-order extrapolation
        if(GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax")
        {
            CE.update_all_pos(GlobalC::ucell);
        }

		// PLEASE design a proper interface to output potentials,
		// not only electrostatic potential but also others
		// mohan add 2021-03-25
		// we need to have a proper
        if(GlobalC::pot.out_pot == 2)
        {
            std::stringstream ssp;
            std::stringstream ssp_ave;
            ssp << GlobalV::global_out_dir << "ElecStaticPot";
            ssp_ave << GlobalV::global_out_dir << "ElecStaticPot_AVE";
            GlobalC::pot.write_elecstat_pot(ssp.str(), ssp_ave.str()); //output 'Hartree + local pseudopot'
        }

        time_t fstart = time(NULL);
        if (GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax")
        {
            stop = this->force_stress(istep, force_step, stress_step, p_esolver);
        }
        time_t fend = time(NULL);

		// PLEASE move the details of CE to other places
		// mohan add 2021-03-25
        //xiaohui add 2014-07-07, for second-order extrapolation
        if(GlobalV::CAL_FORCE)
        {
            CE.save_pos_next(GlobalC::ucell);
        }

        if(GlobalV::OUT_LEVEL=="i")
        {
            double etime_min = difftime(eend, estart)/60.0;
            double ftime_min = difftime(fend, fstart)/60.0;
            std::stringstream ss;
            ss << GlobalV::RELAX_METHOD << istep;

            std::cout << std::setiosflags(ios::scientific)
            << " " << std::setw(7) << ss.str()
            << std::setw(5) << ELEC_scf::iter
            << std::setw(18) << std::setprecision(6) << GlobalC::en.etot * ModuleBase::Ry_to_eV;

            std::cout << std::setprecision(2) << std::setiosflags(ios::scientific)
            << std::setw(10) << IMM.get_ediff() * ModuleBase::Ry_to_eV * 1000
            << std::setw(10) << IMM.get_largest_grad() * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A;
            //<< std::setw(12) << IMM.get_trust_radius();

            std::cout << std::resetiosflags(ios::scientific)
//            << std::setw(8) << IMM.get_update_iter()
            << std::setprecision(2) << std::setw(10) << etime_min + ftime_min;
            std::cout << std::endl;
        }

//#ifdef __MPI
//    MPI_Barrier(MPI_COMM_WORLD);
//    for (int i=0;i<GlobalC::ucell.ntype;i++)
//    {
//        GlobalC::ucell.atoms[i].bcast_atom(); // bcast tau array
//    }
//#endif

        ++istep;
    }

    if(GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax")
    {
        GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
        GlobalV::ofs_running << std::setprecision(16);
        GlobalV::ofs_running << " !FINAL_ETOT_IS " << GlobalC::en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
        GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;

    }

    p_esolver->postprocess();

    ModuleBase::timer::tick("LOOP_ions", "opt_ions");
    return;
}


bool LOOP_ions::force_stress(
	const int &istep,
	int &force_step,
    int& stress_step,
    ModuleESolver::ESolver* p_esolver)
{
    ModuleBase::TITLE("LOOP_ions","force_stress");

    if(!GlobalV::CAL_FORCE && !GlobalV::CAL_STRESS)
    {
        return 1;
    }
    ModuleBase::timer::tick("LOOP_ions","force_stress");

	// set force matrix
	ModuleBase::matrix fcs;
	// set stress matrix
	ModuleBase::matrix scs;

    p_esolver->cal_Force(fcs);
    p_esolver->cal_Stress(scs);

	//--------------------------------------------------
	// only forces are needed, no stresses are needed
	//--------------------------------------------------
    if(GlobalV::CAL_FORCE && !GlobalV::CAL_STRESS)
    {

#ifdef __MPI
        atom_arrange::delete_vector(
			GlobalV::ofs_running,
			GlobalV::SEARCH_PBC,
			GlobalC::GridD,
			GlobalC::ucell,
			GlobalV::SEARCH_RADIUS,
			GlobalV::test_atom_input);
#endif

        if(GlobalV::CALCULATION=="relax")
        {
            IMM.cal_movement(istep, istep, fcs, GlobalC::en.etot);

            if(IMM.get_converged() || (istep==GlobalV::RELAX_NMAX))
            {
                ModuleBase::timer::tick("LOOP_ions","force_stress");
                return 1; // 1 means converged
            }
            else // ions are not converged
            {
                CE.update_istep(istep);
                CE.extrapolate_charge();

                if(GlobalC::pot.chg_extrap=="dm")
                {
                }
                else
                {
                    GlobalC::pot.init_pot( istep, GlobalC::pw.strucFac );
                }
            }
            ModuleBase::timer::tick("LOOP_ions","force_stress");
            return 0;
        }
        else
        {
            ModuleBase::timer::tick("LOOP_ions","force_stress");
            return 1;
        }

        // mohan update 2013-04-11
        // setup the structure factor
        // and do the density extraploation.
        // for both ionic iteration and
        // force calculations.

        //xiaohui modify 2014-08-09
        //GlobalC::pw.setup_structure_factor();

        // charge extrapolation if istep>0.
        //xiaohui modify 2014-08-09
        //CE.extrapolate_charge();

/*xiaohui modify 2014-08-09
        if(GlobalC::pot.chg_extrap==4)
        {
            // done after grid technique.
        }
        else
        {
            GlobalC::pot.init_pot( istep );
        }
xiaohui modify 2014-08-09*/
    }

//    static bool converged_force = false;
    static bool converged_stress = false;

    if(!GlobalV::CAL_FORCE&&GlobalV::CAL_STRESS)
    {

#ifdef __MPI
		atom_arrange::delete_vector(
			GlobalV::ofs_running,
			GlobalV::SEARCH_PBC,
			GlobalC::GridD,
			GlobalC::ucell,
			GlobalV::SEARCH_RADIUS,
			GlobalV::test_atom_input);
#endif
		if(GlobalV::CALCULATION=="cell-relax")
		{
           	LCM.cal_lattice_change(stress_step, scs, GlobalC::en.etot);
           	converged_stress = LCM.get_converged();
           	if(converged_stress)
           	{
                ModuleBase::timer::tick("LOOP_ions","force_stress");
               	return 1;
           	}
           	else
           	{
               	Variable_Cell::init_after_vc();
               	GlobalC::pot.init_pot(stress_step, GlobalC::pw.strucFac);

               	++stress_step;
                ModuleBase::timer::tick("LOOP_ions","force_stress");
               	return 0;
           	}
		}
        else
        {
            ModuleBase::timer::tick("LOOP_ions","force_stress");
            return 1;
        }
	}

    if(GlobalV::CAL_FORCE&&GlobalV::CAL_STRESS)
    {
        atom_arrange::delete_vector(
			GlobalV::ofs_running,
			GlobalV::SEARCH_PBC,
			GlobalC::GridD,
			GlobalC::ucell,
			GlobalV::SEARCH_RADIUS,
			GlobalV::test_atom_input);

        if(GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax")
        {
            IMM.cal_movement(istep, force_step, fcs, GlobalC::en.etot);

            if(IMM.get_converged())
            {
                force_step = 1;


			    if(GlobalV::CALCULATION=="cell-relax")
			    {
            	    LCM.cal_lattice_change(stress_step, scs, GlobalC::en.etot);
            	    converged_stress = LCM.get_converged();
            	    if(converged_stress)
            	    {
                        ModuleBase::timer::tick("LOOP_ions","force_stress");
                	    return 1;
            	    }
            	    else
            	    {
                	    Variable_Cell::init_after_vc();
                	    GlobalC::pot.init_pot(stress_step, GlobalC::pw.strucFac);

                	    ++stress_step;
                        ModuleBase::timer::tick("LOOP_ions","force_stress");
                	    return 0;
                    }
                }
                else
                {
                    ModuleBase::timer::tick("LOOP_ions","force_stress");
                    return 1;
                }

            }
            else
            {
                CE.update_istep(force_step);
                CE.extrapolate_charge();

                if(GlobalC::pot.chg_extrap=="dm")
                {
                }
                else
                {
                    GlobalC::pot.init_pot( istep, GlobalC::pw.strucFac );
                }
                ++force_step;
                ModuleBase::timer::tick("LOOP_ions","force_stress");
                return 0;
            }
        }
        else
        {
            ModuleBase::timer::tick("LOOP_ions","force_stress");
            return 1;
        }
    }
    ModuleBase::timer::tick("LOOP_ions","force_stress");
    return 0;
}


