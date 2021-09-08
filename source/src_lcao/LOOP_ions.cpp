#include "LOOP_ions.h"
#include "../src_pw/global.h"
#include "../src_parallel/parallel_orbitals.h"
#include "../src_pdiag/pdiag_double.h"
#include "LCAO_nnr.h"
#include "FORCE_STRESS.h"
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
#ifdef __DEEPKS
#include "LCAO_descriptor.h"    //caoyu add 2021-07-26
#endif

LOOP_ions::LOOP_ions()
{}

LOOP_ions::~LOOP_ions()
{}

void LOOP_ions::opt_ions(void)
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
    if(GlobalV::FORCE)
    {
        //Ions_Move_Methods
        IMM.allocate();
        //Charge_Extrapolation
        CE.allocate_ions();
    }

    // pengfei Li 2018-05-14
    if(GlobalV::STRESS)
    {
        // allocate arrays related to changes of lattice vectors
        LCM.allocate();
    }



    this->istep = 1;
    int force_step = 1;
    int stress_step = 1;
    bool stop = false;
    while(istep <= GlobalV::NSTEP && !stop)
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
		LOE.solve_elec_stru(this->istep);


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
        if(GlobalC::pot.out_potential == 2)
        {
            std::stringstream ssp;
            std::stringstream ssp_ave;
            ssp << GlobalV::global_out_dir << "ElecStaticPot";
            ssp_ave << GlobalV::global_out_dir << "ElecStaticPot_AVE";
            GlobalC::pot.write_elecstat_pot(ssp.str(), ssp_ave.str()); //output 'Hartree + local pseudopot'
        }

        if(GlobalC::ParaO.out_hsR)
		{
			this->output_HS_R(); //LiuXh add 2019-07-15
		}
        //caoyu add 2021-03-31
#ifdef __DEEPKS
        if (INPUT.out_descriptor)
        {
            //ld.init(ORB.get_lmax_d(), ORB.get_nchimax_d(), ucell.nat* ORB.Alpha[0].getTotal_nchi());
            //ld.build_S_descriptor(0);  //cal overlap, no need dm
            GlobalC::ld.cal_projected_DM(GlobalC::LOC.wfc_dm_2d.dm_gamma[0]);  //need dm
            GlobalC::ld.cal_descriptor();    //final descriptor
            GlobalC::ld.save_npy_d();            //libnpy needed
            if (INPUT.deepks_scf)
            {
                //ld.print_H_V_delta();   //final H_delta
                if (FORCE)
                {
                    GlobalC::ld.build_S_descriptor(1);   //for F_delta calculation
                    GlobalC::ld.cal_f_delta(GlobalC::LOC.wfc_dm_2d.dm_gamma[0]);
                    GlobalC::ld.print_F_delta();
                }
            }
        }
#endif
        time_t fstart = time(NULL);
        if (GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax")
        {
            stop = this->force_stress(istep, force_step, stress_step);
        }
        time_t fend = time(NULL);

		// PLEASE move the details of CE to other places
		// mohan add 2021-03-25
        //xiaohui add 2014-07-07, for second-order extrapolation
        if(GlobalV::FORCE)
        {
            CE.save_pos_next(GlobalC::ucell);
        }

        if(GlobalV::OUT_LEVEL=="i")
        {
            double etime_min = difftime(eend, estart)/60.0;
            double ftime_min = difftime(fend, fstart)/60.0;
            std::stringstream ss;
            ss << GlobalV::MOVE_IONS << istep;

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


    ModuleBase::timer::tick("LOOP_ions","opt_ions"); 
    return;
}


bool LOOP_ions::force_stress(
	const int &istep,
	int &force_step,
	int &stress_step)
{
    ModuleBase::TITLE("LOOP_ions","force_stress");

    if(!GlobalV::FORCE && !GlobalV::STRESS)
    {
        return 1;
    }
    ModuleBase::timer::tick("LOOP_ions","force_stress");

	// set force matrix
	ModuleBase::matrix fcs;
	// set stress matrix
	ModuleBase::matrix scs;
	Force_Stress_LCAO FSL;
	FSL.allocate ();
	FSL.getForceStress(GlobalV::FORCE, GlobalV::STRESS, GlobalV::TEST_FORCE, GlobalV::TEST_STRESS, fcs, scs);

	//--------------------------------------------------
	// only forces are needed, no stresses are needed
	//--------------------------------------------------
    if(GlobalV::FORCE && !GlobalV::STRESS)
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

            if(IMM.get_converged() || (istep==GlobalV::NSTEP))
            {
                return 1; // 1 means converged
            }
            else // ions are not converged
            {
                CE.update_istep(istep);
                CE.extrapolate_charge();

                if(GlobalC::pot.extra_pot=="dm")
                {
                }
                else
                {
                    GlobalC::pot.init_pot( istep, GlobalC::pw.strucFac );
                }
            }
            return 0;
        }
        else
        {
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
        if(GlobalC::pot.extra_pot==4)
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

    if(!GlobalV::FORCE&&GlobalV::STRESS)
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
               	return 1;
           	}
           	else
           	{
               	Variable_Cell::init_after_vc();
               	GlobalC::pot.init_pot(stress_step, GlobalC::pw.strucFac);

               	++stress_step;
               	return 0;
           	}
		}
        else
        {
            return 1;
        }
	}

    if(GlobalV::FORCE&&GlobalV::STRESS)
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
                	    return 1;
            	    }
            	    else
            	    {
                	    Variable_Cell::init_after_vc();
                	    GlobalC::pot.init_pot(stress_step, GlobalC::pw.strucFac);

                	    ++stress_step;
                	    return 0;
                    }
                }
                else
                {
                    return 1;
                }

            }
            else
            {
                CE.update_istep(force_step);
                CE.extrapolate_charge();

                if(GlobalC::pot.extra_pot=="dm")
                {
                }
                else
                {
                    GlobalC::pot.init_pot( istep, GlobalC::pw.strucFac );
                }
                ++force_step;
                return 0;
            }
        }
        else
        {
            return 1;
        }
    }

    return 0;

    ModuleBase::timer::tick("LOOP_ions","force_stress");
}

void LOOP_ions::final_scf(void)
{
    ModuleBase::TITLE("LOOP_ions","final_scf");

    GlobalV::FINAL_SCF = true;

    Variable_Cell::final_calculation_after_vc();

	//------------------------------------------------------------------
	// THIS PART IS THE SAME AS LOOP_elec::set_matrix_grid
    GlobalV::SEARCH_RADIUS = atom_arrange::set_sr_NL(
		GlobalV::ofs_running,
		GlobalV::OUT_LEVEL,
		GlobalC::ORB.get_rcutmax_Phi(),
		GlobalC::ucell.infoNL.get_rcutmax_Beta(),
		GlobalV::GAMMA_ONLY_LOCAL);

    atom_arrange::search(
		GlobalV::SEARCH_PBC,
		GlobalV::ofs_running,
		GlobalC::GridD,
		GlobalC::ucell,
		GlobalV::SEARCH_RADIUS,
		GlobalV::test_atom_input);

    GlobalC::GridT.set_pbc_grid(
        GlobalC::pw.ncx, GlobalC::pw.ncy, GlobalC::pw.ncz,
        GlobalC::pw.bx, GlobalC::pw.by, GlobalC::pw.bz,
        GlobalC::pw.nbx, GlobalC::pw.nby, GlobalC::pw.nbz,
        GlobalC::pw.nbxx, GlobalC::pw.nbzp_start, GlobalC::pw.nbzp);

    // (2) If k point is used here, allocate HlocR after atom_arrange.
    if(!GlobalV::GAMMA_ONLY_LOCAL)
    {
        // For each atom, calculate the adjacent atoms in different cells
        // and allocate the space for H(R) and S(R).
        GlobalC::LNNR.cal_nnr();
        GlobalC::LM.allocate_HS_R(GlobalC::LNNR.nnr);

		// need to first calculae lgd.
        // using GlobalC::GridT.init.
        GlobalC::LNNR.cal_nnrg(GlobalC::GridT);
    }
	//------------------------------------------------------------------




	//------------------------------------------------------------------
	// THIS PART IS THE SAME AS LOOP_elec::before_solver
    // (4) set the augmented orbitals index.
    // after ParaO and GridT,
    // this information is used to calculate
    // the force.
    GlobalC::LOWF.set_trace_aug(GlobalC::GridT);

	GlobalC::LOC.allocate_dm_wfc(GlobalC::GridT);

    GlobalC::UHM.set_lcao_matrices();
	//------------------------------------------------------------------




    if(GlobalC::vdwd2_para.flag_vdwd2)							//Peize Lin add 2014-04-04, update 2021-03-09
    {
        Vdwd2 vdwd2(GlobalC::ucell,GlobalC::vdwd2_para);
        vdwd2.cal_energy();
        GlobalC::en.evdw = vdwd2.get_energy();
    }
	else if(GlobalC::vdwd3_para.flag_vdwd3)							//jiyy add 2019-05-18, update 2021-05-02
    {
        Vdwd3 vdwd3(GlobalC::ucell,GlobalC::vdwd3_para);
        vdwd3.cal_energy();
        GlobalC::en.evdw = vdwd3.get_energy();
    }



	ELEC_scf es;
	es.scf(0);

    if(GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax")
    {
        GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
        GlobalV::ofs_running << std::setprecision(16);
        GlobalV::ofs_running << " !FINAL_ETOT_IS " << GlobalC::en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
        GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;
    }

    return;
}
