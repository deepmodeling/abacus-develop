#include "LOOP_ions.h"
#include "../src_pw/global.h"
#include "../src_parallel/parallel_orbitals.h"
#include "src_pdiag/pdiag_double.h"
#include "LCAO_nnr.h"
#include "FORCE_STRESS.h"
#include "../module_base/global_function.h"
#include "../src_io/write_HS.h"
#include "../src_io/cal_r_overlap_R.h"
#include "../src_ions/variable_cell.h" // mohan add 2021-02-01
#include "src_ri/exx_abfs.h"
#include "src_ri/exx_opt_orb.h"
#include "ELEC_scf.h"
#include "../module_neighbor/sltk_atom_arrange.h"
#include "../src_pw/vdwd2.h"
#include "../src_pw/vdwd3.h"
#include "../src_pw/vdwd2_parameters.h"
#include "../src_pw/vdwd3_parameters.h"
#include "LCAO_descriptor.h"

LOOP_ions::LOOP_ions()
{}

LOOP_ions::~LOOP_ions() 
{}

void LOOP_ions::opt_ions(void)
{
    TITLE("LOOP_ions","opt_ions"); 
    timer::tick("LOOP_ions","opt_ions",'B'); 
		
    if(OUT_LEVEL=="i")
    {
        cout << setprecision(12);
        cout<< " " << setw(7)<< "ISTEP"
        <<setw(5)<< "NE"
        <<setw(18)<< "ETOT(eV)"
        <<setw(10)<< "dE(meV)"
        <<setw(10)<< "F(eV/A)"
        <<setw(10)<< "T(MIN)"
        <<endl;
    }

    // Geometry optimization algorithm setup.
    if(FORCE)
    {
        //Ions_Move_Methods
        IMM.allocate();
        //Charge_Extrapolation
        CE.allocate_ions();
    }

    // pengfei Li 2018-05-14
    if(STRESS)
    {
        // allocate arrays related to changes of lattice vectors
        LCM.allocate();
    }



    this->istep = 1;
    int force_step = 1;
    int stress_step = 1;
    bool stop = false;
    while(istep <= NSTEP && !stop)
    {
        time_t estart = time(NULL);
	
		// xiaohui add "m" option, 2015-09-16
        if(OUT_LEVEL=="ie" || OUT_LEVEL=="m")
        {
            cout << " ---------------------------------------------------------" << endl;
            if(CALCULATION=="relax") 
            {
                cout << " RELAX IONS : " << istep << endl;
            }
            else if(CALCULATION=="cell-relax")
            {
                cout << " RELAX CELL : " << stress_step << endl;
                cout << " RELAX IONS : " << force_step << " (in total: " << istep << ")" << endl;
            }
            else if(CALCULATION=="scf")
            {
                cout << " SELF-CONSISTENT : " << endl;
            }
            cout << " ---------------------------------------------------------" << endl;

            ofs_running << " ---------------------------------------------------------" << endl;
            if(CALCULATION=="relax")
            {
                ofs_running << " RELAX IONS : " << istep << endl;
                ofs_running << " ---------------------------------------------------------" << endl;
            }
            else if(CALCULATION=="cell-relax")
            {
                ofs_running << " RELAX CELL : " << stress_step << endl;
                ofs_running << " RELAX IONS : " << force_step << " (in total: " << istep << ")" << endl;
                ofs_running << " ---------------------------------------------------------" << endl;
            }
            else if(CALCULATION=="scf")
            {
                ofs_running << " SELF-CONSISTENT" << endl;
                ofs_running << " ---------------------------------------------------------" << endl;
            }
        }

    //----------------------------------------------------------
    // about vdw, jiyy add vdwd3 and linpz add vdwd2
    //----------------------------------------------------------	
        if(INPUT.vdw_method=="d2")
        {
            // setup vdwd2 parameters
	        vdwd2_para.initial_parameters(INPUT);
	        vdwd2_para.initset(ucell);
        }
        if(INPUT.vdw_method=="d3_0" || INPUT.vdw_method=="d3_bj")
        {
            vdwd3_para.initial_parameters(INPUT);
        }
        // Peize Lin add 2014.04.04, update 2021.03.09
        if(vdwd2_para.flag_vdwd2)
        {
            Vdwd2 vdwd2(ucell,vdwd2_para);
            vdwd2.cal_energy();
            en.evdw = vdwd2.get_energy();
        }
        // jiyy add 2019-05-18, update 2021.05.02
        else if(vdwd3_para.flag_vdwd3)
        {
            Vdwd3 vdwd3(ucell,vdwd3_para);
            vdwd3.cal_energy();
            en.evdw = vdwd3.get_energy();
        }


		// solve electronic structures in terms of LCAO
		// mohan add 2021-02-09
		LOE.solve_elec_stru(this->istep);

		
		time_t eend = time(NULL);

		//for second-order extrapolation
        if(CALCULATION=="relax" || CALCULATION=="cell-relax")
        {
            CE.update_all_pos(ucell);
        }

		// PLEASE design a proper interface to output potentials,
		// not only electrostatic potential but also others
		// mohan add 2021-03-25
		// we need to have a proper
        if(pot.out_potential == 2)
        {
            stringstream ssp;
            stringstream ssp_ave;
            ssp << global_out_dir << "ElecStaticPot";
            ssp_ave << global_out_dir << "ElecStaticPot_AVE";
            pot.write_elecstat_pot(ssp.str(), ssp_ave.str()); //output 'Hartree + local pseudopot'
        }

        if(ParaO.out_hsR) 
		{
			this->output_HS_R(); //LiuXh add 2019-07-15
		}
        //caoyu add 2021-03-31
        if (INPUT.out_descriptor)
        {
            LCAO_Descriptor ld(ORB.get_lmax_d(), ucell.nat* ORB.Alpha[0].getTotal_nchi());
            ld.build_S_descriptor(0);  //derivation not needed yet
            ld.cal_projected_DM();
            ld.cal_descriptor();
            if (INPUT.deepks_scf)
            {
                ld.build_S_descriptor(1);   //for F_delta calculation
                ld.cal_v_delta("cmodel.pt");
                ld.print_H_V_delta();
                ld.cal_f_delta(LOC.wfc_dm_2d.dm_gamma[0]);
                ld.print_F_delta();
            }
        }

        time_t fstart = time(NULL);
        if (CALCULATION=="scf" || CALCULATION=="relax" || CALCULATION=="cell-relax")
        {
            stop = this->force_stress(istep, force_step, stress_step);
        }            
        time_t fend = time(NULL);

		// PLEASE move the details of CE to other places
		// mohan add 2021-03-25
        //xiaohui add 2014-07-07, for second-order extrapolation
        if(FORCE)
        {
            CE.save_pos_next(ucell);
        }
		
        if(OUT_LEVEL=="i")
        {
            double etime_min = difftime(eend, estart)/60.0;
            double ftime_min = difftime(fend, fstart)/60.0;
            stringstream ss;
            ss << MOVE_IONS << istep;

            cout << setiosflags(ios::scientific)
            << " " << setw(7) << ss.str()
            << setw(5) << ELEC_scf::iter
            << setw(18) << setprecision(6) << en.etot * Ry_to_eV;

            cout << setprecision(2) << setiosflags(ios::scientific)
            << setw(10) << IMM.get_ediff() * Ry_to_eV * 1000
            << setw(10) << IMM.get_largest_grad() * Ry_to_eV / BOHR_TO_A;
            //<< setw(12) << IMM.get_trust_radius();

            cout << resetiosflags(ios::scientific)
//            << setw(8) << IMM.get_update_iter()
            << setprecision(2) << setw(10) << etime_min + ftime_min;
            cout << endl;
        }

//#ifdef __MPI
//    MPI_Barrier(MPI_COMM_WORLD);
//    for (int i=0;i<ucell.ntype;i++)
//    {
//        ucell.atoms[i].bcast_atom(); // bcast tau array
//    }
//#endif

        ++istep;
    }

    if(CALCULATION=="scf" || CALCULATION=="relax" || CALCULATION=="cell-relax")
    {
        ofs_running << "\n\n --------------------------------------------" << endl;
        ofs_running << setprecision(16);
        ofs_running << " !FINAL_ETOT_IS " << en.etot * Ry_to_eV << " eV" << endl; 
        ofs_running << " --------------------------------------------\n\n" << endl;

    }


    timer::tick("LOOP_ions","opt_ions",'B'); 
    return;
}


bool LOOP_ions::force_stress(
	const int &istep, 
	int &force_step, 
	int &stress_step)
{
    TITLE("LOOP_ions","force_stress");

    if(!FORCE && !STRESS)
    {
        return 1;
    }
    timer::tick("LOOP_ions","force_stress",'D');

	// set force matrix
	matrix fcs;
	// set stress matrix
	matrix scs;
	Force_Stress_LCAO FSL;
	FSL.allocate (); 
	FSL.getForceStress(FORCE, STRESS, TEST_FORCE, TEST_STRESS, fcs, scs);

	//--------------------------------------------------
	// only forces are needed, no stresses are needed
	//--------------------------------------------------
    if(FORCE && !STRESS)
    {

#ifdef __MPI
        atom_arrange::delete_vector(
			ofs_running,
			SEARCH_PBC, 
			GridD, 
			ucell, 
			SEARCH_RADIUS, 
			test_atom_input);
#endif

        if(CALCULATION=="relax") 
        {
            IMM.cal_movement(istep, istep, fcs, en.etot);

            if(IMM.get_converged() || (istep==NSTEP))
            {
                return 1; // 1 means converged
            }
            else // ions are not converged
            {
                CE.update_istep(istep); 
                CE.extrapolate_charge();

                if(pot.extra_pot=="dm")
                {
                }
                else
                {
                    pot.init_pot( istep, pw.strucFac );
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
        //pw.setup_structure_factor();

        // charge extrapolation if istep>0.
        //xiaohui modify 2014-08-09
        //CE.extrapolate_charge();

/*xiaohui modify 2014-08-09
        if(pot.extra_pot==4)
        {
            // done after grid technique.
        }
        else
        {
            pot.init_pot( istep );
        }
xiaohui modify 2014-08-09*/
    }

//    static bool converged_force = false;
    static bool converged_stress = false;

    if(!FORCE&&STRESS)
    {

#ifdef __MPI
		atom_arrange::delete_vector(
			ofs_running,
			SEARCH_PBC, 
			GridD, 
			ucell, 
			SEARCH_RADIUS, 
			test_atom_input);
#endif
		if(CALCULATION=="cell-relax")
		{
           	LCM.cal_lattice_change(stress_step, scs, en.etot);
           	converged_stress = LCM.get_converged();
           	if(converged_stress)
           	{
               	return 1;
           	}
           	else
           	{
               	Variable_Cell::init_after_vc();
               	pot.init_pot(stress_step, pw.strucFac);

               	++stress_step;
               	return 0;
           	}
		}
        else
        {
            return 1;
        }
	}

    if(FORCE&&STRESS)
    {
        atom_arrange::delete_vector(
			ofs_running,
			SEARCH_PBC, 
			GridD, 
			ucell, 
			SEARCH_RADIUS, 
			test_atom_input);
        
        if(CALCULATION=="relax" || CALCULATION=="cell-relax")
        {
            IMM.cal_movement(istep, force_step, fcs, en.etot);

            if(IMM.get_converged())
            {
                force_step = 1;


			    if(CALCULATION=="cell-relax")
			    {
            	    LCM.cal_lattice_change(stress_step, scs, en.etot);
            	    converged_stress = LCM.get_converged();
            	    if(converged_stress)
            	    {
                	    return 1;
            	    }
            	    else
            	    {
                	    Variable_Cell::init_after_vc();
                	    pot.init_pot(stress_step, pw.strucFac);

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

                if(pot.extra_pot=="dm")
                {
                }
                else
                {
                    pot.init_pot( istep, pw.strucFac );
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

    timer::tick("LOOP_ions","force_stress",'D');
}

void LOOP_ions::final_scf(void)
{
    TITLE("LOOP_ions","final_scf"); 

    FINAL_SCF = true;

    Variable_Cell::final_calculation_after_vc();

	//------------------------------------------------------------------
	// THIS PART IS THE SAME AS LOOP_elec::set_matrix_grid
    SEARCH_RADIUS = atom_arrange::set_sr_NL(
		ofs_running,
		OUT_LEVEL,
		ORB.get_rcutmax_Phi(), 
		ORB.get_rcutmax_Beta(), 
		GAMMA_ONLY_LOCAL);

    atom_arrange::search(
		SEARCH_PBC,
		ofs_running,
		GridD, 
		ucell, 
		SEARCH_RADIUS, 
		test_atom_input);

    GridT.set_pbc_grid(
        pw.ncx, pw.ncy, pw.ncz,
        pw.bx, pw.by, pw.bz,
        pw.nbx, pw.nby, pw.nbz,
        pw.nbxx, pw.nbzp_start, pw.nbzp);

    // (2) If k point is used here, allocate HlocR after atom_arrange.
    if(!GAMMA_ONLY_LOCAL)
    {
        // For each atom, calculate the adjacent atoms in different cells 
        // and allocate the space for H(R) and S(R).
        LNNR.cal_nnr();
        LM.allocate_HS_R(LNNR.nnr);
        
		// need to first calculae lgd.
        // using GridT.init.
        LNNR.cal_nnrg(GridT);
    }
	//------------------------------------------------------------------




	//------------------------------------------------------------------
	// THIS PART IS THE SAME AS LOOP_elec::before_solver 
    // (4) set the augmented orbitals index.
    // after ParaO and GridT, 
    // this information is used to calculate
    // the force.
    LOWF.set_trace_aug(GridT);

	LOC.allocate_dm_wfc(GridT);
		
    UHM.set_lcao_matrices();
	//------------------------------------------------------------------
	



    if(vdwd2_para.flag_vdwd2)							//Peize Lin add 2014-04-04, update 2021-03-09
    {
        Vdwd2 vdwd2(ucell,vdwd2_para);
        vdwd2.cal_energy();
        en.evdw = vdwd2.get_energy();
    }
	else if(vdwd3_para.flag_vdwd3)							//jiyy add 2019-05-18, update 2021-05-02
    {
        Vdwd3 vdwd3(ucell,vdwd3_para);
        vdwd3.cal_energy();
        en.evdw = vdwd3.get_energy();
    }											  
    


	ELEC_scf es;
	es.scf(0);

    if(CALCULATION=="scf" || CALCULATION=="relax" || CALCULATION=="cell-relax")
    {
        ofs_running << "\n\n --------------------------------------------" << endl;
        ofs_running << setprecision(16);
        ofs_running << " !FINAL_ETOT_IS " << en.etot * Ry_to_eV << " eV" << endl; 
        ofs_running << " --------------------------------------------\n\n" << endl;
    }

    return;
}
