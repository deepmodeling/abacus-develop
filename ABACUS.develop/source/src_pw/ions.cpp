#include "tools.h"
#include "ions.h"
#include "forces.h"
#include "stress_pw.h"
#include "../src_pw/global.h" // use chr.
#include "vdwd2.h"
#include "vdwd3.h"
#include "../src_pw/pw_complement.h"
#include "../src_pw/pw_basis.h"
#include "../src_ions/variable_cell.h" // mohan add 2021-02-01

void Ions::opt_ions_pw(void)
{
	TITLE("Ions","opt_ions_pw");
	timer::tick("Ions","opt_ions_pw",'C');
	
	if(OUT_LEVEL=="i")
	{
		cout << setprecision(12);
    	cout<< " " << setw(7)<< "ISTEP" 
		<<setw(5)<< "NE"
		<<setw(15)<< "ETOT(eV)"
		<<setw(15)<< "EDIFF(eV)"
        <<setw(15)<< "MAX_F(eV/A)"
        <<setw(15)<< "TRADIUS(Bohr)"
		<<setw(8)<< "UPDATE"
		<<setw(11)<< "ETIME(MIN)"
		<<setw(11)<< "FTIME(MIN)"
        <<endl;
	}

	// allocation for ion movement.	
	if(FORCE)
	{
		IMM.allocate();
		CE.allocate_ions();
	}

	if(STRESS)                    // pengfei Li 2018-05-14
	{
		LCM.allocate();
	}

    this->istep = 1;
	int force_step = 1;           // pengfei Li 2018-05-14
	int stress_step = 1;
	bool stop= false;
	
    while(istep <= NSTEP && !stop)
    {
		time_t estart = time(NULL);

		if(OUT_LEVEL=="ie")
		{
	
        	cout << " -------------------------------------------" << endl;
			if(CALCULATION=="relax") //pengfei 2014-10-13
			{
        		cout << " STEP OF ION RELAXATION : " << istep << endl;
			}
            else if(CALCULATION=="cell-relax")
            {
                cout << " RELAX CELL : " << stress_step << endl;
                cout << " RELAX IONS : " << force_step << " (in total: " << istep << ")" << endl;
                cout << " ---------------------------------------------------------" << endl;
            }
			else if(CALCULATION=="scf") //add 4 lines 2015-09-06, xiaohui
			{
        			cout << " SELF-CONSISTENT : " << endl;
			}
			else if(CALCULATION=="md")
			{
        		cout << " STEP OF MOLECULAR DYNAMICS : " << istep << endl;
			}
        	cout << " -------------------------------------------" << endl;

        	ofs_running << " -------------------------------------------" << endl;
			if(CALCULATION=="relax")
			{
        		ofs_running << " STEP OF ION RELAXATION : " << istep << endl;
			}
            else if(CALCULATION=="cell-relax")
            {
                ofs_running << " RELAX CELL : " << stress_step << endl;
                ofs_running << " RELAX IONS : " << force_step << " (in total: " << istep << ")" << endl;
                ofs_running << " ---------------------------------------------------------" << endl;
            }
			else if(CALCULATION=="md")
			{
        		ofs_running << " STEP OF MOLECULAR DYNAMICS : " << istep << endl;
			}
        	ofs_running << " -------------------------------------------" << endl;
		}

			
		if(vdwd2_para.flag_vdwd2)		//Peize Lin add 2014-04-03, update 2021-03-09
		{
			Vdwd2 vdwd2(ucell,vdwd2_para);
			vdwd2.cal_energy();
			en.evdw = vdwd2.get_energy();
		}
		if(vdwd3.vdwD3)		//jiyy add 2019-05-18
		{
			vdwd3.energy();
			en.evdw = vdwd3.energy_result;
		}


		// mohan added eiter to count for the electron iteration number, 2021-01-28
		int eiter=0;		
        if (CALCULATION=="scf" || CALCULATION=="md" || CALCULATION=="relax" || CALCULATION=="cell-relax")  // pengfei 2014-10-13
        {
			if( Exx_Global::Hybrid_Type::No==exx_global.info.hybrid_type  )
			{			
				elec.self_consistent(istep-1);
				eiter = elec.iter;
			}
			else if( Exx_Global::Hybrid_Type::Generate_Matrix == exx_global.info.hybrid_type )
			{
				throw invalid_argument(TO_STRING(__FILE__)+TO_STRING(__LINE__));
			}
			else	// Peize Lin add 2019-03-09
			{
				if( exx_global.info.separate_loop )
				{
					for( size_t hybrid_step=0; hybrid_step!=exx_global.info.hybrid_step; ++hybrid_step )
					{
						elec.self_consistent(istep-1);
						eiter += elec.iter;
						if( elec.iter==1 || hybrid_step==exx_global.info.hybrid_step-1 )		// exx converge
							break;
						exx_global.info.set_xcfunc(xcf);							
						exx_lip.cal_exx();
					}						
				}
				else
				{
					elec.self_consistent(istep-1);	
					eiter += elec.iter;
					exx_global.info.set_xcfunc(xcf);
					elec.self_consistent(istep-1);
					eiter += elec.iter;
				}
			}
        }
        else if(CALCULATION=="nscf")
        {
            elec.non_self_consistent(istep-1);
			eiter = elec.iter;
        }
		// mohan added 2021-01-28, perform stochastic calculations
		else if(CALCULATION=="scf-sto" || CALCULATION=="relax-sto" || CALCULATION=="md-sto")
		{
			elec_sto.scf_stochastic(istep-1);
			eiter = elec_sto.iter;
		}
	

		if(CALCULATION=="relax"|| CALCULATION=="md" || CALCULATION=="cell-relax")
		{
			CE.update_all_pos(ucell);
		}

		if(pot.out_potential == 2)
		{
			stringstream ssp;
			stringstream ssp_ave;
			ssp << global_out_dir << "ElecStaticPot";
			ssp_ave << global_out_dir << "ElecStaticPot_AVE";
			pot.write_elecstat_pot(ssp.str(), ssp_ave.str()); //output 'Hartree + local pseudopot'
		}

		time_t eend = time(NULL);
		time_t fstart = time(NULL);


        if (CALCULATION=="scf" || CALCULATION=="relax" || CALCULATION=="cell-relax")
        {
			stop = this->force_stress(istep, force_step, stress_step);    // pengfei Li 2018-05-14
		}
		time_t fend = time(NULL);


		if(OUT_LEVEL=="i")
		{
			double etime_min = difftime(eend, estart)/60.0; 
			double ftime_min = difftime(fend, fstart)/60.0; 
			stringstream ss;
			ss << MOVE_IONS << istep;
			
			cout << " " << setw(7) << ss.str() 
			<< setw(5) << eiter 
			<< setw(15) << setprecision(6) << en.etot * Ry_to_eV 
			<< setw(15) << IMM.get_ediff() * Ry_to_eV
			<< setprecision(3)
			<< setw(15) << IMM.get_largest_grad() * Ry_to_eV / 0.529177
			<< setw(15) << IMM.get_trust_radius()
			<< setw(8) << IMM.get_update_iter()
			<< setprecision(2) << setw(11) << etime_min
			<< setw(11) << ftime_min << endl;
		}

		++istep;

    }

    if(CALCULATION=="scf" || CALCULATION=="relax" || CALCULATION=="cell-relax")
    {
        ofs_running << "\n\n --------------------------------------------" << endl;
        ofs_running << setprecision(16);
        ofs_running << " !FINAL_ETOT_IS " << en.etot * Ry_to_eV << " eV" << endl; 
        ofs_running << " --------------------------------------------\n\n" << endl;
    }


	if(OUT_LEVEL=="i")
	{
		cout << " ION DYNAMICS FINISHED :)" << endl;
	}

	timer::tick("Ions","opt_ions_pw",'C');
    return;
}


bool Ions::force_stress(const int &istep, int &force_step, int &stress_step)  // pengfei Li 2018-05-14
{
	TITLE("Ions","force_stress");

	if(!FORCE && !STRESS)
	{
		return 1;
	}


	if(FORCE&&!STRESS)
	{
		// (1) calculate the force.
		matrix force;
		Forces fcs;
		fcs.init(force);

		// (2) move the ions.
		bool converged = false;
		//if(CALCULATION=="md") //mohan add
		//{
		//	md.init_md(istep, fcs.force);
		//}
		//else	
		//{
			//IMM.cal_movement(istep, fcs.force, en.etot);
		//}
		if(CALCULATION=="relax")
		{
			IMM.cal_movement(istep, istep, force, en.etot);
			converged = IMM.get_converged();

            if(converged || (istep==NSTEP) ) 
            {
                return 1;
            }
            else
            {
                ofs_running << " Setup the structure factor in plane wave basis." << endl;
                pw.setup_structure_factor();

                ofs_running << " Setup the extrapolated charge." << endl;
                // charge extrapolation if istep>0.
                CE.extrapolate_charge();
			
                ofs_running << " Setup the Vl+Vh+Vxc according to new structure factor and new charge." << endl;
                // calculate the new potential accordint to
                // the new charge density.
                pot.init_pot( istep, pw.strucFac );

                ofs_running << " Setup the new wave functions?" << endl;
                wf.wfcinit();
            }
        }
        else
        {
            return 1;
        }
    }

	static bool converged_force = false;              // pengfe Li  2018-05-14
	static bool converged_stress = false;

	if(!FORCE&&STRESS)
	{
		Stress_PW ss;
		matrix stress;
		ss.cal_stress(stress);

		double unit_transform = 0.0;
		unit_transform = RYDBERG_SI / pow(BOHR_RADIUS_SI,3) * 1.0e-8;
		double external_stress[3] = {PRESS1,PRESS2,PRESS3};
		for(int i=0;i<3;i++)
		{
			stress(i,i) -= external_stress[i]/unit_transform;
		}
		PRESSURE = (stress(0,0)+stress(1,1)+stress(2,2))/3;
		if(CALCULATION=="cell-relax")
		{
			LCM.cal_lattice_change(stress_step, stress, en.etot);
			converged_stress = LCM.get_converged();
            //cout <<"converged_stress = "<<converged_stress<<endl;
            if(converged_stress)
            {
                return 1;
            }
            else
            {
                Variable_Cell::init_after_vc();
                pot.init_pot(stress_step, pw.strucFac); //LiuXh add 20180619
                ofs_running << " Setup the new wave functions?" << endl; //LiuXh add 20180619
                wf.wfcinit(); //LiuXh add 20180619
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
		//cout<<" istep  force_step  stress_step  converged_force  converged_stress = "<<istep<<"  "<<force_step<<"  "<<stress_step<<"  "<<converged_force<<"  "<<converged_stress<<endl;

		matrix force;
		Forces fcs;
		fcs.init(force);
                matrix stress;
                Stress_PW ss;
                ss.cal_stress(stress);
		//IMM.cal_movement(force_step, fcs.force, en.etot);
        if(CALCULATION=="relax" || CALCULATION=="cell-relax")
        {
            IMM.cal_movement(istep, force_step, force, en.etot);
            converged_force = IMM.get_converged();

            //cout<<"converged_force = "<<converged_force<<endl;
            if(converged_force)
            {
                force_step = 1;


                double unit_transform = 0.0;
                unit_transform = RYDBERG_SI / pow(BOHR_RADIUS_SI,3) * 1.0e-8;
                double external_stress[3] = {PRESS1,PRESS2,PRESS3};
                for(int i=0;i<3;i++)
                {
                    stress(i,i) -= external_stress[i]/unit_transform;
                }
				PRESSURE = (stress(0,0)+stress(1,1)+stress(2,2))/3;

                if(CALCULATION=="cell-relax")
                {
                    LCM.cal_lattice_change(stress_step, stress, en.etot);
                    converged_stress = LCM.get_converged();
                    //cout <<"converged_stress = "<<converged_stress<<endl;
                    if(converged_stress)
                    {
                        return 1;
                    }
                    else
                    {
                        Variable_Cell::init_after_vc();
                        pot.init_pot(stress_step, pw.strucFac); //LiuXh add 20180619

                        ofs_running << " Setup the new wave functions?" << endl; //LiuXh add 20180619
                        wf.wfcinit(); //LiuXh add 20180619

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
				ofs_running << " Setup the structure factor in plane wave basis." << endl;
                pw.setup_structure_factor();
				ofs_running << " Setup the extrapolated charge." << endl;
                CE.save_pos_next(ucell);
                CE.update_istep(force_step);
                CE.extrapolate_charge();
				ofs_running << " Setup the Vl+Vh+Vxc according to new structure factor and new charge." << endl;
                pot.init_pot( istep, pw.strucFac );
				ofs_running << " Setup the new wave functions?" << endl;
                wf.wfcinit();
                ++force_step;
                return 0;
            }
        }
        else
        {

            double unit_transform = 0.0;
            unit_transform = RYDBERG_SI / pow(BOHR_RADIUS_SI,3) * 1.0e-8;
            double external_stress[3] = {PRESS1,PRESS2,PRESS3};
            for(int i=0;i<3;i++)
            {
                stress(i,i) -= external_stress[i]/unit_transform;
            }
			PRESSURE = (stress(0,0)+stress(1,1)+stress(2,2))/3;
            return 1;
        }
    }

    return 0;
}


