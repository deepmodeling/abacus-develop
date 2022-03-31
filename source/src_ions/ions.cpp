#include "ions.h"
#include "../src_pw/global.h" // use chr.
#include "../src_io/print_info.h"
#include "variable_cell.h" // mohan add 2021-02-01
#include "src_io/write_wfc_realspace.h"

void Ions::opt_ions_pw(ModuleESolver::ESolver *p_esolver)
{
	ModuleBase::TITLE("Ions","opt_ions_pw");
	ModuleBase::timer::tick("Ions","opt_ions_pw");
	
	if(GlobalV::OUT_LEVEL=="i")
	{
		std::cout << std::setprecision(12);
    	std::cout<< " " << std::setw(7)<< "ISTEP" 
		<<std::setw(5)<< "NE"
		<<std::setw(15)<< "ETOT(eV)"
		<<std::setw(15)<< "EDIFF(eV)"
        <<std::setw(15)<< "MAX_F(eV/A)"
        <<std::setw(15)<< "TRADIUS(Bohr)"
		<<std::setw(8)<< "UPDATE"
		<<std::setw(11)<< "ETIME(MIN)"
		<<std::setw(11)<< "FTIME(MIN)"
        <<std::endl;
	}

	// allocation for ion movement.	
	if(GlobalV::CAL_FORCE)
	{
		IMM.allocate();
		CE.allocate_ions();
	}

	if(GlobalV::CAL_STRESS)                    // pengfei Li 2018-05-14
	{
		LCM.allocate();
	}

    this->istep = 1;
	int force_step = 1;           // pengfei Li 2018-05-14
	int stress_step = 1;
	bool stop= false;
	
    while(istep <= GlobalV::RELAX_NMAX && !stop)
    {
		time_t estart = time(NULL);

		if(GlobalV::OUT_LEVEL=="ie")
		{
			Print_Info::print_screen(stress_step, force_step, istep);
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
		if(GlobalC::vdwd2_para.flag_vdwd2)		//Peize Lin add 2014-04-03, update 2021-03-09
		{
			Vdwd2 vdwd2(GlobalC::ucell,GlobalC::vdwd2_para);
			vdwd2.cal_energy();
			GlobalC::en.evdw = vdwd2.get_energy();
		}
		if(GlobalC::vdwd3_para.flag_vdwd3)		//jiyy add 2019-05-18, update 2021-05-02
		{
			Vdwd3 vdwd3(GlobalC::ucell,GlobalC::vdwd3_para);
			vdwd3.cal_energy();
			GlobalC::en.evdw = vdwd3.get_energy();
		}


		// mohan added eiter to count for the electron iteration number, 2021-01-28
		int eiter=0;		
        if (GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="md" || GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax")  // pengfei 2014-10-13
        {
#ifdef __LCAO
#ifdef __MPI
			if( Exx_Global::Hybrid_Type::No==GlobalC::exx_global.info.hybrid_type  )
			{	
#endif
#endif		
				p_esolver->Run(istep,GlobalC::ucell);
				p_esolver->cal_Energy(GlobalC::en);
				eiter = p_esolver->getiter();
#ifdef __LCAO
#ifdef __MPI
			}
			else if( Exx_Global::Hybrid_Type::Generate_Matrix == GlobalC::exx_global.info.hybrid_type )
			{
				throw std::invalid_argument(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
			}
			else	// Peize Lin add 2019-03-09
			{
				if( GlobalC::exx_global.info.separate_loop )
				{
					for( size_t hybrid_step=0; hybrid_step!=GlobalC::exx_global.info.hybrid_step; ++hybrid_step )
					{
						elec.self_consistent(istep-1);
						eiter += elec.iter;
						if( elec.iter==1 || hybrid_step==GlobalC::exx_global.info.hybrid_step-1 )		// exx converge
							break;
						XC_Functional::set_xc_type(GlobalC::ucell.atoms[0].xc_func);					
						GlobalC::exx_lip.cal_exx();
					}						
				}
				else
				{
					elec.self_consistent(istep-1);	
					eiter += elec.iter;
					XC_Functional::set_xc_type(GlobalC::ucell.atoms[0].xc_func);
					elec.self_consistent(istep-1);
					eiter += elec.iter;
				}
			}
#endif //__MPI
#endif //__LCAO
        }
        else if(GlobalV::CALCULATION=="nscf")
        {
            elec.non_self_consistent(istep-1);
			eiter = elec.iter;
        }

		if(GlobalC::pot.out_pot == 2)
		{
			std::stringstream ssp;
			std::stringstream ssp_ave;
			ssp << GlobalV::global_out_dir << "ElecStaticPot";
			ssp_ave << GlobalV::global_out_dir << "ElecStaticPot_AVE";
			GlobalC::pot.write_elecstat_pot(ssp.str(), ssp_ave.str()); //output 'Hartree + local pseudopot'
		}

		time_t eend = time(NULL);
		time_t fstart = time(NULL);


        if (GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax")
        {
			stop = this->after_scf(p_esolver, istep, force_step, stress_step);    // pengfei Li 2018-05-14
		}
		time_t fend = time(NULL);


		if(GlobalV::OUT_LEVEL=="i")
		{
			double etime_min = difftime(eend, estart)/60.0; 
			double ftime_min = difftime(fend, fstart)/60.0; 
			std::stringstream ss;
			ss << GlobalV::RELAX_METHOD << istep;
			
			std::cout << " " << std::setw(7) << ss.str() 
			<< std::setw(5) << eiter 
			<< std::setw(15) << std::setprecision(6) << GlobalC::en.etot * ModuleBase::Ry_to_eV 
			<< std::setw(15) << IMM.get_ediff() * ModuleBase::Ry_to_eV
			<< std::setprecision(3)
			<< std::setw(15) << IMM.get_largest_grad() * ModuleBase::Ry_to_eV / 0.529177
			<< std::setw(15) << IMM.get_trust_radius()
			<< std::setw(8) << IMM.get_update_iter()
			<< std::setprecision(2) << std::setw(11) << etime_min
			<< std::setw(11) << ftime_min << std::endl;
		}

		++istep;

    }

    if(GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax")
    {
        GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
        GlobalV::ofs_running << std::setprecision(16);
        GlobalV::ofs_running << " !FINAL_ETOT_IS " << GlobalC::en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl; 
        GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;
    }


	if(GlobalV::OUT_LEVEL=="i")
	{
		std::cout << " ION DYNAMICS FINISHED :)" << std::endl;
	}

	if(GlobalC::wf.out_wfc_r == 1)				// Peize Lin add 2021.11.21
	{
		Write_Wfc_Realspace::write_wfc_realspace_1(GlobalC::wf.evc, "wfc_realspace", true);
	}	

	ModuleBase::timer::tick("Ions","opt_ions_pw");
    return;
}

bool Ions::after_scf(ModuleESolver::ESolver *p_esolver, const int &istep, int &force_step, int &stress_step)
{
	ModuleBase::TITLE("Ions","after_scf");
	//calculate and gather all parts of total ionic forces
	ModuleBase::matrix force;
	if(GlobalV::CAL_FORCE)
	{
		this->gather_force_pw(p_esolver, force);
	}
	//calculate and gather all parts of stress
	ModuleBase::matrix stress;
	if(GlobalV::CAL_STRESS)
	{
		this->gather_stress_pw(p_esolver, stress);
	}
	//stop in last step
	if(istep==GlobalV::RELAX_NMAX)
	{
		return 1;
	}
	//choose what to do next
	if(GlobalV::CALCULATION!="cell-relax") force_step = istep;
	if(this->if_do_relax()) 
	{
		//do relax calculation and generate next structure 
		bool converged = 0;
		converged = this->do_relax(istep, force_step, force, GlobalC::en.etot);
		if(!converged) 
		{
			this->reset_after_relax(istep);
			return converged;
		}
		else if(GlobalV::CALCULATION!="cell-relax")
		{
			return converged;
		}
	}
	if(this->if_do_cellrelax())
	{
		//do cell relax calculation and generate next structure
		bool converged = 0;
		converged = this->do_cellrelax(stress_step, stress, GlobalC::en.etot);
		if(!converged) this->reset_after_cellrelax(force_step, stress_step);
		return converged;
	}

    return 1;
}
void Ions::gather_force_pw(ModuleESolver::ESolver *p_esolver, ModuleBase::matrix &force)
{
	ModuleBase::TITLE("Ions","gather_force_pw");
	// Forces fcs;
	// fcs.init(force);
	p_esolver->cal_Force(force);
}

void Ions::gather_stress_pw(ModuleESolver::ESolver *p_esolver, ModuleBase::matrix& stress)
{
	ModuleBase::TITLE("Ions","gather_stress_pw");
	//basic stress
	// Stress_PW ss;
	// ss.cal_stress(stress);
	p_esolver->cal_Stress(stress);
	//external stress
	double unit_transform = 0.0;
	unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI,3) * 1.0e-8;
	double external_stress[3] = {GlobalV::PRESS1,GlobalV::PRESS2,GlobalV::PRESS3};
	for(int i=0;i<3;i++)
	{
		stress(i,i) -= external_stress[i]/unit_transform;
	}
	GlobalV::PRESSURE = (stress(0,0)+stress(1,1)+stress(2,2))/3;
}

bool Ions::if_do_relax()
{
	ModuleBase::TITLE("Ions","if_do_relax");
	if(GlobalV::CALCULATION=="relax"||GlobalV::CALCULATION=="cell-relax")
	{
		if(!GlobalC::ucell.if_atoms_can_move()) 
		{
			ModuleBase::WARNING("Ions","No atom is allowed to move!");
			return 0;
		}
//		if(!IMM.get_converged()) return 1;
		else 
		{
			assert(GlobalV::CAL_FORCE==1);
			return 1;
		}
	}
	else return 0;
}
bool Ions::if_do_cellrelax()
{
	ModuleBase::TITLE("Ions","if_do_cellrelax");
	if(GlobalV::CALCULATION=="cell-relax")
	{
		if(!GlobalC::ucell.if_cell_can_change()) 
		{
			ModuleBase::WARNING("Ions", "Lattice vectors are not allowed to change!");
			return 0;
		}
		else if(GlobalC::ucell.if_atoms_can_move()&&!IMM.get_converged())
		{
			GlobalV::ofs_running<<"Note: Need to wait for atomic relaxation first!";
			return 0;
		}
		else 
		{
			assert(GlobalV::CAL_STRESS==1);
			return 1;
		}
	}
	else return 0;
}
bool Ions::do_relax(const int& istep, int& jstep, const ModuleBase::matrix& ionic_force, const double& total_energy)
{
	ModuleBase::TITLE("Ions","do_relax");
	CE.update_istep(jstep);
	CE.update_all_pos(GlobalC::ucell);
	IMM.cal_movement(istep, jstep, ionic_force, total_energy);
	++jstep;
	return IMM.get_converged();
}
bool Ions::do_cellrelax(const int& istep, const ModuleBase::matrix& stress, const double& total_energy)
{
	ModuleBase::TITLE("Ions","do_cellrelax");
	LCM.cal_lattice_change(istep, stress, total_energy);
    return LCM.get_converged();
}
void Ions::reset_after_relax(const int& istep)
{
	ModuleBase::TITLE("Ions","reset_after_relax");
	GlobalV::ofs_running << " Setup the structure factor in plane wave basis." << std::endl;
	GlobalC::pw.setup_structure_factor();

	GlobalV::ofs_running << " Setup the extrapolated charge." << std::endl;
	// charge extrapolation if istep>0.
	CE.extrapolate_charge();
	CE.save_pos_next(GlobalC::ucell);

	GlobalV::ofs_running << " Setup the Vl+Vh+Vxc according to new structure factor and new charge." << std::endl;
	// calculate the new potential accordint to
	// the new charge density.
	GlobalC::pot.init_pot( istep, GlobalC::pw.strucFac );

	GlobalV::ofs_running << " Setup the new wave functions?" << std::endl;
	GlobalC::wf.wfcinit();
}
void Ions::reset_after_cellrelax(int& f_step, int& s_step)
{
	ModuleBase::TITLE("Ions","reset_after_cellrelax");
	Variable_Cell::init_after_vc();
	GlobalC::pot.init_pot(s_step, GlobalC::pw.strucFac); //LiuXh add 20180619

	GlobalV::ofs_running << " Setup the new wave functions?" << std::endl; //LiuXh add 20180619
	GlobalC::wf.wfcinit(); //LiuXh add 20180619
	f_step = 1;
	++s_step;
}
