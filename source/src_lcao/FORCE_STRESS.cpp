#include "FORCE_STRESS.h"
#include "../src_pw/global.h"
#include "./dftu.h"  //Quxin add for DFT+U on 20201029
// new
#include "../src_pw/vdwd2.h"
#include "../src_pw/vdwd3.h"
#include "../module_base/timer.h"
#include "../module_surchem/efield.h"        // liuyu add 2022-05-18
#ifdef __DEEPKS
#include "../module_deepks/LCAO_deepks.h"	//caoyu add for deepks 2021-06-03
#endif

double Force_Stress_LCAO::force_invalid_threshold_ev = 0.00;
double Force_Stress_LCAO::output_acc = 1.0e-8;

Force_Stress_LCAO::Force_Stress_LCAO(Record_adj& ra) :
    RA(&ra){}
Force_Stress_LCAO::~Force_Stress_LCAO() {}

void Force_Stress_LCAO::getForceStress(
	const bool isforce,
	const bool isstress,
	const bool istestf,
    const bool istests,
    Local_Orbital_Charge& loc,
    const psi::Psi<double>* psid,
	const psi::Psi<std::complex<double>>* psi,
    LCAO_Hamilt &uhm,
    ModuleBase::matrix& fcs,
	ModuleBase::matrix &scs)
{
    ModuleBase::TITLE("Force_Stress_LCAO","getForceStress");
	ModuleBase::timer::tick("Force_Stress_LCAO","getForceStress");
	
	if(!isforce&&!isstress)
	{
		ModuleBase::timer::tick("Force_Stress_LCAO","getForceStress");
		return;
	} 

	const int nat = GlobalC::ucell.nat;

	//total force : ModuleBase::matrix fcs;

	// part of total force
    ModuleBase::matrix foverlap;
    ModuleBase::matrix ftvnl_dphi;
    ModuleBase::matrix fvnl_dbeta;
    ModuleBase::matrix fvl_dphi;
    ModuleBase::matrix fvl_dvl;
    ModuleBase::matrix fewalds;
    ModuleBase::matrix fcc;
	ModuleBase::matrix fscc;

	fvl_dphi.create (nat, 3);//must do it now, update it later, noted by zhengdy

	if(isforce)
	{
		fcs.create (nat, 3);
		foverlap.create (nat, 3);
		ftvnl_dphi.create (nat, 3);
		fvnl_dbeta.create (nat, 3);
		fvl_dvl.create (nat, 3);
		fewalds.create (nat, 3);
		fcc.create (nat, 3);
		fscc.create (nat, 3);
		//calculate basic terms in Force, same method with PW base
		this->calForcePwPart(fvl_dvl, fewalds, fcc, fscc);
	}

	//total stress : ModuleBase::matrix scs
	ModuleBase::matrix sigmacc;
	ModuleBase::matrix sigmadvl;
	ModuleBase::matrix sigmaewa;
	ModuleBase::matrix sigmaxc;
	ModuleBase::matrix sigmahar;
	ModuleBase::matrix soverlap;
	ModuleBase::matrix stvnl_dphi;
	ModuleBase::matrix svnl_dbeta;
	ModuleBase::matrix svl_dphi;	
#ifdef __DEEPKS
	ModuleBase::matrix svnl_dalpha; //deepks
#endif

	if(isstress)
	{
		scs.create(3,3);
		sigmacc.create(3,3);
		sigmadvl.create(3,3);
		sigmaewa.create(3,3);
		sigmaxc.create(3,3);
		sigmahar.create(3,3);

		soverlap.create(3,3);
		stvnl_dphi.create(3,3);
		svnl_dbeta.create(3,3);
		svl_dphi.create(3,3);
#ifdef __DEEPKS
		svnl_dalpha.create(3,3);
#endif
		//calculate basic terms in Stress, similar method with PW base
		this->calStressPwPart(
				sigmadvl,
				sigmahar,
				sigmaewa,
				sigmacc,
				sigmaxc);
		this->stress_mgga(sigmaxc, loc, uhm);
	}
	//--------------------------------------------------------
	// implement four terms which needs integration
	//--------------------------------------------------------
	this->calForceStressIntegralPart(
				GlobalV::GAMMA_ONLY_LOCAL,
				isforce,
				isstress,
                loc,
                psid,
				psi,
                foverlap,
				ftvnl_dphi,
				fvnl_dbeta,
				fvl_dphi,
				soverlap,
				stvnl_dphi,
				svnl_dbeta,
#ifdef __DEEPKS
				svl_dphi,
				svnl_dalpha,
#else
				svl_dphi,
#endif
                uhm);
    //implement vdw force or stress here
	// Peize Lin add 2014-04-04, update 2021-03-09
	ModuleBase::matrix force_vdw;
	ModuleBase::matrix stress_vdw;
	if(GlobalC::vdwd2_para.flag_vdwd2)
	{
		if(isforce)
		{
			force_vdw.create(nat,3);
			Vdwd2 vdwd2(GlobalC::ucell,GlobalC::vdwd2_para);
			vdwd2.cal_force();
			for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
			{
				force_vdw(iat,0) = vdwd2.get_force()[iat].x;
				force_vdw(iat,1) = vdwd2.get_force()[iat].y;
				force_vdw(iat,2) = vdwd2.get_force()[iat].z;
			}
		}
		if(isstress)
		{
			Vdwd2 vdwd2(GlobalC::ucell,GlobalC::vdwd2_para);
			vdwd2.cal_stress();
			stress_vdw = vdwd2.get_stress().to_matrix();
		}
	}
	// jiyy add 2019-05-18, update 2021-05-02
	else if(GlobalC::vdwd3_para.flag_vdwd3)
	{
		if(isforce)
		{
			force_vdw.create(nat,3);
			Vdwd3 vdwd3(GlobalC::ucell,GlobalC::vdwd3_para);
			vdwd3.cal_force();
			for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
			{
				force_vdw(iat,0) = vdwd3.get_force()[iat].x;
				force_vdw(iat,1) = vdwd3.get_force()[iat].y;
				force_vdw(iat,2) = vdwd3.get_force()[iat].z;
			}
		}
		if(isstress)
		{
			Vdwd3 vdwd3(GlobalC::ucell,GlobalC::vdwd3_para);
			vdwd3.cal_stress();
			stress_vdw = vdwd3.get_stress().to_matrix();
		}
	}
	//implement force from E-field
    ModuleBase::matrix fefield;
    if(GlobalV::EFIELD_FLAG&&isforce)
    {
        fefield.create(nat, 3);
        Efield::compute_force(GlobalC::ucell, fefield);
    }
	//Force contribution from DFT+U
	ModuleBase::matrix force_dftu;
	ModuleBase::matrix stress_dftu;
	if (INPUT.dft_plus_u)
	{
		// Quxin add for DFT+U on 20201029
		GlobalC::dftu.force_stress(loc.dm_gamma, loc.dm_k, *uhm.LM);
		
        if (isforce) {
            force_dftu.create(nat, 3);
		}
		if(isstress)
		{
			stress_dftu.create(3, 3);
		}
		for(int i=0; i<3; i++)
		{
			if(isstress)
			{
				for(int j=0; j<3; j++)
				{
					stress_dftu(j,i) = GlobalC::dftu.stress_dftu.at(j).at(i);
				}
			}
			if(isforce)
			{
				for (int iat = 0; iat < nat; iat++)
				{
					force_dftu(iat, i) = GlobalC::dftu.force_dftu.at(iat).at(i);
				}
			}
		}
	}
	//--------------------------------
	//begin calculate and output force
	//--------------------------------
	if(isforce)
	{
		//---------------------------------
		//sum all parts of force!
		//---------------------------------
		for(int i=0; i<3; i++)
		{
			double sum = 0.0;

			for (int iat = 0; iat < nat; iat++)
			{
				fcs(iat, i) += foverlap(iat, i)
					+ ftvnl_dphi(iat, i)
					+ fvnl_dbeta(iat, i)
					+ fvl_dphi(iat, i)
					+ fvl_dvl(iat, i) // derivative of local potential force (pw)
					+ fewalds(iat, i) // ewald force (pw)
					+ fcc(iat, i) //nonlinear core correction force (pw)
					+ fscc(iat, i);//self consistent corretion force (pw)

				// Force contribution from DFT+U, Quxin add on 20201029
				if(INPUT.dft_plus_u)
				{
					fcs(iat, i) += force_dftu(iat, i);
				}
				//VDW force of vdwd2 or vdwd3
				if(GlobalC::vdwd2_para.flag_vdwd2||GlobalC::vdwd3_para.flag_vdwd3)
				{
					fcs(iat,i) += force_vdw(iat,i);
				}
				//E-field force
				if(GlobalV::EFIELD_FLAG)
				{
					fcs(iat, i) += fefield(iat, i);
				}
#ifdef __DEEPKS
				// mohan add 2021-08-04
				if (GlobalV::deepks_scf)
				{
					fcs(iat, i) += GlobalC::ld.F_delta(iat, i);
				}
#endif
				//sum total force for correction
				sum += fcs(iat, i);
			}

			for(int iat=0; iat<nat; ++iat)
			{
				fcs(iat, i) -= sum/nat;
			}

			//xiaohui add "OUT_LEVEL", 2015-09-16
			if(GlobalV::OUT_LEVEL != "m")
			{
				GlobalV::ofs_running << " correction force for each atom along direction "
					<< i+1 << " is " << sum/nat << std::endl;
			}
		}

		// pengfei 2016-12-20
		if(ModuleSymmetry::Symmetry::symm_flag)
		{
			this->forceSymmetry(fcs);
		}

#ifdef __DEEPKS
		//DeePKS force, caoyu add 2021-06-03
		if (GlobalV::deepks_out_labels) //not parallelized yet
        {
            const Parallel_Orbitals* pv = loc.ParaV;
            GlobalC::ld.save_npy_f(fcs, "f_tot.npy", GlobalC::ucell.nat); //Ty/Bohr, F_tot
            if (GlobalV::deepks_scf)
            {
                GlobalC::ld.save_npy_f(fcs - GlobalC::ld.F_delta, "f_base.npy", GlobalC::ucell.nat); //Ry/Bohr, F_base

				if(GlobalV::GAMMA_ONLY_LOCAL)
				{
    				GlobalC::ld.cal_gdmx(loc.dm_gamma[0],
						GlobalC::ucell,
						GlobalC::ORB,
						GlobalC::GridD,
						pv->trace_loc_row,
    					pv->trace_loc_col,
						isstress);
				}
				else
				{			
					GlobalC::ld.cal_gdmx_k(loc.dm_k,
						GlobalC::ucell,
						GlobalC::ORB,
						GlobalC::GridD,
						pv->trace_loc_row,
    					pv->trace_loc_col,
						GlobalC::kv.nks,
						GlobalC::kv.kvec_d,
						isstress);	
				}
				if(GlobalV::deepks_out_unittest) GlobalC::ld.check_gdmx(GlobalC::ucell.nat);
				GlobalC::ld.cal_gvx(GlobalC::ucell.nat);
				
				if(GlobalV::deepks_out_unittest) GlobalC::ld.check_gvx(GlobalC::ucell.nat);
				GlobalC::ld.save_npy_gvx(GlobalC::ucell.nat);//  /Bohr, grad_vx
            }
            else
            {
                GlobalC::ld.save_npy_f(fcs, "f_base.npy", GlobalC::ucell.nat); //no scf, F_base=F_tot
            }

        }
#endif
		// print Rydberg force or not
		bool ry = false;
		if(istestf)
		{
			// test
			//ModuleBase::matrix fvlocal;
			//fvlocal.create(nat,3);
			ModuleBase::matrix ftvnl;
			ftvnl.create(nat, 3);
			for (int iat = 0; iat < nat; iat++)
			{
				for(int i=0; i<3; i++)
				{
					//fvlocal(iat,i) = fvl_dphi(iat,i) + fvl_dvl(iat,i);
					ftvnl(iat,i) = ftvnl_dphi(iat,i) + fvnl_dbeta(iat,i);
				}
			}

			GlobalV::ofs_running << "\n PARTS OF FORCE: " << std::endl;
			GlobalV::ofs_running << std::setiosflags(ios::showpos);
			GlobalV::ofs_running << std::setiosflags(ios::fixed) << std::setprecision(8) << std::endl;
			//-----------------------------
			//regular force terms test.
			//-----------------------------
			//this->print_force("OVERLAP    FORCE",foverlap,1,ry);
			f_pw.print("OVERLAP    FORCE", foverlap,0);
			//  this->print_force("TVNL_DPHI  force",ftvnl_dphi,GlobalV::TEST_FORCE);
			//  this->print_force("VNL_DBETA  force",fvnl_dbeta,GlobalV::TEST_FORCE);
			//this->print_force("T_VNL      FORCE",ftvnl,1,ry);
			f_pw.print("T_VNL      FORCE", ftvnl,0);
			f_pw.print("VL_dPHI    FORCE", fvl_dphi,0);
			//this->print_force("VL_dPHI    FORCE",fvl_dphi,1,ry);
			//this->print_force("VL_dVL     FORCE",fvl_dvl,1,ry);
			f_pw.print("VL_dVL     FORCE", fvl_dvl,0);
			f_pw.print("EWALD      FORCE", fewalds,0);
			// 	this->print_force("VLOCAL     FORCE",fvlocal,GlobalV::TEST_FORCE);
			//this->print_force("EWALD      FORCE",fewalds,1,ry);
			f_pw.print("NLCC       FORCE", fcc,0);
			f_pw.print("SCC        FORCE", fscc,0);
			//this->print_force("NLCC       FORCE",fcc,1,ry);
			//this->print_force("SCC        FORCE",fscc,1,ry);
			//-------------------------------
			//put extra force here for test!
			//-------------------------------
			if(GlobalV::EFIELD_FLAG)
			{
				f_pw.print("EFIELD     FORCE", fefield,0);
				//this->print_force("EFIELD     FORCE",fefield,1,ry);
			}
			if(GlobalC::vdwd2_para.flag_vdwd2||GlobalC::vdwd3_para.flag_vdwd3)
			{
				f_pw.print("VDW        FORCE", force_vdw,0);
				//this->print_force("VDW        FORCE",force_vdw,1,ry);
			}
#ifdef __DEEPKS
			//caoyu add 2021-06-03
			if (GlobalV::deepks_scf)
			{
				this->print_force("DeePKS 	FORCE", GlobalC::ld.F_delta, 1, ry);
			}
#endif
		}

		GlobalV::ofs_running << std::setiosflags(ios::left);

		//this->printforce_total(ry, istestf, fcs);
		f_pw.print("   TOTAL-FORCE (eV/Angstrom)", fcs,0);
		if(istestf)
		{
			GlobalV::ofs_running << "\n FORCE INVALID TABLE." << std::endl;
			GlobalV::ofs_running << " " << std::setw(8) << "atom" << std::setw(5) << "x" << std::setw(5) << "y" << std::setw(5) << "z" << std::endl;
			for(int iat=0; iat<GlobalC::ucell.nat; iat++)
			{
				GlobalV::ofs_running << " " << std::setw(8) << iat;
				for(int i=0; i<3; i++)
				{
					if( abs( fcs(iat,i)*ModuleBase::Ry_to_eV/0.529177 ) < Force_Stress_LCAO::force_invalid_threshold_ev)
					{
						fcs(iat,i) = 0.0;
						GlobalV::ofs_running << std::setw(5) << "1";
					}
					else
					{
						GlobalV::ofs_running << std::setw(5) << "0";
					}
				}
				GlobalV::ofs_running << std::endl;
			}
		}
	}//end of force calculation
	//---------------------------------
	//begin calculate and output stress
	//---------------------------------
	if(isstress)
	{
		for(int i=0; i<3; i++)
		{
			for (int j=0;j<3;j++)
			{
				scs(i,j) += soverlap(i,j)
					+ stvnl_dphi(i,j)
					+ svnl_dbeta(i,j)
					+ svl_dphi(i,j)
					+ sigmadvl(i,j) // derivative of local potential stress (pw)
					+ sigmaewa(i,j) // ewald stress (pw)
					+ sigmacc(i,j) //nonlinear core correction stress (pw)
					+ sigmaxc(i,j)//exchange corretion stress
					+ sigmahar(i,j);// hartree stress

					//VDW stress from linpz and jiyy
				if(GlobalC::vdwd2_para.flag_vdwd2||GlobalC::vdwd3_para.flag_vdwd3)
				{
					scs(i,j) += stress_vdw(i , j);
				}
				//DFT plus U stress from qux
				if(INPUT.dft_plus_u)
				{
					scs(i, j) += stress_dftu(i, j);
				}
			}
		}

#ifdef __DEEPKS
		if (GlobalV::deepks_out_labels) //not parallelized yet
        {
			GlobalC::ld.save_npy_s(scs, "s_base.npy", GlobalC::ucell.omega); //Ry/Bohr^3, S_base; no scf case -- same as S_tot
			// wenfei add 2021/11/2
			if (GlobalV::deepks_scf)
			{
				for (int i=0; i<3; i++)
				{
					for (int j=0; j<3; j++)
					{
						scs(i,j) += svnl_dalpha(i,j);
					}
				}
				GlobalC::ld.save_npy_s(scs, "s_tot.npy", GlobalC::ucell.omega); //Ry/Bohr^3, S_tot
				GlobalC::ld.cal_gvepsl(GlobalC::ucell.nat);
				GlobalC::ld.save_npy_gvepsl(GlobalC::ucell.nat);//  /Bohr^3, grad_vepsl
			}
			else
			{
				GlobalC::ld.save_npy_s(scs, "s_tot.npy", GlobalC::ucell.omega); //Ry/Bohr^3, S_tot; no scf case -- same as S_base
			}
		}	
#endif


		if(ModuleSymmetry::Symmetry::symm_flag)
		{
			GlobalC::symm.stress_symmetry(scs, GlobalC::ucell);
		}//end symmetry

		// print Rydberg stress or not
		bool ry = false;

		//test stress each terms if needed
		if(istests)
		{
			// test
			ModuleBase::matrix svlocal;
			svlocal.create(3,3);
			ModuleBase::matrix stvnl;
			stvnl.create(3,3);
			for (int i = 0; i<3; i++)
			{
				for(int j=0; j<3; j++)
				{
					svlocal(i,j) = svl_dphi(i,j) + sigmadvl(i,j);
					stvnl(i,j) = stvnl_dphi(i,j) + svnl_dbeta(i,j);
				}
			}

			GlobalV::ofs_running << "\n PARTS OF STRESS: " << std::endl;
			GlobalV::ofs_running << std::setiosflags(ios::showpos);
			GlobalV::ofs_running << std::setiosflags(ios::fixed) << std::setprecision(10) << std::endl;
			sc_pw.print_stress("OVERLAP  STRESS",soverlap,GlobalV::TEST_STRESS,ry);
			//test
			sc_pw.print_stress("T        STRESS",stvnl_dphi,GlobalV::TEST_STRESS,ry);
			sc_pw.print_stress("VNL      STRESS",svnl_dbeta,GlobalV::TEST_STRESS,ry);

			sc_pw.print_stress("T_VNL    STRESS",stvnl,GlobalV::TEST_STRESS,ry);

			sc_pw.print_stress("VL_dPHI  STRESS",svl_dphi,GlobalV::TEST_STRESS,ry);
			sc_pw.print_stress("VL_dVL   STRESS",sigmadvl,GlobalV::TEST_STRESS,ry);
			sc_pw.print_stress("HAR      STRESS",sigmahar,GlobalV::TEST_STRESS,ry);

			sc_pw.print_stress("EWALD    STRESS",sigmaewa,GlobalV::TEST_STRESS,ry);
			sc_pw.print_stress("cc       STRESS",sigmacc,GlobalV::TEST_STRESS,ry);
			//		sc_pw.print_stress("NLCC       STRESS",sigmacc,GlobalV::TEST_STRESS,ry);
			sc_pw.print_stress("XC       STRESS",sigmaxc,GlobalV::TEST_STRESS,ry);
			if(GlobalC::vdwd2_para.flag_vdwd2||GlobalC::vdwd3_para.flag_vdwd3)
			{
				sc_pw.print_stress("VDW      STRESS",sigmaxc,GlobalV::TEST_STRESS,ry);
			}
			if(INPUT.dft_plus_u)
			{
				sc_pw.print_stress("DFTU     STRESS",sigmaxc,GlobalV::TEST_STRESS,ry);
			}
			sc_pw.print_stress("TOTAL    STRESS",scs,GlobalV::TEST_STRESS,ry);

		}//end of test
		GlobalV::ofs_running << std::setiosflags(ios::left);
		//print total stress
		sc_pw.printstress_total(scs, ry);

		double unit_transform = 0.0;
		unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI,3) * 1.0e-8;
		double external_stress[3] = {GlobalV::PRESS1,GlobalV::PRESS2,GlobalV::PRESS3};

		for(int i=0;i<3;i++)
		{
			scs(i,i) -= external_stress[i]/unit_transform;
		}
		GlobalV::PRESSURE = (scs(0,0)+scs(1,1)+scs(2,2))/3;
	}//end of stress calculation
	
	ModuleBase::timer::tick("Force_Stress_LCAO","getForceStress");
	return;
}

//print force term for test
void Force_Stress_LCAO::print_force(const std::string &name, ModuleBase::matrix& f, const bool screen, bool ry)const
{
	GlobalV::ofs_running << " --------------------------- " << name << " ----------------------------" << std::endl;
	GlobalV::ofs_running << " " << std::setw(8) << "atom" << std::setw(15) << "x" << std::setw(15) << "y" << std::setw(15) << "z" << std::endl;

	double fac = 1.0;

	if(!ry)
	{
	 	fac = ModuleBase::Ry_to_eV / 0.529177;
	}

	std::cout << std::setprecision(5);
	std::cout << std::setiosflags(ios::showpos);

	if(screen)
	{
		std::cout << " ------------------- " << name << " --------------------" << std::endl;
		std::cout << " " << std::setw(8) << "atom" << std::setw(15) << "x" << std::setw(15) << "y" << std::setw(15) << "z" << std::endl;
	}

    int iat = 0;
    for (int it = 0;it < GlobalC::ucell.ntype;it++)
    {
        for (int ia = 0;ia < GlobalC::ucell.atoms[it].na;ia++)
        {
			std::stringstream ss;
			ss << GlobalC::ucell.atoms[it].label << ia+1;

			GlobalV::ofs_running << " " << std::setw(8) << ss.str();
			if( abs(f(iat,0)) >output_acc) GlobalV::ofs_running << std::setw(15) << f(iat,0)*fac;
			else GlobalV::ofs_running << std::setw(15) << "0";
			if( abs(f(iat,1)) >output_acc) GlobalV::ofs_running << std::setw(15) << f(iat,1)*fac;
			else GlobalV::ofs_running << std::setw(15) << "0";
			if( abs(f(iat,2)) >output_acc) GlobalV::ofs_running << std::setw(15) << f(iat,2)*fac;
			else GlobalV::ofs_running << std::setw(15) << "0";
			GlobalV::ofs_running << std::endl;

			if(screen)
			{
				std::cout << " " << std::setw(8) << ss.str();
				if( abs(f(iat,0)) >output_acc) std::cout << std::setw(15) << f(iat,0)*fac;
				else std::cout << std::setw(15) << "0";
				if( abs(f(iat,1)) >output_acc) std::cout << std::setw(15) << f(iat,1)*fac;
				else std::cout << std::setw(15) << "0";
				if( abs(f(iat,2)) >output_acc) std::cout << std::setw(15) << f(iat,2)*fac;
				else std::cout << std::setw(15) << "0";
				std::cout << std::endl;
			}

            iat++;
        }
    }


	std::cout << std::resetiosflags(ios::showpos);

    return;
}


void Force_Stress_LCAO::printforce_total (const bool ry, const bool istestf, ModuleBase::matrix& fcs)
{
	ModuleBase::TITLE("Force_Stress_LCAO","printforce_total");
	double unit_transform = 1;

	if(!ry)
	{
		unit_transform = ModuleBase::Ry_to_eV / 0.529177;
	}
//	std::cout.setf(ios::fixed);

    int iat=0;

	//GlobalV::ofs_running << std::setiosflags(ios::right);
 	GlobalV::ofs_running << std::setprecision(6) << std::setiosflags(ios::showpos) << std::setiosflags(ios::fixed) << std::endl;
	ModuleBase::GlobalFunc::NEW_PART("TOTAL-FORCE (eV/Angstrom)");

	// print out forces
	if(INPUT.out_force == 1)
	{
		std::ofstream ofs("FORCE.dat");
		if(!ofs)
		{
			std::cout << "open FORCE.dat error !" <<std::endl;
		}

		for(int iat=0; iat<GlobalC::ucell.nat; iat++)
		{
			ofs << "   " << fcs(iat,0)*ModuleBase::Ry_to_eV / 0.529177
				<< "   " << fcs(iat,1)*ModuleBase::Ry_to_eV / 0.529177
				<< "   " << fcs(iat,2)*ModuleBase::Ry_to_eV / 0.529177 << std::endl;
		}
		ofs.close();
	}

 	if(istestf)
	{
		cout << setprecision(6);
		//cout << setiosflags(ios::showpos);
		//cout << setiosflags(ios::fixed) << endl;
		cout << " ------------------- TOTAL      FORCE --------------------" << endl;
    	cout << " " << setw(8) << "Atom" << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;
    	GlobalV::ofs_running << " " << setw(12) << "Atom" << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;
	}

    iat=0;
    for (int it=0; it<GlobalC::ucell.ntype; it++)
    {
        for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
        {
            std::stringstream ss;
            ss << GlobalC::ucell.atoms[it].label << ia+1;

			if(istestf)
			{
            	std::cout << " " << std::setw(8) << ss.str()
					<< std::setw(15) << fcs(iat,0)*unit_transform
					<< std::setw(15) << fcs(iat,1)*unit_transform
					<< std::setw(15) << fcs(iat,2)*unit_transform << std::endl;
			}

            GlobalV::ofs_running << " " << std::setw(12) << ss.str()
				<< std::setw(15) << fcs(iat,0)*unit_transform
				<< std::setw(15) << fcs(iat,1)*unit_transform
				<< std::setw(15) << fcs(iat,2)*unit_transform << std::endl;

            ++iat;
        }
    }
	GlobalV::ofs_running << std::setiosflags(ios::left);
	std::cout << std::resetiosflags(ios::showpos);

    return;
}

//local pseudopotential, ewald, core correction, scc terms in force
void Force_Stress_LCAO::calForcePwPart(
	ModuleBase::matrix &fvl_dvl,
	ModuleBase::matrix &fewalds,
	ModuleBase::matrix &fcc,
	ModuleBase::matrix &fscc)
{
	//--------------------------------------------------------
	// local pseudopotential force:
	// use charge density; plane wave; local pseudopotential;
	//--------------------------------------------------------
	f_pw.cal_force_loc (fvl_dvl, GlobalC::rhopw);
	//--------------------------------------------------------
	// ewald force: use plane wave only.
	//--------------------------------------------------------
	f_pw.cal_force_ew (fewalds,  GlobalC::rhopw); //remain problem
	//--------------------------------------------------------
	// force due to core correlation.
	//--------------------------------------------------------
	f_pw.cal_force_cc(fcc, GlobalC::rhopw);
	//--------------------------------------------------------
	// force due to self-consistent charge.
	//--------------------------------------------------------
	f_pw.cal_force_scc(fscc, GlobalC::rhopw);
	return;
}

//overlap, kinetic, nonlocal pseudopotential, Local potential terms in force and stress
void Force_Stress_LCAO::calForceStressIntegralPart(
	const bool isGammaOnly,
	const bool isforce,
    const bool isstress,
    Local_Orbital_Charge& loc,
    const psi::Psi<double>* psid,
	const psi::Psi<std::complex<double>>* psi,
    ModuleBase::matrix& foverlap,
	ModuleBase::matrix& ftvnl_dphi,
	ModuleBase::matrix& fvnl_dbeta,
	ModuleBase::matrix& fvl_dphi,
	ModuleBase::matrix& soverlap,
	ModuleBase::matrix& stvnl_dphi,
	ModuleBase::matrix& svnl_dbeta,
#if __DEEPKS
	ModuleBase::matrix& svl_dphi,
	ModuleBase::matrix& svnl_dalpha,
#else
	ModuleBase::matrix& svl_dphi,
#endif
    LCAO_Hamilt &uhm)
{
	if(isGammaOnly)
	{
    	flk.ftable_gamma(
				isforce,
				isstress,
                psid,
                loc,
                foverlap,
				ftvnl_dphi,
				fvnl_dbeta,
				fvl_dphi,
				soverlap,
				stvnl_dphi,
				svnl_dbeta,
#if __DEEPKS
				svl_dphi,
				svnl_dalpha,
#else
				svl_dphi,
#endif
                uhm);
    }
	else
	{
		flk.ftable_k(
				isforce,
				isstress,
                *this->RA,
                psi,
                loc,
                foverlap,
				ftvnl_dphi,
				fvnl_dbeta,
				fvl_dphi,
				soverlap,
				stvnl_dphi,
				svnl_dbeta,
#if __DEEPKS
				svl_dphi,
				svnl_dalpha,
#else
				svl_dphi,
#endif
                uhm);
    }
	return;
}

//vlocal, hartree, ewald, core correction, exchange-correlation terms in stress
void Force_Stress_LCAO::calStressPwPart(
	ModuleBase::matrix& sigmadvl,
	ModuleBase::matrix& sigmahar,
	ModuleBase::matrix& sigmaewa,
	ModuleBase::matrix& sigmacc,
	ModuleBase::matrix& sigmaxc
)
{
	//--------------------------------------------------------
	// local pseudopotential stress:
	// use charge density; plane wave; local pseudopotential;
	//--------------------------------------------------------
    sc_pw.stress_loc (sigmadvl, GlobalC::rhopw, 0);

	//--------------------------------------------------------
	//hartree term
	//--------------------------------------------------------
	sc_pw.stress_har (sigmahar, GlobalC::rhopw, 0);

	//--------------------------------------------------------
	// ewald stress: use plane wave only.
	//--------------------------------------------------------
    sc_pw.stress_ewa (sigmaewa,  GlobalC::rhopw, 0); //remain problem


	//--------------------------------------------------------
	// stress due to core correlation.
	//--------------------------------------------------------
	sc_pw.stress_cc(sigmacc,  GlobalC::rhopw, 0);

	//--------------------------------------------------------
	// stress due to self-consistent charge.
	//--------------------------------------------------------
	for(int i=0;i<3;i++)
	{
		sigmaxc(i,i) =  -(GlobalC::en.etxc) / GlobalC::ucell.omega;
	}
	//Exchange-correlation for PBE
	sc_pw.stress_gga(sigmaxc);

	return;
}

void Force_Stress_LCAO::stress_mgga(ModuleBase::matrix &sigmaxc,
		Local_Orbital_Charge& loc,
		LCAO_Hamilt &uhm)
{
	if(XC_Functional::get_func_type()==3)
	{
		double deband0 = 0.0;
		double deband_aux = 0.0;
		int ipol2xy[3][3];
		double sigma_mgga[3][3];
		double ***crosstaus;
		crosstaus = new double**[GlobalV::NSPIN];

		ipol2xy[0][0] = 0;
		ipol2xy[0][1] = 1;
		ipol2xy[0][2] = 2;
		ipol2xy[1][0] = 1;
		ipol2xy[1][1] = 3;
		ipol2xy[1][2] = 4;
		ipol2xy[2][0] = 2;
		ipol2xy[2][1] = 4;
		ipol2xy[2][2] = 5;

		for(int is = 0;is<GlobalV::NSPIN;is++)
		{
			crosstaus[is] = new double*[GlobalC::wfcpw->nrxx];
			for(int ir = 0;ir<GlobalC::wfcpw->nrxx; ir++)
			{
				crosstaus[is][ir] = new double [6];
				ModuleBase::GlobalFunc::ZEROS(crosstaus[is][ir],6);
			}
			Gint_inout inout(loc.DM_R, crosstaus[is], Gint_Tools::job_type::crosstaus);
			uhm.GK.cal_gint(&inout);
		}

		for(int is = 0; is < GlobalV::NSPIN; is++)
		{
			for (int ix = 0; ix < 3; ix++)
			{
				for (int iy = 0; iy < 3; iy++)
				{
					sigma_mgga[ix][iy] = 0.0;
					for(int ir = 0;ir<GlobalC::wfcpw->nrxx;ir++)
					{
						double x = GlobalC::pot.vofk(is,ir) * crosstaus[is][ir][ipol2xy[ix][iy]];
						sigma_mgga[ix][iy] += x;
					}
				}
			}
		}

#ifdef __MPI
		for(int l = 0;l<3;l++)
		{
			for(int m = 0;m<3;m++)
			{
				Parallel_Reduce::reduce_double_all( sigma_mgga[l][m] );
			}
		}
#endif

		for(int is = 0; is < GlobalV::NSPIN; is++)
		{
			for (int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
			{
				deband_aux += GlobalC::CHR.kin_r[is][ir] * GlobalC::pot.vofk(is,ir);
			}
		}
	#ifdef __MPI
		MPI_Allreduce(&deband_aux,&deband0,1,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
	#else
		deband0 = deband_aux;
	#endif

		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				sigmaxc(i,j) += sigma_mgga[i][j] * 2.0 / GlobalC::wfcpw->nxyz;
			}
		}

		for(int i=0;i<3;i++)   
		{
			sigmaxc(i,i) += deband0 / GlobalC::rhopw->nxyz;
		}

	}
}

#include "../module_base/mathzone.h"
//do symmetry for total force
void Force_Stress_LCAO::forceSymmetry(ModuleBase::matrix& fcs)
{
	double *pos;
	double d1,d2,d3;
	pos = new double[GlobalC::ucell.nat*3];
	ModuleBase::GlobalFunc::ZEROS(pos, GlobalC::ucell.nat*3);
	int iat = 0;
	for(int it = 0;it < GlobalC::ucell.ntype;it++)
	{
		for(int ia =0;ia< GlobalC::ucell.atoms[it].na;ia++)
		{
			pos[3*iat  ] = GlobalC::ucell.atoms[it].taud[ia].x ;
			pos[3*iat+1] = GlobalC::ucell.atoms[it].taud[ia].y ;
			pos[3*iat+2] = GlobalC::ucell.atoms[it].taud[ia].z;
			for(int k=0; k<3; ++k)
			{
				GlobalC::symm.check_translation( pos[iat*3+k], -floor(pos[iat*3+k]));
				GlobalC::symm.check_boundary( pos[iat*3+k] );
			}
			iat++;
		}
	}

	for(int iat=0; iat<GlobalC::ucell.nat; iat++)
	{
		ModuleBase::Mathzone::Cartesian_to_Direct(fcs(iat,0),fcs(iat,1),fcs(iat,2),
							GlobalC::ucell.a1.x, GlobalC::ucell.a1.y, GlobalC::ucell.a1.z,
							GlobalC::ucell.a2.x, GlobalC::ucell.a2.y, GlobalC::ucell.a2.z,
							GlobalC::ucell.a3.x, GlobalC::ucell.a3.y, GlobalC::ucell.a3.z,
							d1,d2,d3);

		fcs(iat,0) = d1;fcs(iat,1) = d2;fcs(iat,2) = d3;
	}
	GlobalC::symm.force_symmetry(fcs , pos, GlobalC::ucell);
	for(int iat=0; iat<GlobalC::ucell.nat; iat++)
	{
		ModuleBase::Mathzone::Direct_to_Cartesian(fcs(iat,0),fcs(iat,1),fcs(iat,2),
							GlobalC::ucell.a1.x, GlobalC::ucell.a1.y, GlobalC::ucell.a1.z,
							GlobalC::ucell.a2.x, GlobalC::ucell.a2.y, GlobalC::ucell.a2.z,
							GlobalC::ucell.a3.x, GlobalC::ucell.a3.y, GlobalC::ucell.a3.z,
							d1,d2,d3);

		fcs(iat,0) = d1;fcs(iat,1) = d2;fcs(iat,2) = d3;
	}
	//std::cout << "nrotk =" << GlobalC::symm.nrotk << std::endl;
	delete[] pos;
	return;
}
