#include "FORCE_STRESS.h"
#include "../src_pw/global.h"
#include "../src_pw/potential_libxc.h"
#include "./dftu.h"  //Quxin add for DFT+U on 20201029
// new
#include "../src_pw/H_XC_pw.h"
#include "../src_pw/vdwd2.h"
#include "../src_pw/vdwd3.h"
#include "LCAO_descriptor.h"	//caoyu add for deepks 2021-06-03

double Force_Stress_LCAO::force_invalid_threshold_ev = 0.00;
double Force_Stress_LCAO::output_acc = 1.0e-8;

Force_Stress_LCAO::Force_Stress_LCAO (){}
Force_Stress_LCAO::~Force_Stress_LCAO (){}


void Force_Stress_LCAO::allocate(void)
{
    TITLE("Force_Stress_LCAO","allocate");

    // reduce memory occupy by vlocal
    delete[] ParaO.sender_local_index;
    delete[] ParaO.sender_size_process;
    delete[] ParaO.sender_displacement_process;
    delete[] ParaO.receiver_global_index;
    delete[] ParaO.receiver_size_process;
    delete[] ParaO.receiver_displacement_process;

    ParaO.sender_local_index = new int[1];
    ParaO.sender_size_process = new int[1];
    ParaO.sender_displacement_process = new int[1];
    ParaO.receiver_global_index = new int[1];
    ParaO.receiver_size_process = new int[1];
    ParaO.receiver_displacement_process = new int[1];

    return;
}


#include "../src_pw/efield.h"
void Force_Stress_LCAO::getForceStress(
	const bool isforce, 
	const bool isstress, 
	const bool istestf, 
	const bool istests, 
	matrix &fcs, 
	matrix &scs)
{
    TITLE("Force_Stress_LCAO","getForceStress");
	timer::tick("Force_Stress_LCAO","getForceStress",'E');
	
	if(!isforce&&!isstress) return;

	const int nat = ucell.nat;

	//total force : matrix fcs;

	// part of total force
    matrix foverlap;
    matrix ftvnl_dphi;
    matrix fvnl_dbeta;
    matrix fvl_dphi;
    matrix fvl_dvl;
    matrix fewalds;
    matrix fcc;
	matrix fscc;

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

	//total stress : matrix scs
	matrix sigmacc;
	matrix sigmadvl;
	matrix sigmaewa;
	matrix sigmaxc;
	matrix sigmahar;
	matrix soverlap;
	matrix stvnl_dphi;
	matrix svnl_dbeta;
	matrix svl_dphi;

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
		//calculate basic terms in Stress, similar method with PW base
		this->calStressPwPart(
				sigmadvl,
				sigmahar,
				sigmaewa,
				sigmacc,
				sigmaxc);
	}
	
	//--------------------------------------------------------
	// implement four terms which needs integration
	//--------------------------------------------------------
	this->calForceStressIntegralPart(
				GAMMA_ONLY_LOCAL,
				isforce,
				isstress,
				foverlap,
				ftvnl_dphi,
				fvnl_dbeta,	
				fvl_dphi,
				soverlap,
				stvnl_dphi,
				svnl_dbeta,
				svl_dphi);
	
	//implement vdw force or stress here
	// Peize Lin add 2014-04-04, update 2021-03-09
	matrix force_vdw;
	matrix stress_vdw;
	if(vdwd2_para.flag_vdwd2)
	{
		if(isforce)
		{
			force_vdw.create(nat,3);
			Vdwd2 vdwd2(ucell,vdwd2_para);
			vdwd2.cal_force();
			for(int iat=0; iat<ucell.nat; ++iat)
			{
				force_vdw(iat,0) = vdwd2.get_force()[iat].x;
				force_vdw(iat,1) = vdwd2.get_force()[iat].y;
				force_vdw(iat,2) = vdwd2.get_force()[iat].z;
			}			
		}
		if(isstress)
		{
			Vdwd2 vdwd2(ucell,vdwd2_para);
			vdwd2.cal_stress();
			stress_vdw = vdwd2.get_stress().to_matrix();
		}
	}
	// jiyy add 2019-05-18, update 2021-05-02
	else if(vdwd3_para.flag_vdwd3)
	{
		if(isforce)
		{
			force_vdw.create(nat,3);
			Vdwd3 vdwd3(ucell,vdwd3_para);
			vdwd3.cal_force();
			for(int iat=0; iat<ucell.nat; ++iat)
			{
				force_vdw(iat,0) = vdwd3.get_force()[iat].x;
				force_vdw(iat,1) = vdwd3.get_force()[iat].y;
				force_vdw(iat,2) = vdwd3.get_force()[iat].z;
			}			
		}
		if(isstress)
		{
			Vdwd3 vdwd3(ucell,vdwd3_para);
			vdwd3.cal_stress();
			stress_vdw = vdwd3.get_stress().to_matrix();
		}
	}
	//implement force from E-field
    matrix fefield;
    if(EFIELD&&isforce)
    {
        fefield.create(nat, 3);
        Efield::compute_force(fefield);
    }
	//Force contribution from DFT+U
	matrix force_dftu;
	matrix stress_dftu;
	if (INPUT.dft_plus_u)
	{
		if(isforce)
		{
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
					stress_dftu(j,i) = dftu.stress_dftu.at(j).at(i);
				}
			}
			if(isforce)
			{
				for (int iat = 0; iat < nat; iat++)
				{
					force_dftu(iat, i) = dftu.force_dftu.at(iat).at(i);
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
				if(vdwd2_para.flag_vdwd2||vdwd3_para.flag_vdwd3)
				{
					fcs(iat,i) += force_vdw(iat,i);
				}
				//E-field force
				if(EFIELD)
				{
					fcs(iat, i) += fefield(iat, i);
				}
				//DFT plus U force
				if(INPUT.dft_plus_u)
				{
					fcs(iat, i) += force_dftu(iat, i);
				}
				//sum total force for correction
				sum += fcs(iat, i);
			}

			for(int iat=0; iat<nat; ++iat)
			{
				fcs(iat, i) -= sum/nat;
			}

			//xiaohui add "OUT_LEVEL", 2015-09-16
			if(OUT_LEVEL != "m") 
			{
				ofs_running << " correction force for each atom along direction " 
					<< i+1 << " is " << sum/nat << endl;
			}
		}
		
		// pengfei 2016-12-20
		if(Symmetry::symm_flag)
		{
			this->forceSymmetry(fcs);
		}
		
		//DeePKS force, caoyu add 2021-06-03
		if (INPUT.deepks_scf)
		{
			ld.save_npy_f(fcs);	//save fbase
		}
		
		// print Rydberg force or not
		bool ry = false;
		if(istestf)
		{
			// test
			//matrix fvlocal;
			//fvlocal.create(nat,3);
			matrix ftvnl;
			ftvnl.create(nat, 3);
			for (int iat = 0; iat < nat; iat++)
			{
				for(int i=0; i<3; i++)
				{
					//fvlocal(iat,i) = fvl_dphi(iat,i) + fvl_dvl(iat,i);
					ftvnl(iat,i) = ftvnl_dphi(iat,i) + fvnl_dbeta(iat,i);
				}
			}
			
			ofs_running << "\n PARTS OF FORCE: " << endl;
			ofs_running << setiosflags(ios::showpos);
			ofs_running << setiosflags(ios::fixed) << setprecision(8) << endl;
			//-----------------------------
			//regular force terms test.
			//-----------------------------
			this->print_force("OVERLAP    FORCE",foverlap,1,ry);
			//  this->print_force("TVNL_DPHI  force",ftvnl_dphi,TEST_FORCE);
			//  this->print_force("VNL_DBETA  force",fvnl_dbeta,TEST_FORCE);
			this->print_force("T_VNL      FORCE",ftvnl,1,ry);
			this->print_force("VL_dPHI    FORCE",fvl_dphi,1,ry);
			this->print_force("VL_dVL     FORCE",fvl_dvl,1,ry);
			// 	this->print_force("VLOCAL     FORCE",fvlocal,TEST_FORCE);
			this->print_force("EWALD      FORCE",fewalds,1,ry);
			this->print_force("NLCC       FORCE",fcc,1,ry);
			this->print_force("SCC        FORCE",fscc,1,ry);
			//-------------------------------
			//put extra force here for test! 
			//-------------------------------
			if(EFIELD)
			{
				this->print_force("EFIELD     FORCE",fefield,1,ry);
			}
			if(vdwd2_para.flag_vdwd2||vdwd3_para.flag_vdwd3)
			{
				this->print_force("VDW        FORCE",force_vdw,1,ry);
			}
			//caoyu add 2021-06-03
			if (INPUT.deepks_scf)
			{
				this->print_force("DeePKS 	FORCE", ld.F_delta, 1, ry);
			}
		}
		
		ofs_running << setiosflags(ios::left);  
			
		this->printforce_total(ry, istestf, fcs);
		if(istestf)
		{
			ofs_running << "\n FORCE INVALID TABLE." << endl;
			ofs_running << " " << setw(8) << "atom" << setw(5) << "x" << setw(5) << "y" << setw(5) << "z" << endl;
			for(int iat=0; iat<ucell.nat; iat++)
			{
				ofs_running << " " << setw(8) << iat;
				for(int i=0; i<3; i++)
				{
					if( abs( fcs(iat,i)*Ry_to_eV/0.529177 ) < Force_Stress_LCAO::force_invalid_threshold_ev)
					{
						fcs(iat,i) = 0.0;
						ofs_running << setw(5) << "1";
					}
					else
					{
						ofs_running << setw(5) << "0";
					}
				}
				ofs_running << endl;
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
				if(vdwd2_para.flag_vdwd2||vdwd3_para.flag_vdwd3)
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
		

		if(Symmetry::symm_flag)
		{
			symm.stress_symmetry(scs);
		}//end symmetry

		// print Rydberg stress or not
		bool ry = false;

		//test stress each terms if needed 
		if(istests)
		{
			// test
			matrix svlocal;
			svlocal.create(3,3);
			matrix stvnl;
			stvnl.create(3,3);
			for (int i = 0; i<3; i++)
			{
				for(int j=0; j<3; j++)
				{
					svlocal(i,j) = svl_dphi(i,j) + sigmadvl(i,j);
					stvnl(i,j) = stvnl_dphi(i,j) + svnl_dbeta(i,j);
				}
			}

			ofs_running << "\n PARTS OF STRESS: " << endl;
			ofs_running << setiosflags(ios::showpos);
			ofs_running << setiosflags(ios::fixed) << setprecision(10) << endl;
			sc_pw.print_stress("OVERLAP  STRESS",soverlap,TEST_STRESS,ry);
			//test
			sc_pw.print_stress("T        STRESS",stvnl_dphi,TEST_STRESS,ry);
			sc_pw.print_stress("VNL      STRESS",svnl_dbeta,TEST_STRESS,ry);	

			sc_pw.print_stress("T_VNL    STRESS",stvnl,TEST_STRESS,ry);

			sc_pw.print_stress("VL_dPHI  STRESS",svl_dphi,TEST_STRESS,ry);
			sc_pw.print_stress("VL_dVL   STRESS",sigmadvl,TEST_STRESS,ry);
			sc_pw.print_stress("HAR      STRESS",sigmahar,TEST_STRESS,ry);

			sc_pw.print_stress("EWALD    STRESS",sigmaewa,TEST_STRESS,ry);
			sc_pw.print_stress("cc       STRESS",sigmacc,TEST_STRESS,ry);
			//		sc_pw.print_stress("NLCC       STRESS",sigmacc,TEST_STRESS,ry);
			sc_pw.print_stress("XC       STRESS",sigmaxc,TEST_STRESS,ry);
			if(vdwd2_para.flag_vdwd2||vdwd3_para.flag_vdwd3)
			{
				sc_pw.print_stress("VDW      STRESS",sigmaxc,TEST_STRESS,ry);
			}
			if(INPUT.dft_plus_u)
			{
				sc_pw.print_stress("DFTU     STRESS",sigmaxc,TEST_STRESS,ry);
			}
			sc_pw.print_stress("TOTAL    STRESS",scs,TEST_STRESS,ry);
			
		}//end of test
		ofs_running << setiosflags(ios::left);  
		//print total stress 
		sc_pw.printstress_total(scs, ry);

		double unit_transform = 0.0;
		unit_transform = RYDBERG_SI / pow(BOHR_RADIUS_SI,3) * 1.0e-8;
		double external_stress[3] = {PRESS1,PRESS2,PRESS3};

		for(int i=0;i<3;i++)
		{
			scs(i,i) -= external_stress[i]/unit_transform;
		}
		PRESSURE = (scs(0,0)+scs(1,1)+scs(2,2))/3;
	}//end of stress calculation
	
	timer::tick("Force_LCAO","start_force",'E');
	return;
}

//print force term for test
void Force_Stress_LCAO::print_force(const string &name, matrix& f, const bool screen, bool ry)const
{
	ofs_running << " --------------------------- " << name << " ----------------------------" << endl;
	ofs_running << " " << setw(8) << "atom" << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;
	
	double fac = 1.0;
	
	if(!ry)
	{
	 	fac = Ry_to_eV / 0.529177;
	}

	cout << setprecision(5);
	cout << setiosflags(ios::showpos);

	if(screen)
	{
		cout << " ------------------- " << name << " --------------------" << endl;
		cout << " " << setw(8) << "atom" << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;
	}

    int iat = 0;
    for (int it = 0;it < ucell.ntype;it++)
    {
        for (int ia = 0;ia < ucell.atoms[it].na;ia++)
        {
			stringstream ss;
			ss << ucell.atoms[it].label << ia+1;

			ofs_running << " " << setw(8) << ss.str();
			if( abs(f(iat,0)) >output_acc) ofs_running << setw(15) << f(iat,0)*fac;
			else ofs_running << setw(15) << "0";
			if( abs(f(iat,1)) >output_acc) ofs_running << setw(15) << f(iat,1)*fac;
			else ofs_running << setw(15) << "0";
			if( abs(f(iat,2)) >output_acc) ofs_running << setw(15) << f(iat,2)*fac;
			else ofs_running << setw(15) << "0";
			ofs_running << endl;

			if(screen)
			{
				cout << " " << setw(8) << ss.str();
				if( abs(f(iat,0)) >output_acc) cout << setw(15) << f(iat,0)*fac;
				else cout << setw(15) << "0";
				if( abs(f(iat,1)) >output_acc) cout << setw(15) << f(iat,1)*fac;
				else cout << setw(15) << "0";
				if( abs(f(iat,2)) >output_acc) cout << setw(15) << f(iat,2)*fac;
				else cout << setw(15) << "0";
				cout << endl;
			}	
				
            iat++;
        }
    }


	cout << resetiosflags(ios::showpos);

    return;
}


void Force_Stress_LCAO::printforce_total (const bool ry, const bool istestf, matrix& fcs)
{
	TITLE("Force_Stress_LCAO","printforce_total");
	double unit_transform = 1;

	if(!ry)
	{
		unit_transform = Ry_to_eV / 0.529177;
	}
//	cout.setf(ios::fixed);

    int iat=0;

	//ofs_running << setiosflags(ios::right);
 	ofs_running << setprecision(6) << setiosflags(ios::showpos) << setiosflags(ios::fixed) << endl;
	NEW_PART("TOTAL-FORCE (eV/Angstrom)");

	// print out forces
	if(INPUT.force_set == 1)
	{
		ofstream ofs("FORCE.dat");
		if(!ofs)
		{
			cout << "open FORCE.dat error !" <<endl;
		}

		for(int iat=0; iat<ucell.nat; iat++)
		{
			ofs << "   " << fcs(iat,0)*Ry_to_eV / 0.529177 
				<< "   " << fcs(iat,1)*Ry_to_eV / 0.529177 
				<< "   " << fcs(iat,2)*Ry_to_eV / 0.529177 << endl;
		}
		ofs.close();
	}

 	if(istestf) 
	{
		cout << setprecision(6) << setiosflags(ios::showpos) << setiosflags(ios::fixed) << endl;
		cout << " ------------------- TOTAL      FORCE --------------------" << endl;
    	cout << " " << setw(8) << "Atom" << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;
    	ofs_running << " " << setw(12) << "Atom" << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;
	}

    iat=0;
    for (int it=0; it<ucell.ntype; it++)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            stringstream ss;
            ss << ucell.atoms[it].label << ia+1;

			if(istestf)
			{
            	cout << " " << setw(8) << ss.str() 
					<< setw(15) << fcs(iat,0)*unit_transform 
					<< setw(15) << fcs(iat,1)*unit_transform 
					<< setw(15) << fcs(iat,2)*unit_transform << endl;
			}

            ofs_running << " " << setw(12) << ss.str() 
				<< setw(15) << fcs(iat,0)*unit_transform 
				<< setw(15) << fcs(iat,1)*unit_transform 
				<< setw(15) << fcs(iat,2)*unit_transform << endl;

            ++iat;
        }
    }
	ofs_running << setiosflags(ios::left);
	cout << resetiosflags(ios::showpos);

    return;
}

//local pseudopotential, ewald, core correction, scc terms in force
void Force_Stress_LCAO::calForcePwPart(
	matrix &fvl_dvl, 
	matrix &fewalds, 
	matrix &fcc, 
	matrix &fscc)
{
	//--------------------------------------------------------
	// local pseudopotential force: 
	// use charge density; plane wave; local pseudopotential;
	//--------------------------------------------------------
	f_pw.cal_force_loc (fvl_dvl);
	//--------------------------------------------------------
	// ewald force: use plane wave only.
	//--------------------------------------------------------
	f_pw.cal_force_ew (fewalds); //remain problem
	//--------------------------------------------------------
	// force due to core correlation. 
	//--------------------------------------------------------
	f_pw.cal_force_cc(fcc);
	//--------------------------------------------------------
	// force due to self-consistent charge.
	//--------------------------------------------------------
	f_pw.cal_force_scc(fscc);	
	return;
}

//overlap, kinetic, nonlocal pseudopotential, Local potential terms in force and stress
void Force_Stress_LCAO::calForceStressIntegralPart(
	const bool isGammaOnly,
	const bool isforce,
	const bool isstress,
	matrix& foverlap,
	matrix& ftvnl_dphi,
	matrix& fvnl_dbeta,	
	matrix& fvl_dphi,
	matrix& soverlap,
	matrix& stvnl_dphi,
	matrix& svnl_dbeta,
	matrix& svl_dphi)
{
	if(isGammaOnly)
	{
    	flk.ftable_gamma(
				isforce,
				isstress,
				foverlap,
				ftvnl_dphi,
				fvnl_dbeta,	
				fvl_dphi,
				soverlap,
				stvnl_dphi,
				svnl_dbeta,
				svl_dphi);
	}
	else
	{
		flk.ftable_k(
				isforce,
				isstress,
				foverlap,
				ftvnl_dphi,
				fvnl_dbeta,	
				fvl_dphi,
				soverlap,
				stvnl_dphi,
				svnl_dbeta,
				svl_dphi);
	}
	return;
}

//vlocal, hartree, ewald, core correction, exchange-correlation terms in stress
void Force_Stress_LCAO::calStressPwPart(
	matrix& sigmadvl,
	matrix& sigmahar,
	matrix& sigmaewa,
	matrix& sigmacc,
	matrix& sigmaxc
)
{
	//--------------------------------------------------------
	// local pseudopotential stress: 
	// use charge density; plane wave; local pseudopotential;
	//--------------------------------------------------------
    sc_pw.stress_loc (sigmadvl, 0);
 
	//--------------------------------------------------------
	//hartree term
	//--------------------------------------------------------
	sc_pw.stress_har (sigmahar, 0);
    
	//--------------------------------------------------------
	// ewald stress: use plane wave only.
	//--------------------------------------------------------
    sc_pw.stress_ewa (sigmaewa, 0); //remain problem


	//--------------------------------------------------------
	// stress due to core correlation. 
	//--------------------------------------------------------
	sc_pw.stress_cc(sigmacc, 0);

	//--------------------------------------------------------
	// stress due to self-consistent charge.
	//--------------------------------------------------------
	for(int i=0;i<3;i++)
	{
		sigmaxc(i,i) =  -(H_XC_pw::etxc) / ucell.omega;
	}
	//Exchange-correlation for PBE
	sc_pw.stress_gga(sigmaxc);
	return;
}

//do symmetry for total force
void Force_Stress_LCAO::forceSymmetry(matrix& fcs)
{
	double *pos;
	double d1,d2,d3;
	pos = new double[ucell.nat*3];
	ZEROS(pos, ucell.nat*3);
	int iat = 0;
	for(int it = 0;it < ucell.ntype;it++)
	{
		for(int ia =0;ia< ucell.atoms[it].na;ia++)
		{
			pos[3*iat  ] = ucell.atoms[it].taud[ia].x ;
			pos[3*iat+1] = ucell.atoms[it].taud[ia].y ;
			pos[3*iat+2] = ucell.atoms[it].taud[ia].z;
			for(int k=0; k<3; ++k)
			{
				symm.check_translation( pos[iat*3+k], -floor(pos[iat*3+k]));
				symm.check_boundary( pos[iat*3+k] );
			}
			iat++;
		}
	}
		
	for(int iat=0; iat<ucell.nat; iat++)
	{
		Mathzone::Cartesian_to_Direct(fcs(iat,0),fcs(iat,1),fcs(iat,2),
							ucell.a1.x, ucell.a1.y, ucell.a1.z,
							ucell.a2.x, ucell.a2.y, ucell.a2.z,
							ucell.a3.x, ucell.a3.y, ucell.a3.z,
							d1,d2,d3);
									
		fcs(iat,0) = d1;fcs(iat,1) = d2;fcs(iat,2) = d3;
	}
	symm.force_symmetry(fcs , pos);
	for(int iat=0; iat<ucell.nat; iat++)
	{
		Mathzone::Direct_to_Cartesian(fcs(iat,0),fcs(iat,1),fcs(iat,2),
							ucell.a1.x, ucell.a1.y, ucell.a1.z,
							ucell.a2.x, ucell.a2.y, ucell.a2.z,
							ucell.a3.x, ucell.a3.y, ucell.a3.z,
							d1,d2,d3);
										
		fcs(iat,0) = d1;fcs(iat,1) = d2;fcs(iat,2) = d3;
	}
	//cout << "nrotk =" << symm.nrotk << endl;
	delete[] pos;
	return;
}
