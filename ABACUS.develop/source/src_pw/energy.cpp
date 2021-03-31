//==========================================================
// AUTHOR : Lixin He, mohan
// DATE : 2008-11-21
//==========================================================
#include "tools.h"
#include "global.h"
#include "energy.h"
#include "src_global/mymath.h"
#include <vector>
#include <mpi.h>
#include <sys/time.h>
#include "../src_pw/tools.h"
#include "../src_lcao/dftu.h"  //Quxin adds for DFT+U on 20201029
#include "src_pw/myfunc.h"
//new
#include "H_Ewald_pw.h"
#include "H_Hartree_pw.h"
#include "H_XC_pw.h"


energy::energy()
{
    // the maximum number of R vectors included in r
    // the square of the electron charge (Ry atomic units)
    this->etot   = 0;          // the total energy of the solid
	this->etot_harris = 0;	   // total energy of harris functional
    this->eband  = 0;          // the band energy
    this->deband = 0;          // correction for variational energy
	this->deband_harris = 0;   // harris energy
    this->etxcc  = 0;          // the nlcc exchange and correlation
	this->exx    = 0;          // the exact exchange energy.

    this->demet  = 0;          // correction for metals
    this->ef     = 0;          // the fermi energy
}

energy::~energy()
{
}

#include "efield.h"

void energy::calculate_harris(const int &flag)
{
//	TITLE("energy","calculate_harris");
	
	if(flag==1)
	{
		this->deband_harris = this->delta_e();
	}
	else if(flag==2)
	{
		this->etot_harris = eband + deband_harris 
		+ (H_XC_pw::etxc - etxcc) 
		+ H_Ewald_pw::ewald_energy 
		+ H_Hartree_pw::hartree_energy 
		+ demet
		+ exx
		+ Efield::etotefield
		+ evdw;							// Peize Lin add evdw 2021.03.09

        if(INPUT.dft_plus_u) 
		{
			this->etot_harris += dftu.EU;  //Energy correction from DFT+U; Quxin adds on 20201029
		}
	}
	
	return;
}

void energy::calculate_etot(void)
{
	TITLE("energy","calculate_etot");
	//cout << "\n demet in etot = " << demet << endl;
	this->etot = eband + deband 
	+ (H_XC_pw::etxc - etxcc) 
	+ H_Ewald_pw::ewald_energy 
	+ H_Hartree_pw::hartree_energy 
	+ demet
	+ descf
	+ exx
	+ Efield::etotefield
	+ evdw;							// Peize Lin add evdw 2021.03.09

    //Quxin adds for DFT+U energy correction on 20201029

	/*cout << resetiosflags(ios::scientific) << endl;
	cout << setprecision(6) << endl;
	cout << " eband=" << eband << endl;
	cout << " deband=" << deband << endl;
	cout << " etxc-etxcc=" << etxc-etxcc << endl;
	cout << " ewld=" << ewld << endl;
	cout << " ehart=" << ehart << endl;
	cout << " demet=" << demet << endl;
	cout << " descf=" << descf << endl;
	cout << " efiled=" << Efield::etotefield << endl;
	cout << " fermienergy= "<<ef<<endl;*/
    if(INPUT.dft_plus_u) 
	{
		this->etot += dftu.EU;																	  
	}
	return;	
}

void energy::print_etot(
	const bool converged, 
	const int &istep, 
	const int &iter_in, 
	const double &dr2, 
	const double &duration, 
	const double &ethr, 
	const double &avg_iter,
	bool print)
{
	TITLE("energy","print_etot");
	this->iter = iter_in;

	ofs_running << setprecision(12);
	ofs_running << setiosflags(ios::left);

	ofs_running << "\n Density error is " << dr2 << endl;

	if(OUT_LEVEL != "m") //xiaohui add "OUT_LEVEL", 2015-09-16
	{
		if(BASIS_TYPE=="pw")OUT(ofs_running,"Error Threshold",ethr); //xiaohui add 2013-09-02

		if( this->printe>0 && ( (iter+1) % this->printe == 0 || converged || iter == NITER) )	
		{
			ofs_running << "\n " << setw(12) << "Energy" << setw(30) << "Rydberg" << setw(30) << "eV" << endl;
			this->print_format("E_KohnSham",etot);
			this->print_format("E_Harris",etot_harris);
			this->print_format("E_band",eband);
			this->print_format("E_one_elec",eband+deband);
			this->print_format("E_Hartree",H_Hartree_pw::hartree_energy);
			this->print_format("E_xc",H_XC_pw::etxc-etxcc);
			this->print_format("E_Ewald",H_Ewald_pw::ewald_energy);
			this->print_format("E_demet",demet); //mohan add 2011-12-02
			this->print_format("E_descf",descf);
			this->print_format("E_efield",Efield::etotefield);
			if(vdwd2_para.flag_vdwd2)					//Peize Lin add 2014-04, update 2021-03-09
			{
				this->print_format("E_vdwD2",evdw);
			}
			if(vdwd3.vdwD3)					//jiyy add 2019-05
			{
				this->print_format("E_vdwD3",vdwd3.energy_result);
			}
			this->print_format("E_exx",exx);
		}
		else
		{
			ofs_running << "\n " << setw(12) << "Energy" << setw(30) << "Rydberg" << setw(30) << "eV" << endl;
			this->print_format("E_KohnSham",etot);
			this->print_format("E_Harris",etot_harris);
		}

		if(TWO_EFERMI)
		{
			this->print_format("E_Fermi_up",ef_up);
			this->print_format("E_Fermi_dw",ef_dw);
		}
		else
		{
			this->print_format("E_Fermi",this->ef);
		}
	}//xiaohui add "OUT_LEVEL", 2015-09-16

	if (iter_in == 1)   // pengfei Li added 2015-1-31
	{
		this->etot_old = this->etot;
	}
	
	// mohan update 2011-02-26
	stringstream ss;

	//xiaohui add 2013-09-02, Peize Lin update 2020.11.14
    string label;
	if(KS_SOLVER=="cg")
	{
		label = "CG";
	}
	else if (KS_SOLVER=="lapack")
	{
		label = "LA";
	}
    else if(KS_SOLVER=="hpseps")
	{
		label = "HP";
	}
    else if(KS_SOLVER=="genelpa")
	{
        label = "GE";
	}
	else if(KS_SOLVER=="dav")
	{
		label = "DA";
	}
    else if(KS_SOLVER=="scalapack_gvx")
	{
        label = "GV";
	}
	else
	{
		WARNING_QUIT("Energy","print_etot");
	}
    ss << label << iter;
	//xiaohui add 2013-09-02

	bool scientific=true;
	int prec = 6;


	if(!print) return;

	if(OUT_LEVEL=="ie" || OUT_LEVEL=="m") //xiaohui add 'm' option, 2015-09-16
	{
		cout << " " << setw(7) << ss.str();
		//cout << setiosflags(ios::fixed);
		//cout << setiosflags(ios::showpos);
		if(scientific)
		{
			cout << setiosflags(ios::scientific);
		}

		if(COLOUR)
		{
			if(MY_RANK==0)
			{
				//printf( "\e[36m%-15f\e[0m", en.etot);	
				printf( "[36m%-15f[0m", en.etot);	
				if(NSPIN==2)
				{
					cout << setprecision(2);
					cout<<setw(10)<<mag.tot_magnetization;
					cout<<setw(10)<<mag.abs_magnetization;
				}
				else if(NSPIN==4 && NONCOLIN)
				{
					cout << setprecision(2);
					cout<<setw(10)<<mag.tot_magnetization_nc[0]
					<<setw(10)<<mag.tot_magnetization_nc[1]
					<<setw(10)<<mag.tot_magnetization_nc[2];
					cout<<setw(10)<<mag.abs_magnetization;
				}
				if(dr2>1.0)
				{
					// 31 is red
					//printf( "\e[31m%-14e\e[0m", dr2);
					printf( "[31m%-14e[0m", dr2);
				}
				else
				{
					// 32 is green
					//printf( "\e[32m%-14e\e[0m", dr2);
					printf( "[32m%-14e[0m", dr2);
				}
				// 34 is blue
				//printf( "\e[36m%-15f\e[0m", en.etot*Ry_to_eV);	
				printf( "[36m%-15f[0m", en.etot*Ry_to_eV);	
				cout << setprecision(3);
	//			cout << setw(11) << en.eband;
	//			cout << setw(11) << H_Hartree_pw::hartree_energy;
	//			cout << setw(11) << en.etxc - en.etxcc;
				cout << resetiosflags(ios::scientific);
				//if(DIAGO_TYPE=="cg") xiaohui modify 2013-09-02
				if(KS_SOLVER=="cg") //xiaohui add 2013-09-02
				{
					cout << setw(11) << avg_iter;
				}
				//xiaohui modified 2013-03-23
				//else if(DIAGO_TYPE=="selinv")
				//{
					// because Selinv::iter starts from 0.
				//	cout << setw(11) << Selinv::iter;
				//}
				cout << setw(11) << duration;
				cout << endl;
			}
		}
		else
		{
			cout << setprecision(prec);
			//cout << setw(15) << en.etot;
			if(NSPIN==2)
			{
				cout << setprecision(2);
				cout<<setw(10)<<mag.tot_magnetization;
				cout<<setw(10)<<mag.abs_magnetization;
			}
			cout << setprecision(6);
			cout << setw(15) << en.etot*Ry_to_eV;
                        cout << setw(15) << (en.etot - en.etot_old) *Ry_to_eV;  //pengfei Li added 2015-1-31
                        cout << setprecision(3);
                        cout << setw(11) << dr2;
			cout << setprecision(3);
	//		cout << setw(11) << en.eband;
	//		cout << setw(11) << H_Hartree_pw::hartree_energy;
	//		cout << setw(11) << en.etxc - en.etxcc;
			//if(DIAGO_TYPE=="cg") xiaohui modify 2013-09-02
			if(KS_SOLVER=="cg") //xiaohui add 2013-09-02
			{
				cout << setw(11) << avg_iter;
			}
			//xiaohui modified 2013-03-23
			//else if(DIAGO_TYPE=="selinv")
			//{
				// because Selinv::iter starts from 0.
			//	cout << setw(11) << Selinv::iter+1;
			//}
			cout << setw(11) << duration;
			cout << endl;
		}

	}
	else
	{
	}

    this->etot_old = this->etot;
	return;
}

void energy::print_format(const string &name, const double &value)
{
	ofs_running << setiosflags(ios::showpos);
	stringstream name2;
	name2 << name;
	ofs_running << " " << setw(12) << name2.str() << setw(30) <<  value 
	<< setw(30) << value * Ry_to_eV << endl;
	ofs_running << resetiosflags(ios::showpos);
	return;
}


// from ddelta_e.f90
double energy::delta_e(void)
{
    // out potentials from potential mixing
    // total energy and band energy corrections
	double deband0 = 0.0;

    double deband_aux = 0.0;

    for (int ir=0; ir<pw.nrxx; ir++) deband_aux -= CHR.rho[0][ir] * pot.vr(0, ir);

    if (NSPIN == 2)
    {
        for (int ir=0; ir<pw.nrxx; ir++)
        {
            deband_aux -= CHR.rho[1][ir] * pot.vr(1, ir);
        }

    }
    else if(NSPIN == 4)
    {
        for (int ir=0; ir<pw.nrxx; ir++)
        {
            deband_aux -= CHR.rho[1][ir] * pot.vr(1, ir);
            deband_aux -= CHR.rho[2][ir] * pot.vr(2, ir);
            deband_aux -= CHR.rho[3][ir] * pot.vr(3, ir);
        }
    }

#ifdef __MPI
    MPI_Allreduce(&deband_aux,&deband0,1,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
#else
    deband0 = deband_aux;
#endif

    deband0 *= ucell.omega / pw.ncxyz;
	
	// \int rho(r) v_{exx}(r) dr = 2 E_{exx}[rho]
	deband0 -= 2*exx;				// Peize Lin add 2017-10-16
	
    return deband0;
} // end subroutine delta_e



void energy::delta_escf(void)
{
	TITLE("energy","delta_escf");
    this->descf = 0.0;

	// now rho1 is "mixed" charge density
	// and rho1_save is "output" charge density
	// because in "deband" the energy is calculated from "output" charge density,
	// so here is the correction.
    for (int ir=0; ir<pw.nrxx; ir++) 
	{
		this->descf -= ( CHR.rho[0][ir]- CHR.rho_save[0][ir] ) * pot.vr(0,ir);
	}

    if (NSPIN==2)
    {
        for (int ir=0; ir<pw.nrxx; ir++)
        {
            this->descf -= ( CHR.rho[1][ir] - CHR.rho_save[1][ir] ) * pot.vr(1, ir);
        }
    }
    if (NSPIN==4)
    {
        for(int ir=0; ir<pw.nrxx; ir++)
        {
            this->descf -= ( CHR.rho[1][ir] - CHR.rho_save[1][ir] ) * pot.vr(1, ir);
            this->descf -= ( CHR.rho[2][ir] - CHR.rho_save[2][ir] ) * pot.vr(2, ir);
            this->descf -= ( CHR.rho[3][ir] - CHR.rho_save[3][ir] ) * pot.vr(3, ir);
        }
    }

    Parallel_Reduce::reduce_double_pool( descf );

    this->descf *= ucell.omega / pw.ncxyz;
    return;
}


void energy::print_band(const int &ik)
{
	//check the band energy.
    bool wrong = false;
	for(int ib=0; ib<NBANDS; ++ib)
	{
		if( abs( wf.ekb[ik][ib] ) > 1.0e10)
		{
			ofs_warning << " ik=" << ik+1 << " ib=" << ib+1 << " " << wf.ekb[ik][ib] << " Ry" << endl;
			wrong = true;
		}
	}
	if(wrong)
    {
        WARNING_QUIT("Threshold_Elec::print_eigenvalue","Eigenvalues are too large!");
    }



	if(MY_RANK==0)
	{
		//if( DIAGO_TYPE == "selinv" ) xiaohui modify 2013-09-02
		if(KS_SOLVER=="selinv") //xiaohui add 2013-09-02
		{
			ofs_running << " No eigenvalues are available for selected inversion methods." << endl;	
		}
		else
		{
			if( printe>0 && ((this->iter+1) % this->printe == 0))
			{
				//	NEW_PART("ENERGY BANDS (Rydberg), (eV)");
				ofs_running << setprecision(6);
				ofs_running << " Energy (eV) & Occupations  for spin=" << CURRENT_SPIN+1 << " K-point=" << ik+1 << endl;
				ofs_running << setiosflags(ios::showpoint);
				for(int ib=0;ib<NBANDS;ib++)
				{
					ofs_running << " "<< setw(6) << ib+1  
						<< setw(15) << wf.ekb[ik][ib] * Ry_to_eV;
					// for the first electron iteration, we don't have the energy
					// spectrum, so we can't get the occupations. 
					ofs_running << setw(15) << wf.wg(ik,ib);
					ofs_running << endl;
				}
			}
		}
	}
	return;
}

// Peize Lin add 2016-12-03
void energy::set_exx()
{
	TITLE("energy", "set_exx");

	auto exx_energy = []() -> double
	{
		if("lcao_in_pw"==BASIS_TYPE)
		{
			return exx_lip.get_exx_energy();
		}
		else if("lcao"==BASIS_TYPE)
		{
			return exx_lcao.get_energy();
		}
		else
		{
			throw invalid_argument(TO_STRING(__FILE__)+TO_STRING(__LINE__));
		}
	};

	if( 5==xcf.iexch_now && 0==xcf.igcx_now )				// HF
	{
		this->exx = exx_energy();
	}
	else if( 6==xcf.iexch_now && 8==xcf.igcx_now )			// PBE0
	{
		this->exx = exx_global.info.hybrid_alpha * exx_energy();
	}
	else if( 9==xcf.iexch_now && 12==xcf.igcx_now )			// HSE
	{
		this->exx = exx_global.info.hybrid_alpha * exx_energy();
	}

	return;
}
