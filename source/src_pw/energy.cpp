#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../src_parallel/parallel_reduce.h"
#include "global.h"
#include "energy.h"
#include "../module_base/mymath.h"
#include <vector>
#ifdef __MPI
#include "mpi.h"
#endif
#include <sys/time.h>
#ifdef __LCAO
#include "../src_lcao/dftu.h"  //Quxin adds for DFT+U on 20201029
#endif
#include "myfunc.h"
//new
#include "H_Ewald_pw.h"
#include "H_Hartree_pw.h"
#include "../module_surchem/efield.h"    // liuyu add 2022-05-06
#include "../module_surchem/gatefield.h"    // liuyu add 2022-09-13
#include "../module_surchem/surchem.h"
#ifdef __DEEPKS
#include "../module_deepks/LCAO_deepks.h"
#endif


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
	this->esol_el = 0;		   // the implicit solvation energy Ael
	this->esol_cav = 0;		   // the implicit solvation energy Acav
}

energy::~energy()
{
}

void energy::calculate_harris(const int &flag)
{
//	ModuleBase::TITLE("energy","calculate_harris");
	
	if(flag==1)
	{
		this->deband_harris = this->delta_e();
	}
	else if(flag==2)
	{
		this->etot_harris = eband + deband_harris 
		+ (etxc - etxcc) 
		+ H_Ewald_pw::ewald_energy 
		+ H_Hartree_pw::hartree_energy 
		+ demet
		+ exx
		+ Efield::etotefield
        + Gatefield::etotgatefield
		+ evdw;  						// Peize Lin add evdw 2021.03.09

#ifdef __LCAO
        if(INPUT.dft_plus_u) 
		{
			this->etot_harris += GlobalC::dftu.EU;  //Energy correction from DFT+U; Quxin adds on 20201029
		}
#endif
#ifdef __DEEPKS
        if (GlobalV::deepks_scf)
        {
            this->etot_harris += GlobalC::ld.E_delta - GlobalC::ld.e_delta_band;
        }
#endif
    }
	
	return;
}

void energy::calculate_etot(void)
{
	ModuleBase::TITLE("energy","calculate_etot");
	//std::cout << "\n demet in etot = " << demet << std::endl;
	this->etot = eband + deband 
	+ (etxc - etxcc) 
	+ H_Ewald_pw::ewald_energy 
	+ H_Hartree_pw::hartree_energy 
	+ demet
	+ descf
	+ exx
	+ Efield::etotefield
    + Gatefield::etotgatefield
	+ evdw;							// Peize Lin add evdw 2021.03.09
	if (GlobalV::imp_sol)
    {
	this->etot += GlobalC::solvent_model.cal_Ael(GlobalC::ucell, GlobalC::rhopw)
				 + GlobalC::solvent_model.cal_Acav(GlobalC::ucell, GlobalC::rhopw);
	}
    //Quxin adds for DFT+U energy correction on 20201029

	// std::cout << std::resetiosflags(ios::scientific) << std::endl;
	// std::cout << std::setprecision(16) << std::endl;
	// std::cout << " eband=" << eband << std::endl;
	// std::cout << " deband=" << deband << std::endl;
	// std::cout << " etxc-etxcc=" <<etxc-etxcc << std::endl;
	// std::cout << " ewld=" << H_Ewald_pw::ewald_energy << std::endl;
	// std::cout << " ehart=" << H_Hartree_pw::hartree_energy << std::endl;
	// std::cout << " demet=" << demet << std::endl;
	// std::cout << " descf=" << descf << std::endl;
	// std::cout << " exx=" << exx << std::endl;
	// std::cout << " efiled=" << Efield::etotefield << std::endl;
	// std::cout << " total= "<<etot<<std::endl;
	// std::cout << " fermienergy= "<<ef<<std::endl;

#ifdef __LCAO
    if(INPUT.dft_plus_u) 
	{
		this->etot += GlobalC::dftu.EU;																	  
	}
#endif
#ifdef __DEEPKS
	if (GlobalV::deepks_scf)
	{
		this->etot += GlobalC::ld.E_delta - GlobalC::ld.e_delta_band;
	}
#endif
	return;
}

void energy::print_etot(
	const bool converged, 
	const int &iter_in, 
	const double &scf_thr, 
	const double &duration, 
	const double &pw_diag_thr, 
	const double &avg_iter,
	bool print)
{
	ModuleBase::TITLE("energy","print_etot");
	this->iter = iter_in;

	GlobalV::ofs_running << std::setprecision(12);
	GlobalV::ofs_running << std::setiosflags(ios::left);

	GlobalV::ofs_running << "\n Density error is " << scf_thr << std::endl;

	if(GlobalV::OUT_LEVEL != "m") //xiaohui add "OUT_LEVEL", 2015-09-16
	{
		if(GlobalV::BASIS_TYPE=="pw")ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Error Threshold",pw_diag_thr); //xiaohui add 2013-09-02

		if (this->printe > 0 && ((iter + 1) % this->printe == 0 || converged || iter == GlobalV::SCF_NMAX))
		{
			GlobalV::ofs_running << "\n " << std::setw(12) << "Energy" << std::setw(30) << "Rydberg" << std::setw(30) << "eV" << std::endl;
			this->print_format("E_KohnSham", etot);
			this->print_format("E_Harris", etot_harris);
			this->print_format("E_band", eband);
			this->print_format("E_one_elec", eband + deband);
			this->print_format("E_Hartree", H_Hartree_pw::hartree_energy);
			this->print_format("E_xc", etxc - etxcc);
			this->print_format("E_Ewald", H_Ewald_pw::ewald_energy);
			this->print_format("E_demet", demet); //mohan add 2011-12-02
			this->print_format("E_descf", descf);
			if (GlobalC::vdwd2_para.flag_vdwd2)					//Peize Lin add 2014-04, update 2021-03-09
			{
				this->print_format("E_vdwD2", evdw);
			}
			if(GlobalC::vdwd3_para.flag_vdwd3)					//jiyy add 2019-05, update 2021-05-02
			{
				this->print_format("E_vdwD3", evdw);
			}
			this->print_format("E_exx", exx);	
					
			if (GlobalV::imp_sol)
        	{
				esol_el = GlobalC::solvent_model.cal_Ael(GlobalC::ucell, GlobalC::rhopw);
            	esol_cav = GlobalC::solvent_model.cal_Acav(GlobalC::ucell, GlobalC::rhopw);
				this->print_format("E_sol_el", esol_el);
				this->print_format("E_sol_cav", esol_cav);
			}
            if(GlobalV::EFIELD_FLAG)
            {
                this->print_format("E_efield", Efield::etotefield);
            }
            if(GlobalV::GATE_FLAG)
            {
                this->print_format("E_gatefield", Gatefield::etotgatefield);
            }

#ifdef __DEEPKS
			if (GlobalV::deepks_scf)	//caoyu add 2021-08-10
			{
				this->print_format("E_DeePKS", GlobalC::ld.E_delta);
			}
#endif
		}
			else
		{
			GlobalV::ofs_running << "\n " << std::setw(12) << "Energy" << std::setw(30) << "Rydberg" << std::setw(30) << "eV" << std::endl;
			this->print_format("E_KohnSham",etot);
			this->print_format("E_Harris",etot_harris);
		}

		if(GlobalV::TWO_EFERMI)
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
	std::stringstream ss;

	//xiaohui add 2013-09-02, Peize Lin update 2020.11.14
    std::string label;
	if(GlobalV::KS_SOLVER=="cg")
	{
		label = "CG";
	}
	else if (GlobalV::KS_SOLVER=="lapack")
	{
		label = "LA";
	}
    else if(GlobalV::KS_SOLVER=="hpseps")
	{
		label = "HP";
	}
    else if(GlobalV::KS_SOLVER=="genelpa")
	{
        label = "GE";
	}
	else if(GlobalV::KS_SOLVER=="dav")
	{
		label = "DA";
	}
    else if(GlobalV::KS_SOLVER=="scalapack_gvx")
	{
        label = "GV";
	}
	else if(GlobalV::KS_SOLVER=="cusolver")
	{
        label = "CU";
	}
	else
	{
		ModuleBase::WARNING_QUIT("Energy","print_etot");
	}
    ss << label << iter;
	//xiaohui add 2013-09-02

	bool scientific=true;
	int prec = 6;


	if(!print) return;

	if(GlobalV::OUT_LEVEL=="ie" || GlobalV::OUT_LEVEL=="m") //xiaohui add 'm' option, 2015-09-16
	{
		std::cout << " " << std::setw(7) << ss.str();
		//std::cout << std::setiosflags(ios::fixed);
		//std::cout << std::setiosflags(ios::showpos);
		if(scientific)
		{
			std::cout << std::setiosflags(ios::scientific);
		}

		if(GlobalV::COLOUR)
		{
			if(GlobalV::MY_RANK==0)
			{
				printf( "\e[36m%-15f\e[0m", GlobalC::en.etot);	
				//printf( "[36m%-15f[0m", GlobalC::en.etot);	
				if(GlobalV::NSPIN==2)
				{
					std::cout << std::setprecision(2);
					std::cout<<std::setw(10)<<GlobalC::ucell.magnet.tot_magnetization;
					std::cout<<std::setw(10)<<GlobalC::ucell.magnet.abs_magnetization;
				}
				else if(GlobalV::NSPIN==4 && GlobalV::NONCOLIN)
				{
					std::cout << std::setprecision(2);
					std::cout<<std::setw(10)<<GlobalC::ucell.magnet.tot_magnetization_nc[0]
					<<std::setw(10)<<GlobalC::ucell.magnet.tot_magnetization_nc[1]
					<<std::setw(10)<<GlobalC::ucell.magnet.tot_magnetization_nc[2];
					std::cout<<std::setw(10)<<GlobalC::ucell.magnet.abs_magnetization;
				}
				if(scf_thr>1.0)
				{
					// 31 is red
					printf( "\e[31m%-14e\e[0m", scf_thr);
					//printf( "[31m%-14e[0m", scf_thr);
				}
				else
				{
					// 32 is green
					printf( "\e[32m%-14e\e[0m", scf_thr);
					//printf( "[32m%-14e[0m", scf_thr);
				}
				// 34 is blue
				printf( "\e[36m%-15f\e[0m", GlobalC::en.etot*ModuleBase::Ry_to_eV);	
				//printf( "[36m%-15f[0m", GlobalC::en.etot*ModuleBase::Ry_to_eV);	
				std::cout << std::setprecision(3);
	//			std::cout << std::setw(11) << GlobalC::en.eband;
	//			std::cout << std::setw(11) << H_Hartree_pw::hartree_energy;
	//			std::cout << std::setw(11) << GlobalC::en.etxc - GlobalC::en.etxcc;
				std::cout << std::resetiosflags(ios::scientific);
				//if(GlobalV::DIAGO_TYPE=="cg") xiaohui modify 2013-09-02
				// if(GlobalV::KS_SOLVER=="cg") //xiaohui add 2013-09-02
				// {
				// 	std::cout << std::setw(11) << avg_iter;
				// }
				//xiaohui modified 2013-03-23
				//else if(GlobalV::DIAGO_TYPE=="selinv")
				//{
					// because Selinv::iter starts from 0.
				//	std::cout << std::setw(11) << Selinv::iter;
				//}
				std::cout << std::setw(11) << duration;
				std::cout << std::endl;
			}
		}
		else
		{
			std::cout << std::setprecision(prec);
			//std::cout << std::setw(15) << GlobalC::en.etot;
			if(GlobalV::NSPIN==2)
			{
				std::cout << std::setprecision(2);
				std::cout<<std::setw(10)<<GlobalC::ucell.magnet.tot_magnetization;
				std::cout<<std::setw(10)<<GlobalC::ucell.magnet.abs_magnetization;
			}
			std::cout << std::setprecision(6);
			std::cout << std::setw(15) << GlobalC::en.etot*ModuleBase::Ry_to_eV;
                        std::cout << std::setw(15) << (GlobalC::en.etot - GlobalC::en.etot_old) *ModuleBase::Ry_to_eV;  //pengfei Li added 2015-1-31
                        std::cout << std::setprecision(3);
                        std::cout << std::setw(11) << scf_thr;
			std::cout << std::setprecision(3);
	//		std::cout << std::setw(11) << GlobalC::en.eband;
	//		std::cout << std::setw(11) << H_Hartree_pw::hartree_energy;
	//		std::cout << std::setw(11) << GlobalC::en.etxc - GlobalC::en.etxcc;
			//if(GlobalV::DIAGO_TYPE=="cg") xiaohui modify 2013-09-02
			// if(GlobalV::KS_SOLVER=="cg") //xiaohui add 2013-09-02
			// {
			// 	std::cout << std::setw(11) << avg_iter;
			// }
			//xiaohui modified 2013-03-23
			//else if(GlobalV::DIAGO_TYPE=="selinv")
			//{
				// because Selinv::iter starts from 0.
			//	std::cout << std::setw(11) << Selinv::iter+1;
			//}
			std::cout << std::setw(11) << duration;
			std::cout << std::endl;
		}

	}
	else
	{
	}

    this->etot_old = this->etot;
	return;
}

void energy::print_format(const std::string &name, const double &value)
{
	GlobalV::ofs_running << std::setiosflags(ios::showpos);
	std::stringstream name2;
	name2 << name;
	GlobalV::ofs_running << " " << std::setw(12) << name2.str() << std::setw(30) <<  value 
	<< std::setw(30) << value * ModuleBase::Ry_to_eV << std::endl;
	GlobalV::ofs_running << std::resetiosflags(ios::showpos);
	return;
}


// from ddelta_e.f90
double energy::delta_e(void)
{
    // out potentials from potential mixing
    // total energy and band energy corrections
	double deband0 = 0.0;

    double deband_aux = 0.0;

    for (int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
    {
    	deband_aux -= GlobalC::CHR.rho[0][ir] * GlobalC::pot.vr(0, ir);
		if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
		{
			deband_aux -= GlobalC::CHR.kin_r[0][ir] * GlobalC::pot.vofk(0,ir);
		}
	}

    if (GlobalV::NSPIN == 2)
    {
    	for (int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
    	{
    		deband_aux -= GlobalC::CHR.rho[1][ir] * GlobalC::pot.vr(1, ir);
			if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
			{
				deband_aux -= GlobalC::CHR.kin_r[1][ir] * GlobalC::pot.vofk(1,ir);
			}
		}
    }
    else if(GlobalV::NSPIN == 4)
    {
        for (int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
        {
            deband_aux -= GlobalC::CHR.rho[1][ir] * GlobalC::pot.vr(1, ir);
            deband_aux -= GlobalC::CHR.rho[2][ir] * GlobalC::pot.vr(2, ir);
            deband_aux -= GlobalC::CHR.rho[3][ir] * GlobalC::pot.vr(3, ir);
        }
    }

#ifdef __MPI
    MPI_Allreduce(&deband_aux,&deband0,1,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
#else
    deband0 = deband_aux;
#endif

    deband0 *= GlobalC::ucell.omega / GlobalC::rhopw->nxyz;
	
	// \int rho(r) v_{exx}(r) dr = 2 E_{exx}[rho]
	deband0 -= 2*exx;				// Peize Lin add 2017-10-16
	
    return deband0;
} // end subroutine delta_e



void energy::delta_escf(void)
{
	ModuleBase::TITLE("energy","delta_escf");
    this->descf = 0.0;

	// now rho1 is "mixed" charge density
	// and rho1_save is "output" charge density
	// because in "deband" the energy is calculated from "output" charge density,
	// so here is the correction.

    for (int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
    {
		this->descf -= ( GlobalC::CHR.rho[0][ir] - GlobalC::CHR.rho_save[0][ir] ) * GlobalC::pot.vr(0, ir);
		if(XC_Functional::get_func_type() == 3)
		{
         	this->descf -= ( GlobalC::CHR.kin_r[0][ir] - GlobalC::CHR.kin_r_save[0][ir] ) * GlobalC::pot.vofk(0, ir);
		}
    }

    if (GlobalV::NSPIN==2)
    {
       	for (int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
       	{
           	this->descf -= ( GlobalC::CHR.rho[1][ir] - GlobalC::CHR.rho_save[1][ir] ) * GlobalC::pot.vr(1, ir);
			if(XC_Functional::get_func_type() == 3)
			{
           		this->descf -= ( GlobalC::CHR.kin_r[1][ir] - GlobalC::CHR.kin_r_save[1][ir] ) * GlobalC::pot.vofk(1, ir);
			}
       	}
    }
    if (GlobalV::NSPIN==4)
    {
        for(int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
        {
            this->descf -= ( GlobalC::CHR.rho[1][ir] - GlobalC::CHR.rho_save[1][ir] ) * GlobalC::pot.vr(1, ir);
            this->descf -= ( GlobalC::CHR.rho[2][ir] - GlobalC::CHR.rho_save[2][ir] ) * GlobalC::pot.vr(2, ir);
            this->descf -= ( GlobalC::CHR.rho[3][ir] - GlobalC::CHR.rho_save[3][ir] ) * GlobalC::pot.vr(3, ir);
        }
    }

    Parallel_Reduce::reduce_double_pool( descf );

    this->descf *= GlobalC::ucell.omega / GlobalC::rhopw->nxyz;
    return;
}


void energy::print_band(const int &ik)
{
	//check the band energy.
    bool wrong = false;
	for(int ib=0; ib<GlobalV::NBANDS; ++ib)
	{
		if( abs( GlobalC::wf.ekb[ik][ib] ) > 1.0e10)
		{
			GlobalV::ofs_warning << " ik=" << ik+1 << " ib=" << ib+1 << " " << GlobalC::wf.ekb[ik][ib] << " Ry" << std::endl;
			wrong = true;
		}
	}
	if(wrong)
    {
        ModuleBase::WARNING_QUIT("Threshold_Elec::print_eigenvalue","Eigenvalues are too large!");
    }



	if(GlobalV::MY_RANK==0)
	{
		//if( GlobalV::DIAGO_TYPE == "selinv" ) xiaohui modify 2013-09-02
		if(GlobalV::KS_SOLVER=="selinv") //xiaohui add 2013-09-02
		{
			GlobalV::ofs_running << " No eigenvalues are available for selected inversion methods." << std::endl;	
		}
		else
		{
			if( printe>0 && ((this->iter+1) % this->printe == 0))
			{
				//	NEW_PART("ENERGY BANDS (Rydberg), (eV)");
				GlobalV::ofs_running << std::setprecision(6);
				GlobalV::ofs_running << " Energy (eV) & Occupations  for spin=" << GlobalV::CURRENT_SPIN+1 << " K-point=" << ik+1 << std::endl;
				GlobalV::ofs_running << std::setiosflags(ios::showpoint);
				for(int ib=0;ib<GlobalV::NBANDS;ib++)
				{
					GlobalV::ofs_running << " "<< std::setw(6) << ib+1  
						<< std::setw(15) << GlobalC::wf.ekb[ik][ib] * ModuleBase::Ry_to_eV;
					// for the first electron iteration, we don't have the energy
					// spectrum, so we can't get the occupations. 
					GlobalV::ofs_running << std::setw(15) << GlobalC::wf.wg(ik,ib);
					GlobalV::ofs_running << std::endl;
				}
			}
		}
	}
	return;
}

// Peize Lin add 2016-12-03
#ifdef __LCAO
#ifdef __MPI
void energy::set_exx()
{
	ModuleBase::TITLE("energy", "set_exx");

	auto exx_energy = []() -> double
	{
		if("lcao_in_pw"==GlobalV::BASIS_TYPE)
		{
			return GlobalC::exx_lip.get_exx_energy();
		}
		else if("lcao"==GlobalV::BASIS_TYPE)
		{
			return GlobalC::exx_lcao.get_energy();
		}
		else
		{
			throw std::invalid_argument(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}
	};
	if( Exx_Global::Hybrid_Type::HF   == GlobalC::exx_lcao.info.hybrid_type )				// HF
	{
		this->exx = exx_energy();
	}
	else if( Exx_Global::Hybrid_Type::PBE0 == GlobalC::exx_lcao.info.hybrid_type ||
			Exx_Global::Hybrid_Type::SCAN0 == GlobalC::exx_lcao.info.hybrid_type || 
			 Exx_Global::Hybrid_Type::HSE  == GlobalC::exx_lcao.info.hybrid_type )			// PBE0 or HSE
	{
		this->exx = GlobalC::exx_global.info.hybrid_alpha * exx_energy();
	}

	return;
}
#endif //__MPI
#endif //_LCAO
