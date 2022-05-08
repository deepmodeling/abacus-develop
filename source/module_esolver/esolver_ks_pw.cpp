#include "esolver_ks_pw.h"
#include <iostream>
#include "../src_io/wf_io.h"

//--------------temporary----------------------------
#include "../src_pw/global.h"
#include "../module_base/global_function.h"
#include "../module_symmetry/symmetry.h"
#include "../src_pw/vdwd2.h"
#include "../src_pw/vdwd3.h"
#include "../src_pw/vdwd2_parameters.h"
#include "../src_pw/vdwd3_parameters.h"
#include "../src_pw/pw_complement.h"
#include "../src_pw/pw_basis.h"
#include "../src_pw/symmetry_rho.h"
#include "../src_io/print_info.h"
#include "../src_pw/H_Ewald_pw.h"
#include "../src_pw/electrons.h"
#include "../src_pw/occupy.h"
#include "../src_io/chi0_standard.h"
#include "../src_io/chi0_hilbert.h"
#include "../src_io/epsilon0_pwscf.h"
#include "../src_io/epsilon0_vasp.h"
//-----force-------------------
#include "../src_pw/forces.h"
//-----stress------------------
#include "../src_pw/stress_pw.h"
//---------------------------------------------------

namespace ModuleESolver
{

ESolver_KS_PW::ESolver_KS_PW()
{
    classname = "ESolver_KS_PW";
    basisname = "PW";
}

void ESolver_KS_PW::Init(Input &inp, UnitCell_pseudo &ucell)
{
    // setup GlobalV::NBANDS 
	// Yu Liu add 2021-07-03
	GlobalC::CHR.cal_nelec();

	if(GlobalC::ucell.atoms[0].xc_func=="HSE"||GlobalC::ucell.atoms[0].xc_func=="PBE0")
	{
		XC_Functional::set_xc_type("pbe");
	}
	else
	{
		XC_Functional::set_xc_type(GlobalC::ucell.atoms[0].xc_func);
	}

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SETUP UNITCELL");

    //-----------------------------------------------------

    // symmetry analysis should be performed every time the cell is changed
    if (ModuleSymmetry::Symmetry::symm_flag)
    {
        GlobalC::symm.analy_sys(GlobalC::ucell, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
    }

    // Setup the k points according to symmetry.
    GlobalC::kv.set( GlobalC::symm, GlobalV::global_kpoint_card, GlobalV::NSPIN, GlobalC::ucell.G, GlobalC::ucell.latvec );
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"INIT K-POINTS");

    // print information
    // mohan add 2021-01-30
    Print_Info::setup_parameters(GlobalC::ucell, GlobalC::kv);

    // Initalize the plane wave basis set
    GlobalC::pw.gen_pw(GlobalV::ofs_running, GlobalC::ucell, GlobalC::kv);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"INIT PLANEWAVE");
    std::cout << " UNIFORM GRID DIM     : " << GlobalC::pw.nx <<" * " << GlobalC::pw.ny <<" * "<< GlobalC::pw.nz << std::endl;
    std::cout << " UNIFORM GRID DIM(BIG): " << GlobalC::pw.nbx <<" * " << GlobalC::pw.nby <<" * "<< GlobalC::pw.nbz << std::endl;

    // mohan add 2010-09-13
    // initialize the real-space uniform grid for FFT and parallel
    // distribution of plane waves
    GlobalC::Pgrid.init(GlobalC::pw.ncx, GlobalC::pw.ncy, GlobalC::pw.ncz, GlobalC::pw.nczp,
    GlobalC::pw.nrxx, GlobalC::pw.nbz, GlobalC::pw.bz); // mohan add 2010-07-22, update 2011-05-04
        

    // Calculate Structure factor
    GlobalC::pw.setup_structure_factor();
    // cout<<"after pgrid init nrxx = "<<GlobalC::pw.nrxx<<endl;
    
    //----------------------------------------------------------
    // 1 read in initial data:
    //   a lattice structure:atom_species,atom_positions,lattice vector
    //   b k_points
    //   c pseudopotential
    // 2 setup planeware basis, FFT,structure factor, ...
    // 3 initialize local and nonlocal pseudopotential in G_space
    // 4 initialize charge desity and warefunctios in G_space
    //----------------------------------------------------------

    //=====================================
    // init charge/potential/wave functions
    //=====================================
    GlobalC::CHR.allocate(GlobalV::NSPIN, GlobalC::pw.nrxx, GlobalC::pw.ngmc);
    GlobalC::pot.allocate(GlobalC::pw.nrxx);

    GlobalC::wf.allocate(GlobalC::kv.nks);

    // cout<<GlobalC::pw.nrxx<<endl;
    // cout<<"before ufft allocate"<<endl;
    GlobalC::UFFT.allocate();

    // cout<<"after ufft allocate"<<endl;

    //=======================
    // init pseudopotential
    //=======================
    GlobalC::ppcell.init(GlobalC::ucell.ntype);

    //=====================
    // init hamiltonian
    // only allocate in the beginning of ELEC LOOP!
    //=====================
    GlobalC::hm.hpw.allocate(GlobalC::wf.npwx, GlobalV::NPOL, GlobalC::ppcell.nkb, GlobalC::pw.nrxx);

    //=================================
    // initalize local pseudopotential
    //=================================
    GlobalC::ppcell.init_vloc(GlobalC::pw.nggm, GlobalC::ppcell.vloc);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LOCAL POTENTIAL");

    //======================================
    // Initalize non local pseudopotential
    //======================================
    GlobalC::ppcell.init_vnl(GlobalC::ucell);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "NON-LOCAL POTENTIAL");

    //=========================================================
    // calculate the total local pseudopotential in real space
    //=========================================================
    GlobalC::pot.init_pot(0, GlobalC::pw.strucFac); //atomic_rho, v_of_rho, set_vrs

    GlobalC::pot.newd();

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT POTENTIAL");

    //==================================================
    // create GlobalC::ppcell.tab_at , for trial wave functions.
    //==================================================
    GlobalC::wf.init_at_1();

    //================================
    // Initial start wave functions
    //================================
    if (GlobalV::NBANDS != 0 || (GlobalV::CALCULATION != "scf-sto" && GlobalV::CALCULATION != "relax-sto" && GlobalV::CALCULATION != "md-sto")) //qianrui add
    {
        GlobalC::wf.wfcinit();
    }

#ifdef __LCAO
#ifdef __MPI
    switch (GlobalC::exx_global.info.hybrid_type) // Peize Lin add 2019-03-09
    {
    case Exx_Global::Hybrid_Type::HF:
    case Exx_Global::Hybrid_Type::PBE0:
    case Exx_Global::Hybrid_Type::HSE:
        GlobalC::exx_lip.init(&GlobalC::kv, &GlobalC::wf, &GlobalC::pw, &GlobalC::UFFT, &GlobalC::ucell);
        break;
    case Exx_Global::Hybrid_Type::No:
        break;
    case Exx_Global::Hybrid_Type::Generate_Matrix:
    default:
        throw std::invalid_argument(ModuleBase::GlobalFunc::TO_STRING(__FILE__) + ModuleBase::GlobalFunc::TO_STRING(__LINE__));
    }
#endif
#endif

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT BASIS");



}

void ESolver_KS_PW:: beforescf(int istep)
{
    //calculate ewald energy
    H_Ewald_pw::compute_ewald(GlobalC::ucell, GlobalC::pw);
    //Symmetry_rho should be moved to Init()
    Symmetry_rho srho;
    for(int is=0; is<GlobalV::NSPIN; is++)
    {
        srho.begin(is, GlobalC::CHR,GlobalC::pw, GlobalC::Pgrid, GlobalC::symm);
    }
} 

void ESolver_KS_PW:: eachiterinit(const int istep, const int iter)
{
    // mohan add 2010-07-16
    if(iter==1) GlobalC::CHR.set_new_e_iteration(true);
    else GlobalC::CHR.set_new_e_iteration(false);

    if(GlobalV::FINAL_SCF && iter==1) 
    {
        GlobalC::CHR.irstep=0;
        GlobalC::CHR.idstep=0;
        GlobalC::CHR.totstep=0;
    }

    // mohan move harris functional to here, 2012-06-05
	// use 'rho(in)' and 'v_h and v_xc'(in)
	GlobalC::en.calculate_harris(1);

    //(2) save change density as previous charge,
    // prepared fox mixing.
    GlobalC::CHR.save_rho_before_sum_band();
}

//Temporary, it should be replaced by hsolver later.
void ESolver_KS_PW:: hamilt2density(const int istep, const int iter, const double ethr)
{
    GlobalV::PW_DIAG_THR = ethr;
    this->c_bands(istep,iter);
    
    GlobalC::en.eband  = 0.0;
    GlobalC::en.demet  = 0.0;
    GlobalC::en.ef     = 0.0;
    GlobalC::en.ef_up  = 0.0;
    GlobalC::en.ef_dw  = 0.0;

    // calculate weights of each band.
    Occupy::calculate_weights();

    // calculate new charge density according to
    // new wave functions.
    // calculate the new eband here.
    GlobalC::CHR.sum_band();

	// add exx
#ifdef __LCAO
#ifdef __MPI
	GlobalC::en.set_exx();		// Peize Lin add 2019-03-09
#endif
#endif
    // calculate the delta_harris energy
	// according to new charge density.
	// mohan add 2009-01-23
	GlobalC::en.calculate_harris(2);
	Symmetry_rho srho;
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		srho.begin(is, GlobalC::CHR,GlobalC::pw, GlobalC::Pgrid, GlobalC::symm);
	}

    // compute magnetization, only for LSDA(spin==2)
    GlobalC::ucell.magnet.compute_magnetization();
    // deband is calculated from "output" charge density calculated
    // in sum_band
    // need 'rho(out)' and 'vr (v_h(in) and v_xc(in))'

    GlobalC::en.deband = GlobalC::en.delta_e();
    //if (LOCAL_BASIS) xiaohui modify 2013-09-02
}

//Temporary:
void ESolver_KS_PW:: c_bands(const int istep, const int iter)
{
    Electrons elec;
    elec.iter = iter;
    elec.c_bands(istep);
}

//Temporary, it should be rewritten with Hamilt class. 
void ESolver_KS_PW:: updatepot(const int istep, const int iter, const bool conv_elec)
{
    if (!conv_elec)
    {
        // not converged yet, calculate new potential from mixed charge density
        GlobalC::pot.vr = GlobalC::pot.v_of_rho(GlobalC::CHR.rho, GlobalC::CHR.rho_core);
        // because <T+V(ionic)> = <eband+deband> are calculated after sum
        // band, using output charge density.
        // but E_Hartree and Exc(GlobalC::en.etxc) are calculated in v_of_rho above,
        // using the mixed charge density.
        // so delta_escf corrects for this difference at first order.
        GlobalC::en.delta_escf();
    }
    else
    {
        for(int is=0; is<GlobalV::NSPIN; ++is)
        {
            for(int ir=0; ir<GlobalC::pw.nrxx; ++ir)
            {
                GlobalC::pot.vnew(is,ir) = GlobalC::pot.vr(is,ir);
            }
        }
        // the new potential V(PL)+V(H)+V(xc)
        GlobalC::pot.vr = GlobalC::pot.v_of_rho(GlobalC::CHR.rho, GlobalC::CHR.rho_core);
        //std::cout<<"Exc = "<<GlobalC::en.etxc<<std::endl;
        //( vnew used later for scf correction to the forces )
        GlobalC::pot.vnew = GlobalC::pot.vr - GlobalC::pot.vnew;
        GlobalC::en.descf = 0.0;
    }
    GlobalC::pot.set_vr_eff();
}

void ESolver_KS_PW:: eachiterfinish(const int iter, const bool conv_elec)
{
    //print_eigenvalue(GlobalV::ofs_running);
    GlobalC::en.calculate_etot();
    //We output it for restarting the scf.
    bool print = false;
    if(this->out_freq_elec == 0 )
    {
        if(conv_elec) print = true;
    }
    else
    {
        if(iter % this->out_freq_elec == 0 || conv_elec) print = true;
    }
    
    if(print)
    {
        if(GlobalC::CHR.out_chg > 0)
        {
            for(int is=0; is<GlobalV::NSPIN; is++)
            {
                std::stringstream ssc;
                std::stringstream ss1;
                ssc << GlobalV::global_out_dir << "tmp" << "_SPIN" << is + 1 << "_CHG";
                GlobalC::CHR.write_rho(GlobalC::CHR.rho_save[is], is, iter, ssc.str(), 3);//mohan add 2007-10-17
	        	ss1 << GlobalV::global_out_dir << "tmp" << "_SPIN" << is + 1 << "_CHG.cube";
	        	GlobalC::CHR.write_rho_cube(GlobalC::CHR.rho_save[is], is, ssc.str(), 3);
            }
        }
        //output wavefunctions
        if(GlobalC::wf.out_wfc_pw == 1 || GlobalC::wf.out_wfc_pw == 2)
        {
            std::stringstream ssw;
            ssw << GlobalV::global_out_dir << "WAVEFUNC";
            //WF_io::write_wfc( ssw.str(), GlobalC::wf.evc );
            // mohan update 2011-02-21
	    	//qianrui update 2020-10-17
            WF_io::write_wfc2( ssw.str(), GlobalC::wf.evc, GlobalC::pw.gcar);
            //ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"write wave functions into file WAVEFUNC.dat");
        }
        
    }
    
    
}


void ESolver_KS_PW::afterscf(const int iter, const bool conv_elec)
{
#ifdef __LCAO
    if(GlobalC::chi0_hilbert.epsilon)                 // pengfei 2016-11-23
    {
        std::cout <<"eta = "<<GlobalC::chi0_hilbert.eta<<std::endl;
        std::cout <<"domega = "<<GlobalC::chi0_hilbert.domega<<std::endl;
        std::cout <<"nomega = "<<GlobalC::chi0_hilbert.nomega<<std::endl;
        std::cout <<"dim = "<<GlobalC::chi0_hilbert.dim<<std::endl;
        //std::cout <<"oband = "<<GlobalC::chi0_hilbert.oband<<std::endl;
        GlobalC::chi0_hilbert.Chi();
    }
#endif

    if(GlobalC::chi0_standard.epsilon)
    {
        std::cout <<"eta = "<<GlobalC::chi0_standard.eta<<std::endl;
        std::cout <<"domega = "<<GlobalC::chi0_standard.domega<<std::endl;
        std::cout <<"nomega = "<<GlobalC::chi0_standard.nomega<<std::endl;
        std::cout <<"dim = "<<GlobalC::chi0_standard.dim<<std::endl;
        //std::cout <<"oband = "<<GlobalC::chi0_standard.oband<<std::endl;
        GlobalC::chi0_standard.Chi();
    }
    if(GlobalC::epsilon0_pwscf.epsilon)
    {
        GlobalC::epsilon0_pwscf.Cal_epsilon0();
    }
    if(GlobalC::epsilon0_vasp.epsilon)
    {
        GlobalC::epsilon0_vasp.cal_epsilon0();
    }

    for(int is=0; is<GlobalV::NSPIN; is++)
    {
        std::stringstream ssc;
        std::stringstream ss1;
        ssc << GlobalV::global_out_dir << "SPIN" << is + 1 << "_CHG";
		ss1 << GlobalV::global_out_dir << "SPIN" << is + 1 << "_CHG.cube";
        GlobalC::CHR.write_rho(GlobalC::CHR.rho_save[is], is, 0, ssc.str() );//mohan add 2007-10-17
	    GlobalC::CHR.write_rho_cube(GlobalC::CHR.rho_save[is], is, ss1.str(), 3);
    }
    if(conv_elec)
    {
        //GlobalV::ofs_running << " convergence is achieved" << std::endl;			
        //GlobalV::ofs_running << " !FINAL_ETOT_IS " << GlobalC::en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl; 
        GlobalV::ofs_running << "\n charge density convergence is achieved" << std::endl;
        GlobalV::ofs_running << " final etot is " << GlobalC::en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
    }
    else
    {
        GlobalV::ofs_running << " convergence has NOT been achieved!" << std::endl;
    }

    if(GlobalV::OUT_LEVEL != "m")
	{
		this->print_eigenvalue(GlobalV::ofs_running);
	}
}

void ESolver_KS_PW::print_eigenvalue(std::ofstream &ofs)
{
	bool wrong = false;
	for(int ik=0; ik<GlobalC::kv.nks; ++ik)
	{
		for(int ib=0; ib<GlobalV::NBANDS; ++ib)
		{
			if( abs( GlobalC::wf.ekb[ik][ib] ) > 1.0e10)
			{
				GlobalV::ofs_warning << " ik=" << ik+1 << " ib=" << ib+1 << " " << GlobalC::wf.ekb[ik][ib] << " Ry" << std::endl;
				wrong = true;
			}
		}
	}
	if(wrong)
	{
		ModuleBase::WARNING_QUIT("Threshold_Elec::print_eigenvalue","Eigenvalues are too large!");
	}


	if(GlobalV::MY_RANK!=0)
	{
		return;
	}

	ModuleBase::TITLE("Threshold_Elec","print_eigenvalue");

    ofs << "\n STATE ENERGY(eV) AND OCCUPATIONS ";
	ofs << std::setprecision(5);
    for (int ik = 0;ik < GlobalC::kv.nks;ik++)
    {
        if(ik==0)
        {
            ofs << "   NSPIN == " << GlobalV::NSPIN << std::endl;
            if(GlobalV::NSPIN == 2)
            {
                ofs << "SPIN UP : " << std::endl;
            }
        }
        else if(ik == GlobalC::kv.nks/2)
        {
            if(GlobalV::NSPIN == 2)
            {
                ofs << "SPIN DOWN : " << std::endl;
            }
        }

        if (GlobalV::NSPIN==2)
        {
            if (GlobalC::kv.isk[ik] == 0)
            {
                ofs << " " << ik+1 << "/" << GlobalC::kv.nks/2 << " kpoint (Cartesian) = "
                << GlobalC::kv.kvec_c[ik].x << " " << GlobalC::kv.kvec_c[ik].y << " " << GlobalC::kv.kvec_c[ik].z
                << " (" << GlobalC::kv.ngk[ik] << " pws)" << std::endl;

                ofs << std::setprecision(6);

            }
            if (GlobalC::kv.isk[ik] == 1)
            {
                ofs << " " << ik+1-GlobalC::kv.nks/2 << "/" << GlobalC::kv.nks/2 << " kpoint (Cartesian) = "
                << GlobalC::kv.kvec_c[ik].x << " " << GlobalC::kv.kvec_c[ik].y << " " << GlobalC::kv.kvec_c[ik].z
                << " (" << GlobalC::kv.ngk[ik] << " pws)" << std::endl;

                ofs << std::setprecision(6);

			}
		}       // Pengfei Li  added  14-9-9
		else
		{
			ofs << " " << ik+1 << "/" << GlobalC::kv.nks << " kpoint (Cartesian) = "
				<< GlobalC::kv.kvec_c[ik].x << " " << GlobalC::kv.kvec_c[ik].y << " " << GlobalC::kv.kvec_c[ik].z
				<< " (" << GlobalC::kv.ngk[ik] << " pws)" << std::endl;

			ofs << std::setprecision(6);
		}

		//----------------------
		// no energy to output
		//----------------------
		if(GlobalV::KS_SOLVER=="selinv")
		{
			ofs << " USING SELINV, NO BAND ENERGY IS AVAILABLE." << std::endl;
		}
		//----------------------
		// output energy
		//----------------------
		else
		{
			GlobalV::ofs_running << std::setprecision(6);
			GlobalV::ofs_running << std::setiosflags(ios::showpoint);
			for (int ib = 0; ib < GlobalV::NBANDS; ib++)
			{
				ofs << std::setw(8) << ib+1 
				    << std::setw(15) << GlobalC::wf.ekb[ik][ib] * ModuleBase::Ry_to_eV 
                    << std::setw(15) << GlobalC::wf.wg(ik, ib) << std::endl;
			}
			ofs << std::endl;
		}
    }//end ik
    return;
}



void ESolver_KS_PW::cal_Energy(energy &en)
{

}

void ESolver_KS_PW::cal_Force(ModuleBase::matrix &force)
{
	Forces ff;
	ff.init(force);
}
void ESolver_KS_PW::cal_Stress(ModuleBase::matrix &stress)
{
	Stress_PW ss;
	ss.cal_stress(stress);
}
void ESolver_KS_PW::postprocess()
{
	// compute density of states
    GlobalC::en.perform_dos_pw();
}

}