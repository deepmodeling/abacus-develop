#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "global.h"
#include "electrons.h"
#include "../src_pw/symmetry_rho.h"
#include "../src_io/print_info.h"
#include "../src_io/wf_io.h"
#include "../src_io/chi0_hilbert.h"
#include "../src_io/chi0_standard.h"
#include "../src_io/epsilon0_pwscf.h"
#include "../src_io/epsilon0_vasp.h"
#include "../src_io/to_wannier90.h"
#include "../src_io/berryphase.h"
// new
#include "H_Ewald_pw.h"

double Electrons::avg_iter = 0;

Electrons::Electrons()
{
    iter = 0;
    test = 0;
    delta_total_energy = 0.0;
}

Electrons::~Electrons()
{
}

void Electrons::non_self_consistent(const int &istep)
{
    ModuleBase::TITLE("Electrons","non_self_consistent");
    ModuleBase::timer::tick("Electrons","non_self_consistent");

    //========================================
    // diagonalization of the KS hamiltonian
    // =======================================
    Electrons::c_bands(istep);

    GlobalV::ofs_running << "\n End of Band Structure Calculation \n" << std::endl;


    for (int ik = 0; ik < GlobalC::kv.nks; ik++)
    {
        if (GlobalV::NSPIN==2)
        {
            if (ik == 0) GlobalV::ofs_running << " spin up :" << std::endl;
            if (ik == ( GlobalC::kv.nks / 2)) GlobalV::ofs_running << " spin down :" << std::endl;
        }
        //out.printV3(GlobalV::ofs_running, GlobalC::kv.kvec_c[ik]);

        GlobalV::ofs_running << " k-points" << ik+1
        << "(" << GlobalC::kv.nkstot << "): "
        << GlobalC::kv.kvec_c[ik].x
        << " " << GlobalC::kv.kvec_c[ik].y
        << " " << GlobalC::kv.kvec_c[ik].z << std::endl;

        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            GlobalV::ofs_running << " spin" << GlobalC::kv.isk[ik]+1
            << "_final_band " << ib+1
            << " " << GlobalC::wf.ekb[ik][ib] * ModuleBase::Ry_to_eV
            << " " << GlobalC::wf.wg(ik, ib)*GlobalC::kv.nks << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }

    // add by jingan in 2018.11.7
    if(GlobalV::CALCULATION == "nscf" && INPUT.towannier90)
    {
        toWannier90 myWannier(GlobalC::kv.nkstot,GlobalC::ucell.G);
        myWannier.init_wannier();
    }

    //=======================================================
    // Do a Berry phase polarization calculation if required
    //=======================================================

    if (berryphase::berry_phase_flag && ModuleSymmetry::Symmetry::symm_flag == 0)
    {
        berryphase bp;
        bp.Macroscopic_polarization();
    }

    ModuleBase::timer::tick("Electrons","non_self_consistent");
    return;
}

/*
#include "occupy.h"
void Electrons::self_consistent(const int &istep)
{
    ModuleBase::timer::tick("Electrons","self_consistent");

	// mohan update 2021-02-25
	H_Ewald_pw::compute_ewald(GlobalC::ucell, GlobalC::pw);

    set_pw_diag_thr();

    if(GlobalV::OUT_LEVEL=="ie" || GlobalV::OUT_LEVEL=="m")
    {
        if(GlobalV::COLOUR && GlobalV::MY_RANK==0)
        {
            printf( "\e[33m%-7s\e[0m", "ITER");
            printf( "\e[33m%-15s\e[0m", "ETOT(Ry)");
            if(GlobalV::NSPIN==2)
            {
                printf( "\e[33m%-10s\e[0m", "TMAG");
                printf( "\e[33m%-10s\e[0m", "AMAG");
            }
            printf( "\e[33m%-14s\e[0m", "SCF_THR");
            printf( "\e[33m%-15s\e[0m", "ETOT(eV)");
            if(GlobalV::KS_SOLVER=="cg")
            {
                printf( "\e[33m%-15s\e[0m", "CG_ITER");
            }
            printf( "\e[33m%-11s\e[0m\n", "TIME(s)");
        }
        else
        {
            std::cout << " " << std::setw(7) << "ITER";

            if(GlobalV::NSPIN==2)
            {
                std::cout << std::setw(10) << "TMAG";
                std::cout << std::setw(10) << "AMAG";
            }

            std::cout << std::setw(15) << "ETOT(eV)";
            std::cout << std::setw(15) << "EDIFF(eV)";
            std::cout << std::setw(11) << "SCF_THR";
            if(GlobalV::KS_SOLVER == "cg")
            {
                std::cout << std::setw(11) << "CG_ITER";
            }
			std::cout << std::setw(11) << "TIME(s)" << std::endl;
		}
	}

    Symmetry_rho srho;
    for(int is=0; is<GlobalV::NSPIN; is++)
    {
        srho.begin(is, GlobalC::CHR,GlobalC::pw, GlobalC::Pgrid, GlobalC::symm);
    }

    // conv_elec is a member of Threshold_Elec
    this->conv_elec = false;//mohan add 2008-05-25


    // mohan add 2010/3/25
    // output the charge mixing data :
    // iteration && scf_thr.
    // std::stringstream ss;
    // ss << GlobalV::global_out_dir << "ChargeMixing.dat";
    // std::ofstream ofs_mix;

    // if(GlobalV::MY_RANK==0)
    // {
    //     ofs_mix.open(ss.str().c_str());
    // }
    // ###

    for (this->iter = 1;iter <= GlobalV::SCF_NMAX;iter++)
    {
        GlobalV::ofs_running
        << "\n PW ALGORITHM --------------- ION=" << std::setw(4) << istep + 1
        << "  ELEC=" << std::setw(4) << iter
        << "--------------------------------\n";

        // mohan add 2010-07-16
        if(iter==1) GlobalC::CHR.set_new_e_iteration(true);
        else GlobalC::CHR.set_new_e_iteration(false);

        // record the start time.
		// the clock is not accurate, needs to be fixed 2021-03-15 mohan
        clock_t start=std::clock();

        //(1) set converged threshold,
        // automatically updated during self consistency.
        //this->update_pw_diag_thr(iter);
        if(GlobalV::FINAL_SCF && iter==1) GlobalV::PW_DIAG_THR = 1.0e-2;
        else this->update_pw_diag_thr(iter);
        if(GlobalV::FINAL_SCF && iter==1)
        {
            init_mixstep_final_scf();
            //GlobalC::CHR.irstep=0;
            //GlobalC::CHR.idstep=0;
            //GlobalC::CHR.totstep=0;
        }

		// mohan move harris functional to here, 2012-06-05
		// use 'rho(in)' and 'v_h and v_xc'(in)
		GlobalC::en.calculate_harris(1);

		// first_iter_again:					// Peize Lin delete 2019-05-01

		// calculate exact-exchange
#ifdef __MPI //liyuanbo 2022/2/23
#ifdef __LCAO
		if( Exx_Global::Hybrid_Type::HF   == GlobalC::exx_lcao.info.hybrid_type || 
			Exx_Global::Hybrid_Type::PBE0 == GlobalC::exx_lcao.info.hybrid_type || 
			Exx_Global::Hybrid_Type::HSE  == GlobalC::exx_lcao.info.hybrid_type )
        {
            if( !GlobalC::exx_global.info.separate_loop )
            {
                GlobalC::exx_lip.cal_exx();
            }
            break;
		}
#endif
#endif
        //(2) save change density as previous charge,
        // prepared fox mixing.
        GlobalC::CHR.save_rho_before_sum_band();

		bool onescf = false;
    scf_step:
		//(3) calculate band energy using cg or davidson method.
		// output the new eigenvalues and wave functions.
        this->c_bands(istep);

        if (check_stop_now())
        {
            ModuleBase::timer::tick("Electrons","self_consistent");
            return;
        } 

        GlobalC::en.eband  = 0.0;
        GlobalC::en.demet  = 0.0;
        GlobalC::en.ef     = 0.0;
        GlobalC::en.ef_up  = 0.0;
        GlobalC::en.ef_dw  = 0.0;

        //(4) calculate weights of each band.
        Occupy::calculate_weights();

        //(5) calculate new charge density according to
        // new wave functions.

        // calculate the new eband here.
        GlobalC::CHR.sum_band();


		// add exx
#ifdef __LCAO
#ifdef __MPI
		GlobalC::en.set_exx();		// Peize Lin add 2019-03-09
#endif
#endif

		//(6) calculate the delta_harris energy
		// according to new charge density.
		// mohan add 2009-01-23
		GlobalC::en.calculate_harris(2);

		Symmetry_rho srho;
		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			srho.begin(is, GlobalC::CHR,GlobalC::pw, GlobalC::Pgrid, GlobalC::symm);
		}

        //(7) compute magnetization, only for LSDA(spin==2)
        GlobalC::ucell.magnet.compute_magnetization();

        //(8) deband is calculated from "output" charge density calculated
        // in sum_band
        // need 'rho(out)' and 'vr (v_h(in) and v_xc(in))'
        GlobalC::en.deband = GlobalC::en.delta_e();

        //if (LOCAL_BASIS) xiaohui modify 2013-09-02
		if(GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw") //xiaohui add 2013-09-02
        {
            GlobalC::CHR.tmp_mixrho(scf_thr,0,GlobalV::SCF_THR,iter,conv_elec);
        }
        else
        {
            // tr2_min used only in first scf iteraton
            double diago_error = 0.0;
            if(iter==1)
            {
                // if 'scf_thr < GlobalV::PW_DIAG_THR * nelec' happen,
                // in other word, 'scf_thr < diago_error'
                // we update GlobalV::PW_DIAG_THR.
                diago_error = GlobalV::PW_DIAG_THR*std::max(1.0, GlobalC::CHR.nelec);
            }

            // if converged is achieved, or the self-consistent error(scf_thr)
            // is samller than the estimated error due to diagonalization(diago_error)
            // rhoin and rho are unchanged:
            // rhoin contain the input charge density and
            // rho contain the output charge density.
            // in other cases rhoin contains the mixed charge density
            // (the new input density) while rho is unchanged.
            GlobalC::CHR.tmp_mixrho(scf_thr,diago_error,GlobalV::SCF_THR,iter,conv_elec);

            //if(GlobalV::MY_RANK==0)
            //{
            //    ofs_mix << std::setw(5) << iter << std::setw(20) << scf_thr << std::endl;
            //}

            if (iter==1 && !onescf)
            {
                onescf = true;
                if (scf_thr < diago_error)
                {
                    GlobalV::ofs_running << " Notice: Threshold on eigenvalues was too large.\n";

                    ModuleBase::WARNING("scf","Threshold on eigenvalues was too large.");
                    GlobalV::ofs_running << " scf_thr=" << scf_thr << " < diago_error=" << diago_error << std::endl;

                    // update GlobalV::PW_DIAG_THR.
                    GlobalV::ofs_running << " Origin PW_DIAG_THR = " << GlobalV::PW_DIAG_THR << std::endl;
                    GlobalV::PW_DIAG_THR = 0.1 * scf_thr / GlobalC::CHR.nelec;
                    GlobalV::ofs_running << " New    PW_DIAG_THR = " << GlobalV::PW_DIAG_THR << std::endl;
                    //goto first_iter_again;
                    goto scf_step;
                }
            }
        }

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
            // mohan add 2012-06-05
            for(int is=0; is<GlobalV::NSPIN; ++is)
            {
                for(int ir=0; ir<GlobalC::pw.nrxx; ++ir)
                {
                    GlobalC::pot.vnew(is,ir) = GlobalC::pot.vr(is,ir);
                }
            }

            // mohan fix bug 2012-06-05,
            // the new potential V(PL)+V(H)+V(xc)
            GlobalC::pot.vr = GlobalC::pot.v_of_rho(GlobalC::CHR.rho, GlobalC::CHR.rho_core);
            //std::cout<<"Exc = "<<GlobalC::en.etxc<<std::endl;
            //( vnew used later for scf correction to the forces )
            GlobalC::pot.vnew = GlobalC::pot.vr - GlobalC::pot.vnew;
            GlobalC::en.descf = 0.0;
        }


        // output for tmp.
        for(int is=0; is<GlobalV::NSPIN; is++)
        {
            std::stringstream ssc;
            std::stringstream ss1;
            ssc << GlobalV::global_out_dir << "tmp" << "_SPIN" << is + 1 << "_CHG";
            GlobalC::CHR.write_rho(GlobalC::CHR.rho_save[is], is, iter, ssc.str(), 3);//mohan add 2007-10-17
			ss1 << GlobalV::global_out_dir << "tmp" << "_SPIN" << is + 1 << "_CHG.cube";
			GlobalC::CHR.write_rho_cube(GlobalC::CHR.rho_save[is], is, ssc.str(), 3);
        }

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

		GlobalC::pot.set_vr_eff();

        //print_eigenvalue(GlobalV::ofs_running);
        GlobalC::en.calculate_etot();

		// the clock is not accurate, needs to be fixed 2021-03-15 mohan
        clock_t finish=clock();
        double duration = (double)(finish - start) / CLOCKS_PER_SEC;

		GlobalC::en.print_etot(conv_elec, iter, scf_thr, duration, GlobalV::PW_DIAG_THR, avg_iter);

        if (conv_elec || iter==GlobalV::SCF_NMAX)
        {

            //--------------------------------------
            // output charge density for converged,
            // 0 means don't need to consider iter,
            //--------------------------------------
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
				print_eigenvalue(GlobalV::ofs_running);
			}
            ModuleBase::timer::tick("Electrons","self_consistent");
            return;
        }

        //if ( imix >= 0 )  GlobalC::CHR.rho = GlobalC::CHR.rho_save;
        //GlobalV::ofs_running << "\n start next iterate for idum ";
    } //END DO

    ModuleBase::timer::tick("Electrons","self_consistent");
    return;
} // end Electrons
*/

bool Electrons::check_stop_now(void)
{
    bool check_stop_now = false;

    if (check_stop_now)
    {
        conv_elec = false;
    }

    return check_stop_now;
} // END FUNCTION check_stop_now


void Electrons::c_bands(const int &istep)
{
    if (GlobalV::test_elec) ModuleBase::TITLE("Electrons","c_bands");
    ModuleBase::timer::tick("Electrons", "c_bands"
    );

    int precondition_type = 2;

    double *h_diag = new double[GlobalC::wf.npwx * GlobalV::NPOL];//added by zhengdy-soc
    ModuleBase::GlobalFunc::ZEROS(h_diag, GlobalC::wf.npwx * GlobalV::NPOL);

    avg_iter = 0.0;

    GlobalV::ofs_running << " "  <<std::setw(8) << "K-point" << std::setw(15) << "CG iter num" << std::setw(15) << "Time(Sec)"<< std::endl;
    GlobalV::ofs_running << std::setprecision(6) << std::setiosflags(ios::fixed) << std::setiosflags(ios::showpoint);
    for (int ik = 0;ik < GlobalC::kv.nks;ik++)
    {
        GlobalC::hm.hpw.init_k(ik);
        // GlobalC::hm.hpw_gpu.init_k(ik);

        //===========================================
        // Conjugate-Gradient diagonalization
        // h_diag is the precondition matrix
        // h_diag(1:npw) = MAX( 1.0, g2kin(1:npw) );
        //===========================================
        if (precondition_type==1)
        {
            for (int ig = 0;ig < GlobalC::wf.npw; ig++)
            {
                h_diag[ig] = std::max(1.0, GlobalC::wf.g2kin[ig]);
                if(GlobalV::NPOL==2) h_diag[ig+GlobalC::wf.npwx] = h_diag[ig];
            }
        }
        else if (precondition_type==2)
        {
            for (int ig = 0;ig < GlobalC::wf.npw; ig++)
            {
                h_diag[ig] = 1 + GlobalC::wf.g2kin[ig] + sqrt( 1 + (GlobalC::wf.g2kin[ig] - 1) * (GlobalC::wf.g2kin[ig] - 1));
                if(GlobalV::NPOL==2) h_diag[ig+GlobalC::wf.npwx] = h_diag[ig];
            }
        }
        //h_diag can't be zero!  //zhengdy-soc
		if(GlobalV::NPOL==2)
        {
            for(int ig = GlobalC::wf.npw;ig < GlobalC::wf.npwx; ig++)
            {
                h_diag[ig] = 1.0;
                h_diag[ig+ GlobalC::wf.npwx] = 1.0;
            }
        }

        clock_t start=clock();

        //============================================================
        // diago the hamiltonian!!
        // In plane wave method, firstly using diagH_subspace to diagnolize,
        // then using cg method.
        //
        // In localized orbital presented in plane wave case,
        // only use diagH_subspace.
        //=============================================================
        double avg_iter_k = 0.0;
        GlobalC::hm.diagH_pw(istep, this->iter, ik, h_diag, avg_iter_k);

        avg_iter += avg_iter_k;

        GlobalC::en.print_band(ik); //mohan add 2012-04-16

        clock_t finish=clock();
        const double duration = static_cast<double>(finish - start) / CLOCKS_PER_SEC;

        GlobalV::ofs_running << " " << std::setw(8)
        << ik+1 << std::setw(15)
        << avg_iter_k << std::setw(15) << duration << std::endl;
    }//End K Loop

    //if (!LOCAL_BASIS) xiaohui modify 2013-09-02
    if(GlobalV::BASIS_TYPE=="pw") //xiaohui add 2013-09-02
    {
        //GlobalV::ofs_running << " avg_iteri " << avg_iter << std::endl;
        Parallel_Reduce::reduce_double_allpool(avg_iter); //mohan fix bug 2012-06-05
        //GlobalV::ofs_running << " avg_iter_after " << avg_iter << std::endl;
        avg_iter /= static_cast<double>(GlobalC::kv.nkstot);
    }
    delete [] h_diag;
    ModuleBase::timer::tick("Electrons","c_bands");
    return;
} // END SUBROUTINE c_bands_k


void Electrons::init_mixstep_final_scf(void)
{
    ModuleBase::TITLE("electrons","init_mixstep_final_scf");

    GlobalC::CHR.irstep=0;
    GlobalC::CHR.idstep=0;
    GlobalC::CHR.totstep=0;

    return;
}
