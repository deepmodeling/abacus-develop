#include "tools.h"
#include "global.h"
#include "electrons.h"
#include "symmetry_rho.h"
#include "src_io/wf_io.h"
#include "src_io/chi0_hilbert.h"
#include "src_io/chi0_standard.h"
#include "src_io/epsilon0_pwscf.h"
#include "src_io/epsilon0_vasp.h"
#include "src_io/to_wannier90.h"
#include "src_io/berryphase.h"
// new
#include "H_Ewald_pw.h"

double Electrons::avg_iter = 0;

Electrons::Electrons()
{
}

Electrons::~Electrons()
{
}

void Electrons::non_self_consistent(const int &istep)
{
    TITLE("Electrons","non_self_consistent");
    timer::tick("Electrons","non_self_consistent",'D');

    //========================================
    // diagonalization of the KS hamiltonian
    // =======================================
    Electrons::c_bands(istep);

    ofs_running << " End of Band Structure Calculation " << endl;
    ofs_running << " Band eigenvalue (eV) :" << endl;

    for (int ik = 0; ik < kv.nks; ik++)
    {
        if (NSPIN==2)
        {
            if (ik == 0) ofs_running << " spin up :" << endl;
            if (ik == ( kv.nks / 2)) ofs_running << " spin down :" << endl;
        }
        //out.printV3(ofs_running, kv.kvec_c[ik]);

        ofs_running << " k-points" << ik+1 
        << "(" << kv.nkstot << "): " 
        << kv.kvec_c[ik].x 
        << " " << kv.kvec_c[ik].y 
        << " " << kv.kvec_c[ik].z << endl;

        for (int ib = 0; ib < NBANDS; ib++)
        {			
            ofs_running << " spin" << kv.isk[ik]+1 
            << "_final_band " << ib+1 
            << " " << wf.ekb[ik][ib] * Ry_to_eV 
            << " " << wf.wg(ik, ib)*kv.nks << endl;
        }
        ofs_running << endl;
    }

    // add by jingan in 2018.11.7
    if(CALCULATION == "nscf" && INPUT.towannier90)
    {
        toWannier90 myWannier(kv.nkstot,ucell.G);
        myWannier.init_wannier();
    }

    //=======================================================
    // Do a Berry phase polarization calculation if required
    //=======================================================

    if (berryphase::berry_phase_flag && Symmetry::symm_flag == 0)
    {
        berryphase bp;
        bp.Macroscopic_polarization();
    }

    timer::tick("Electrons","non_self_consistent",'D');
    return;
}


#include "occupy.h"
void Electrons::self_consistent(const int &istep)
{
    timer::tick("Electrons","self_consistent",'D');

	// mohan update 2021-02-25
	H_Ewald_pw::compute_ewald(ucell,pw); 

    set_ethr();

    this->unit = 0;

    if(OUT_LEVEL=="ie")
    {
        cout << setprecision(12);
        cout<< " " << setw(7)<< "ITER"; // pengfei Li added 2015-1-31 

        if(NSPIN==2)
        {
            cout<<setw(10)<<"TMAG";
            cout<<setw(10)<<"AMAG";
        }

        cout<<setw(15)<< "ETOT(eV)"<<setw(15)<< "EDIFF(eV)"<<setw(11)<< "DRHO2"; // pengfei Li added 2015-1-31
        //if(DIAGO_TYPE=="cg") xiaohui modify 2013-09-02
        if(KS_SOLVER=="cg") //xiaohui add 2013-09-02
        {
            cout<<setw(11)<<"CG_ITER";
        }

        cout<<setw(11)<< "TIME(S)";
        cout<<endl;
    }
    else
    {

    }

    Symmetry_rho srho;
    for(int is=0; is<NSPIN; is++)
    {
        srho.begin(is);
    }

    // conv_elec is a member of Threshold_Elec
    this->conv_elec = false;//mohan add 2008-05-25


    // mohan add 2010/3/25
    // output the charge mixing data :
    // iteration && dr2.
    // stringstream ss;
    // ss << global_out_dir << "ChargeMixing.dat"; 
    // ofstream ofs_mix;

    // if(MY_RANK==0)
    // {
    //     ofs_mix.open(ss.str().c_str());
    // }
    // ###

    for (this->iter = 1;iter <= NITER;iter++)
    {
        ofs_running 
        << "\n PW ALGORITHM --------------- ION=" << setw(4) << istep + 1
        << "  ELEC=" << setw(4) << iter 
        << "--------------------------------\n";
        // mohan add 2010-07-16
        if(iter==1) CHR.set_new_e_iteration(true);
        else CHR.set_new_e_iteration(false);

        // record the start time.
		// the clock is not accurate, needs to be fixed 2021-03-15 mohan
        clock_t start=std::clock();

        //(1) set converged threshold, 
        // automatically updated during self consistency.
        //this->update_ethr(iter);
        if(FINAL_SCF && iter==1) ETHR = 1.0e-2;
        else this->update_ethr(iter);
        if(FINAL_SCF && iter==1)
        {
            init_mixstep_final_scf();
            //CHR.irstep=0;
            //CHR.idstep=0;
            //CHR.totstep=0;
        }

		// mohan move harris functional to here, 2012-06-05
		// use 'rho(in)' and 'v_h and v_xc'(in)
		en.calculate_harris(1);
	
		// first_iter_again:					// Peize Lin delete 2019-05-01
		
		// calculate exact-exchange
		switch(xcf.iexch_now)						// Peize Lin add 2019-03-09
		{
			case 5:    case 6:   case 9:
				if( !exx_global.info.separate_loop )				
				{
					exx_lip.cal_exx();			
				}
				break;
		}
		
		//(2) calculate band energy using cg or davidson method.
		// output the new eigenvalues and wave functions.
        this->c_bands(istep);

        if (check_stop_now()) return;

        en.eband  = 0.0;
        en.demet  = 0.0;
        en.ef     = 0.0;
        en.ef_up  = 0.0;
        en.ef_dw  = 0.0;

        //(3) calculate weights of each band.
        Occupy::calculate_weights();	

        //(4) save change density as previous charge,
        // prepared fox mixing.
        CHR.save_rho_before_sum_band();

        //(5) calculate new charge density according to
        // new wave functions.
		
        // calculate the new eband here.
        CHR.sum_band();

		// add exx
		en.set_exx();		// Peize Lin add 2019-03-09
		
		//(6) calculate the delta_harris energy 
		// according to new charge density.
		// mohan add 2009-01-23
		en.calculate_harris(2);

		Symmetry_rho srho;
		for(int is=0; is<NSPIN; is++)
		{
			srho.begin(is);
		}

        //(7) compute magnetization, only for LSDA(spin==2)
        mag.compute_magnetization();

        //(8) deband is calculated from "output" charge density calculated 
        // in sum_band
        // need 'rho(out)' and 'vr (v_h(in) and v_xc(in))'
        en.deband = en.delta_e();

        //if (LOCAL_BASIS) xiaohui modify 2013-09-02
		if(BASIS_TYPE=="lcao" || BASIS_TYPE=="lcao_in_pw") //xiaohui add 2013-09-02
        {
            CHR.mix_rho(dr2,0,DRHO2,iter,conv_elec);
        }
        else
        {
            // tr2_min used only in first scf iteraton
            double diago_error = 0.0;
            if(iter==1) 
            {
                // if 'dr2 < ETHR * nelec' happen,
                // in other word, 'dr2 < diago_error'
                // we update ETHR.
                diago_error = ETHR*std::max(1.0, ucell.nelec);
            }

            // if converged is achieved, or the self-consistent error(dr2)
            // is samller than the estimated error due to diagonalization(diago_error)
            // rhoin and rho are unchanged:
            // rhoin contain the input charge density and 
            // rho contain the output charge density.
            // in other cases rhoin contains the mixed charge density
            // (the new input density) while rho is unchanged.
            CHR.mix_rho(dr2,diago_error,DRHO2,iter,conv_elec);

            //if(MY_RANK==0)
            //{
            //    ofs_mix << setw(5) << iter << setw(20) << dr2 << endl; 
            //}

            if (iter==1)
            {
                if (dr2 < diago_error)
                {
                    ofs_running << " Notice: Threshold on eigenvalues was too large.\n";

                    WARNING("scf","Threshold on eigenvalues was too large.");
                    ofs_running << " dr2=" << dr2 << " < diago_error=" << diago_error << endl;

                    // update ETHR.
                    ofs_running << " Origin ETHR = " << ETHR << endl;
                    ETHR = dr2 / ucell.nelec;
                    ofs_running << " New    ETHR = " << ETHR << endl;
                    //goto first_iter_again;
                }
            }
        }

        if (!conv_elec)
        {
            // not converged yet, calculate new potential from mixed charge density
            pot.v_of_rho(CHR.rho, pot.vr);

            // because <T+V(ionic)> = <eband+deband> are calculated after sum
            // band, using output charge density.
            // but E_Hartree and Exc(en.etxc) are calculated in v_of_rho above,
            // using the mixed charge density.
            // so delta_escf corrects for this difference at first order. 
            en.delta_escf();
        }
        else
        {
            // mohan add 2012-06-05
            for(int is=0; is<NSPIN; ++is)
            {
                for(int ir=0; ir<pw.nrxx; ++ir)
                {
                    pot.vnew(is,ir) = pot.vr(is,ir);
                }
            }

            // mohan fix bug 2012-06-05,
            // the new potential V(PL)+V(H)+V(xc)
            pot.v_of_rho(CHR.rho, pot.vr);
            //cout<<"Exc = "<<en.etxc<<endl;
            //( vnew used later for scf correction to the forces )
            pot.vnew = pot.vr - pot.vnew;
            en.descf = 0.0;
        }

        stringstream ssw;
        ssw << global_out_dir << "WAVEFUNC.dat";

		//qianrui add 2020-10-12
		stringstream ssgk;
		ssgk << global_out_dir << "GKK.dat";

        // output for tmp.
        for(int is=0; is<NSPIN; is++)
        {
            stringstream ssc;
            ssc << global_out_dir << "tmp" << "_SPIN" << is + 1 << "_CHG";
            CHR.write_rho(CHR.rho_save[is], is, iter, ssc.str(), 3);//mohan add 2007-10-17
        }

        if(wf.out_wf)
        {
            //WF_io::write_wfc( ssw.str(), wf.evc );
            // mohan update 2011-02-21
			//qianrui update 2020-10-17
            WF_io::write_wfc2( ssw.str(), wf.evc, pw.gcar);
            //DONE(ofs_running,"write wave functions into file WAVEFUNC.dat");
        }

        if(vext == 0) 
		{
			pot.set_vr_eff();
		}
        else 
		{
			pot.set_vrs_tddft(istep);
		}

        //print_eigenvalue(ofs_running);
        en.calculate_etot();

		// the clock is not accurate, needs to be fixed 2021-03-15 mohan
        clock_t finish=clock();
        double duration = (double)(finish - start) / CLOCKS_PER_SEC;

		en.print_etot(conv_elec, istep, iter, dr2, duration, ETHR, avg_iter);

        if (conv_elec || iter==NITER)
        {
            //--------------------------------------
            // output charge density for converged,
            // 0 means don't need to consider iter,
            //--------------------------------------
            if(chi0_hilbert.epsilon)                 // pengfei 2016-11-23
            {
                cout <<"eta = "<<chi0_hilbert.eta<<endl;
                cout <<"domega = "<<chi0_hilbert.domega<<endl;
                cout <<"nomega = "<<chi0_hilbert.nomega<<endl;
                cout <<"dim = "<<chi0_hilbert.dim<<endl;
                //cout <<"oband = "<<chi0_hilbert.oband<<endl;
                chi0_hilbert.Chi();
            }

            if(chi0_standard.epsilon)
            {
                cout <<"eta = "<<chi0_standard.eta<<endl;
                cout <<"domega = "<<chi0_standard.domega<<endl;
                cout <<"nomega = "<<chi0_standard.nomega<<endl;
                cout <<"dim = "<<chi0_standard.dim<<endl;
                //cout <<"oband = "<<chi0_standard.oband<<endl;
                chi0_standard.Chi();
            }

            if(epsilon0_pwscf.epsilon)
            {
                epsilon0_pwscf.Cal_epsilon0();
            }

            if(epsilon0_vasp.epsilon)
            {
                epsilon0_vasp.cal_epsilon0();
            }

            for(int is=0; is<NSPIN; is++)
            {
                stringstream ssc;
                ssc << global_out_dir << "SPIN" << is + 1 << "_CHG";
                CHR.write_rho(CHR.rho_save[is], is, 0, ssc.str() );//mohan add 2007-10-17
            }

            if(conv_elec)
            {
                //ofs_running << " convergence is achieved" << endl;			
                //ofs_running << " !FINAL_ETOT_IS " << en.etot * Ry_to_eV << " eV" << endl; 
                ofs_running << " charge density convergence is achieved" << endl;
                ofs_running << " final etot is " << en.etot * Ry_to_eV << " eV" << endl;
            }
            else
            {
                ofs_running << " convergence has NOT been achieved!" << endl;			
            }

            iter_end(ofs_running);
            timer::tick("Electrons","self_consistent",'D');
            return;
        }

        //if ( imix >= 0 )  CHR.rho = CHR.rho_save;
        //ofs_running << "\n start next iterate for idum ";
    } //END DO

    timer::tick("Electrons","self_consistent",'D');
    return;
} // end Electrons


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
    if (test_elec) TITLE("Electrons","c_bands");
    timer::tick("Electrons","c_bands",'E');

    int precondition_type = 2;

    double *h_diag = new double[wf.npwx * NPOL];//added by zhengdy-soc
    ZEROS(h_diag, wf.npwx * NPOL);

    avg_iter = 0.0;

    ofs_running << " "  <<setw(8) << "K-point" << setw(15) << "CG iter num" << setw(15) << "Time(Sec)"<< endl;
    ofs_running << setprecision(6) << setiosflags(ios::fixed) << setiosflags(ios::showpoint);
    for (int ik = 0;ik < kv.nks;ik++)
    {
        hm.hpw.init_k(ik);

        //===========================================
        // Conjugate-Gradient diagonalization
        // h_diag is the precondition matrix
        // h_diag(1:npw) = MAX( 1.0, g2kin(1:npw) );
        //===========================================
        if (precondition_type==1)
        {
            for (int ig = 0;ig < wf.npw; ig++)
            {
                h_diag[ig] = max(1.0, wf.g2kin[ig]);
                if(NPOL==2) h_diag[ig+wf.npwx] = h_diag[ig];
            }
        }
        else if (precondition_type==2)
        {
            for (int ig = 0;ig < wf.npw; ig++)
            {
                h_diag[ig] = 1 + wf.g2kin[ig] + sqrt( 1 + (wf.g2kin[ig] - 1) * (wf.g2kin[ig] - 1));
                if(NPOL==2) h_diag[ig+wf.npwx] = h_diag[ig];
            }
        }
        //h_diag can't be zero!  //zhengdy-soc
		if(NPOL==2)
        {
            for(int ig = wf.npw;ig < wf.npwx; ig++)
            {
                h_diag[ig] = 1.0;
                h_diag[ig+ wf.npwx] = 1.0;
            }
        }

        clock_t start=clock();

        //============================================================
        // diago the hamiltonian!!
        // In plane wave method, firstly using cinitcgg to diagnolize,
        // then using cg method.
        //
        // In localized orbital presented in plane wave case,
        // only using cinitcgg.
        //=============================================================
        double avg_iter_k = 0.0;
        hm.diago(istep, this->iter, ik, h_diag, avg_iter_k);

        avg_iter += avg_iter_k;

        en.print_band(ik); //mohan add 2012-04-16

        clock_t finish=clock();
        const double duration = static_cast<double>(finish - start) / CLOCKS_PER_SEC;

        ofs_running << " " << setw(8) 
        << ik+1 << setw(15) 
        << avg_iter_k << setw(15) << duration << endl;
    }//End K Loop
	
    //if (!LOCAL_BASIS) xiaohui modify 2013-09-02
    if(BASIS_TYPE=="pw") //xiaohui add 2013-09-02
    {
        //ofs_running << " avg_iteri " << avg_iter << endl;
        Parallel_Reduce::reduce_double_allpool(avg_iter); //mohan fix bug 2012-06-05
        //ofs_running << " avg_iter_after " << avg_iter << endl;
        avg_iter /= static_cast<double>(kv.nkstot);
    }
    delete [] h_diag;
    timer::tick("electrons","c_bands",'E');
    return;
} // END SUBROUTINE c_bands_k


void Electrons::init_mixstep_final_scf(void)
{
    TITLE("electrons","init_mixstep_final_scf");

    CHR.irstep=0;
    CHR.idstep=0;
    CHR.totstep=0;

    return;
}
