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
#include "../src_pw/structure_factor.h"
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
#include "module_hsolver/hsolver_pw.h"
#include "module_elecstate/elecstate_pw.h"
#include "module_hamilt/hamilt_pw.h"
#include "module_hsolver/diago_iter_assist.h"

#include "src_io/write_wfc_realspace.h"
#include "src_io/winput.h"
#include "src_io/numerical_descriptor.h"
#include "src_io/numerical_basis.h"
#include "src_io/to_wannier90.h"
#include "src_io/berryphase.h"

namespace ModuleESolver
{

    ESolver_KS_PW::ESolver_KS_PW()
    {
        classname = "ESolver_KS_PW";
        basisname = "PW";
    }
    ESolver_KS_PW::~ESolver_KS_PW()
    {
        //delete HSolver and ElecState
        if(this->phsol != nullptr)
        {
            delete (hsolver::HSolverPW*)this->phsol;
            this->phsol = nullptr;
        }
        if(this->pelec != nullptr)
        {
            delete (elecstate::ElecStatePW*)this->pelec;
            this->pelec = nullptr;
        }
        //delete Hamilt
        if(this->phami != nullptr)
        {
            delete (hamilt::HamiltPW*)this->phami;
            this->phami = nullptr;
        }
    }

    void ESolver_KS_PW::Init_GlobalC(Input& inp, UnitCell_pseudo& cell)
    {
        this->psi = GlobalC::wf.allocate(GlobalC::kv.nks);

        // cout<<GlobalC::rhopw->nrxx<<endl;
        // cout<<"before ufft allocate"<<endl;

        // cout<<"after ufft allocate"<<endl;

        //=======================
        // init pseudopotential
        //=======================
        GlobalC::ppcell.init(GlobalC::ucell.ntype);

        //=====================
        // init hamiltonian
        // only allocate in the beginning of ELEC LOOP!
        //=====================
        GlobalC::hm.hpw.allocate(GlobalC::wf.npwx, GlobalV::NPOL, GlobalC::ppcell.nkb, GlobalC::rhopw->nrxx);

        //=================================
        // initalize local pseudopotential
        //=================================
        GlobalC::ppcell.init_vloc(GlobalC::ppcell.vloc,GlobalC::rhopw);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LOCAL POTENTIAL");

        //======================================
        // Initalize non local pseudopotential
        //======================================
        GlobalC::ppcell.init_vnl(GlobalC::ucell);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "NON-LOCAL POTENTIAL");

        //=========================================================
        // calculate the total local pseudopotential in real space
        //=========================================================
        GlobalC::pot.init_pot(0, GlobalC::sf.strucFac); //atomic_rho, v_of_rho, set_vrs

        GlobalC::pot.newd();

        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT POTENTIAL");

        //==================================================
        // create GlobalC::ppcell.tab_at , for trial wave functions.
        //==================================================
        GlobalC::wf.init_at_1();

        //================================
        // Initial start wave functions
        //================================
        if (GlobalV::NBANDS != 0 || GlobalV::CALCULATION.substr(0,3) != "sto")
        // qianrui add temporarily. In the future, wfcinit() should be compatible with cases when NBANDS=0
        {
            GlobalC::wf.wfcinit(this->psi);
        }

#ifdef __LCAO
#ifdef __MPI
        switch (GlobalC::exx_global.info.hybrid_type) // Peize Lin add 2019-03-09
        {
        case Exx_Global::Hybrid_Type::HF:
        case Exx_Global::Hybrid_Type::PBE0:
        case Exx_Global::Hybrid_Type::HSE:
            GlobalC::exx_lip.init(&GlobalC::kv, &GlobalC::wf, GlobalC::wfcpw, GlobalC::rhopw, &GlobalC::ucell);
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

    void ESolver_KS_PW::Init(Input& inp, UnitCell_pseudo& ucell)
    {
        ESolver_KS::Init(inp,ucell);

        //temporary
        this->Init_GlobalC(inp,ucell);

        //init ElecState,
        if(this->pelec == nullptr)
        {
            this->pelec = new elecstate::ElecStatePW( GlobalC::wfcpw, (Charge*)(&(GlobalC::CHR)), (K_Vectors*)(&(GlobalC::kv)), GlobalV::NBANDS);
        }
        //init HSolver
        if(this->phsol == nullptr)
        {
            this->phsol = new hsolver::HSolverPW(GlobalC::wfcpw);
        }
    }

    void ESolver_KS_PW::beforescf(int istep)
    {
        //init Hamilt, this should be allocated before each scf loop
        //Operators in HamiltPW should be reallocated once cell changed
        //delete Hamilt if not first scf
        if(this->phami != nullptr)
        {
            delete (hamilt::HamiltPW*)this->phami;
            this->phami = nullptr;
        }
        //allocate HamiltPW
        if(this->phami == nullptr)
        {
            this->phami = new hamilt::HamiltPW();
        }

        //calculate ewald energy
        H_Ewald_pw::compute_ewald(GlobalC::ucell, GlobalC::rhopw);
        //Symmetry_rho should be moved to Init()
        Symmetry_rho srho;
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            srho.begin(is, GlobalC::CHR, GlobalC::rhopw, GlobalC::Pgrid, GlobalC::symm);
        }
    } 

    void ESolver_KS_PW::eachiterinit(const int istep, const int iter)
    {
        // mohan add 2010-07-16
        if (iter == 1) GlobalC::CHR.set_new_e_iteration(true);
        else GlobalC::CHR.set_new_e_iteration(false);

        if (GlobalV::FINAL_SCF && iter == 1)
        {
            GlobalC::CHR.irstep = 0;
            GlobalC::CHR.idstep = 0;
            GlobalC::CHR.totstep = 0;
        }

        // mohan move harris functional to here, 2012-06-05
        // use 'rho(in)' and 'v_h and v_xc'(in)
        GlobalC::en.calculate_harris(1);

        //(2) save change density as previous charge,
        // prepared fox mixing.
        if(GlobalV::MY_STOGROUP == 0)
	    {
            GlobalC::CHR.save_rho_before_sum_band();
        }
    }

    //Temporary, it should be replaced by hsolver later.
    void ESolver_KS_PW:: hamilt2density(const int istep, const int iter, const double ethr)
    {
        if(this->phsol != nullptr)
        {
            // reset energy 
            this->pelec->eband  = 0.0;
            this->pelec->demet  = 0.0;
            this->pelec->ef     = 0.0;
            GlobalC::en.ef_up  = 0.0;
            GlobalC::en.ef_dw  = 0.0;
            // choose if psi should be diag in subspace
            // be careful that istep start from 0 and iter start from 1
            if((istep==0||istep==1)&&iter==1) 
            {
                hsolver::DiagoIterAssist::need_subspace = false;
            }
            else 
            {
                hsolver::DiagoIterAssist::need_subspace = true;
            }

            hsolver::DiagoIterAssist::PW_DIAG_THR = ethr; 
            hsolver::DiagoIterAssist::PW_DIAG_NMAX = GlobalV::PW_DIAG_NMAX;
            this->phsol->solve(this->phami, this->psi[0], this->pelec, GlobalV::KS_SOLVER);

            // transform energy for print
            GlobalC::en.eband = this->pelec->eband;
            GlobalC::en.demet = this->pelec->demet;
            GlobalC::en.ef = this->pelec->ef;
        }
        else
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_PW", "HSolver has not been initialed!");
        }

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
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            srho.begin(is, GlobalC::CHR, GlobalC::rhopw, GlobalC::Pgrid, GlobalC::symm);
        }

        // compute magnetization, only for LSDA(spin==2)
        GlobalC::ucell.magnet.compute_magnetization();
        // deband is calculated from "output" charge density calculated
        // in sum_band
        // need 'rho(out)' and 'vr (v_h(in) and v_xc(in))'

        GlobalC::en.deband = GlobalC::en.delta_e();
        //if (LOCAL_BASIS) xiaohui modify 2013-09-02
    }

    //Temporary, it should be rewritten with Hamilt class. 
    void ESolver_KS_PW::updatepot(const int istep, const int iter)
    {
        if (!this->conv_elec)
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
            for (int is = 0; is < GlobalV::NSPIN; ++is)
            {
                for (int ir = 0; ir < GlobalC::rhopw->nrxx; ++ir)
                {
                    GlobalC::pot.vnew(is, ir) = GlobalC::pot.vr(is, ir);
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

    void ESolver_KS_PW::eachiterfinish(const int iter)
    {
        //print_eigenvalue(GlobalV::ofs_running);
        GlobalC::en.calculate_etot();
        //We output it for restarting the scf.
        bool print = false;
        if (this->out_freq_elec == 0)
        {
            if (this->conv_elec) print = true;
        }
        else
        {
            if (iter % this->out_freq_elec == 0 || this->conv_elec) print = true;
        }

        if (print)
        {
            if (GlobalC::CHR.out_chg > 0)
            {
                for (int is = 0; is < GlobalV::NSPIN; is++)
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
            if (GlobalC::wf.out_wfc_pw == 1 || GlobalC::wf.out_wfc_pw == 2)
            {
                std::stringstream ssw;
                ssw << GlobalV::global_out_dir << "WAVEFUNC";
                // mohan update 2011-02-21
                //qianrui update 2020-10-17
                WF_io::write_wfc(ssw.str(), this->psi[0], &GlobalC::kv, GlobalC::wfcpw);
                //ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"write wave functions into file WAVEFUNC.dat");
            }

        }


    }


    void ESolver_KS_PW::afterscf()
    {
        for(int ik=0; ik<this->pelec->ekb.nr; ++ik)
        {
            for(int ib=0; ib<this->pelec->ekb.nc; ++ib)
            {
                GlobalC::wf.ekb[ik][ib] = this->pelec->ekb(ik, ib);
                GlobalC::wf.wg(ik, ib) = this->pelec->wg(ik, ib);
            }
        }
#ifdef __LCAO
        if (GlobalC::chi0_hilbert.epsilon)                 // pengfei 2016-11-23
        {
            std::cout << "eta = " << GlobalC::chi0_hilbert.eta << std::endl;
            std::cout << "domega = " << GlobalC::chi0_hilbert.domega << std::endl;
            std::cout << "nomega = " << GlobalC::chi0_hilbert.nomega << std::endl;
            std::cout << "dim = " << GlobalC::chi0_hilbert.dim << std::endl;
            //std::cout <<"oband = "<<GlobalC::chi0_hilbert.oband<<std::endl;
            GlobalC::chi0_hilbert.Chi();
        }
#endif

        if (GlobalC::chi0_standard.epsilon)
        {
            std::cout << "eta = " << GlobalC::chi0_standard.eta << std::endl;
            std::cout << "domega = " << GlobalC::chi0_standard.domega << std::endl;
            std::cout << "nomega = " << GlobalC::chi0_standard.nomega << std::endl;
            std::cout << "dim = " << GlobalC::chi0_standard.dim << std::endl;
            //std::cout <<"oband = "<<GlobalC::chi0_standard.oband<<std::endl;
            GlobalC::chi0_standard.Chi();
        }
        if (GlobalC::epsilon0_pwscf.epsilon)
        {
            GlobalC::epsilon0_pwscf.Cal_epsilon0();
        }
        if (GlobalC::epsilon0_vasp.epsilon)
        {
            GlobalC::epsilon0_vasp.cal_epsilon0();
        }

        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            std::stringstream ssc;
            std::stringstream ss1;
            ssc << GlobalV::global_out_dir << "SPIN" << is + 1 << "_CHG";
            ss1 << GlobalV::global_out_dir << "SPIN" << is + 1 << "_CHG.cube";
            GlobalC::CHR.write_rho(GlobalC::CHR.rho_save[is], is, 0, ssc.str());//mohan add 2007-10-17
            GlobalC::CHR.write_rho_cube(GlobalC::CHR.rho_save[is], is, ss1.str(), 3);
        }
        if (this->conv_elec)
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

        if (GlobalV::OUT_LEVEL != "m")
        {
            this->print_eigenvalue(GlobalV::ofs_running);
        }
    }

    void ESolver_KS_PW::print_eigenvalue(std::ofstream& ofs)
    {
        bool wrong = false;
        for (int ik = 0; ik < GlobalC::kv.nks; ++ik)
        {
            for (int ib = 0; ib < GlobalV::NBANDS; ++ib)
            {
                if (abs(GlobalC::wf.ekb[ik][ib]) > 1.0e10)
                {
                    GlobalV::ofs_warning << " ik=" << ik + 1 << " ib=" << ib + 1 << " " << GlobalC::wf.ekb[ik][ib] << " Ry" << std::endl;
                    wrong = true;
                }
            }
        }
        if (wrong)
        {
            ModuleBase::WARNING_QUIT("Threshold_Elec::print_eigenvalue", "Eigenvalues are too large!");
        }


        if (GlobalV::MY_RANK != 0)
        {
            return;
        }

        ModuleBase::TITLE("Threshold_Elec", "print_eigenvalue");

        ofs << "\n STATE ENERGY(eV) AND OCCUPATIONS ";
        ofs << std::setprecision(5);
        for (int ik = 0;ik < GlobalC::kv.nks;ik++)
        {
            if (ik == 0)
            {
                ofs << "   NSPIN == " << GlobalV::NSPIN << std::endl;
                if (GlobalV::NSPIN == 2)
                {
                    ofs << "SPIN UP : " << std::endl;
                }
            }
            else if (ik == GlobalC::kv.nks / 2)
            {
                if (GlobalV::NSPIN == 2)
                {
                    ofs << "SPIN DOWN : " << std::endl;
                }
            }

            if (GlobalV::NSPIN == 2)
            {
                if (GlobalC::kv.isk[ik] == 0)
                {
                    ofs << " " << ik + 1 << "/" << GlobalC::kv.nks / 2 << " kpoint (Cartesian) = "
                        << GlobalC::kv.kvec_c[ik].x << " " << GlobalC::kv.kvec_c[ik].y << " " << GlobalC::kv.kvec_c[ik].z
                        << " (" << GlobalC::kv.ngk[ik] << " pws)" << std::endl;

                    ofs << std::setprecision(6);

                }
                if (GlobalC::kv.isk[ik] == 1)
                {
                    ofs << " " << ik + 1 - GlobalC::kv.nks / 2 << "/" << GlobalC::kv.nks / 2 << " kpoint (Cartesian) = "
                        << GlobalC::kv.kvec_c[ik].x << " " << GlobalC::kv.kvec_c[ik].y << " " << GlobalC::kv.kvec_c[ik].z
                        << " (" << GlobalC::kv.ngk[ik] << " pws)" << std::endl;

                    ofs << std::setprecision(6);

                }
            }       // Pengfei Li  added  14-9-9
            else
            {
                ofs << " " << ik + 1 << "/" << GlobalC::kv.nks << " kpoint (Cartesian) = "
                    << GlobalC::kv.kvec_c[ik].x << " " << GlobalC::kv.kvec_c[ik].y << " " << GlobalC::kv.kvec_c[ik].z
                    << " (" << GlobalC::kv.ngk[ik] << " pws)" << std::endl;

                ofs << std::setprecision(6);
            }

            //----------------------
            // no energy to output
            //----------------------
            if (GlobalV::KS_SOLVER == "selinv")
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
                    ofs << std::setw(8) << ib + 1
                        << std::setw(15) << GlobalC::wf.ekb[ik][ib] * ModuleBase::Ry_to_eV
                        << std::setw(15) << GlobalC::wf.wg(ik, ib) << std::endl;
                }
                ofs << std::endl;
            }
        }//end ik
        return;
    }



    void ESolver_KS_PW::cal_Energy(energy& en)
    {

    }

    void ESolver_KS_PW::cal_Force(ModuleBase::matrix& force)
    {
        Forces ff;
        ff.init(force, this->psi);
    }

    void ESolver_KS_PW::cal_Stress(ModuleBase::matrix& stress)
    {
        Stress_PW ss;
        ss.cal_stress(stress, this->psi);
    }

    void ESolver_KS_PW::postprocess()
    {
        // compute density of states
        GlobalC::en.perform_dos_pw();

        // caoyu add 2020-11-24, mohan updat 2021-01-03
        if(GlobalV::BASIS_TYPE=="pw" && GlobalV::deepks_out_labels)
        {
            Numerical_Descriptor nc;
            nc.output_descriptor(this->psi[0], INPUT.deepks_descriptor_lmax);
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"GENERATE DESCRIPTOR FOR DEEPKS");
        }

        if(GlobalV::BASIS_TYPE=="pw" && winput::out_spillage) //xiaohui add 2013-09-01
        {
            //std::cout << "\n Output Spillage Information : " << std::endl;
            // calculate spillage value.
#ifdef __LCAO
            if ( winput::out_spillage == 3)
            {
                GlobalV::BASIS_TYPE="pw"; 
                std::cout << " NLOCAL = " << GlobalV::NLOCAL << std::endl;

                for (int ik=0; ik<GlobalC::kv.nks; ik++)
                {
                    GlobalC::wf.wanf2[ik].create(GlobalV::NLOCAL, GlobalC::wf.npwx);
                    if(GlobalV::BASIS_TYPE=="pw")
                    {
                        std::cout << " ik=" << ik + 1 << std::endl;

                        GlobalV::BASIS_TYPE="lcao_in_pw";
                        GlobalC::wf.LCAO_in_pw_k(ik, GlobalC::wf.wanf2[ik]);
                        GlobalV::BASIS_TYPE="pw";
                    }
                }

                //Spillage sp;
                //sp.get_both(GlobalV::NBANDS, GlobalV::NLOCAL, GlobalC::wf.wanf2, GlobalC::wf.evc);
            }
#endif

            // output overlap
            if ( winput::out_spillage <= 2 )
            {
                Numerical_Basis numerical_basis;
                numerical_basis.output_overlap(this->psi[0]);
                ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"BASIS OVERLAP (Q and S) GENERATION.");
            }
        }

        if(GlobalC::wf.out_wfc_r == 1)				// Peize Lin add 2021.11.21
        {
            Write_Wfc_Realspace::write_wfc_realspace_1(this->psi[0], "wfc_realspace", true);
        }	
    }

    void ESolver_KS_PW::hamilt2estates(const double ethr)
    {
        if(this->phsol != nullptr)
        {
            hsolver::DiagoIterAssist::need_subspace = false;
            hsolver::DiagoIterAssist::PW_DIAG_THR = ethr; 
            this->phsol->solve(this->phami, this->psi[0], this->pelec, GlobalV::KS_SOLVER, true);
        }
        else
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_PW", "HSolver has not been initialed!");
        }
    }

    void ESolver_KS_PW::nscf()
    {
        ModuleBase::TITLE("ESolver_KS_PW","nscf");
        ModuleBase::timer::tick("ESolver_KS_PW","nscf");

        this->beforescf(1);
        //========================================
        // diagonalization of the KS hamiltonian
        // =======================================
        double diag_ethr = this->phsol->set_diagethr(1, 1, drho);

        this->hamilt2estates(diag_ethr);

        for(int ik=0; ik<this->pelec->ekb.nr; ++ik)
        {
            for(int ib=0; ib<this->pelec->ekb.nc; ++ib)
            {
                GlobalC::wf.ekb[ik][ib] = this->pelec->ekb(ik, ib);
            }
        }

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
                << " " << this->pelec->ekb(ik, ib) * ModuleBase::Ry_to_eV
                << " " << GlobalC::wf.wg(ik, ib)*GlobalC::kv.nks << std::endl;
            }
            GlobalV::ofs_running << std::endl;
        }

        // add by jingan in 2018.11.7
        if(INPUT.towannier90)
        {
            toWannier90 myWannier(GlobalC::kv.nkstot,GlobalC::ucell.G);
            myWannier.init_wannier(this->psi);
        }

        //=======================================================
        // Do a Berry phase polarization calculation if required
        //=======================================================

        if (berryphase::berry_phase_flag && ModuleSymmetry::Symmetry::symm_flag == 0)
        {
            berryphase bp;
            bp.Macroscopic_polarization(this->psi);
        }

        ModuleBase::timer::tick("ESolver_KS_PW","nscf");
        return;
    }

}
