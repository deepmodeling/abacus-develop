#include "esolver_ks_pw.h"
#include <iostream>
#include "../src_io/wf_io.h"

//--------------temporary----------------------------
#include "../src_pw/global.h"
#include "../src_pw/symmetry_rho.h"
#include "../src_io/print_info.h"
#include "../src_pw/H_Ewald_pw.h"
#include "../src_pw/occupy.h"
#include "../module_relax/relax_old/variable_cell.h"    // liuyu 2022-11-07
//-----force-------------------
#include "../src_pw/forces.h"
//-----stress------------------
#include "../src_pw/stress_pw.h"
//---------------------------------------------------
#include "module_hsolver/hsolver_pw.h"
#include "module_elecstate/elecstate_pw.h"
#include "module_hamilt/hamilt_pw.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_vdw/vdw.h"

#include "src_io/write_wfc_realspace.h"
#include "src_io/winput.h"
#include "src_io/numerical_descriptor.h"
#include "src_io/numerical_basis.h"
#include "src_io/to_wannier90.h"
#include "src_io/berryphase.h"
#include "module_psi/include/device.h"
#include "module_hsolver/include/math_kernel.h"

namespace ModuleESolver
{

    template<typename FPTYPE, typename Device>
    ESolver_KS_PW<FPTYPE, Device>::ESolver_KS_PW()
    {
        this->classname = "ESolver_KS_PW";
        this->basisname = "PW";
        this->device = psi::device::get_device_type<Device>(this->ctx);
    #if ((defined __CUDA) || (defined __ROCM))
        if (this->device == psi::GpuDevice) {
            hsolver::createBLAShandle();
        }
    #endif
    }

    template<typename FPTYPE, typename Device>
    ESolver_KS_PW<FPTYPE, Device>::~ESolver_KS_PW()
    {
        //delete HSolver and ElecState
        if(this->phsol != nullptr)
        {
            delete reinterpret_cast<hsolver::HSolverPW<FPTYPE, Device>*>(this->phsol);
            this->phsol = nullptr;
        }
        if(this->pelec != nullptr)
        {
            delete reinterpret_cast<elecstate::ElecStatePW<FPTYPE, Device>*>(this->pelec);
            this->pelec = nullptr;
        }
        //delete Hamilt
        if(this->p_hamilt != nullptr)
        {
            delete reinterpret_cast<hamilt::HamiltPW<FPTYPE, Device>*>(this->p_hamilt);
            this->p_hamilt = nullptr;
        }
    #if ((defined __CUDA) || (defined __ROCM))
        if (this->device == psi::GpuDevice) {
            delete reinterpret_cast<psi::Psi<std::complex<FPTYPE>, Device>*>(this->kspw_psi);
            hsolver::destoryBLAShandle();
        }
    #endif
    }

    template<typename FPTYPE, typename Device>
    void ESolver_KS_PW<FPTYPE, Device>::Init_GlobalC(Input& inp, UnitCell& cell)
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
        //not used anymore
        //GlobalC::hm.hpw.allocate(GlobalC::wf.npwx, GlobalV::NPOL, GlobalC::ppcell.nkb, GlobalC::rhopw->nrxx);

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

        GlobalC::ppcell.cal_effective_D();


        //==================================================
        // create GlobalC::ppcell.tab_at , for trial wave functions.
        //==================================================
        GlobalC::wf.init_at_1();

        //================================
        // Initial start wave functions
        //================================
        if (GlobalV::NBANDS != 0 || GlobalV::ESOLVER_TYPE != "sdft")
        // qianrui add temporarily. In the future, wfcinit() should be compatible with cases when NBANDS=0
        {
            GlobalC::wf.wfcinit(this->psi);
        }

        // denghui added 20221116
        this->kspw_psi = this->device == psi::GpuDevice ?
                         new psi::Psi<std::complex<FPTYPE>, Device>(this->psi[0]) :
                         reinterpret_cast<psi::Psi<std::complex<FPTYPE>, Device>*> (this->psi);

        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT BASIS");
    }


    template<typename FPTYPE, typename Device>
    void ESolver_KS_PW<FPTYPE, Device>::Init(Input& inp, UnitCell& ucell)
    {
        ESolver_KS<FPTYPE, Device>::Init(inp,ucell);

        //init HSolver
        if(this->phsol == nullptr)
        {
            this->phsol = new hsolver::HSolverPW<FPTYPE, Device>(GlobalC::wfcpw);
        }

        //init ElecState,
        if(this->pelec == nullptr)
        {
            this->pelec = new elecstate::ElecStatePW<FPTYPE, Device>( GlobalC::wfcpw, &(this->chr), (K_Vectors*)(&(GlobalC::kv)));
        }

        // Inititlize the charge density.
        this->pelec->charge->allocate(GlobalV::NSPIN, GlobalC::rhopw->nrxx, GlobalC::rhopw->npw);

        // Initialize the potential.
        if(this->pelec->pot == nullptr)
        {
            this->pelec->pot = new elecstate::Potential(
                GlobalC::rhopw,
                &GlobalC::ucell,
                &(GlobalC::ppcell.vloc),
                &(GlobalC::sf.strucFac),
                &(GlobalC::en.etxc),
                &(GlobalC::en.vtxc)
            );
        }
        
        //temporary
        this->Init_GlobalC(inp,ucell);

        //Fix pelec->wg by ocp_kb
        if(GlobalV::ocp)
        {
            this->pelec->fixed_weights(GlobalV::ocp_kb.data());
        }
    }

    template<typename FPTYPE, typename Device>
    void ESolver_KS_PW<FPTYPE, Device>::beforescf(int istep)
    {
        ModuleBase::TITLE("ESolver_KS_PW", "beforescf");

        // Temporary, md and relax will merge later   liuyu add 2022-11-07
        if(GlobalV::CALCULATION == "md" && istep)
        {
            this->CE.update_istep();
            this->CE.save_pos_next(GlobalC::ucell);
            this->CE.extrapolate_charge(this->pelec->charge);

            if(GlobalC::ucell.cell_parameter_updated)
            {
                Variable_Cell::init_after_vc();
            }

            //this->pelec->init_scf(istep, GlobalC::sf.strucFac);
        }

        if(GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax")
        {
            if(GlobalC::ucell.ionic_position_updated)
            {
                GlobalV::ofs_running << " Setup the extrapolated charge." << std::endl;
                // charge extrapolation if istep>0.
                this->CE.update_istep();
                this->CE.update_all_pos(GlobalC::ucell);
                this->CE.extrapolate_charge(this->pelec->charge);
                this->CE.save_pos_next(GlobalC::ucell);

                GlobalV::ofs_running << " Setup the Vl+Vh+Vxc according to new structure factor and new charge." << std::endl;
                // calculate the new potential accordint to
                // the new charge density.
                //this->pelec->init_scf( istep, GlobalC::sf.strucFac );
            }
        }
        if(GlobalC::ucell.cell_parameter_updated)
        {
            GlobalC::wfcpw->initgrids(GlobalC::ucell.lat0, GlobalC::ucell.latvec, GlobalC::wfcpw->nx, GlobalC::wfcpw->ny, GlobalC::wfcpw->nz);
            GlobalC::wfcpw->initparameters(false, INPUT.ecutwfc, GlobalC::kv.nks, GlobalC::kv.kvec_d.data());
            GlobalC::wfcpw->collect_local_pw(); 
            GlobalC::wf.init_after_vc(GlobalC::kv.nks);
            GlobalC::wf.init_at_1();
        }
        //init Hamilt, this should be allocated before each scf loop
        //Operators in HamiltPW should be reallocated once cell changed
        //delete Hamilt if not first scf
        if(this->p_hamilt != nullptr)
        {
            delete reinterpret_cast<hamilt::HamiltPW<FPTYPE, Device>*>(this->p_hamilt);
            this->p_hamilt = nullptr;
        }
        //allocate HamiltPW
        if(this->p_hamilt == nullptr)
        {
            this->p_hamilt = new hamilt::HamiltPW<FPTYPE, Device>(this->pelec->pot);
        }

        //----------------------------------------------------------
        // about vdw, jiyy add vdwd3 and linpz add vdwd2
        //----------------------------------------------------------
        auto vdw_solver = vdw::make_vdw(GlobalC::ucell, INPUT);
        if (vdw_solver != nullptr)
        {
            GlobalC::en.evdw = vdw_solver->get_energy();
        }

        //calculate ewald energy
        if(!GlobalV::test_skip_ewald)
        {
            H_Ewald_pw::compute_ewald(GlobalC::ucell, GlobalC::rhopw);
        }

        //=========================================================
        // calculate the total local pseudopotential in real space
        //=========================================================
        this->pelec->init_scf(istep, GlobalC::sf.strucFac);
        //Symmetry_rho should behind init_scf, because charge should be initialized first.
        Symmetry_rho srho;
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            srho.begin(is, *(this->pelec->charge), GlobalC::rhopw, GlobalC::Pgrid, GlobalC::symm);
        }

    } 

    template<typename FPTYPE, typename Device>
    void ESolver_KS_PW<FPTYPE, Device>::othercalculation(const int istep)
    {
        ModuleBase::TITLE("ESolver_KS_PW", "othercalculation");
        ModuleBase::timer::tick("ESolver_KS_PW", "othercalculation");
        if(GlobalV::CALCULATION == "test_memory")
        {
            Cal_Test::test_memory();
            return;
        }

        if (GlobalV::CALCULATION == "gen_bessel")
        {
            // caoyu add 2020-11-24, mohan updat 2021-01-03
            Numerical_Descriptor nc;
            nc.output_descriptor(this->psi[0], INPUT.bessel_lmax, INPUT.bessel_rcut, INPUT.bessel_tol);
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"GENERATE DESCRIPTOR FOR DEEPKS");
            return;
        }

        // self consistent calculations for electronic ground state
        if (GlobalV::CALCULATION == "nscf")
        {
            this->nscf();
        }        
        else
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::othercalculation", "CALCULATION type not supported");
        }

        ModuleBase::timer::tick("ESolver_KS_PW", "othercalculation");
        return;
    }

    template<typename FPTYPE, typename Device>
    void ESolver_KS_PW<FPTYPE, Device>::eachiterinit(const int istep, const int iter)
    {
        // mohan add 2010-07-16
        if (iter == 1) GlobalC::CHR_MIX.reset();

        // mohan move harris functional to here, 2012-06-05
        // use 'rho(in)' and 'v_h and v_xc'(in)
        GlobalC::en.deband_harris = GlobalC::en.delta_e(this->pelec);

        //(2) save change density as previous charge,
        // prepared fox mixing.
        if(GlobalV::MY_STOGROUP == 0)
	    {
            this->pelec->charge->save_rho_before_sum_band();
        }
    }

    //Temporary, it should be replaced by hsolver later.
    template<typename FPTYPE, typename Device>
    void ESolver_KS_PW<FPTYPE, Device>::hamilt2density(const int istep, const int iter, const FPTYPE ethr)
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
                hsolver::DiagoIterAssist<FPTYPE>::need_subspace = false;
            }
            else 
            {
                hsolver::DiagoIterAssist<FPTYPE>::need_subspace = true;
            }

            hsolver::DiagoIterAssist<FPTYPE>::PW_DIAG_THR = ethr; 
            hsolver::DiagoIterAssist<FPTYPE>::PW_DIAG_NMAX = GlobalV::PW_DIAG_NMAX;
            this->phsol->solve(this->p_hamilt, this->kspw_psi[0], this->pelec, GlobalV::KS_SOLVER);
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
#ifdef __EXX
        GlobalC::en.set_exx();		// Peize Lin add 2019-03-09
#endif
#endif
    // calculate the delta_harris energy
    // according to new charge density.
    // mohan add 2009-01-23
        GlobalC::en.calculate_harris();
        Symmetry_rho srho;
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            srho.begin(is, *(this->pelec->charge), GlobalC::rhopw, GlobalC::Pgrid, GlobalC::symm);
        }

        // compute magnetization, only for LSDA(spin==2)
        GlobalC::ucell.magnet.compute_magnetization(this->pelec->charge, this->pelec->nelec_spin.data());
        // deband is calculated from "output" charge density calculated
        // in sum_band
        // need 'rho(out)' and 'vr (v_h(in) and v_xc(in))'

        GlobalC::en.deband = GlobalC::en.delta_e(this->pelec);
        //if (LOCAL_BASIS) xiaohui modify 2013-09-02
    }

    //Temporary, it should be rewritten with Hamilt class. 
    template<typename FPTYPE, typename Device>
    void ESolver_KS_PW<FPTYPE, Device>::updatepot(const int istep, const int iter)
    {
        if (!this->conv_elec)
        {
            this->pelec->pot->update_from_charge(this->pelec->charge, &GlobalC::ucell);
            GlobalC::en.delta_escf(this->pelec);
        }
        else
        {
            GlobalC::en.cal_converged(this->pelec);
        }
    }

    template<typename FPTYPE, typename Device>
    void ESolver_KS_PW<FPTYPE, Device>::eachiterfinish(const int iter)
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
            if (GlobalV::out_chg > 0)
            {
                for (int is = 0; is < GlobalV::NSPIN; is++)
                {
                    std::stringstream ssc;
                    std::stringstream ss1;
                    ssc << GlobalV::global_out_dir << "tmp" << "_SPIN" << is + 1 << "_CHG";
                    this->pelec->charge->write_rho(this->pelec->charge->rho_save[is], is, iter, ssc.str(), 3);//mohan add 2007-10-17
                    ss1 << GlobalV::global_out_dir << "tmp" << "_SPIN" << is + 1 << "_CHG.cube";
                    this->pelec->charge->write_rho_cube(this->pelec->charge->rho_save[is], is, ss1.str(), 3);
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

    template<typename FPTYPE, typename Device>
    void ESolver_KS_PW<FPTYPE, Device>::afterscf(const int istep)
    {
        // Temporary liuyu add 2022-11-07
        this->CE.update_all_pos(GlobalC::ucell);

        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            std::stringstream ssc;
            std::stringstream ss1;
            ssc << GlobalV::global_out_dir << "SPIN" << is + 1 << "_CHG";
            ss1 << GlobalV::global_out_dir << "SPIN" << is + 1 << "_CHG.cube";
            this->pelec->charge->write_rho(this->pelec->charge->rho_save[is], is, 0, ssc.str());//mohan add 2007-10-17
            this->pelec->charge->write_rho_cube(this->pelec->charge->rho_save[is], is, ss1.str(), 3);
        }
        if (this->conv_elec)
        {
            GlobalV::ofs_running << "\n charge density convergence is achieved" << std::endl;
            GlobalV::ofs_running << " final etot is " << GlobalC::en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
        }
        else
        {
            GlobalV::ofs_running << " convergence has NOT been achieved!" << std::endl;
        }

		if(GlobalV::out_pot == 2)
		{
			std::stringstream ssp;
			std::stringstream ssp_ave;
			ssp << GlobalV::global_out_dir << "ElecStaticPot";
			ssp_ave << GlobalV::global_out_dir << "ElecStaticPot_AVE";
			this->pelec->pot->write_elecstat_pot(ssp.str(), ssp_ave.str(), GlobalC::rhopw, this->pelec->charge); //output 'Hartree + local pseudopot'
		}

        if (GlobalV::OUT_LEVEL != "m")
        {
            this->print_eigenvalue(GlobalV::ofs_running);
        }
        if (this->device == psi::GpuDevice) {
            syncmem_complex_d2h_op()(
                this->psi[0].get_device(),
                this->kspw_psi[0].get_device(),
                this->psi[0].get_pointer() - this->psi[0].get_psi_bias(),
                this->kspw_psi[0].get_pointer() - this->kspw_psi[0].get_psi_bias(),
                this->psi[0].size());
        }
    }

    template<typename FPTYPE, typename Device>
    void ESolver_KS_PW<FPTYPE, Device>::print_eigenvalue(std::ofstream& ofs)
    {
        bool wrong = false;
        for (int ik = 0; ik < GlobalC::kv.nks; ++ik)
        {
            for (int ib = 0; ib < GlobalV::NBANDS; ++ib)
            {
                if (abs(this->pelec->ekb(ik, ib)) > 1.0e10)
                {
                    GlobalV::ofs_warning << " ik=" << ik + 1 << " ib=" << ib + 1 << " " << this->pelec->ekb(ik, ib) << " Ry" << std::endl;
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

            GlobalV::ofs_running << std::setprecision(6);
            GlobalV::ofs_running << std::setiosflags(ios::showpoint);
            for (int ib = 0; ib < GlobalV::NBANDS; ib++)
            {
                ofs << std::setw(8) << ib + 1
                    << std::setw(15) << this->pelec->ekb(ik, ib) * ModuleBase::Ry_to_eV
                    << std::setw(15) << this->pelec->wg(ik, ib) << std::endl;
            }
            ofs << std::endl;
        }//end ik
        return;
    }



    template<typename FPTYPE, typename Device>
    void ESolver_KS_PW<FPTYPE, Device>::cal_Energy(FPTYPE& etot)
    {
        etot = GlobalC::en.etot;
    }

    template<typename FPTYPE, typename Device>
    void ESolver_KS_PW<FPTYPE, Device>::cal_Force(ModuleBase::matrix& force)
    {
        Forces<FPTYPE, Device> ff;
        ff.init(force, this->pelec->wg, this->pelec->charge, this->kspw_psi);
    }

    template<typename FPTYPE, typename Device>
    void ESolver_KS_PW<FPTYPE, Device>::cal_Stress(ModuleBase::matrix& stress)
    {
        Stress_PW<FPTYPE, Device> ss(this->pelec);
        ss.cal_stress(stress, this->psi, this->kspw_psi);

        //external stress
        FPTYPE unit_transform = 0.0;
        unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI,3) * 1.0e-8;
        FPTYPE external_stress[3] = {GlobalV::PRESS1,GlobalV::PRESS2,GlobalV::PRESS3};
        for(int i=0;i<3;i++)
        {
            stress(i,i) -= external_stress[i]/unit_transform;
        }
        GlobalV::PRESSURE = (stress(0,0)+stress(1,1)+stress(2,2))/3;
    }

    template<typename FPTYPE, typename Device>
    void ESolver_KS_PW<FPTYPE, Device>::postprocess()
    {

        GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
        GlobalV::ofs_running << std::setprecision(16);
        GlobalV::ofs_running << " !FINAL_ETOT_IS " << GlobalC::en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
        GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;
        
        //print occupation in istate.info

	    GlobalC::en.print_occ(this->pelec);
        // compute density of states
        GlobalC::en.perform_dos_pw(this->pelec);

        if(GlobalV::BASIS_TYPE=="pw" && winput::out_spillage) //xiaohui add 2013-09-01
        {
            //std::cout << "\n Output Spillage Information : " << std::endl;
            // calculate spillage value.
#ifdef __LCAO
//We are not goint to support lcao_in_paw until
//the obsolete GlobalC::hm is replaced by the 
//refactored moeules (psi, hamilt, etc.)
/*
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
*/
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

        if(INPUT.cal_cond)
	    {
            this->KG(INPUT.cond_nche,INPUT.cond_fwhm,INPUT.cond_wcut,INPUT.cond_dw,INPUT.cond_wenlarge, this->pelec->wg);
        }
    }

    template<typename FPTYPE, typename Device>
    void ESolver_KS_PW<FPTYPE, Device>::hamilt2estates(const FPTYPE ethr)
    {
        if(this->phsol != nullptr)
        {
            hsolver::DiagoIterAssist<FPTYPE, Device>::need_subspace = false;
            hsolver::DiagoIterAssist<FPTYPE, Device>::PW_DIAG_THR = ethr; 
            this->phsol->solve(this->p_hamilt, this->kspw_psi[0], this->pelec, GlobalV::KS_SOLVER, true);
        }
        else
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_PW", "HSolver has not been initialed!");
        }
    }

    template<typename FPTYPE, typename Device>
    void ESolver_KS_PW<FPTYPE, Device>::nscf()
    {
        ModuleBase::TITLE("ESolver_KS_PW","nscf");
        ModuleBase::timer::tick("ESolver_KS_PW","nscf");

        this->beforescf(0);
        //========================================
        // diagonalization of the KS hamiltonian
        // =======================================
        FPTYPE diag_ethr = GlobalV::PW_DIAG_THR;
        if(diag_ethr - 1e-2 > -1e-5)   
            diag_ethr = std::max(1e-13, 0.1*std::min(1e-2,GlobalV::SCF_THR / GlobalV::nelec));
        GlobalV::ofs_running << " PW_DIAG_THR  = "<< diag_ethr << std::endl;

        this->hamilt2estates(diag_ethr);
        this->pelec->calculate_weights();

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
                << " " << this->pelec->wg(ik, ib)*GlobalC::kv.nks << std::endl;
            }
            GlobalV::ofs_running << std::endl;
        }

        // add by jingan in 2018.11.7
        if(INPUT.towannier90)
        {
            toWannier90 myWannier(GlobalC::kv.nkstot,GlobalC::ucell.G);
            myWannier.init_wannier(this->pelec->ekb, this->psi);
        }

        //=======================================================
        // Do a Berry phase polarization calculation if required
        //=======================================================

        if (berryphase::berry_phase_flag && ModuleSymmetry::Symmetry::symm_flag != 1)
        {
            berryphase bp;
            bp.Macroscopic_polarization(this->psi);
        }

        ModuleBase::timer::tick("ESolver_KS_PW","nscf");
        return;
    }

template class ESolver_KS_PW<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class ESolver_KS_PW<double, psi::DEVICE_GPU>;
#endif
}
