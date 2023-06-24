#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_esolver/esolver_ks_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
//
#include "module_base/timer.h"
#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_io/berryphase.h"
#include "module_io/istate_charge.h"
#include "module_io/istate_envelope.h"
#include "module_io/to_wannier90.h"
#include "module_io/write_HS_R.h"
#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#endif
#include "module_hamilt_general/module_ewald/H_Ewald_pw.h"
#include "module_hamilt_general/module_vdw/vdw.h"

#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/op_exx_lcao.h"
#include "module_io/dm_io.h"

namespace ModuleESolver
{

    void ESolver_KS_LCAO::set_matrix_grid(Record_adj& ra)
    {
        ModuleBase::TITLE("ESolver_KS_LCAO", "set_matrix_grid");
        ModuleBase::timer::tick("ESolver_KS_LCAO", "set_matrix_grid");

        // (1) Find adjacent atoms for each atom.
        GlobalV::SEARCH_RADIUS = atom_arrange::set_sr_NL(
            GlobalV::ofs_running,
            GlobalV::OUT_LEVEL,
            GlobalC::ORB.get_rcutmax_Phi(),
            GlobalC::ucell.infoNL.get_rcutmax_Beta(),
            GlobalV::GAMMA_ONLY_LOCAL);

        atom_arrange::search(
            GlobalV::SEARCH_PBC,
            GlobalV::ofs_running,
            GlobalC::GridD,
            GlobalC::ucell,
            GlobalV::SEARCH_RADIUS,
            GlobalV::test_atom_input);

        //ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"SEARCH ADJACENT ATOMS");

        // (3) Periodic condition search for each grid.
        this->GridT.set_pbc_grid(pw_rho->nx,
                                 pw_rho->ny,
                                 pw_rho->nz,
                                 this->pw_big->bx,
                                 this->pw_big->by,
                                 this->pw_big->bz,
                                 this->pw_big->nbx,
                                 this->pw_big->nby,
                                 this->pw_big->nbz,
                                 this->pw_big->nbxx,
                                 this->pw_big->nbzp_start,
                                 this->pw_big->nbzp,
                                 pw_rho->ny,
                                 pw_rho->nplane,
                                 pw_rho->startz_current);

        // (2)For each atom, calculate the adjacent atoms in different cells
        // and allocate the space for H(R) and S(R).
        // If k point is used here, allocate HlocR after atom_arrange.
        Parallel_Orbitals* pv = this->UHM.LM->ParaV;
        ra.for_2d(*pv, GlobalV::GAMMA_ONLY_LOCAL);
        if (!GlobalV::GAMMA_ONLY_LOCAL)
        {
            this->UHM.LM->allocate_HS_R(pv->nnr);
#ifdef __DEEPKS
            GlobalC::ld.allocate_V_deltaR(pv->nnr);
#endif

            // need to first calculae lgd.
            // using GridT.init.
            this->GridT.cal_nnrg(pv);
        }

        ModuleBase::timer::tick("ESolver_KS_LCAO", "set_matrix_grid");
        return;
    }


    void ESolver_KS_LCAO::beforesolver(const int istep)
    {
        ModuleBase::TITLE("ESolver_KS_LCAO", "beforesolver");
        ModuleBase::timer::tick("ESolver_KS_LCAO", "beforesolver");

        // 1. prepare HS matrices, prepare grid integral
        this->set_matrix_grid(this->RA);

        // 2. density matrix extrapolation 

        // set the augmented orbitals index.
        // after ParaO and GridT, 
        // this information is used to calculate
        // the force.

        // init psi
        if (GlobalV::GAMMA_ONLY_LOCAL)
        {
            if(this->psid==nullptr)
            {
                int ncol = this->LOWF.ParaV->ncol_bands;
                if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="lapack_gvx"
#ifdef __CUSOLVER_LCAO
                ||GlobalV::KS_SOLVER=="cusolver"
#endif
                )
                {
                    ncol = this->LOWF.ParaV->ncol;
                }
                this->psid = new psi::Psi<double>(GlobalV::NSPIN, ncol, this->LOWF.ParaV->nrow, nullptr);
            }
        }
        else
        {   
            if(this->psi == nullptr)
            {
                #ifdef __MPI
                int ncol = this->LOWF.ParaV->ncol_bands;
                #else
                int ncol = GlobalV::NBANDS;
                #endif
#ifdef __CUSOLVER_LCAO
                if(GlobalV::KS_SOLVER=="cusolver")
                {
                    ncol = this->LOWF.paraV->ncol;
                }
#endif
                this->psi = new psi::Psi<std::complex<double>>(this->kv.nks, ncol, this->LOWF.ParaV->nrow, nullptr);
            }
        }
        
        // prepare grid in Gint
        this->UHM.grid_prepare(this->GridT, *pw_rho, *this->pw_big);

        // init Hamiltonian
        if (this->p_hamilt != nullptr)
        {
            delete this->p_hamilt;
            this->p_hamilt = nullptr;
        }
        if(this->p_hamilt == nullptr)
        {
            // two cases for hamilt class
            // Gamma_only case
            if (GlobalV::GAMMA_ONLY_LOCAL)
            {
                this->p_hamilt = new hamilt::HamiltLCAO<double>(&(this->UHM.GG),
                                                                &(this->UHM.genH),
                                                                &(this->LM),
                                                                &(this->LOC),
                                                                this->pelec->pot,
                                                                this->kv);
            }
            // multi_k case
            else
            {
                this->p_hamilt = new hamilt::HamiltLCAO<std::complex<double>>(&(this->UHM.GK),
                                                                              &(this->UHM.genH),
                                                                              &(this->LM),
                                                                              &(this->LOC),
                                                                              this->pelec->pot,
                                                                              this->kv);
            }
        }

        // init density kernel and wave functions.
        this->LOC.allocate_dm_wfc(this->GridT, this->pelec, this->LOWF, this->psid, this->psi, this->kv);

        //======================================
        // do the charge extrapolation before the density matrix is regenerated.
        // mohan add 2011-04-08
        // because once atoms are moving out of this processor,
        // the density matrix will not std::map the new atomic configuration,
        //======================================
        // THIS IS A BUG, BECAUSE THE INDEX GlobalC::GridT.trace_lo
        // HAS BEEN REGENERATED, SO WE NEED TO
        // REALLOCATE DENSITY MATRIX FIRST, THEN READ IN DENSITY MATRIX,
        // AND USE DENSITY MATRIX TO DO RHO GlobalV::CALCULATION.-- mohan 2013-03-31
        //======================================
        if (GlobalV::chg_extrap == "dm" && istep > 1)//xiaohui modify 2015-02-01
        {
            for (int is = 0; is < GlobalV::NSPIN; is++)
            {
                ModuleBase::GlobalFunc::ZEROS(pelec->charge->rho[is], pw_rho->nrxx);
                std::stringstream ssd;
                ssd << GlobalV::global_out_dir << "SPIN" << is + 1 << "_DM";
                // reading density matrix,
                double& ef_tmp = this->pelec->eferm.get_ef(is);
                ModuleIO::read_dm(
#ifdef __MPI
		            this->GridT.nnrg,
		            this->GridT.trace_lo,
#endif
		            is,
		            ssd.str(),
		            this->LOC.DM,
		            this->LOC.DM_R,
		            ef_tmp,
		            &(GlobalC::ucell));
            }

            // calculate the charge density
            if (GlobalV::GAMMA_ONLY_LOCAL)
            {
                Gint_inout inout(this->LOC.DM, pelec->charge->rho, Gint_Tools::job_type::rho);
                this->UHM.GG.cal_gint(&inout);
                if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type()==5)
                {
                    for(int is=0; is<GlobalV::NSPIN; is++)
                    {
                        ModuleBase::GlobalFunc::ZEROS(pelec->charge->kin_r[0], pw_rho->nrxx);
                    }
                    Gint_inout inout1(this->LOC.DM, pelec->charge->kin_r, Gint_Tools::job_type::tau);
                    this->UHM.GG.cal_gint(&inout1);
                }
            }
            else
            {
                Gint_inout inout(this->LOC.DM_R, pelec->charge->rho, Gint_Tools::job_type::rho);
                this->UHM.GK.cal_gint(&inout);
                if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type()==5)
                {
                    for(int is=0; is<GlobalV::NSPIN; is++)
                    {
                        ModuleBase::GlobalFunc::ZEROS(pelec->charge->kin_r[0], pw_rho->nrxx);
                    }
                    Gint_inout inout1(this->LOC.DM_R, pelec->charge->kin_r, Gint_Tools::job_type::tau);
                    this->UHM.GK.cal_gint(&inout1);
                }
            }

            // renormalize the charge density
            pelec->charge->renormalize_rho();
        }

#ifdef __DEEPKS
        //for each ionic step, the overlap <psi|alpha> must be rebuilt
        //since it depends on ionic positions
        if (GlobalV::deepks_setorb)
        {
            const Parallel_Orbitals* pv = this->UHM.LM->ParaV;
            //build and save <psi(0)|alpha(R)> at beginning
            GlobalC::ld.build_psialpha(GlobalV::CAL_FORCE,
                GlobalC::ucell,
                GlobalC::ORB,
                GlobalC::GridD,
                pv->trace_loc_row,
                pv->trace_loc_col,
                GlobalC::UOT);

            if (GlobalV::deepks_out_unittest)
            {
                GlobalC::ld.check_psialpha(GlobalV::CAL_FORCE,
                    GlobalC::ucell,
                    GlobalC::ORB,
                    GlobalC::GridD,
                    pv->trace_loc_row,
                    pv->trace_loc_col,
                    GlobalC::UOT);
            }
        }
#endif
        ModuleBase::timer::tick("ESolver_KS_LCAO", "beforesolver");

    }

    void ESolver_KS_LCAO::beforescf(int istep)
    {
        ModuleBase::TITLE("ESolver_KS_LCAO", "beforescf");
        ModuleBase::timer::tick("ESolver_KS_LCAO", "beforescf");

        if (GlobalC::ucell.cell_parameter_updated)
        {
            this->init_after_vc(INPUT, GlobalC::ucell);
        }
        if (GlobalC::ucell.ionic_position_updated && GlobalV::md_prec_level != 2)
        {
            CE.update_all_dis(GlobalC::ucell);
            CE.extrapolate_charge(
#ifdef __MPI
                &(GlobalC::Pgrid),
#endif
                GlobalC::ucell,
                pelec->charge,
                &(sf));
        }

        //----------------------------------------------------------
        // about vdw, jiyy add vdwd3 and linpz add vdwd2
        //----------------------------------------------------------
        auto vdw_solver = vdw::make_vdw(GlobalC::ucell, INPUT);
        if (vdw_solver != nullptr)
        {
            this->pelec->f_en.evdw = vdw_solver->get_energy();
        }
        
        this->beforesolver(istep);
        this->pelec->init_scf(istep, sf.strucFac);
        // the electron charge density should be symmetrized,
        // here is the initialization
        Symmetry_rho srho;
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            srho.begin(is, *(pelec->charge), pw_rho, GlobalC::Pgrid, this->symm);
        }
//Peize Lin add 2016-12-03
#ifdef __EXX
        if (GlobalC::exx_info.info_ri.real_number)
            this->exd->exx_beforescf(kv, *this->p_chgmix);
        else
            this->exc->exx_beforescf(kv, *this->p_chgmix);
#endif // __EXX
        // 1. calculate ewald energy.
        // mohan update 2021-02-25
        if(!GlobalV::test_skip_ewald)
        {
            this->pelec->f_en.ewald_energy = H_Ewald_pw::compute_ewald(GlobalC::ucell, pw_rho, sf.strucFac);
        }

        p_hamilt->non_first_scf = istep;

        ModuleBase::timer::tick("ESolver_KS_LCAO", "beforescf");
        return;
    }

    void ESolver_KS_LCAO::othercalculation(const int istep)
    {
        ModuleBase::TITLE("ESolver_KS_LCAO", "othercalculation");
        ModuleBase::timer::tick("ESolver_KS_LCAO", "othercalculation");
        if(GlobalV::CALCULATION == "get_S")
        {
            this->get_S();
            ModuleBase::QUIT();
        }
        
        if(GlobalV::CALCULATION == "test_memory")
        {
            Cal_Test::test_memory(this->pw_rho,
                                  this->pw_wfc,
                                  this->p_chgmix->get_mixing_mode(),
                                  this->p_chgmix->get_mixing_ndim());
            return;
        }

        if(GlobalV::CALCULATION == "test_neighbour")
        {
            //test_search_neighbor();
            GlobalV::SEARCH_RADIUS = atom_arrange::set_sr_NL(
                GlobalV::ofs_running,
                GlobalV::OUT_LEVEL,
                GlobalC::ORB.get_rcutmax_Phi(),
                GlobalC::ucell.infoNL.get_rcutmax_Beta(),
                GlobalV::GAMMA_ONLY_LOCAL);

            atom_arrange::search(
                GlobalV::SEARCH_PBC,
                GlobalV::ofs_running,
                GlobalC::GridD,
                GlobalC::ucell,
                GlobalV::SEARCH_RADIUS,
                GlobalV::test_atom_input,
                1);
            return;
        }

        this->beforesolver(istep);
        //pelec should be initialized before these calculations
        this->pelec->init_scf(istep, sf.strucFac);
        // self consistent calculations for electronic ground state
        if (GlobalV::CALCULATION == "nscf")
        {
            this->nscf();
        }
        else if (GlobalV::CALCULATION == "get_pchg")
        {
            IState_Charge ISC(this->psid, this->LOC);
            ISC.begin(this->UHM.GG, this->pelec, this->pw_rho, this->pw_big);
        }
        else if (GlobalV::CALCULATION == "get_wf")
        {
            IState_Envelope IEP(this->pelec);
            if (GlobalV::GAMMA_ONLY_LOCAL)
                IEP.begin(this->psid,
                          this->pw_rho,
                          this->pw_wfc,
                          this->pw_big,
                          this->LOWF,
                          this->UHM.GG,
                          INPUT.out_wfc_pw,
                          wf.out_wfc_r,
                          this->kv);
            else
                IEP.begin(this->psi,
                          this->pw_rho,
                          this->pw_wfc,
                          this->pw_big,
                          this->LOWF,
                          this->UHM.GK,
                          INPUT.out_wfc_pw,
                          wf.out_wfc_r,
                          this->kv);
        }
        else
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::othercalculation", "CALCULATION type not supported");
        }

        ModuleBase::timer::tick("ESolver_KS_LCAO", "othercalculation");
        return;
    }

    void ESolver_KS_LCAO::get_S()
    {
        ModuleBase::TITLE("ESolver_KS_LCAO", "get_S");
        if(GlobalV::GAMMA_ONLY_LOCAL)
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::get_S", "not implemented for");
        }
        else
        {
            // (1) Find adjacent atoms for each atom.
            GlobalV::SEARCH_RADIUS = atom_arrange::set_sr_NL(
                GlobalV::ofs_running,
                GlobalV::OUT_LEVEL,
                GlobalC::ORB.get_rcutmax_Phi(),
                GlobalC::ucell.infoNL.get_rcutmax_Beta(),
                GlobalV::GAMMA_ONLY_LOCAL);

            atom_arrange::search(
                GlobalV::SEARCH_PBC,
                GlobalV::ofs_running,
                GlobalC::GridD,
                GlobalC::ucell,
                GlobalV::SEARCH_RADIUS,
                GlobalV::test_atom_input);

            this->RA.for_2d(this->orb_con.ParaV, GlobalV::GAMMA_ONLY_LOCAL);
            this->UHM.genH.LM->ParaV = &this->orb_con.ParaV;
            this->LM.allocate_HS_R(this->orb_con.ParaV.nnr);
            this->LM.zeros_HSR('S');
            this->UHM.genH.calculate_S_no(this->LM.SlocR.data());
            ModuleIO::output_S_R(this->UHM,"SR.csr");
        }
    }

    void ESolver_KS_LCAO::nscf()
    {
        ModuleBase::TITLE("ESolver_KS_LCAO", "nscf");

        std::cout << " NON-SELF CONSISTENT CALCULATIONS" << std::endl;

        time_t time_start = std::time(NULL);

#ifdef __EXX
#ifdef __MPI
        // Peize Lin add 2018-08-14
        if ( GlobalC::exx_info.info_global.cal_exx )
        {
            //GlobalC::exx_lcao.cal_exx_elec_nscf(this->LOWF.ParaV[0]);
			const std::string file_name_exx = GlobalV::global_out_dir + "HexxR_" + std::to_string(GlobalV::MY_RANK);
			if(GlobalC::exx_info.info_ri.real_number)
				this->exd->read_Hexxs(file_name_exx);
			else
				this->exc->read_Hexxs(file_name_exx);

            // This is a temporary fix
            if(GlobalV::GAMMA_ONLY_LOCAL)
            {
                hamilt::Operator<double>* exx
                    = new hamilt::OperatorEXX<hamilt::OperatorLCAO<double>>(&LM,
                                                                            nullptr, // no explicit call yet
                                                                            &(LM.Hloc),
                                                                            this->kv);
                p_hamilt->opsd->add(exx);
            }
            else
            {
                hamilt::Operator<std::complex<double>>* exx
                    = new hamilt::OperatorEXX<hamilt::OperatorLCAO<std::complex<double>>>(
                        &LM,
                        nullptr, // no explicit call yet
                        &(LM.Hloc2),
                        this->kv);
                p_hamilt->ops->add(exx);
            }
        }
#endif // __MPI
#endif // __EXX

        // mohan add 2021-02-09
        // in ions, istep starts from 1,
        // then when the istep is a variable of scf or nscf,
        // istep becomes istep-1, this should be fixed in future
        int istep = 0;
        if(this->phsol != nullptr)
        {
            if(this->psi != nullptr)
            {
                this->phsol->solve(this->p_hamilt, this->psi[0], this->pelec, GlobalV::KS_SOLVER, true);
            }
            else if(this->psid != nullptr)
            {
                this->phsol->solve(this->p_hamilt, this->psid[0], this->pelec, GlobalV::KS_SOLVER, true);
            }
        }
        else
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_PW", "HSolver has not been initialed!");
        }

        time_t time_finish = std::time(NULL);
        ModuleBase::GlobalFunc::OUT_TIME("cal_bands", time_start, time_finish);

        GlobalV::ofs_running << " end of band structure calculation " << std::endl;
        GlobalV::ofs_running << " band eigenvalue in this processor (eV) :" << std::endl;

        for (int ik = 0; ik < this->kv.nks; ik++)
        {
            if (GlobalV::NSPIN == 2)
            {
                if (ik == 0)
                {
                    GlobalV::ofs_running << " spin up :" << std::endl;
                }
                if (ik == (this->kv.nks / 2))
                {
                    GlobalV::ofs_running << " spin down :" << std::endl;
                }
            }

            GlobalV::ofs_running << " k-points" << ik + 1 << "(" << this->kv.nkstot << "): " << this->kv.kvec_c[ik].x << " "
                                 << this->kv.kvec_c[ik].y << " " << this->kv.kvec_c[ik].z << std::endl;

            for (int ib = 0; ib < GlobalV::NBANDS; ib++)
            {
                GlobalV::ofs_running << " spin" << this->kv.isk[ik] + 1 << "final_state " << ib + 1 << " "
                                     << this->pelec->ekb(ik, ib) * ModuleBase::Ry_to_eV << " "
                                     << this->pelec->wg(ik, ib) * this->kv.nks << std::endl;
            }
            GlobalV::ofs_running << std::endl;
        }
        if (GlobalV::out_bandgap)
        {
            if (!GlobalV::TWO_EFERMI)
            {
                this->pelec->cal_bandgap();
                GlobalV::ofs_running << " E_bandgap " << this->pelec->bandgap * ModuleBase::Ry_to_eV << " eV"
                                     << std::endl;
            }
            else
            {
                this->pelec->cal_bandgap_updw();
                GlobalV::ofs_running << " E_bandgap_up " << this->pelec->bandgap_up * ModuleBase::Ry_to_eV << " eV"
                                     << std::endl;
                GlobalV::ofs_running << " E_bandgap_dw " << this->pelec->bandgap_dw * ModuleBase::Ry_to_eV << " eV"
                                     << std::endl;
            }
        
        }

        // add by jingan in 2018.11.7
        if (GlobalV::CALCULATION == "nscf" && INPUT.towannier90)
        {
            toWannier90 myWannier(this->kv.nkstot, GlobalC::ucell.G, this->LOWF.wfc_k_grid);
            myWannier.init_wannier_lcao(this->GridT,
                                        this->pelec->ekb,
                                        this->pw_rho,
                                        this->pw_wfc,
                                        this->pw_big,
                                        this->sf,
                                        this->kv,
                                        nullptr);
        }

        // add by jingan
        if (berryphase::berry_phase_flag && ModuleSymmetry::Symmetry::symm_flag != 1)
        {
            berryphase bp(this->LOWF);
            bp.Macroscopic_polarization(this->pw_wfc->npwk_max, this->psi, this->pw_rho, this->pw_wfc, this->kv);
        }

        //below is for DeePKS NSCF calculation
#ifdef __DEEPKS
        const Parallel_Orbitals* pv = this->LOWF.ParaV;
        if (GlobalV::deepks_out_labels || GlobalV::deepks_scf)
        {
            if (GlobalV::GAMMA_ONLY_LOCAL)
            {
                GlobalC::ld.cal_projected_DM(this->LOC.dm_gamma[0],
                                         GlobalC::ucell,
                                         GlobalC::ORB,
                                         GlobalC::GridD,
                                         pv->trace_loc_row,
                                         pv->trace_loc_col);
            }
            else
            {
                GlobalC::ld.cal_projected_DM_k(this->LOC.dm_k,
                                               GlobalC::ucell,
                                               GlobalC::ORB,
                                               GlobalC::GridD,
                                               pv->trace_loc_row,
                                               pv->trace_loc_col,
                                               this->kv.nks,
                                               this->kv.kvec_d);
            }
            GlobalC::ld.cal_descriptor(); // final descriptor
            GlobalC::ld.cal_gedm(GlobalC::ucell.nat);
            if (GlobalV::GAMMA_ONLY_LOCAL)
            {
                GlobalC::ld.add_v_delta(GlobalC::ucell,
                                        GlobalC::ORB,
                                        GlobalC::GridD,
                                        pv->trace_loc_row,
                                        pv->trace_loc_col,
                                        pv->nrow,
                                        pv->ncol);
            }
            else
            {
                GlobalC::ld.add_v_delta_k(GlobalC::ucell, 
                                          GlobalC::ORB, 
                                          GlobalC::GridD, 
                                          pv->trace_loc_row,
                                          pv->trace_loc_col,
                                          pv->nnr);
            }
        }
#endif
        return;
    }

}
