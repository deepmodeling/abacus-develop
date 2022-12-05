#include "module_esolver/esolver_ks_lcao.h"
#include "../src_pw/global.h"
#include "../src_pw/symmetry_rho.h"
#include "src_lcao/LCAO_evolve.h"
#include "src_lcao/dftu.h"
//
#include "../module_neighbor/sltk_atom_arrange.h"
#include "../src_io/istate_charge.h"
#include "../src_io/istate_envelope.h"
#include "src_lcao/ELEC_evolve.h"
//
#include "../src_ri/exx_abfs-jle.h"
#include "../src_ri/exx_opt_orb.h"
#include "../src_io/berryphase.h"
#include "../src_io/to_wannier90.h"
#include "../module_base/timer.h"
#ifdef __DEEPKS
#include "../module_deepks/LCAO_deepks.h"
#endif
#include "../src_pw/H_Ewald_pw.h"
#include "module_vdw/vdw.h"
#include "../module_relax/relax_old/variable_cell.h"    // liuyu 2022-11-07

#include "module_hamilt/ks_lcao/op_exx_lcao.h"

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
        GlobalC::GridT.set_pbc_grid(
            GlobalC::rhopw->nx, GlobalC::rhopw->ny, GlobalC::rhopw->nz,
            GlobalC::bigpw->bx, GlobalC::bigpw->by, GlobalC::bigpw->bz,
            GlobalC::bigpw->nbx, GlobalC::bigpw->nby, GlobalC::bigpw->nbz,
            GlobalC::bigpw->nbxx, GlobalC::bigpw->nbzp_start, GlobalC::bigpw->nbzp);

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
            // using GlobalC::GridT.init.
            GlobalC::GridT.cal_nnrg(pv);
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
                this->psi = new psi::Psi<std::complex<double>>(GlobalC::kv.nks, ncol, this->LOWF.ParaV->nrow, nullptr);
            }
        }

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
                                                            this->pelec->pot);
            }
            // multi_k case
            else
            {
                this->p_hamilt = new hamilt::HamiltLCAO<std::complex<double>>(&(this->UHM.GK),
                                                                        &(this->UHM.genH),
                                                                        &(this->LM),
                                                                        &(this->LOC),
                                                                        this->pelec->pot);
            }
        }

        // prepare grid in Gint
        this->UHM.grid_prepare();

        // init density kernel and wave functions.
        this->LOC.allocate_dm_wfc(GlobalC::GridT.lgd, this->pelec, this->LOWF, this->psid, this->psi);

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
                ModuleBase::GlobalFunc::ZEROS(pelec->charge->rho[is], GlobalC::rhopw->nrxx);
                std::stringstream ssd;
                ssd << GlobalV::global_out_dir << "SPIN" << is + 1 << "_DM";
                // reading density matrix,
                this->LOC.read_dm(is, ssd.str());
            }

            // calculate the charge density
            if (GlobalV::GAMMA_ONLY_LOCAL)
            {
                Gint_inout inout(this->LOC.DM, pelec->charge, Gint_Tools::job_type::rho);
                this->UHM.GG.cal_gint(&inout);
                if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type()==5)
                {
                    for(int is=0; is<GlobalV::NSPIN; is++)
                    {
                        ModuleBase::GlobalFunc::ZEROS(pelec->charge->kin_r[0], GlobalC::rhopw->nrxx);
                    }
                    Gint_inout inout1(this->LOC.DM, pelec->charge, Gint_Tools::job_type::tau);
                    this->UHM.GG.cal_gint(&inout1);
                }
            }
            else
            {
                Gint_inout inout(this->LOC.DM_R, pelec->charge, Gint_Tools::job_type::rho);
                this->UHM.GK.cal_gint(&inout);
                if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type()==5)
                {
                    for(int is=0; is<GlobalV::NSPIN; is++)
                    {
                        ModuleBase::GlobalFunc::ZEROS(pelec->charge->kin_r[0], GlobalC::rhopw->nrxx);
                    }
                    Gint_inout inout1(this->LOC.DM_R, pelec->charge, Gint_Tools::job_type::tau);
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

        // Temporary, md and relax will merge later   liuyu add 2022-11-07
        if(GlobalV::CALCULATION == "md" && istep)
        {
            CE.update_istep();
            CE.save_pos_next(GlobalC::ucell);
            CE.extrapolate_charge(pelec->charge);

            if(GlobalC::ucell.cell_parameter_updated)
            {
                Variable_Cell::init_after_vc();
            }
        }

        if(GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax")
        {
            if(GlobalC::ucell.ionic_position_updated)
            {
                GlobalV::ofs_running << " Setup the extrapolated charge." << std::endl;
                // charge extrapolation if istep>0.
                CE.update_istep();
                CE.update_all_pos(GlobalC::ucell);
                CE.extrapolate_charge(pelec->charge);
                CE.save_pos_next(GlobalC::ucell);

                GlobalV::ofs_running << " Setup the Vl+Vh+Vxc according to new structure factor and new charge." << std::endl;
                // calculate the new potential accordint to
                // the new charge density.
            }
        }

        //----------------------------------------------------------
        // about vdw, jiyy add vdwd3 and linpz add vdwd2
        //----------------------------------------------------------
        auto vdw_solver = vdw::make_vdw(GlobalC::ucell, INPUT);
        if (vdw_solver != nullptr)
        {
            GlobalC::en.evdw = vdw_solver->get_energy();
        }
        
        this->beforesolver(istep);
        this->pelec->init_scf( istep, GlobalC::sf.strucFac );
        // the electron charge density should be symmetrized,
        // here is the initialization
        Symmetry_rho srho;
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            srho.begin(is, *(pelec->charge), GlobalC::rhopw, GlobalC::Pgrid, GlobalC::symm);
        }
//Peize Lin add 2016-12-03
#ifdef __EXX
#ifdef __MPI
		if ( GlobalC::exx_info.info_global.cal_exx )
		{
            if (GlobalC::ucell.atoms[0].ncpp.xc_func == "HSE" || GlobalC::ucell.atoms[0].ncpp.xc_func == "PBE0")
            {
                XC_Functional::set_xc_type("pbe");
            }
            else if (GlobalC::ucell.atoms[0].ncpp.xc_func == "SCAN0")
            {
                XC_Functional::set_xc_type("scan");
            }

			//GlobalC::exx_lcao.cal_exx_ions(*this->LOWF.ParaV);
			if(GlobalV::GAMMA_ONLY_LOCAL)
				GlobalC::exx_lri_double.cal_exx_ions();
			else
				GlobalC::exx_lri_complex.cal_exx_ions();
		}

		if (Exx_Abfs::Jle::generate_matrix)
		{
			//program should be stopped after this judgement
			Exx_Opt_Orb exx_opt_orb;
			exx_opt_orb.generate_matrix();
			ModuleBase::timer::tick("ESolver_KS_LCAO", "beforescf");
			return;
		}
#endif // __MPI
#endif // __EXX
        // 1. calculate ewald energy.
        // mohan update 2021-02-25
        if(!GlobalV::test_skip_ewald)
        {
            H_Ewald_pw::compute_ewald(GlobalC::ucell, GlobalC::rhopw);
        }

        p_hamilt->non_first_scf = istep;

        // for exx two_level scf
        this->two_level_step = 0;

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
            Cal_Test::test_memory();
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
        this->pelec->init_scf( istep, GlobalC::sf.strucFac );
        // self consistent calculations for electronic ground state
        if (GlobalV::CALCULATION == "nscf")
        {
            this->nscf();
        }
        else if (GlobalV::CALCULATION == "istate")
        {
            IState_Charge ISC(this->psid, this->LOC);
            ISC.begin(this->UHM.GG, this->pelec);
        }
        else if (GlobalV::CALCULATION == "ienvelope")
        {
            IState_Envelope IEP(this->pelec);
            if (GlobalV::GAMMA_ONLY_LOCAL)
                IEP.begin(this->psid, this->LOWF, this->UHM.GG, INPUT.out_wfc_pw, GlobalC::wf.out_wfc_r);
            else
                IEP.begin(this->psi, this->LOWF, this->UHM.GK, INPUT.out_wfc_pw, GlobalC::wf.out_wfc_r);
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
            this->output_SR("SR.csr");

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
			if(GlobalV::GAMMA_ONLY_LOCAL)
				GlobalC::exx_lri_double.read_Hexxs(file_name_exx);
			else
				GlobalC::exx_lri_complex.read_Hexxs(file_name_exx);

            // This is a temporary fix
            if(GlobalV::GAMMA_ONLY_LOCAL)
            {
                hamilt::Operator<double>* exx
                    = new hamilt::OperatorEXX<hamilt::OperatorLCAO<double>>(
                        &LM,
                        nullptr, //no explicit call yet
                        &(LM.Hloc)
                    );
                p_hamilt->opsd->add(exx);
            }
            else
            {
                hamilt::Operator<std::complex<double>>* exx
                    = new hamilt::OperatorEXX<hamilt::OperatorLCAO<std::complex<double>>>(
                        &LM,
                        nullptr, //no explicit call yet
                        &(LM.Hloc2)
                    );
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

        for (int ik = 0; ik < GlobalC::kv.nks; ik++)
        {
            if (GlobalV::NSPIN == 2)
            {
                if (ik == 0)
                {
                    GlobalV::ofs_running << " spin up :" << std::endl;
                }
                if (ik == (GlobalC::kv.nks / 2))
                {
                    GlobalV::ofs_running << " spin down :" << std::endl;
                }
            }

            GlobalV::ofs_running << " k-points"
                << ik + 1 << "(" << GlobalC::kv.nkstot << "): "
                << GlobalC::kv.kvec_c[ik].x << " " << GlobalC::kv.kvec_c[ik].y << " " << GlobalC::kv.kvec_c[ik].z << std::endl;

            for (int ib = 0; ib < GlobalV::NBANDS; ib++)
            {
                GlobalV::ofs_running << " spin" << GlobalC::kv.isk[ik] + 1
                    << "final_state " << ib + 1 << " "
                    << this->pelec->ekb(ik, ib) * ModuleBase::Ry_to_eV
                    << " " << this->pelec->wg(ik, ib) * GlobalC::kv.nks << std::endl;
            }
            GlobalV::ofs_running << std::endl;
        }

        // add by jingan in 2018.11.7
        if (GlobalV::CALCULATION == "nscf" && INPUT.towannier90)
        {
            toWannier90 myWannier(GlobalC::kv.nkstot, GlobalC::ucell.G, this->LOWF.wfc_k_grid);
            myWannier.init_wannier(this->pelec->ekb, nullptr);
        }

        // add by jingan
        if (berryphase::berry_phase_flag && ModuleSymmetry::Symmetry::symm_flag != 1)
        {
            berryphase bp(this->LOWF);
            bp.Macroscopic_polarization(this->psi);
        }

        return;
    }

}
