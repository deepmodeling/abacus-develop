#include "module_esolver/esolver_ks_lcao.h"
#include "src_lcao/LCAO_diago.h"
#include "../src_pw/global.h"
#include "../src_pw/symmetry_rho.h"
#include "input_update.h"
#include "../src_io/chi0_hilbert.h"
#include "src_lcao/LCAO_evolve.h"
#include "src_lcao/dftu.h"
//
#include "../module_neighbor/sltk_atom_arrange.h"
#include "../src_io/istate_charge.h"
#include "../src_io/istate_envelope.h"
#include "src_lcao/ELEC_evolve.h"
//
#include "../src_ri/exx_abfs.h"
#include "../src_ri/exx_opt_orb.h"
#include "../src_io/berryphase.h"
#include "../src_io/to_wannier90.h"
#include "../src_pw/vdwd2.h"
#include "../src_pw/vdwd3.h"
#include "../module_base/timer.h"
#ifdef __DEEPKS
#include "../module_deepks/LCAO_deepks.h"
#endif
#include "../src_pw/H_Ewald_pw.h"

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

        // 2. density matrix extrapolation and prepare S,T,VNL matrices 

        // set the augmented orbitals index.
        // after ParaO and GridT, 
        // this information is used to calculate
        // the force.

        // init psi
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
                int ncol = this->LOWF.ParaV->ncol_bands;
#ifdef __CUSOLVER_LCAO
                if(GlobalV::KS_SOLVER=="cusolver")
                {
                    ncol = this->LOWF.paraV->ncol;
                }
#endif
                this->psi = new psi::Psi<std::complex<double>>(GlobalC::kv.nks, ncol, this->LOWF.ParaV->nrow, nullptr);
            }
        }
        // init density kernel and wave functions.
        this->LOC.allocate_dm_wfc(GlobalC::GridT.lgd, this->LOWF, this->psid, this->psi);

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
        if (GlobalC::pot.chg_extrap == "dm" && istep > 1)//xiaohui modify 2015-02-01
        {
            for (int is = 0; is < GlobalV::NSPIN; is++)
            {
                ModuleBase::GlobalFunc::ZEROS(GlobalC::CHR.rho[is], GlobalC::rhopw->nrxx);
                std::stringstream ssd;
                ssd << GlobalV::global_out_dir << "SPIN" << is + 1 << "_DM";
                // reading density matrix,
                this->LOC.read_dm(is, ssd.str());
            }

            // calculate the charge density
            if (GlobalV::GAMMA_ONLY_LOCAL)
            {
                Gint_inout inout(this->LOC.DM, (Charge*)(&GlobalC::CHR), Gint_Tools::job_type::rho);
                this->UHM.GG.cal_gint(&inout);
                if (XC_Functional::get_func_type() == 3)
                {
                    for(int is=0; is<GlobalV::NSPIN; is++)
                    {
                        ModuleBase::GlobalFunc::ZEROS(GlobalC::CHR.kin_r[0], GlobalC::rhopw->nrxx);
                    }
                    Gint_inout inout1(this->LOC.DM, (Charge*)(&GlobalC::CHR), Gint_Tools::job_type::tau);
                    this->UHM.GG.cal_gint(&inout1);
                }
            }
            else
            {
                Gint_inout inout(this->LOC.DM_R, (Charge*)(&GlobalC::CHR), Gint_Tools::job_type::rho);
                this->UHM.GK.cal_gint(&inout);
                if (XC_Functional::get_func_type() == 3)
                {
                    for(int is=0; is<GlobalV::NSPIN; is++)
                    {
                        ModuleBase::GlobalFunc::ZEROS(GlobalC::CHR.kin_r[0], GlobalC::rhopw->nrxx);
                    }
                    Gint_inout inout1(this->LOC.DM_R, (Charge*)(&GlobalC::CHR), Gint_Tools::job_type::tau);
                    this->UHM.GK.cal_gint(&inout1);
                }
            }

            // renormalize the charge density
            GlobalC::CHR.renormalize_rho();

            // initialize the potential
            GlobalC::pot.init_pot(istep - 1, GlobalC::sf.strucFac);
        }


        // 3. compute S, T, Vnl, Vna matrix.
        this->UHM.set_lcao_matrices();

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
        this->beforesolver(istep);
//Peize Lin add 2016-12-03
#ifdef __MPI
        if(Exx_Global::Hybrid_Type::No != GlobalC::exx_global.info.hybrid_type)
        {
            if (Exx_Global::Hybrid_Type::HF == GlobalC::exx_lcao.info.hybrid_type
                || Exx_Global::Hybrid_Type::PBE0 == GlobalC::exx_lcao.info.hybrid_type
                || Exx_Global::Hybrid_Type::HSE == GlobalC::exx_lcao.info.hybrid_type)
            {
                GlobalC::exx_lcao.cal_exx_ions(*this->LOWF.ParaV);
            }

            if (Exx_Global::Hybrid_Type::Generate_Matrix == GlobalC::exx_global.info.hybrid_type)
            {
                //program should be stopped after this judgement
                Exx_Opt_Orb exx_opt_orb;
                exx_opt_orb.generate_matrix();
                ModuleBase::timer::tick("ESolver_KS_LCAO", "beforescf");
                return;
            }
        }
#endif
        // 1. calculate ewald energy.
        // mohan update 2021-02-25
        H_Ewald_pw::compute_ewald(GlobalC::ucell, GlobalC::rhopw);

        //2. the electron charge density should be symmetrized,
        // here is the initialization
        Symmetry_rho srho;
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            srho.begin(is, GlobalC::CHR, GlobalC::rhopw, GlobalC::Pgrid, GlobalC::symm);
        }

        phami->non_first_scf = istep;

        // for exx two_level scf
        this->two_level_step = 0;

        ModuleBase::timer::tick("ESolver_KS_LCAO", "beforescf");
        return;
    }

    void ESolver_KS_LCAO::othercalculation(const int istep)
    {
        ModuleBase::TITLE("ESolver_KS_LCAO", "othercalculation");
        ModuleBase::timer::tick("ESolver_KS_LCAO", "othercalculation");
        this->beforesolver(istep);
        // self consistent calculations for electronic ground state
        if (GlobalV::CALCULATION == "nscf")
        {
            this->nscf();
        }
        else if (GlobalV::CALCULATION == "istate")
        {
            IState_Charge ISC(this->psid, this->LOC);
            ISC.begin(this->UHM.GG);
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

    void ESolver_KS_LCAO::nscf()
    {
        ModuleBase::TITLE("ESolver_KS_LCAO", "nscf");

        std::cout << " NON-SELF CONSISTENT CALCULATIONS" << std::endl;

        time_t time_start = std::time(NULL);

    #ifdef __MPI
        // Peize Lin add 2018-08-14
        switch (GlobalC::exx_lcao.info.hybrid_type)
        {
        case Exx_Global::Hybrid_Type::HF:
        case Exx_Global::Hybrid_Type::PBE0:
        case Exx_Global::Hybrid_Type::HSE:
            GlobalC::exx_lcao.cal_exx_elec_nscf(this->LOWF.ParaV[0]);
            break;
        }
    #endif

        // mohan add 2021-02-09
        // in LOOP_ions, istep starts from 1,
        // then when the istep is a variable of scf or nscf,
        // istep becomes istep-1, this should be fixed in future
        int istep = 0;
        if(this->phsol != nullptr)
        {
            if(this->psi != nullptr)
            {
                this->phsol->solve(this->phami, this->psi[0], this->pelec, GlobalV::KS_SOLVER, true);
            }
            else if(this->psid != nullptr)
            {
                this->phsol->solve(this->phami, this->psid[0], this->pelec, GlobalV::KS_SOLVER, true);
            }
            for(int ik=0; ik<this->pelec->ekb.nr; ++ik)
            {
                for(int ib=0; ib<this->pelec->ekb.nc; ++ib)
                {
                    GlobalC::wf.ekb[ik][ib] = this->pelec->ekb(ik, ib);
                }
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
                    << GlobalC::wf.ekb[ik][ib] * ModuleBase::Ry_to_eV
                    << " " << GlobalC::wf.wg(ik, ib) * GlobalC::kv.nks << std::endl;
            }
            GlobalV::ofs_running << std::endl;
        }

        // add by jingan in 2018.11.7
        if (GlobalV::CALCULATION == "nscf" && INPUT.towannier90)
        {
            toWannier90 myWannier(GlobalC::kv.nkstot, GlobalC::ucell.G, this->LOWF.wfc_k_grid);
            myWannier.init_wannier(nullptr);
        }

        // add by jingan
        if (berryphase::berry_phase_flag && ModuleSymmetry::Symmetry::symm_flag == 0)
        {
            berryphase bp(this->LOWF);
            bp.Macroscopic_polarization(nullptr);
        }

        return;
    }

}
