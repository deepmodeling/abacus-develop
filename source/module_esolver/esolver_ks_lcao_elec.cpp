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
#include "src_lcao/ELEC_scf.h"
#include "src_lcao/ELEC_nscf.h"
#include "src_lcao/ELEC_cbands_gamma.h"
#include "src_lcao/ELEC_cbands_k.h"
#include "src_lcao/ELEC_evolve.h"
//
#include "../src_ri/exx_abfs.h"
#include "../src_ri/exx_opt_orb.h"
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
            GlobalC::pw.ncx, GlobalC::pw.ncy, GlobalC::pw.ncz,
            GlobalC::pw.bx, GlobalC::pw.by, GlobalC::pw.bz,
            GlobalC::pw.nbx, GlobalC::pw.nby, GlobalC::pw.nbz,
            GlobalC::pw.nbxx, GlobalC::pw.nbzp_start, GlobalC::pw.nbzp);

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

        // init density kernel and wave functions.
        this->LOC.allocate_dm_wfc(GlobalC::GridT.lgd, this->LOWF);

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
                ModuleBase::GlobalFunc::ZEROS(GlobalC::CHR.rho[is], GlobalC::pw.nrxx);
                std::stringstream ssd;
                ssd << GlobalV::global_out_dir << "SPIN" << is + 1 << "_DM";
                // reading density matrix,
                this->LOC.read_dm(is, ssd.str());
            }

            // calculate the charge density
            if (GlobalV::GAMMA_ONLY_LOCAL)
            {
                this->UHM.GG.cal_rho(this->LOC.DM, (Charge*)(&GlobalC::CHR));
            }
            else
            {
                this->UHM.GK.cal_rho_k(this->LOC.DM_R, (Charge*)(&GlobalC::CHR));
            }

            // renormalize the charge density
            GlobalC::CHR.renormalize_rho();

            // initialize the potential
            GlobalC::pot.init_pot(istep - 1, GlobalC::pw.strucFac);
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
    }

    void ESolver_KS_LCAO::beforescf(int istep)
    {
        ModuleBase::TITLE("ESolver_KS_LCAO", "beforescf");
        ModuleBase::timer::tick("ESolver_KS_LCAO", "beforescf");
        this->beforesolver(istep);
        // 1. calculate ewald energy.
        // mohan update 2021-02-25
        H_Ewald_pw::compute_ewald(GlobalC::ucell, GlobalC::pw);

        //2. the electron charge density should be symmetrized,
        // here is the initialization
        Symmetry_rho srho;
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            srho.begin(is, GlobalC::CHR, GlobalC::pw, GlobalC::Pgrid, GlobalC::symm);
        }

        ModuleBase::timer::tick("ESolver_KS_LCAO", "beforescf");
        return;
    }

    void ESolver_KS_LCAO::othercalculation(const int istep)
    {
        ModuleBase::TITLE("ESolver_KS_LCAO", "othercalculation");
        ModuleBase::timer::tick("ESolver_KS_LCAO", "othercalculation");
        this->beforesolver(istep);
        // self consistent calculations for electronic ground state
        if (GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "md"
            || GlobalV::CALCULATION == "relax" || GlobalV::CALCULATION == "cell-relax") //pengfei 2014-10-13
        {
#ifdef __MPI
            //Peize Lin add 2016-12-03
            if (Exx_Global::Hybrid_Type::HF == GlobalC::exx_lcao.info.hybrid_type
                || Exx_Global::Hybrid_Type::PBE0 == GlobalC::exx_lcao.info.hybrid_type
                || Exx_Global::Hybrid_Type::HSE == GlobalC::exx_lcao.info.hybrid_type)
            {
                GlobalC::exx_lcao.cal_exx_ions(*this->LOWF.ParaV);
            }
            if (Exx_Global::Hybrid_Type::Generate_Matrix == GlobalC::exx_global.info.hybrid_type)
            {
                Exx_Opt_Orb exx_opt_orb;
                exx_opt_orb.generate_matrix();
            }
            else    // Peize Lin add 2016-12-03
            {
#endif // __MPI
                ELEC_scf es;
                es.scf(istep, this->LOC, this->LOWF, this->UHM);
#ifdef __MPI
                if (GlobalC::exx_global.info.separate_loop)
                {
                    for (size_t hybrid_step = 0; hybrid_step != GlobalC::exx_global.info.hybrid_step; ++hybrid_step)
                    {
                        XC_Functional::set_xc_type(GlobalC::ucell.atoms[0].xc_func);
                        GlobalC::exx_lcao.cal_exx_elec(this->LOC, this->LOWF.wfc_k_grid);

                        ELEC_scf es;
                        es.scf(istep, this->LOC, this->LOWF, this->UHM);
                        if (ELEC_scf::iter == 1)     // exx converge
                        {
                            break;
                        }
                    }
                }
                else
                {
                    XC_Functional::set_xc_type(GlobalC::ucell.atoms[0].xc_func);
                    ELEC_scf es;
                    es.scf(istep, this->LOC, this->LOWF, this->UHM);
                }
            }
#endif // __MPI
        }
        else if (GlobalV::CALCULATION == "nscf")
        {
            ELEC_nscf::nscf(this->UHM, this->LOC.dm_gamma, this->LOC.dm_k, this->LOWF);
        }
        else if (GlobalV::CALCULATION == "istate")
        {
            IState_Charge ISC(this->LOWF.wfc_gamma, this->LOC);
            ISC.begin(this->UHM.GG);
        }
        else if (GlobalV::CALCULATION == "ienvelope")
        {
            IState_Envelope IEP;
            if (GlobalV::GAMMA_ONLY_LOCAL)
                IEP.begin(this->LOWF, this->UHM.GG, INPUT.out_wfc_pw, GlobalC::wf.out_wfc_r);
            else
                IEP.begin(this->LOWF, this->UHM.GK, INPUT.out_wfc_pw, GlobalC::wf.out_wfc_r);
        }
        else
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::othercalculation", "CALCULATION type not supported");
        }

        ModuleBase::timer::tick("ESolver_KS_LCAO", "othercalculation");
        return;
    }

}
