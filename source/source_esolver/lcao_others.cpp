#include "source_lcao/hamilt_lcao.h"
#include "source_lcao/module_dftu/dftu.h"
#include "source_esolver/esolver_ks_lcao.h"
#include "source_estate/cal_ux.h"
#include "source_estate/module_charge/symmetry_rho.h"
#include "source_pw/module_pwdft/global.h"
//
#include "source_io/berryphase.h"
#include "source_io/get_pchg_lcao.h"
#include "source_io/get_wf_lcao.h"
#include "source_io/write_HS_R.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/timer.h"
#include "source_cell/module_neighbor/sltk_atom_arrange.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_lcao/LCAO_domain.h"
#include "source_lcao/module_operator_lcao/op_exx_lcao.h"
#include "source_lcao/module_operator_lcao/operator_lcao.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"
#include "source_io/read_wfc_nao.h"
#include "source_io/write_elecstat_pot.h"
#include "source_base/formatter.h"
#include "source_estate/elecstate_lcao.h"
#include "source_estate/module_dm/cal_dm_psi.h"

#ifdef __EXX
#include "source_io/restart_exx_csr.h"
#endif

// mohan add 2025-03-06
#include "source_io/cal_test.h"

namespace ModuleESolver
{

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::others(UnitCell& ucell, const int istep)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "others");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "others");

    const std::string cal_type = PARAM.inp.calculation;

    if (cal_type == "test_memory")
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "testing memory");
        Cal_Test::test_memory(ucell.nat,
                              ucell.ntype,
                              ucell.GGT,
                              this->pw_rho,
                              this->pw_wfc,
                              this->p_chgmix->get_mixing_mode(),
                              this->p_chgmix->get_mixing_ndim());
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "testing memory");
        return;
    }
    else if (cal_type == "test_neighbour")
    {
        // test_search_neighbor();
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "testing neighbour");
        double search_radius = PARAM.inp.search_radius;
        atom_arrange::search(PARAM.globalv.search_pbc,
                             GlobalV::ofs_running,
                             this->gd,
                             ucell,
                             search_radius,
                             PARAM.inp.test_atom_input,
                             true);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "testing neighbour");
        return;
    }
    else if (cal_type == "gen_opt_abfs")
    {
        return;
    }

    // 1. prepare HS matrices, prepare grid integral
    // (1) Find adjacent atoms for each atom.
    double search_radius = atom_arrange::set_sr_NL(GlobalV::ofs_running,
                                                   PARAM.inp.out_level,
                                                   orb_.get_rcutmax_Phi(),
                                                   ucell.infoNL.get_rcutmax_Beta(),
                                                   PARAM.globalv.gamma_only_local);

    atom_arrange::search(PARAM.globalv.search_pbc,
                         GlobalV::ofs_running,
                         this->gd,
                         ucell,
                         search_radius,
                         PARAM.inp.test_atom_input);

    // (3) Periodic condition search for each grid.
    gint_info_.reset(
        new ModuleGint::GintInfo(
        this->pw_big->nbx,
        this->pw_big->nby,
        this->pw_big->nbz,
        this->pw_rho->nx,
        this->pw_rho->ny,
        this->pw_rho->nz,
        0,
        0,
        this->pw_big->nbzp_start,
        this->pw_big->nbx,
        this->pw_big->nby,
        this->pw_big->nbzp,
        orb_.Phi,
        ucell,
        this->gd));
    ModuleGint::Gint::set_gint_info(gint_info_.get());

    // (2)For each atom, calculate the adjacent atoms in different cells
    // and allocate the space for H(R) and S(R).
    // If k point is used here, allocate HlocR after atom_arrange.
    this->RA.for_2d(ucell, this->gd, this->pv, PARAM.globalv.gamma_only_local, orb_.cutoffs());

    // 2. density matrix extrapolation

    // set the augmented orbitals index.
    // after ParaO and GridT,
    // this information is used to calculate
    // the force.

    // init psi
    if (this->psi == nullptr)
    {
        int nsk = 0;
        int ncol = 0;
        if (PARAM.globalv.gamma_only_local)
        {
            nsk = PARAM.inp.nspin;
            ncol = this->pv.ncol_bands;
            if (PARAM.inp.ks_solver == "genelpa" || PARAM.inp.ks_solver == "elpa" || PARAM.inp.ks_solver == "lapack"
                || PARAM.inp.ks_solver == "pexsi" || PARAM.inp.ks_solver == "cusolver"
                || PARAM.inp.ks_solver == "cusolvermp")
            {
                ncol = this->pv.ncol;
            }
        }
        else
        {
            nsk = this->kv.get_nks();
#ifdef __MPI
            ncol = this->pv.ncol_bands;
#else
            ncol = PARAM.inp.nbands;
#endif
        }
        this->psi = new psi::Psi<TK>(nsk, ncol, this->pv.nrow, this->kv.ngk, true);
    }

    // init wfc from file
    if (istep == 0 && PARAM.inp.init_wfc == "file")
    {
        if (!ModuleIO::read_wfc_nao(PARAM.globalv.global_readin_dir,
                                    this->pv,
                                    *(this->psi),
									this->pelec->ekb,
									this->pelec->wg,
                                    this->pelec->klist->ik2iktot,
                                    this->pelec->klist->get_nkstot(),
                                    PARAM.inp.nspin))
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::others", "read wfc nao failed");
        }
    }

    // init Hamiltonian
    if (this->p_hamilt != nullptr)
    {
        delete this->p_hamilt;
        this->p_hamilt = nullptr;
    }
    if (this->p_hamilt == nullptr)
    {
        elecstate::DensityMatrix<TK, double>* DM = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();
        this->p_hamilt = new hamilt::HamiltLCAO<TK, TR>(
            ucell,
            this->gd,
            &this->pv,
            this->pelec->pot,
            this->kv,
            two_center_bundle_,
            orb_,
            DM,
			this->deepks,
			istep,
			this->exx_nao);
	}


    // for each ionic step, the overlap <phi|alpha> must be rebuilt
    // since it depends on ionic positions
    this->deepks.build_overlap(ucell, orb_, pv, gd, *(two_center_bundle_.overlap_orb_alpha), PARAM.inp);

    if (PARAM.inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
        sc.init_sc(PARAM.inp.sc_thr,
                   PARAM.inp.nsc,
                   PARAM.inp.nsc_min,
                   PARAM.inp.alpha_trial,
                   PARAM.inp.sccut,
                   PARAM.inp.sc_drop_thr,
                   ucell,
                   &(this->pv),
                   PARAM.inp.nspin,
                   this->kv,
                   this->p_hamilt,
                   this->psi,
                   this->pelec);
    }

    //=========================================================
    // cal_ux should be called before init_scf because
    // the direction of ux is used in noncoline_rho
    //=========================================================
    elecstate::cal_ux(ucell);

    // pelec should be initialized before these calculations
    this->pelec->init_scf(istep, ucell, this->Pgrid, this->sf.strucFac, this->locpp.numeric, ucell.symm);

    // self consistent calculations for electronic ground state
    if (cal_type == "get_pchg")
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "getting partial charge");
        Get_pchg_lcao get_pchg(this->psi, &(this->pv));
        if (PARAM.globalv.gamma_only_local)
        {
            get_pchg.begin(this->chr.rho,
                           this->pelec->wg,
                           this->pelec->eferm.get_all_ef(),
                           this->pw_rhod->nrxx,
                           PARAM.inp.out_pchg,
                           PARAM.inp.nbands,
                           PARAM.inp.nelec,
                           PARAM.inp.nspin,
                           &ucell,
                           this->Pgrid,
                           &this->gd,
                           this->kv,
                           PARAM.globalv.global_out_dir,
                           GlobalV::ofs_running);
        }
        else
        {
            get_pchg.begin(this->chr.rho,
                           this->chr.rhog,
                           this->pelec->wg,
                           this->pelec->eferm.get_all_ef(),
                           this->pw_rhod,
                           this->pw_rhod->nrxx,
                           PARAM.inp.out_pchg,
                           PARAM.inp.nbands,
                           PARAM.inp.nelec,
                           PARAM.inp.nspin,
                           &ucell,
                           this->Pgrid,
                           &this->gd,
                           this->kv,
                           PARAM.globalv.global_out_dir,
                           GlobalV::ofs_running,
                           PARAM.inp.if_separate_k,
                           this->chr.ngmc);
        }
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "getting partial charge");
    }
    else if (cal_type == "get_wf")
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "getting wave function");
        Get_wf_lcao get_wf(this->pelec);
        if (PARAM.globalv.gamma_only_local)
        {
            get_wf.begin(ucell,
                         this->psi,
                         this->pw_wfc,
                         this->Pgrid,
                         this->pv,
                         PARAM.inp.out_wfc_pw,
                         this->kv,
                         PARAM.inp.nelec,
                         PARAM.inp.out_wfc_norm,
                         PARAM.inp.out_wfc_re_im,
                         PARAM.inp.nbands,
                         PARAM.inp.nspin,
                         PARAM.globalv.nlocal,
                         PARAM.globalv.global_out_dir,
                         GlobalV::ofs_running);
        }
        else
        {
            get_wf.begin(ucell,
                         this->psi,
                         this->pw_wfc,
                         this->Pgrid,
                         this->pv,
                         PARAM.inp.out_wfc_pw,
                         this->kv,
                         PARAM.inp.nelec,
                         PARAM.inp.out_wfc_norm,
                         PARAM.inp.out_wfc_re_im,
                         PARAM.inp.nbands,
                         PARAM.inp.nspin,
                         PARAM.globalv.nlocal,
                         PARAM.globalv.global_out_dir,
                         GlobalV::ofs_running);
        }
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "getting wave function");
    }
    else
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::others", "CALCULATION type not supported");
    }

    ModuleBase::timer::tick("ESolver_KS_LCAO", "others");
    return;
}

template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
} // namespace ModuleESolver
