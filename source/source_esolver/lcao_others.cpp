#include "source_esolver/esolver_ks_lcao.h"
#include "source_cell/cal_ux.h"
#include "source_estate/module_charge/symm_rho.h"
#include "source_lcao/hamilt_lcao.h"
#include "source_lcao/module_dftu/dftu.h"
#include "source_hamilt/module_gint/gint.h"
#include "source_base/formatter.h"
#include "source_base/timer.h"
#include "source_cell/module_neighbor/sltk_atom_arrange.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_estate/elecstate_lcao.h"
#include "source_estate/module_dm/cal_dm_psi.h"
#include "source_io/module_chgpot/get_pchg_lcao.h"
#include "source_io/module_wf/get_wf_lcao.h"
#include "source_io/module_parameter/parameter.h"
#include "source_io/module_hs/write_hs_r.h"
#include "source_lcao/lcao_domain.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"
#include "source_lcao/module_operator_lcao/op_exx_lcao.h"
#include "source_lcao/module_operator_lcao/operator_lcao.h"

#ifdef __EXX
#endif

// mohan add 2025-03-06
#include "source_io/module_output/cal_test.h"

namespace ModuleESolver
{

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::others(BaseCell& basecell, const int istep)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    ModuleBase::TITLE("ESolver_KS_LCAO", "others");
    ModuleBase::timer::start("ESolver_KS_LCAO", "others");

    const std::string cal_type = this->inp_->calculation;
    const std::string global_out_dir = PARAM.globalv.global_out_dir;
    const bool gamma_only_local = PARAM.globalv.gamma_only_local;

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
        double search_radius = this->inp_->search_radius;
        atom_arrange::search(PARAM.globalv.search_pbc,
                             GlobalV::ofs_running,
                             this->gd,
                             ucell,
                             search_radius,
                             this->inp_->test_atom_input,
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
                                                   this->inp_->out_level,
                                                   orb_.get_rcutmax_Phi(),
                                                   ucell.infoNL->get_rcutmax_Beta(),
                                                   gamma_only_local);

    atom_arrange::search(PARAM.globalv.search_pbc,
                         GlobalV::ofs_running,
                         this->gd,
                         ucell,
                         search_radius,
                         this->inp_->test_atom_input);

    // (3) Periodic condition search for each grid.
    gint_info_.reset(new ModuleGint::GintInfo(this->pw_big->nbx,
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
    this->RA.for_2d(ucell, this->gd, this->pv, gamma_only_local, orb_.cutoffs());

    // 2. density matrix extrapolation

    // set the augmented orbitals index.
    // after ParaO and GridT,
    // this information is used to calculate
    // the force.

    // init psi deleted by taoni 2026-01-23
    // don't need to since initialized in LCAO_domain::set_psi_occ_dm_chg in before_all_runners

    // init Hamiltonian
    if (this->p_hamilt != nullptr)
    {
        delete this->p_hamilt;
        this->p_hamilt = nullptr;
    }
    if (this->p_hamilt == nullptr)
    {
        this->p_hamilt = new hamilt::HamiltLCAO<TK, TR>(ucell,
                                                        this->gd,
                                                        &this->pv,
                                                        this->pelec->pot,
                                                        this->kv,
                                                        two_center_bundle_,
                                                        orb_,
                                                        this->dmat.dm,
                                                        &this->dftu,
                                                        this->deepks,
                                                        istep,
                                                        this->exx_nao,
                                                        this->exx_info_);
    }

    // for each ionic step, the overlap <phi|alpha> must be rebuilt
    // since it depends on ionic positions
    this->deepks.build_overlap(ucell, orb_, pv, gd, *(two_center_bundle_.overlap_orb_alpha), *this->inp_);

    if (this->inp_->sc_mag_switch)
    {
        spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
        sc.init_sc(this->inp_->sc_thr,
                   this->inp_->nsc,
                   this->inp_->nsc_min,
                   this->inp_->alpha_trial,
                   this->inp_->sccut,
                   this->inp_->sc_drop_thr,
                   ucell,
                   this->inp_->sc_direction_only,
                   &(this->pv),
                   this->inp_->nspin,
                   this->kv,
                   this->p_hamilt,
                   this->psi,
                   this->dmat.dm,
                   this->pelec);
    }

    //=========================================================
    // cal_ux should be called before init_scf because
    // the direction of ux is used in noncoline_rho
    //=========================================================
    unitcell::cal_ux(ucell, this->inp_->nspin);

    // pelec should be initialized before these calculations
    elecstate::init_scf(ucell, this->Pgrid, this->sf.strucFac, this->locpp.numeric, 
		    istep, global_out_dir, *this->inp_, this->pelec);

    // self consistent calculations for electronic ground state
    if (cal_type == "get_pchg")
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "getting partial charge");
        Get_pchg_lcao get_pchg(this->psi, &(this->pv));
        if (gamma_only_local)
        {
            get_pchg.begin(this->chr.rho,
                           this->pelec->wg,
                           this->pelec->eferm.get_all_ef(),
                           this->pw_rhod->nrxx,
                           this->inp_->out_pchg,
                           this->inp_->nbands,
                           this->inp_->nelec,
                           this->inp_->nspin,
                           &ucell,
                           this->Pgrid,
                           &this->gd,
                           this->kv,
                           global_out_dir,
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
                           this->inp_->out_pchg,
                           this->inp_->nbands,
                           this->inp_->nelec,
                           this->inp_->nspin,
                           &ucell,
                           this->Pgrid,
                           &this->gd,
                           this->kv,
                           global_out_dir,
                           GlobalV::ofs_running,
                           this->inp_->if_separate_k,
                           this->chr.ngmc);
        }
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "getting partial charge");
    }
    else if (cal_type == "get_wf")
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "getting wave function");
        Get_wf_lcao get_wf(this->pelec);
        if (gamma_only_local)
        {
            get_wf.begin(ucell,
                         this->psi,
                         this->pw_wfc,
                         this->Pgrid,
                         this->pv,
                         this->inp_->out_wfc_pw,
                         this->kv,
                         this->inp_->nelec,
                         this->inp_->out_wfc_norm,
                         this->inp_->out_wfc_re_im,
                         this->inp_->nbands,
                         this->inp_->nspin,
                         PARAM.globalv.nlocal,
                         global_out_dir,
                         GlobalV::ofs_running);
        }
        else
        {
            get_wf.begin(ucell,
                         this->psi,
                         this->pw_wfc,
                         this->Pgrid,
                         this->pv,
                         this->inp_->out_wfc_pw,
                         this->kv,
                         this->inp_->nelec,
                         this->inp_->out_wfc_norm,
                         this->inp_->out_wfc_re_im,
                         this->inp_->nbands,
                         this->inp_->nspin,
                         PARAM.globalv.nlocal,
                         global_out_dir,
                         GlobalV::ofs_running);
        }
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "getting wave function");
    }
    else
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::others", "CALCULATION type not supported");
    }

    ModuleBase::timer::end("ESolver_KS_LCAO", "others");
    return;
}

template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
} // namespace ModuleESolver
