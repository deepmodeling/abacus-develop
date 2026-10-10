// Initialize LR from an embedded KS solver or standalone ground-state files.
#include "esolver_lr_lcao_tddft.h"
#include "source_basis/module_pw/pw_basis_big.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/utils/lr_util_print.h"
#include "source_hamilt/module_gint/gint.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_cell/module_neighbor/sltk_atom_arrange.h"
#include "source_io/module_output/print_info.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_io/module_parameter/parameter.h"
#ifdef __JSON
#include "source_io/module_json/output_info.h"
#endif
#ifdef __EXX
#include "source_lcao/module_ri/exx_lri_interface.h"
#include "source_hamilt/module_xc/exx_info.h"
#endif
#include <algorithm>

#ifdef __EXX
namespace
{

    /// One `Exx_LRI` carries ONE Coulomb operator, and it may be needed for two different
    /// reasons: the LR kernel (when `xc_kernel` is a hybrid) and the ground-state force (when
    /// `dft_functional` is a hybrid). That operator is NOT chosen here -- `Exx_LRI` reads
    /// `info_ri.coulomb_param`, which `input_conv` builds from `dft_functional` alone. So when
    /// the two disagree, the kernel silently gets the ground state's screening, not its own.
    ///
    /// This used to be decided by a local `exx_ccp_type()` returning Erfc for "hse" and bare
    /// Hf otherwise, written onto `info_global.ccp_type`. With general range-separated hybrids
    /// that two-way split is not even expressible ($\alpha/r+\beta\,\mathrm{erfc}(\mu r)/r$
    /// is both at once), and the write was in any case dead for the RI path: only `Exx_LRI`
    /// runs here, and it never looks at `ccp_type`.
    void warn_if_kernel_differs_from_gs(const std::string& xc_kernel, const std::string& dft_functional,
        std::ofstream& ofs_running)
    {
        const bool k = LR::exx_kernel_list().count(xc_kernel) > 0;
        const bool g = LR::exx_kernel_list().count(dft_functional) > 0;
        if (k && xc_kernel != dft_functional)
        {
            ofs_running << " WARNING: xc_kernel (" << xc_kernel << ") and dft_functional ("
                << dft_functional << ") are not the same functional. The Coulomb operator of"
                " Exx_LRI follows dft_functional" << (g ? "" : ", which is not a hybrid at all"
                " (no coulomb_param, so the LR exchange kernel vanishes)") << "; set both to the"
                " same hybrid to get the kernel you asked for." << std::endl;
        }
    }
}
#endif

#ifdef __EXX
template<>
void ModuleESolver::ESolver_LR<double>::share_exx_lri(std::shared_ptr<Exx_LRI<double>>& exx_ks)
{
    ModuleBase::TITLE("ESolver_LR<double>", "share_exx_lri");
    this->exx_lri = exx_ks->make_lr_workspace(*this->ucell_, this->kv);
}
template<>
void ModuleESolver::ESolver_LR<std::complex<double>>::share_exx_lri(std::shared_ptr<Exx_LRI<std::complex<double>>>& exx_ks)
{
    ModuleBase::TITLE("ESolver_LR<complex>", "share_exx_lri");
    this->exx_lri = exx_ks->make_lr_workspace(*this->ucell_, this->kv);
}
template<>
void ModuleESolver::ESolver_LR<std::complex<double>>::share_exx_lri(std::shared_ptr<Exx_LRI<double>>& exx_ks)
{
    throw std::runtime_error("ESolver_LR<std::complex<double>>::share_exx_lri: cannot share double to std::complex<double>");
}
template<>
void ModuleESolver::ESolver_LR<double>::share_exx_lri(std::shared_ptr<Exx_LRI<std::complex<double>>>& exx_ks)
{
    throw std::runtime_error("ESolver_LR<double>::share_exx_lri: cannot share std::complex<double> to double");
}
#endif

using namespace LR;

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::before_all_runners(BaseCell& basecell, const Input_para& inp)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);
    this->ucell_ = &ucell;
    this->inp_ = &inp;
    if (inp.esolver_type == "ks-lr")
    {
#ifdef __JSON
        // The embedded KS run happens before Relax_Driver starts its first step.
        Json::init_output_array_obj();
#endif
        // the ground-state solver is a member, not a temporary: `runner` re-runs its SCF on every
        // ionic step, and the objects aliased from it have to stay alive as long as it does.
        // Its SCF is deliberately NOT run here -- `before_all_runners` is called once, outside the
        // relaxation loop, so the ground state has to be recomputed from `runner(istep)` instead.
        this->ks_ = LR_Util::make_unique<ModuleESolver::ESolver_KS_LCAO<T, TR>>();
        this->ks_->before_all_runners(basecell, inp);
        LR_Util::check_force_pp(ucell, inp.cal_force, inp.calculation);
    }
    else
    {
        this->initialize_from_unitcell_(ucell, inp);
    }
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::bind_ground_state_aliases_()
{
    if (this->ks_)
    {
        // `ESolver_FP::before_all_runners` was never called on this object, so its own `pw_rhod`,
        // `Pgrid`, `sf` and `locpp` are empty -- the excited-state force needs all four. Point at
        // the ground-state solver's, which `before_scf` refreshes for the current geometry every
        // ionic step. `pw_rho_flag` stays false so the destructor does not free what it borrowed.
        this->pw_rho = this->ks_->pw_rho;
        this->pw_rhod = this->ks_->pw_rhod;
        this->pw_big = this->ks_->pw_big;
        this->pw_rho_flag = false;
        this->pgrid_ptr_ = &this->ks_->Pgrid;
        this->sf_ptr_ = &this->ks_->sf;
        this->locpp_ptr_ = &this->ks_->locpp;
        this->gd_ptr_ = &this->ks_->gd;
        // the two-center tables depend on the orbitals only, so the ground-state solver's single
        // build serves every geometry. This used to be moved over only for `abs_gauge velocity`,
        // leaving `LR_Force` and `cal_hs_grad` with an empty bundle in the length gauge.
        this->tcb_ptr_ = &this->ks_->two_center_bundle_;
    }
    else
    {
        this->pgrid_ptr_ = &this->Pgrid;
        this->sf_ptr_ = &this->sf;
        this->locpp_ptr_ = &this->locpp;
        this->gd_ptr_ = &this->gd_own_;
        this->tcb_ptr_ = &this->two_center_bundle_own_;
    }
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::initialize_from_ks_(UnitCell& ucell, const Input_para& inp)
{
    ModuleBase::TITLE("ESolver_LR", "ESolver_LR(KS)");
    ModuleESolver::ESolver_KS_LCAO<T, TR>& ks_sol = *this->ks_;
    this->bind_ground_state_aliases_();

	if (this->inp_->lr_solver == "spectrum")
	{
		throw std::invalid_argument("when lr_solver==spectrum, esolver_type must be `lr` to skip KS calculation.");
	}

    // xc kernel
    this->xc_kernel = LR_Util::tolower(inp.xc_kernel);
    //kv
    this->kv = ks_sol.kv;   // copy: cheap, and the KS solver keeps using its own

    this->parameter_check();

    this->set_dimension();

    // KS and LR must use the same AO distribution, including its BLACS grid.
    const int ks_block_size = ks_sol.pv.get_block_size();
    LR_Util::setup_2d_division(this->paraMat_, ks_block_size, this->nbasis, this->nbasis
#ifdef __MPI
        , ks_sol.pv.blacs_ctxt
#endif
    );
    this->set_parallel_orbitals_band(this->paraMat_, this->nbands);
    if (this->inp_->cal_force)
    {
        LR_Util::setup_2d_division(this->paraMat_all_, ks_block_size, this->nbasis, this->nbasis
#ifdef __MPI
            , ks_sol.pv.blacs_ctxt
#endif
        );
        this->set_parallel_orbitals_band(this->paraMat_all_, this->inp_->nbands);
    }

    this->paraMat_.atom_begin_row = ks_sol.pv.atom_begin_row;
    this->paraMat_.atom_begin_col = ks_sol.pv.atom_begin_col;
    this->paraMat_.iat2iwt_ = ucell.get_iat2iwt();

    LR_Util::setup_2d_division(this->paraC_, ks_block_size, this->nbasis, this->nbands
#ifdef __MPI
        , this->paraMat_.blacs_ctxt
#endif
    );

    // allocate psi_ks and eig_ks in the [nocc, nvirt] window
#ifdef __MPI
    this->psi_ks.reset(new psi::Psi<T>(this->kv.get_nks(),
        this->paraC_.get_col_size(),
        this->paraC_.get_row_size(),
        this->kv.ngk,
        true));
#else
    this->psi_ks.reset(new psi::Psi<T>(this->kv.get_nks(), this->nbands, this->nbasis, this->kv.ngk, true));
#endif
    this->eig_ks.create(this->kv.get_nks(), this->nbands);
    this->pelec = new elecstate::ElecStateLCAO<T>();
    orb_cutoff_ = ks_sol.orb_.cutoffs();

#ifdef __EXX
    // Two independent reasons to need an Exx_LRI: the LR exchange kernel, and the ground-state
    // EXX terms of the gradient (the H_gs[T+Z] and W multipliers, gated on gs_is_hybrid).
    // `initialize_from_unitcell_` has always covered both; this path used to test only the first,
    // so a local kernel on top of a hybrid ground state had no exx_lri at all when asked for forces.
    if (exx_kernel_list().count(xc_kernel) || (this->inp_->cal_force && gs_is_hybrid(this->inp_->dft_functional)))
    {
        std::string dft_functional = LR_Util::tolower(this->inp_->dft_functional);
        // Either object would be built from the same `info_ri.coulomb_param`, which `input_conv`
        // derives from dft_functional alone -- so whenever the ground-state solver has one of the
        // right type its geometry tensors are already up to date. Share those tensors in a
        // separate LR electronic workspace, skipping a `cal_exx_ions` per ionic step.
        const bool share = (ks_sol.exx_nao.exd && std::is_same<T, double>::value)
                        || (ks_sol.exx_nao.exc && std::is_same<T, std::complex<double>>::value);
        warn_if_kernel_differs_from_gs(xc_kernel, dft_functional, this->ofs_running_);
        if (share) { this->exx_owned_ = false; }   // `refresh_from_ks_` refreshes the geometry snapshot
        else    // construct C, V from scratch
        {
            // `input_conv` already filled `info_ri.coulomb_param` from INPUT.
            exx_info.sync_from_global();
            // populate ABFs/JLE file lists from UnitCell; keep in sync with Exx_NAO::init
            exx_info.info_ri.files_abfs = ucell.abfs_orbital_files;
            exx_info.info_opt_abfs.files_abfs = ucell.abfs_orbital_files;
            exx_info.info_opt_abfs.files_jles = ucell.jle_orbital_files;
            this->exx_lri = std::make_shared<Exx_LRI<T>>(exx_info.info_ri);
            this->exx_lri->init(MPI_COMM_WORLD, ucell, this->kv, ks_sol.orb_);
            this->exx_owned_ = true;   // the position-dependent `cal_exx_ions` is left to `refresh_from_ks_`
        }
    }
#endif

    refresh_from_ks_(ucell);
    this->ks_initialized_ = true;
}

template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::refresh_from_ks_(UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_LR", "refresh_from_ks_");
    ModuleESolver::ESolver_KS_LCAO<T, TR>& ks_sol = *this->ks_;
    // The ground-state solver owns `psi` and reuses it across ionic steps, so it cannot be stolen.
    // `eig_ks_all` / `wg_ks_all` are nspin x nbands matrices -- a few kB, copied rather than aliased
    // so that they survive the KS solver overwriting `pelec` on the next step.
    this->psi_ks_all_ = ks_sol.psi;
    this->eig_ks_all = ks_sol.pelec->ekb;
    this->wg_ks_all = ks_sol.pelec->wg;
    const int start_band = this->nocc_max - *std::max_element(nocc.begin(), nocc.end());

    for (int ik = 0;ik < this->kv.get_nks();++ik)
    {
        // copy the KS orbitals in the [nocc, nvirt] window
#ifdef __MPI
        Cpxgemr2d(this->nbasis, this->nbands, &(*this->psi_ks_all_)(ik, 0, 0), 1, start_band + 1, ks_sol.pv.desc_wfc,
            &(*this->psi_ks)(ik, 0, 0), 1, 1, this->paraC_.desc, this->paraC_.blacs_ctxt);
#else
        // serial: each band is `nbasis` contiguous coefficients, so the window is a plain
        // band-by-band copy (this loop used to compute the two pointers and copy nothing,
        // leaving `psi_ks` uninitialized in every non-MPI build)
        for (int ib = 0;ib < this->nbands;++ib)
        {
            const auto* start = &(*this->psi_ks_all_)(ik, start_band + ib, 0);
            auto* to = &(*this->psi_ks)(ik, ib, 0);
            std::copy(start, start + this->nbasis, to);
        }
#endif
        // copy the KS bands in the [nocc, nvirt] window
        for (int ib = 0;ib < this->nbands;++ib) { this->eig_ks(ik, ib) = this->eig_ks_all(ik, start_band + ib); }
    }

    if (nspin == 2)
    {
        const int nupdown_now = cal_nupdown_form_occ(ks_sol.pelec->wg);
        if (!this->ks_initialized_)
        {
            this->nupdown = nupdown_now;
            reset_dim_spin2();   // shifts nocc/nvirt between the spin channels: must run exactly once
        }
        else if (nupdown_now != this->nupdown)
        {
            ModuleBase::WARNING_QUIT("ESolver_LR::refresh_from_ks_",
                "the ground-state spin population changed between ionic steps, but nocc/nvirt and"
                " the distributions built from them were fixed at the first step.");
        }
    }
    // `gint_info_` stays owned by the KS solver: its `before_scf` rebuilds and re-publishes it
    //init potential and calculate kernels using ground state charge
    init_pot(*ks_sol.pelec->charge);

#ifdef __EXX
    if (exx_kernel_list().count(xc_kernel) || (this->inp_->cal_force && gs_is_hybrid(this->inp_->dft_functional)))
    {
        if (this->exx_owned_)
        {   // Cs/Vs follow the atoms, so they are rebuilt for every geometry
            this->exx_lri->cal_exx_ions(ucell, this->inp_->out_ri_cv);
        }
        else if (ks_sol.exx_nao.exd) { this->share_exx_lri(ks_sol.exx_nao.exd->exx_ptr); }
        else if (ks_sol.exx_nao.exc) { this->share_exx_lri(ks_sol.exx_nao.exc->exx_ptr); }
    }
#endif
    // the grid-integration tables hang off a static pointer that the ground-state solver
    // re-publishes in its `before_scf`; make sure it names the object we integrate on
    ModuleGint::Gint::set_gint_info(this->ks_->gint_info_.get());

    // the Z-vector window: after `reset_dim_spin2`, so nocc/nvirt/openshell are final
#ifdef __MPI
    this->fill_z_window_(ks_sol.pv.desc_wfc);
#else
    this->fill_z_window_(nullptr);
#endif
}


template <typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::initialize_from_unitcell_(UnitCell& ucell, const Input_para& inp)
{
    ModuleBase::TITLE("ESolver_LR", "ESolver_LR(from scratch)");
    this->bind_ground_state_aliases_();
    // xc kernel
    this->xc_kernel = LR_Util::tolower(inp.xc_kernel);

    // necessary steps in ESolver_FP
    ESolver_FP::before_all_runners(ucell, inp);
    LR_Util::check_force_pp(ucell, inp.cal_force, inp.calculation);
    this->pelec = new elecstate::ElecStateLCAO<T>();

    // necessary steps in ESolver_KS::before_all_runners : symmetry and k-points
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        const int cal_symm_repr[2] = {this->inp_->cal_symm_repr[0], this->inp_->cal_symm_repr[1]};
        ucell.symm.analy_sys(ucell.lat, ucell.atoms,
                             ucell.nat, ucell.ntype, ucell.iat2it, ucell.iat2ia, ucell.itia2iat,
                             this->ofs_running_,
                             this->inp_->symmetry_prec, this->inp_->nspin, this->inp_->calculation, cal_symm_repr);
        ModuleBase::GlobalFunc::DONE(this->ofs_running_, "SYMMETRY");
    }
    const bool use_ibz = false;
    const bool gamma_only_local = PARAM.globalv.gamma_only_local;
    const double kspacing[3] = {this->inp_->kspacing[0], this->inp_->kspacing[1], this->inp_->kspacing[2]};
    const double koffset[3] = {this->inp_->koffset[0], this->inp_->koffset[1], this->inp_->koffset[2]};
    this->kv.set(ucell, ucell.symm, this->inp_->kpoint_file, this->inp_->nspin, ucell.G, ucell.latvec, this->ofs_running_, this->ofs_warning_, use_ibz, this->out_dir, gamma_only_local, kspacing, this->inp_->kmesh_type, koffset);
    ModuleBase::GlobalFunc::DONE(this->ofs_running_, "INIT K-POINTS");
    ModuleIO::print_parameters(ucell, this->kv, inp);

    this->parameter_check();

    /// read orbitals and build the interpolation table
    two_center_bundle_own_.build_orb(ucell.ntype, ucell.orbital_fn.data(), inp.orbital_dir);

    LCAO_Orbitals orb;
    two_center_bundle_own_.to_LCAO_Orbitals(orb, inp.lcao_ecut, inp.lcao_dk, inp.lcao_dr, inp.lcao_rmax,
                                        inp.out_element_info, inp.cal_force);
    orb_cutoff_ = orb.cutoffs();
    if (LR_Util::tolower(this->inp_->abs_gauge) == "velocity")
    {
        setup_2center_table(this->two_center_bundle_own_, orb, ucell);
    }

    this->set_dimension();
    //  setup 2d-block distribution for AO-matrix and KS wfc
    LR_Util::setup_2d_division(this->paraMat_, 1, this->nbasis, this->nbasis);
    this->set_parallel_orbitals_band(this->paraMat_, this->nbands);
    if (this->inp_->cal_force)
    {
        LR_Util::setup_2d_division(this->paraMat_all_, 1, this->nbasis, this->nbasis);
        this->set_parallel_orbitals_band(this->paraMat_all_, this->inp_->nbands);
    }

    // read the ground state info
    // now ModuleIO::read_wfc_nao needs `Parallel_Orbitals` and can only read all the bands
    // it need improvement to read only the bands needed
    this->psi_ks.reset(new psi::Psi<T>(this->kv.get_nks(),
                                       this->paraMat_.ncol_bands,
                                       this->paraMat_.get_row_size(),
                                       this->kv.ngk,
                                       true));
    this->read_ks_wfc();


    if (nspin == 2)
    {
        // Complete file populations were read before selecting the occupied window.
        reset_dim_spin2();
    }
    // the Z-vector window: after `reset_dim_spin2`, so nocc/nvirt/openshell are final
#ifdef __MPI
    this->fill_z_window_(paraMat_all_.desc_wfc);
#else
    this->fill_z_window_(nullptr);
#endif

    LR_Util::setup_2d_division(this->paraC_, 1, this->nbasis, this->nbands
#ifdef __MPI
        , paraMat_.blacs_ctxt
#endif
    );

    // clear ks info, new elecstate for excition
    delete this->pelec;   // the ElecStateLCAO allocated above, only needed while reading the KS data
    this->pelec = new elecstate::ElecState();

    // read the ground state charge density and calculate xc kernel
    pgrid().init(this->pw_rho->nx,
        this->pw_rho->ny,
        this->pw_rho->nz,
        this->pw_rho->nplane,
        this->pw_rho->nrxx,
        pw_big->nbz,
        pw_big->bz,
        GlobalV::NPROC);
    Charge chg_gs;
    if (this->inp_->ri_hartree_benchmark == "none") { this->read_ks_chg(chg_gs); }
    this->init_pot(chg_gs);

    // search adjacent atoms and init Gint
    double search_radius = -1.0;
    search_radius = atom_arrange::set_sr_NL(this->ofs_running_,
        this->inp_->out_level,
        orb.get_rcutmax_Phi(),
        ucell.infoNL->get_rcutmax_Beta(),
        PARAM.globalv.gamma_only_local);
    atom_arrange::search(PARAM.globalv.search_pbc,
                         this->ofs_running_,
                         this->gd(),
                         *this->ucell_,
                         search_radius,
                         this->inp_->test_atom_input);
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
        orb.Phi,
        ucell,
        this->gd(),
        this->inp_->nspin,
        PARAM.globalv.gamma_only_local,
        PARAM.globalv.domag,
        this->inp_->device == "gpu",
        this->inp_->nstream));
    ModuleGint::Gint::set_gint_info(gint_info_.get());
    // if EXX from scratch, init 2-center integral and calculate Cs, Vs
    // when:
    // 1. EXX xc_kernel
    // 2. cal_force with ground state with EXX functional
#ifdef __EXX
    if (((exx_kernel_list().count(xc_kernel)) && this->inp_->lr_solver != "spectrum")
        || (this->inp_->cal_force && gs_is_hybrid(this->inp_->dft_functional)))
    {
        warn_if_kernel_differs_from_gs(xc_kernel, LR_Util::tolower(this->inp_->dft_functional), this->ofs_running_);
        // `input_conv` already filled `info_ri.coulomb_param` from INPUT.
        exx_info.sync_from_global();
        // populate ABFs/JLE file lists from UnitCell; keep in sync with Exx_NAO::init
        exx_info.info_ri.files_abfs = ucell.abfs_orbital_files;
        exx_info.info_opt_abfs.files_abfs = ucell.abfs_orbital_files;
        exx_info.info_opt_abfs.files_jles = ucell.jle_orbital_files;
        this->exx_lri = std::make_shared<Exx_LRI<T>>(exx_info.info_ri);
        this->exx_lri->init(MPI_COMM_WORLD, ucell,this->kv, orb);
        this->exx_lri->cal_exx_ions(ucell,this->inp_->out_ri_cv);
    }
    // else
#endif
        // ModuleBase::Ylm::set_coefficients() is deprecated
}

template class ModuleESolver::ESolver_LR<double, double>;
template class ModuleESolver::ESolver_LR<std::complex<double>, double>;
