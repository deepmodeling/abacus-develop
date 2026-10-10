#include "source_esolver/esolver_lr_lcao_tddft.h"
#include "source_lcao/module_lr/zeq_solver.h"
#include "source_lcao/module_lr/cal_edm.h"
#include "source_lcao/module_lr/lr_force.h"
#include "source_lcao/module_lr/gradient_inputs.h"
#include "source_lcao/module_lr/gradient_output.h"
#include "source_lcao/module_lr/lr_amp.h"
#include "source_lcao/module_lr/grad_degen.h"
#include "source_base/parallel_reduce.h"
#include <algorithm>
#include <complex>
#include <numeric>
#include "source_estate/module_dm/dm_from_psi.h"
#include "source_io/module_output/output_log.h"

using namespace LR;


///========================= excited-state geometry relaxation =========================

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::init_pot_groundstate(const Charge& chg_gs)
{
    ModuleBase::TITLE("ESolver_LR", "init_pot_gs");
    std::vector<std::string> pot_register;
    if (this->inp_->vl_in_h)
    {
        if (!this->ks_)
        {   // on the `ks-lr` path `sfac()`/`vloc()` alias the ground-state solver's, which
            // `ESolver_FP::before_scf` already refreshed for the current geometry
            //! 11) calculate the structure factor
            // the has_float_data flag below is read the same way core ABACUS reads it at the
            // analogous call site (source_esolver/esolver_fp.cpp's `this->sf.setup(...)`), so this
            // mirrors the established convention rather than introducing a new global dependency.
            this->sfac().setup(&(*this->ucell_), pgrid(), this->pw_rhod, PARAM.globalv.has_float_data);
            this->vloc().init_vloc((*this->ucell_), this->pw_rho);
        }
        pot_register.push_back("local");
    }
    if(this->inp_->vh_in_h)
    {
        pot_register.push_back("hartree");
    }
    pot_register.push_back("xc");
    
    // initialize the ground state potential
    this->pot_gs = LR_Util::make_unique<elecstate::Potential>(this->pw_rhod, this->pw_rho,
        &(*this->ucell_), &this->vloc().vloc, &this->sfac(), &this->solvent,
       &this->etxc_gs, &this->vtxc_gs);
    this->pot_gs.get()->pot_register(pot_register);
    XC_Functional::set_xc_type((*this->ucell_).atoms[0].ncpp.xc_func);    // set XC type of the ground state
    this->pot_gs->init_pot(&chg_gs);  // call update_from_charge inside
    if (LR_Util::has_local_xc(this->xc_kernel))
    {
        XC_Functional::set_xc_type(this->xc_kernel);    // recover the excited state xc kernel type
    }
    if (this->inp_->test_force)
    {
        this->pot_gs_hartree = LR_Util::make_unique<elecstate::Potential>(this->pw_rhod, this->pw_rho,
            &(*this->ucell_), &this->vloc().vloc, &this->sfac(), &this->solvent,
            &this->etxc_gs, &this->vtxc_gs);
        this->pot_gs_hartree->pot_register({ "hartree" });
        this->pot_gs_hartree->init_pot(&chg_gs);  // call update_from_charge inside
    }
}

template<typename T, typename TR>
ct::Tensor ModuleESolver::ESolver_LR<T, TR>::pad_X_to_z_(const int ispin, const int istate_begin, const int nst) const
{
    return LR::pad_amplitudes<T>(this->X, this->openshell, this->paraX_, this->paraX_z_,
        this->nocc, this->nvirt, this->nk, this->nloc_per_state, this->nloc_per_state_z_,
        ispin, istate_begin, nst);
}

template<typename T, typename TR>
ct::Tensor ModuleESolver::ESolver_LR<T, TR>::solve_zvector_eqation(const int ispin, const int nst, const ct::Tensor& Xz)
{
    ModuleBase::TITLE("ESolver_LR", "cal_force");
    ModuleBase::timer::start("ESolver_LR", "solve_zvector_eqation");
    // `Z_vector_equation` treats X and Z as `nst` independent blocks of `nloc_per_state_z_`,
    // and the blocks need not be eigenvectors of the Casida equation in the order the
    // diagonalizer returned them -- `cal_grad_matrix_degenerate` feeds it linear combinations.
    // X arrives already widened into the Z window (`Xz`), which spans every virtual band the
    // ground state produced rather than the `nvirt` window X was solved in -- the Brillouin
    // condition the Z-vector enforces holds in EVERY occupied-virtual rotation. The padded
    // entries are zero, so every X-derived quantity (D^X, T) is unchanged; only the space Z is
    // solved in grows. Both the closed- and the open-shell path go through here.
    ct::Tensor Z = LR_Util::newTensor<T>({ nst, this->nloc_per_state_z_ });
    // construct and solve the Z-vector equation
    const LR::GradientInputs<T> inputs = this->gradient_inputs_();
    const T* const x_data = Xz.template data<T>();
    T* const z_data = Z.template data<T>();
    std::weak_ptr<PotHxcLR> pot_weak = inputs.pot[ispin];
    Z_vector_equation(inputs, x_data, z_data, nst, pot_weak,
        this->spin_types[ispin], this->in_dir, this->openshell, this->inp_->lr_grad_solver);
    ModuleBase::timer::end("ESolver_LR", "solve_zvector_eqation");
    return Z;
}

template<typename T, typename TR>
std::vector<ModuleBase::matrix> ModuleESolver::ESolver_LR<T, TR>::cal_force(const int ispin, const int istate_only)
{
    if (this->inp_->test_force && ispin == 0) { this->test_force(); }
    if (this->openshell) { return this->cal_force_openshell(istate_only); }

    const int ist_begin = (istate_only < 0) ? 0 : istate_only;
    const int nst = (istate_only < 0) ? this->nstates : 1;
    // The whole closed-shell gradient runs in the Z window (every virtual band the ground state
    // produced), not the `nvirt` window X was solved in: only there does the Z-vector enforce
    // the Brillouin condition in every occupied-virtual rotation. X is zero-padded into it, so
    // D^X and T come out bit-identical -- and Omega is untouched, so an existing finite-difference
    // reference stays valid. See `fill_z_window_`.
    const ct::Tensor Xz = this->pad_X_to_z_(ispin, ist_begin, nst);
    std::vector<double> omega(nst);
    for (int i = 0; i < nst; ++i)
    {
        omega[i] = this->pelec->ekb.c[ispin * this->nstates + ist_begin + i];
    }
    return this->cal_force_Xz(ispin, Xz, omega, ist_begin);
}

template<typename T, typename TR>
LR::GradientInputs<T> ModuleESolver::ESolver_LR<T, TR>::gradient_inputs_() const
{
    return { *this->ucell_, this->kv, this->gd(), this->orb_cutoff_, this->paraMat_,
        this->paraC_z_, this->paraX_z_, *this->psi_ks_z_, this->eig_ks_z_, this->nocc,
        this->nvirt_z_, this->nspin, this->nk, this->nbasis, this->nloc_per_state_z_,
        this->xc_kernel, this->inp_->dft_functional, this->inp_->ks_solver,
        this->inp_->test_force, this->excited_relax_, this->out_dir, this->my_rank_,
        this->spin_types, this->pot, this->pot_hxc_gs, this->ofs_running_
#ifdef __EXX
        , this->exx_lri, this->exx_info.info_global.hybrid_alpha
#endif
    };
}

template<typename T, typename TR>
std::vector<ModuleBase::matrix> ModuleESolver::ESolver_LR<T, TR>::cal_force_Xz(const int ispin,
    const ct::Tensor& Xz, const std::vector<double>& omega, const int label_begin)
{
    ModuleBase::timer::start("ESolver_LR", "cal_force_Xz");
    const int nst = static_cast<int>(omega.size());
    const int channel = ispin;
    const ct::Tensor Z = this->solve_zvector_eqation(channel, nst, Xz);
    const auto dm_gs = this->cal_dm_gs();
    const LR::GradientInputs<T> inputs = this->gradient_inputs_();
    LR_Force<T> force_terms(*this->ucell_, this->kv.kvec_d, this->paraMat_,
        *this->pw_rhod, *this->pw_rho, this->vloc(), this->sfac(), this->gd(), this->tcb()
#ifdef __EXX
        , inputs.exx_lri, inputs.hybrid_alpha
#endif
    );
    const auto forces = LR::evaluate_closed_shell_force(inputs, force_terms, dm_gs, Xz, Z, omega, label_begin, ispin);
    ModuleBase::timer::end("ESolver_LR", "cal_force_Xz");
    return forces;
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::cal_force_and_grad_matrix_(const int ispin, std::ofstream& ofs)
{
    const std::vector<ModuleBase::matrix> forces = this->cal_force(ispin);
    if (this->inp_->lr_degen_thr <= 0.0) { return; }
    // Rebuilt from scratch for this channel: a relaxation calls this once per ionic step, and a
    // stale multiplet from the previous geometry must not survive into the next one.
    this->multiplet_lvc_.erase(
        std::remove_if(this->multiplet_lvc_.begin(), this->multiplet_lvc_.end(),
            [ispin](const MultipletLVC& m) { return m.ispin == ispin; }),
        this->multiplet_lvc_.end());
    // The per-state gradients above are the diagonal of the degenerate-subspace gradient matrix in
    // whatever basis the eigensolver returned, so inside a multiplet only their trace means
    // anything. `lr_degen_thr` asks for the off-diagonal part too, which is the rest of the
    // first-order information.
    const int ekb_off = this->openshell ? 0 : ispin * this->nstates;
    std::vector<double> omega(this->nstates);
    for (int ist = 0; ist < this->nstates; ++ist) { omega[ist] = this->pelec->ekb.c[ekb_off + ist]; }
    const std::vector<std::vector<int>> groups
        = LR::group_degenerate_states(omega, this->inp_->lr_degen_thr);
    for (const std::vector<int>& group : groups)
    {
        if (group.size() < 2) { continue; }
        std::vector<ModuleBase::matrix> diag;
        for (const int ist : group) { diag.push_back(forces[ist]); }
        MultipletLVC lvc;
        lvc.ispin = ispin;
        lvc.states = group;
        double lo = omega[group.front()];
        double hi = omega[group.front()];
        double sum = 0.0;
        for (const int ist : group)
        {
            lo = std::min(lo, omega[ist]);
            hi = std::max(hi, omega[ist]);
            sum += omega[ist];
        }
        lvc.omega0 = sum / static_cast<double>(group.size());
        lvc.omega_spread = hi - lo;
        lvc.g = this->cal_grad_matrix_degenerate(ispin, group, diag, ofs);
        this->multiplet_lvc_.push_back(lvc);
    }
}

template<typename T, typename TR>
std::vector<std::vector<ModuleBase::matrix>>
ModuleESolver::ESolver_LR<T, TR>::cal_grad_matrix_degenerate(const int ispin,
    const std::vector<int>& group, const std::vector<ModuleBase::matrix>& diag, std::ofstream& ofs)
{
    ModuleBase::TITLE("ESolver_LR", "cal_grad_matrix_degenerate");
    ModuleBase::timer::start("ESolver_LR", "cal_grad_matrix_degenerate");
    const int d = static_cast<int>(group.size());
    assert(d >= 2);
    assert(static_cast<int>(diag.size()) == d);
    const std::vector<std::pair<int, int>> pairs = LR::degenerate_pairs(d);
    const int nloc_g = this->nloc_per_state_z_;
    const int ekb_off = this->openshell ? 0 : ispin * this->nstates;

    // 1. the multiplet's members, each widened into the Z window
    const ct::Tensor Xz = this->pad_group_to_z_(ispin, group);
    std::vector<double> omega_member(d);
    for (int k = 0; k < d; ++k) { omega_member[k] = this->pelec->ekb.c[ekb_off + group[k]]; }

    // 2. the two preconditions of route (A2), reported rather than enforced: the combinations
    //    $X_\pm$ are normalized eigenvectors only if the members are orthonormal, and the gradient
    //    is a single quadratic form only if they share one $\Omega$. An accidental near-degeneracy
    //    passes the grouping threshold but fails the second, and its `omega_spread` says so.
    double max_ovlp_err = 0.0;
    for (int k = 0; k < d; ++k)
    {
        for (int l = k; l < d; ++l)
        {
            const T* const xk = Xz.template data<T>() + static_cast<size_t>(k) * nloc_g;
            const T* const xl = Xz.template data<T>() + static_cast<size_t>(l) * nloc_g;
            T loc = static_cast<T>(0);
            for (int i = 0; i < nloc_g; ++i) { loc += LR_Util::get_conj(xk[i]) * xl[i]; }
            Parallel_Reduce::reduce_all(loc);
            const double ref = (k == l) ? 1.0 : 0.0;
            max_ovlp_err = std::max(max_ovlp_err, std::abs(loc - ref));
        }
    }
    const double omega_spread = *std::max_element(omega_member.begin(), omega_member.end())
        - *std::min_element(omega_member.begin(), omega_member.end());
    // one $\Omega$ for every combination: the members share it up to `omega_spread`, and the mean
    // is the neutral choice for a vector that belongs to no single member
    const double omega_mean
        = std::accumulate(omega_member.begin(), omega_member.end(), 0.0) / static_cast<double>(d);
    ofs << " Degenerate multiplet of " << d << " states, Omega = " << omega_mean
        << " Ry, spread = " << omega_spread << " Ry, max |<X_k|X_l> - delta_kl| = " << max_ovlp_err
        << std::endl;
    // Dimensionless diagnostic tolerance for the orthonormality needed by X_plus.
    // This warning does not change the multiplet or the computed gradient.
    const double orthonormality_warn_tol = 1e-6;
    if (max_ovlp_err > orthonormality_warn_tol)
    {
        ofs << " WARNING: this multiplet's eigenvectors are not orthonormal to"
            << " tolerance " << orthonormality_warn_tol
            << " (dimensionless overlap error); (X_k+X_l)/sqrt(2) may not be normalized"
            << " and the assembled gradient matrix may consequently be inaccurate."
            << std::endl;
    }

    // 3. the $d(d-1)/2$ combinations $X_+=(X_k+X_l)/\sqrt2$, all in one Z-vector solve
    const int npair = static_cast<int>(pairs.size());
    ct::Tensor Xp = LR_Util::newTensor<T>({ npair, nloc_g });
    Xp.zero();
    for (int ip = 0; ip < npair; ++ip)
    {
        LR::combine_normalized(
            Xz.template data<T>() + static_cast<size_t>(pairs[ip].first) * nloc_g,
            Xz.template data<T>() + static_cast<size_t>(pairs[ip].second) * nloc_g,
            static_cast<size_t>(nloc_g),
            Xp.template data<T>() + static_cast<size_t>(ip) * nloc_g);
    }
    const std::vector<double> omega_pair(npair, omega_mean);
    const std::vector<ModuleBase::matrix> plus = this->openshell
        ? this->cal_force_openshell_Xz(Xp, omega_pair, group[0])
        : this->cal_force_Xz(ispin, Xp, omega_pair, group[0]);

    // 4. $G_{kl}=\mathcal F[X_+]-\tfrac12(G_{kk}+G_{ll})$
    const std::vector<std::vector<ModuleBase::matrix>> g
        = LR::assemble_grad_matrix(diag, plus, pairs);
    print_grad_matrix(g, group, (*this->ucell_), ofs);
    ModuleBase::timer::end("ESolver_LR", "cal_grad_matrix_degenerate");
    return g;
}

template<typename T, typename TR>
std::vector<ModuleBase::matrix> ModuleESolver::ESolver_LR<T, TR>::cal_force_openshell(const int istate_only)
{
    ModuleBase::TITLE("ESolver_LR", "cal_force_openshell");

    // Open shell: there is a single eigenproblem whose vector is the concatenation
    // [up-block | down-block], and every density matrix has two independent channels.
    // The spin-orbital formulas apply verbatim -- unlike the closed-shell singlet/triplet
    // algorithm, X here is normalized over BOTH channels, so it carries no implicit sqrt(2)
    // and none of the collapsed 2/4 factors are needed.
    const int ist_begin_ = (istate_only < 0) ? 0 : istate_only;
    const int nst_ = (istate_only < 0) ? this->nstates : 1;
    // Like the closed-shell path, the whole gradient runs in the Z window: X is zero-padded
    // into it (both channels, each re-based -- see `pad_X_to_z_`), so D^X and T are unchanged
    // while the Z-vector gets every occupied-virtual rotation the AO basis supports.
    const ct::Tensor Xz = this->pad_X_to_z_(0, ist_begin_, nst_);
    std::vector<double> omega_(nst_);
    for (int i = 0; i < nst_; ++i) { omega_[i] = this->pelec->ekb.c[ist_begin_ + i]; }
    return this->cal_force_openshell_Xz(Xz, omega_, ist_begin_);
}

template<typename T, typename TR>
std::vector<ModuleBase::matrix> ModuleESolver::ESolver_LR<T, TR>::cal_force_openshell_Xz(
    const ct::Tensor& Xz, const std::vector<double>& omega, const int label_begin)
{
    ModuleBase::timer::start("ESolver_LR", "cal_force_openshell_Xz");
    const int nst = static_cast<int>(omega.size());
    const int channel = 0;
    const ct::Tensor Z = this->solve_zvector_eqation(channel, nst, Xz);
    const auto dm_gs = this->cal_dm_gs();
    const LR::GradientInputs<T> inputs = this->gradient_inputs_();
    LR_Force<T> force_terms(*this->ucell_, this->kv.kvec_d, this->paraMat_,
        *this->pw_rhod, *this->pw_rho, this->vloc(), this->sfac(), this->gd(), this->tcb()
#ifdef __EXX
        , inputs.exx_lri, inputs.hybrid_alpha
#endif
    );
    const auto forces = LR::evaluate_open_shell_force(inputs, force_terms, dm_gs, Xz, Z, omega, label_begin);
    ModuleBase::timer::end("ESolver_LR", "cal_force_openshell_Xz");
    return forces;
}

template<typename T, typename TR>
module_dm::DensityMatrix<T, double> ModuleESolver::ESolver_LR<T, TR>::cal_dm_gs()
{
    module_dm::DensityMatrix<T, double> dm_gs(&this->paraMat_, this->nspin, this->kv.kvec_d, this->nk);
    // `psi_ks_all_` is distributed per `this->ks_->pv` on the ks-lr path (it is aliased straight
    // from the ground-state solver's own `psi`), but per `this->paraMat_all_` on the
    // read-from-file path (where it was allocated with that descriptor's own sizes). Using the
    // wrong one here reads `psi_ks_all_`'s local block with the wrong block size/process-grid
    // assumption -- invisible at nprocs=1, where every descriptor degenerates to one block, but
    // silently wrong at nprocs>1. `refresh_from_ks_`/`fill_z_window_` already make this same
    // distinction for their `Cpxgemr2d` source descriptor; this mirrors it.
    const Parallel_Orbitals* const pv_all = this->ks_ ? &this->ks_->pv : &this->paraMat_all_;
    module_dm::dm_from_psi(pv_all, this->wg_ks_all, *this->psi_ks_all_, dm_gs);   // nbands is important here
    LR_Util::initialize_DMR(dm_gs, this->paraMat_, (*this->ucell_), this->gd(), this->orb_cutoff_);   // nbands is not important here
    dm_gs.cal_dmr(-1);
    return dm_gs;
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::test_force()
{
    LR_Force<T> lr_force((*this->ucell_), this->kv.kvec_d, this->paraMat_, *this->pw_rhod, *this->pw_rho,
        this->vloc(), this->sfac(), this->gd(), this->tcb()
#ifdef __EXX
        , std::weak_ptr<Exx_LRI<T>>(this->exx_lri), this->exx_info.info_global.hybrid_alpha
#endif
    );

    const module_dm::DensityMatrix<T, double>& dm_gs = this->cal_dm_gs();
    // LR_Util::print_DMR(dm_gs, "DM(R) of ground state");
    ///========================== test 1: reproduce the force of ground state =========================
    // energy density matrix of the ground state
    module_dm::DensityMatrix<T, double> edm_gs(&this->paraMat_, this->nspin, this->kv.kvec_d, this->nk);   //DX
    ModuleBase::matrix wg_ekb_ks_all(nspin, this->inp_->nbands);
    std::transform(this->wg_ks_all.c, this->wg_ks_all.c + nspin * this->inp_->nbands,
        this->eig_ks_all.c, wg_ekb_ks_all.c, std::multiplies<double>());
    // see `cal_dm_gs()`: `psi_ks_all_` is distributed per `this->ks_->pv` on the ks-lr path,
    // not `this->paraMat_all_`.
    const Parallel_Orbitals* const pv_all_edm = this->ks_ ? &this->ks_->pv : &this->paraMat_all_;
    module_dm::dm_from_psi(pv_all_edm, wg_ekb_ks_all, *this->psi_ks_all_, edm_gs);
    LR_Util::initialize_DMR(edm_gs, this->paraMat_, (*this->ucell_), this->gd(), this->orb_cutoff_);
    edm_gs.cal_dmr(-1);
    // ground-state force
    ModuleBase::matrix force_gs = lr_force.reproduce_force_gs(kv, dm_gs, edm_gs);
    ModuleIO::print_force(this->ofs_running_, (*this->ucell_), "Ground State FORCE (eV/Angstrom)", force_gs, false);
    /// ======================================= END test 1 =========================================
    ///========================== test 2: reproduce the DX Hartree term =========================
    ModuleBase::matrix f_hxc_potgs = lr_force.reproduce_force_gs_loc(dm_gs, *this->pot_gs_hartree);
    ModuleIO::print_force(this->ofs_running_, (*this->ucell_), "GS Hartree force calculated by 'cal_pulay_fs' from potential (eV/Angstrom)", f_hxc_potgs, false);
    ModuleBase::matrix f_hxc_potlr = lr_force.cal_force_hxc_dmtrans(dm_gs, *this->pot[0]);
    ModuleIO::print_force(this->ofs_running_, (*this->ucell_), "GS Hxc force calculated by 'LR_Force' from kernel (eV/Angstrom)", f_hxc_potlr, false);
    // `cal_force_hxc_dmtrans` now includes the Pulay -> Pulay+Hellmann-Feynman factor 2 itself,
    // so this must match the ground-state Hartree force directly (dm_gs is already symmetric).
    /// ======================================= END test 2 =========================================

}

template class ModuleESolver::ESolver_LR<double, double>;
template class ModuleESolver::ESolver_LR<std::complex<double>, double>;
