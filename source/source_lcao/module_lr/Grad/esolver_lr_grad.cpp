#include "source_esolver/esolver_lr_lcao_tddft.h"
#include "source_lcao/module_lr/Grad/multipliers/zeq_solver.h"
#include "source_lcao/module_lr/Grad/multipliers/cal_edm_from_multipliers.h"
#include "source_lcao/module_lr/Grad/force/lr_force.h"
#include "source_lcao/module_lr/Grad/degenerate/grad_matrix_degenerate.h"
#include "source_base/parallel_reduce.h"
#include <algorithm>
#include <complex>
#include <numeric>
#include "source_estate/module_dm/cal_dm_psi.h"
#include "source_io/module_output/output_log.h"

using namespace LR;


template <typename Tstream>
inline void print_force(const std::vector<ModuleBase::matrix>& force, Tstream& ofs, const int istate_begin = 0)
{
    const int nstate = force.size();
    ofs << "Forces (-gradients) of each excited state: (eV/Angstrom)" << std::endl;
    ofs << std::setprecision(6) << std::setw(6) << "state" << std::setw(6) << "atom"
        << std::setw(15) << "x" << std::setw(15) << "y" << std::setw(15) << "z" << std::endl;
    const double fac = ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A;
    for (int i = 0;i < nstate;++i)
    {
        for (int iat = 0;iat < force[i].nr;++iat)
        {
            std::string istate = iat == 0 ? std::to_string(istate_begin + i) : " ";
            ofs << std::setw(6) << istate << std::setw(6) << iat << std::setw(6) << "force";
            for (int ixyz = 0;ixyz < 3;++ixyz) { ofs << std::setw(15) << force[i](iat, ixyz) * fac; }
            ofs << std::endl;
        }
    }
}

/// @brief The mean of a set of force matrices.
///
/// Used for the multiplet average $\operatorname{Tr}G/d$, the one smooth, basis-independent $3N$
/// vector field a degenerate multiplet has. Factored out because three callers need it from
/// different inputs: the printout and the stored LVC hold the whole gradient matrix, while the
/// relaxation only ever computes its diagonal.
inline ModuleBase::matrix average_forces(const std::vector<ModuleBase::matrix>& f)
{
    assert(!f.empty());
    ModuleBase::matrix avg(f[0].nr, f[0].nc);
    for (size_t k = 0; k < f.size(); ++k) { avg += f[k]; }
    avg *= 1.0 / static_cast<double>(f.size());
    return avg;
}

/// @brief Print the degenerate-subspace gradient matrix $G_{kl}$, one $d\times d$ block per
///        nuclear coordinate, plus the multiplet average on its diagonal.
///
/// The individual diagonal entries are basis-dependent: only the eigenvalues of
/// $M(u)=\sum_a u_aG^{(a)}$ are branch slopes, and only $\operatorname{Tr}G$ is invariant. The
/// average $\operatorname{Tr}G/d$ is printed because it IS a smooth, basis-independent $3N$ vector
/// field -- the one a symmetry-constrained relaxation can follow.
inline void print_grad_matrix(const std::vector<std::vector<ModuleBase::matrix>>& g,
    const std::vector<int>& group, const UnitCell& ucell, std::ofstream& ofs)
{
    const double fac = ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A;
    const int d = static_cast<int>(g.size());
    const int nat = g[0][0].nr;
    ofs << std::endl << " DEGENERATE-SUBSPACE GRADIENT MATRIX G_kl (eV/Angstrom), states";
    for (int k = 0; k < d; ++k) { ofs << " " << group[k]; }
    ofs << std::endl
        << " Branch slopes along a displacement u are the EIGENVALUES of sum_a u_a G^(a); the"
        << std::endl
        << " diagonal entries alone are basis-dependent and only their trace is invariant."
        << std::endl;
    ofs << " For each coordinate the matrix is followed by its eigenvalues, which ARE the branch"
        << std::endl
        << " slopes for displacing that one atom along that one axis. They must NOT be combined"
        << std::endl
        << " across axes: G^(x), G^(y), G^(z) do not commute in general, so eigenvalues are not"
        << std::endl
        << " additive and picking one per axis describes no adiabatic state at all." << std::endl;
    ofs << std::setprecision(6);
    for (int iat = 0; iat < nat; ++iat)
    {
        for (int ixyz = 0; ixyz < 3; ++ixyz)
        {
            ofs << " atom " << std::setw(5) << iat << " dir " << std::setw(2) << ixyz << std::endl;
            std::vector<double> block(static_cast<size_t>(d) * d);
            for (int k = 0; k < d; ++k)
            {
                ofs << "     ";
                for (int l = 0; l < d; ++l)
                {
                    const double v = g[k][l](iat, ixyz) * fac;
                    ofs << std::setw(15) << v;
                    block[static_cast<size_t>(k) * d + l] = v;
                }
                ofs << std::endl;
            }
            // `diag_lapack` overwrites its input with the eigenvectors, hence the scratch copy
            std::vector<double> eig(d, 0.0);
            LR_Util::diag_lapack(d, block.data(), eig.data());
            ofs << "       eig";
            for (int k = 0; k < d; ++k) { ofs << std::setw(15) << eig[k]; }
            ofs << std::endl;
        }
    }
    std::vector<ModuleBase::matrix> diag;
    for (int k = 0; k < d; ++k) { diag.push_back(g[k][k]); }
    ModuleIO::print_force(ofs, ucell, "MULTIPLET-AVERAGE FORCE Tr(G)/d (eV/Angstrom)",
        average_forces(diag), false);
}

// check C_uaC_va-C_uiC_vi of lumo-homo, nocc=1,  nk=1
template<typename T>
inline void test_dm_diff_H2(const T* dm, const psi::Psi<T>& c, const int nbasis)
{
    std::cout << "difference dm cal: " << std::endl;
    LR_Util::print_value(dm, nbasis, nbasis);
    std::cout << "difference dm ref: " << std::endl;
    for (int i = 0;i < nbasis;++i)
    {
        for (int j = 0;j < nbasis;++j)
        {
            std::cout << c(0, 1, i) * c(0, 1, j) - c(0, 0, i) * c(0, 0, j) << " ";
        }
        std::cout << std::endl;
    }
}

// check e_aC_uaC_va-e_iC_uiC_vi of lumo-homo, nocc=1,  nk=1
template<typename T>
inline void test_edm_H2(const T* const edm, const double* const eig_ks, const psi::Psi<T>& c, const int nbasis)
{
    std::cout << "edm cal: " << std::endl;
    LR_Util::print_value(edm, nbasis, nbasis);
    std::cout << "edm ref: " << std::endl;
    for (int i = 0;i < nbasis;++i)
    {
        for (int j = 0;j < nbasis;++j)
        {
            std::cout << eig_ks[1] * c(0, 1, i) * c(0, 1, j) - eig_ks[0] * c(0, 0, i) * c(0, 0, j) << " ";
        }
        std::cout << std::endl;
    }
}

///========================= excited-state geometry relaxation =========================

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::setup_relax_target_()
{
    this->excited_relax_ = (PARAM.inp.calculation == "relax");
    if (!this->excited_relax_) { return; }

    // The Z-vector equation has no complex solver (see Grad/multipliers/zeq_solver.hpp), so an
    // excited-state gradient only exists at gamma. Fail here rather than after the SCF.
    if (!std::is_same<T, double>::value)
    {
        ModuleBase::WARNING_QUIT("ESolver_LR",
            "excited-state relaxation currently requires gamma-only sampling: the complex "
            "Z-vector solver is not implemented.");
    }
    if (this->inp_->lr_target_state >= this->nstates)
    {
        ModuleBase::WARNING_QUIT("ESolver_LR",
            "lr_target_state is beyond the states actually solved (lr_nstates <= 0 expands to "
            "all particle-hole pairs, which may be fewer than requested).");
    }

    // `openshell` is only settled after the ground-state occupations have been read, which is why
    // this cannot live in the input-file checks
    const std::string& spin = this->inp_->lr_target_spin;
    if (this->openshell)
    {
        // An open-shell calculation solves one spin-conserving channel, so there is nothing to
        // choose: whatever lr_target_spin says, this is the state that gets relaxed. Only an
        // explicit `triplet` is worth mentioning -- `singlet` is the default and expresses no
        // intent, and `updown` is already the right name for this channel.
        if (spin == "triplet")
        {
            GlobalV::ofs_running << " WARNING: lr_target_spin=triplet is ignored. This is an"
                " open-shell calculation with a single spin-conserving channel (updown), which is"
                " what the relaxation will follow." << std::endl;
        }
        this->target_is_ = 0;
    }
    else
    {
        // Closed shell is the opposite case: singlet and triplet are genuinely different states
        // with different gradients, so `updown` here is ambiguous rather than redundant -- it
        // usually means lr_unrestricted was meant to be set.
        if (spin == "updown")
        {
            ModuleBase::WARNING_QUIT("ESolver_LR",
                "lr_target_spin=updown, but this is a closed-shell calculation, where singlet and "
                "triplet are separate states with separate gradients. Pick one of them, or set "
                "lr_unrestricted to run spin-unrestricted.");
        }
        this->target_is_ = (spin == "triplet") ? 1 : 0;
        if (this->target_is_ >= this->nspin)
        {
            ModuleBase::WARNING_QUIT("ESolver_LR",
                "lr_target_spin=triplet requires nspin=2: the triplet channel is not built here.");
        }
    }
    this->force_gs_.create(this->ucell_->nat, 3);
    this->lr_force_.create(this->ucell_->nat, 3);
    GlobalV::ofs_running << " Excited-state relaxation follows state " << this->inp_->lr_target_state
        << " of the " << (this->openshell ? "updown" : (this->target_is_ == 1 ? "triplet" : "singlet"))
        << " channel." << std::endl;
}

template<typename T, typename TR>
double ModuleESolver::ESolver_LR<T, TR>::cal_energy()
{
    // Outside a relaxation nothing consumes this, and returning a non-zero value would change
    // what the existing single-point outputs report.
    if (!this->excited_relax_) { return 0.0; }
    // `target_omega_()` is the multiplet average under `lr_relax_degen_mode = average` and the
    // single state otherwise, matching whatever `cal_lr_force_relax_` produced the gradient of.
    // The energy-based optimisers line-search on this, so the two must not describe different
    // surfaces.
    return this->etot_gs_ + this->target_omega_();
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::cal_force(BaseCell& basecell, ModuleBase::matrix& force)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    const UnitCell& ucell = static_cast<const UnitCell&>(basecell);
    if (!this->excited_relax_)
    {   // single-point runs print the gradients of every state from `after_all_runners` instead
        return;
    }
    if (this->lr_force_.nr != ucell.nat)
    {
        ModuleBase::WARNING_QUIT("ESolver_LR::cal_force",
            "the excited-state gradient has not been computed for this geometry.");
    }
    // Both halves already follow the ABACUS force convention F = -dE/dR (Ry/Bohr), so they add.
    // The "Gradients of each excited state" heading that `cal_force(int)` prints under is a
    // misnomer: the finite-difference reference it was validated against (`abacus-fd lr-custom`)
    // computes (E(-h) - E(+h))/h, which is -d(Omega)/dR, and the two agree in sign.
    force.create(ucell.nat, 3);
    force = this->force_gs_ + this->lr_force_;
    ModuleIO::print_force(GlobalV::ofs_running, ucell, "EXCITED-STATE TOTAL-FORCE (eV/Angstrom)", force, false);
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::cal_stress(BaseCell& basecell, ModuleBase::matrix& stress)
{
    static_cast<void>(stress);
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    ModuleBase::WARNING_QUIT("ESolver_LR::cal_stress",
        "the excited-state stress is not implemented (every stress term under module_lr/Grad is a "
        "dummy passed with isstress=false).");
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::init_pot_groundstate(const Charge& chg_gs)
{
    ModuleBase::TITLE("ESolver_LR", "init_pot_gs");
    std::vector<std::string> pot_register;
    if (PARAM.inp.vl_in_h)
    {
        if (!this->ks_)
        {   // on the `ks-lr` path `sfac()`/`vloc()` alias the ground-state solver's, which
            // `ESolver_FP::before_scf` already refreshed for the current geometry
            //! 11) calculate the structure factor
            this->sfac().setup(&(*this->ucell_), pgrid(), this->pw_rhod);
            this->vloc().init_vloc((*this->ucell_), this->pw_rho);
        }
        pot_register.push_back("local");
    }
    if(PARAM.inp.vh_in_h)
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
    if (PARAM.inp.test_force)
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
    ct::Tensor Xz = LR_Util::newTensor<T>({ nst, this->nloc_per_state_z_ });
    Xz.zero();
    // Closed shell: one channel, and `ispin` selects the spin COMBINATION (singlet/triplet)
    // whose X is being widened -- `nocc`/`nvirt`/`paraX_` are the same for both.
    // Open shell: one eigenvector holding both channels back to back, so both sub-blocks are
    // widened and the second one is re-based, because widening the first moves where it starts.
    const std::vector<int> chan = this->openshell ? std::vector<int>{ 0, 1 }
                                                  : std::vector<int>{ ispin };
    const int ix = this->openshell ? 0 : ispin;   // which entry of `X`
    int soff_ch = 0, zoff_ch = 0;                 // channel offsets inside one state's block
    for (const int is : chan)
    {
        const Parallel_2D& pxs = this->paraX_[is];
        const Parallel_2D& pxz = this->paraX_z_[is];
        for (int ist = 0;ist < nst;++ist)
        {
            const T* const src = this->X[ix].template data<T>()
                + (istate_begin + ist) * this->nloc_per_state + soff_ch;
            T* const dst = Xz.data<T>() + ist * this->nloc_per_state_z_ + zoff_ch;
            for (int ik = 0;ik < this->nk;++ik)
            {
                const int soff = ik * pxs.get_local_size();
                const int zoff = ik * pxz.get_local_size();
                // Both windows are block-cyclic with nb2d = 1 on the same process grid and share
                // the occupied dimension, so a global (virt, occ) element lives on the same
                // process in both -- only its local index differs. The widening is therefore a
                // purely local copy.
                for (int o = 0;o < this->nocc[is];++o)
                {
                    for (int v = 0;v < this->nvirt[is];++v)
                    {
                        if (!pxs.in_this_processor(v, o)) { continue; }
                        dst[zoff + pxz.global2local_col(o) * pxz.get_row_size() + pxz.global2local_row(v)]
                            = src[soff + pxs.global2local_col(o) * pxs.get_row_size() + pxs.global2local_row(v)];
                    }
                }
            }
        }
        soff_ch += this->nk * pxs.get_local_size();
        zoff_ch += this->nk * pxz.get_local_size();
    }
    return Xz;
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
    Z_vector_equation(Xz.template data<T>(), Z.template data<T>(),
        this->xc_kernel, nst, this->nspin, this->nbasis, this->nocc, this->nvirt_z_,
        (*this->ucell_), orb_cutoff_, this->gd(), *this->psi_ks_z_, this->eig_ks_z_,
#ifdef __EXX    
        std::weak_ptr<Exx_LRI<T>>(this->exx_lri), this->exx_info.info_global.hybrid_alpha,
#endif
        std::weak_ptr<PotHxcLR>(this->pot[ispin]), std::weak_ptr<PotHxcLR>(this->pot_hxc_gs),
        this->kv, this->paraX_z_, this->paraC_z_,
        this->paraMat_, this->spin_types[ispin], this->openshell);
    ModuleBase::timer::end("ESolver_LR", "solve_zvector_eqation");
    return Z;
}

template<typename T, typename TR>
std::vector<ModuleBase::matrix> ModuleESolver::ESolver_LR<T, TR>::cal_force(const int ispin, const int istate_only)
{
    if (PARAM.inp.test_force && ispin == 0) { this->test_force(); }
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
std::vector<ModuleBase::matrix> ModuleESolver::ESolver_LR<T, TR>::cal_force_Xz(const int ispin,
    const ct::Tensor& Xz, const std::vector<double>& omega, const int label_begin)
{
    // for each block, calculate dm_trans, dm_relaxed_diff, edm and force
    const int nst = static_cast<int>(omega.size());
    assert(static_cast<int>(Xz.shape().dim_size(0)) == nst);
    // `ist_begin`/`ist_end` label the blocks in the output only; the gradient itself never looks
    // up a state, it only uses `omega[i]`. That is what lets a caller pass excitation vectors
    // that are not the stored eigenvectors (see `cal_grad_matrix_degenerate`).
    const int ist_begin = label_begin;
    const int ist_end = label_begin + nst;
    const std::vector<int>& nvirt_g = this->nvirt_z_;
    const std::vector<Parallel_2D>& paraX_g = this->paraX_z_;
    const int nloc_g = this->nloc_per_state_z_;

    const ct::Tensor& Z = this->solve_zvector_eqation(ispin, nst, Xz);

    ModuleBase::TITLE("ESolver_LR", "cal_force");
    ModuleBase::timer::start("ESolver_LR", "cal_force");
    // Spin channel 0, NOT `ispin`. `ispin` indexes `spin_types` = {singlet, triplet}.
    // Closed shell always use spin-up channel of psi_ks, i.e. psi_ks(0).
    const auto& c = LR_Util::get_psi_spin(*this->psi_ks_z_, 0, this->nk);   // wavefunction coefficients of ground state

    // calculate the force (the partial gradient of Lagrangian)
    LR_Force<T> lr_force((*this->ucell_), this->kv.kvec_d, this->paraMat_,
        *this->pw_rhod, *this->pw_rho, this->vloc(), this->sfac(), this->gd(), this->tcb()
#ifdef __EXX
        , std::weak_ptr<Exx_LRI<T>>(this->exx_lri), this->exx_info.info_global.hybrid_alpha
#endif
    );
    GlobalV::ofs_running << "Start to calculate excited-state force of " << this->spin_types[ispin] << std::endl;
    // ground state dm for currrent spin (only for test the correctness of the force)
    // elecstate::DensityMatrix<T, T> dm_gs(this->paraMat_, 1, this->kv.kvec_d, this->nk);

    std::vector<ModuleBase::matrix> forces(ist_end - ist_begin);
    for (int istate = ist_begin;istate < ist_end;++istate)
    {
        const int offset = (istate - ist_begin) * nloc_g;   // block of X (widened into the Z window)
        const int zoffset = offset;                         // block of Z
        // The imag part will be cancelled in the force calculation, so we use double DM(R) to calculate force. 
        // But complex transition DM(R) is still used in energy density matrix calculation.
        const auto& dm_trans_k = cal_dm_trans_pblas(Xz.data<T>() + offset, paraX_g[ispin], c, this->paraC_z_, this->nbasis, this->nocc[ispin], nvirt_g[ispin], this->paraMat_);
        // D(X) complex, for the EXX (LibRI) force. Built FIRST and left UN-symmetrized:
        // the exchange kernel (mu kappa | nu lambda) puts the two indices of one D^X into
        // different electron coordinates, so Tr[D^X D^X K_exx] = (aa|ii) requires the full
        // non-symmetric D^X. Symmetrizing would give 1/2[(aa|ii)+(ai|ia)], which is wrong.
        // (`cal_force_exx_dm_trans` feeds the same tensor to both slots, so it is consistent.)
        auto dm_trans =   // D(X) complex
            LR_Util::build_dm_from_dmk<T, T>(dm_trans_k,
                this->paraMat_, this->nk, this->kv.kvec_d, (*this->ucell_), this->gd(), this->orb_cutoff_);
        LR_Util::transpose_DMR(dm_trans, (*this->ucell_).nat);
        // D(X) real, for the grid Hxc force. The Coulomb kernel (mu nu | kappa lambda) is
        // symmetric within each index pair, so it only ever sees the symmetric part of D^X.
        // In `PulayForceStress::cal_pulay_fs`, `cal_gint_rho` (which builds v) symmetrizes
        // implicitly, while `cal_gint_fvl`'s internal factor 2 assumes D_{mu nu} = D_{nu mu}.
        // Passing an un-symmetrized D^X makes the two slots of the bilinear form disagree.
        // NOTE: `build_dm_from_dmk` symmetrizes `dm_trans_k` IN PLACE, hence the ordering.
        auto dm_trans_real =   // D(X), double (FIXME: not enough for periodic system!)
            LR_Util::build_dm_from_dmk<T, double>(dm_trans_k,
                this->paraMat_, this->nk, this->kv.kvec_d, (*this->ucell_), this->gd(), this->orb_cutoff_,
                /*symmetrize=*/true);
        LR_Util::transpose_DMR(dm_trans_real, (*this->ucell_).nat);
        // LR_Util::print_DMR(dm_trans, "dm_trans of istate " + std::to_string(istate));
        // difference density matrix 
        std::vector<ct::Tensor> dm_diff_k = cal_dm_diff_pblas(Xz.data<T>() + offset, paraX_g[ispin], c, this->paraC_z_, this->nbasis, this->nocc[ispin], nvirt_g[ispin], this->paraMat_);
        std::cout << "dm_diff_k T(k) before symmetrization, istate " + std::to_string(istate) << std::endl;
        LR_Util::print_value(dm_diff_k[0].data<T>(), this->paraMat_.get_col_size(), this->paraMat_.get_row_size());
        // for (auto& d : dm_diff_k) { LR_Util::matsym(d.data<T>(), this->nbasis, this->paraMat_); }   // symmetrize
        // std::cout << "dm_diff_k T(k) after symmetrization, istate " + std::to_string(istate) << std::endl;
        // LR_Util::print_value(dm_diff_k[0].data<T>(), this->paraMat_.get_col_size(), this->paraMat_.get_row_size());

        const std::vector<ct::Tensor>& dm_relaxed_k = cal_dm_trans_pblas(Z.template data<T>() + zoffset, paraX_g[ispin], c, this->paraC_z_, this->nbasis, this->nocc[ispin], nvirt_g[ispin], this->paraMat_);
        std::cout << "dm_relaxed_k Z(k) before symmetrization, istate " + std::to_string(istate) << std::endl;
        LR_Util::print_value(dm_relaxed_k[0].data<T>(), this->paraMat_.get_col_size(), this->paraMat_.get_row_size());
        for (auto& d : dm_relaxed_k) { LR_Util::matsym(d.data<T>(), this->nbasis, this->paraMat_); }    // symmetrize
        std::cout << "dm_relaxed_k Z(k) after symmetrization, istate " + std::to_string(istate) << std::endl;
        LR_Util::print_value(dm_relaxed_k[0].data<T>(), this->paraMat_.get_col_size(), this->paraMat_.get_row_size());
        // relaxed difference density matrix
        const std::vector<ct::Tensor>& relaxed_diff_dm_k = dm_diff_k + dm_relaxed_k;
        const elecstate::DensityMatrix<T, T>& diff_dm =
            LR_Util::build_dm_from_dmk<T, T>(dm_diff_k,
                this->paraMat_, this->nk, this->kv.kvec_d, (*this->ucell_), this->gd(), this->orb_cutoff_);
        const elecstate::DensityMatrix<T, T>& relaxed_diff_dm =
            LR_Util::build_dm_from_dmk<T, T>(relaxed_diff_dm_k,
                this->paraMat_, this->nk, this->kv.kvec_d, (*this->ucell_), this->gd(), this->orb_cutoff_);
        // LR_Util::print_DMR(relaxed_diff_dm, "relaxed_diff_dm T+Z (Z symmetrized) of istate " + std::to_string(istate));

        // elecstate::DensityMatrix<T, T> relaxed_diff_dm =    // T+D(Z), (R) can be complex
        //     LR_Util::build_dm_from_dmk<T, T>(
        //         // LR_Util::operator+(
        //             cal_dm_diff_pb las(Xz.data<T>() + offset, paraX_g[ispin], c, this->paraC_z_, this->nbasis, this->nocc[ispin], nvirt_g[ispin], this->paraMat_)
        //             + cal_dm_trans_pblas(Z.template data<T>() + offset, paraX_g[ispin], c, this->paraC_z_, this->nbasis, this->nocc[ispin], nvirt_g[ispin], this->paraMat_)
        //             ,// ),
        //         this->paraMat_, this->nk, this->kv.kvec_d, (*this->ucell_), this->gd(), this->orb_cutoff_);
        // LR_Util::print_DMR(relaxed_diff_dm, "relaxed_diff_dm of istate " + std::to_string(istate));
        elecstate::DensityMatrix<T, double> relaxed_diff_dm_real(&this->paraMat_, 1, this->kv.kvec_d, this->nk);
        LR_Util::initialize_DMR(relaxed_diff_dm_real, this->paraMat_, (*this->ucell_), this->gd(), this->orb_cutoff_);
        LR_Util::get_DMR_real_imag_part(relaxed_diff_dm, relaxed_diff_dm_real, 'R');

        // get edm of type DensityMatrix
        // weak_ptr here is to avoid "could not match 'weak_ptr' against 'shared_ptr'"
        // but why there're no bug in the previous code (HamiltLR and HamiltULF)? 
        // seems because those two are classes having constructors 
        // but `cal_edm_from_XZ_istate` here is a functions
        std::weak_ptr<PotHxcLR> pot_weak = this->pot[ispin];
        std::weak_ptr<PotHxcLR> pot_hxc_gs_weak = this->pot_hxc_gs;
#ifdef __EXX
        std::weak_ptr<Exx_LRI<T>> exx_lri_weak = this->exx_lri;
#endif
        const std::vector<ct::Tensor>& edm_k =
            cal_edm_from_XZ_istate(Xz.data<T>() + offset,
                Z.template data<T>() + zoffset,
                omega[istate - ist_begin],
                // pack the following as a struct or use parameter package
                this->eig_ks_z_.c, dm_trans,
                c, this->nspin, this->nbasis, this->nocc, nvirt_g, (*this->ucell_), this->orb_cutoff_,
#ifdef __EXX
                exx_lri_weak, this->exx_info.info_global.hybrid_alpha,
#endif
                pot_weak, pot_hxc_gs_weak,
                this->kv, this->gd(), paraX_g, this->paraC_z_, this->paraMat_,
                this->xc_kernel, this->spin_types[ispin]);
        if (PARAM.inp.test_force && nocc[0] == 1 && nvirt_g[0] == 1)
        {
            const std::vector<ct::Tensor>& dm_diff = cal_dm_diff_pblas(Xz.data<T>() + offset, paraX_g[0], c, this->paraC_z_, this->nbasis, this->nocc[0], nvirt_g[0], this->paraMat_);
            // test_dm_diff_H2<T>(relaxed_diff_dm.get_DMK_pointer(0), c, this->nbasis);
            test_dm_diff_H2<T>(dm_diff[0].data<T>(), c, this->nbasis);
            test_edm_H2<T>(edm_k[0].data<T>(), this->eig_ks_z_.c, c, this->nbasis);
        }
        elecstate::DensityMatrix<T, double> edm_real = LR_Util::build_dm_from_dmk<T, double>(edm_k,
            this->paraMat_, this->nk, this->kv.kvec_d, (*this->ucell_), this->gd(), this->orb_cutoff_, /*symmetrize=*/true);
        // print edm_real (R)
        if (PARAM.inp.test_force)
        {
            LR_Util::save_DMR(edm_real, "data-EDMR-sparse" + std::string(this->excited_relax_ ? "_state" + std::to_string(istate) : ""), this->paraMat_);
            // LR_Util::print_DMR(edm_real, "edm_real (R) of istate " + std::to_string(istate));
        }

        ModuleBase::matrix force_hxc_dmtrans = lr_force.cal_force_hxc_dmtrans(dm_trans_real, *this->pot[ispin]);
        if (PARAM.inp.test_force)
            ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "HXC DMTRANS FORCE (eV/Angstrom)", force_hxc_dmtrans, false);

        const elecstate::DensityMatrix<T, double>& dm_gs = this->cal_dm_gs();

        // the $g^{xc}$ half of $\partial_x K[D^X]D^X$, i.e. the derivative of the xc kernel through
        // the ground-state density (see `cal_force_gxc_dmtrans`). Only for local kernels.
        if (LR_Util::has_local_xc(this->xc_kernel))
        {
            PotGradXCLR pot_grad(this->pot_hxc_gs->xc_kernel_components(), this->pot_hxc_gs->get_rho_basis(),
                (*this->ucell_), this->pot_hxc_gs->nrxx, this->spin_types[ispin] == "triplet");
            ModuleBase::matrix force_gxc_dmtrans = lr_force.cal_force_gxc_dmtrans(dm_trans_real, dm_gs, pot_grad);
            if (PARAM.inp.test_force)
                ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "GXC DMTRANS FORCE (eV/Angstrom)", force_gxc_dmtrans, false);
            force_hxc_dmtrans += force_gxc_dmtrans;
        }
        ModuleBase::matrix force_hamiltgs_relaxed_diff = lr_force.cal_force_hamilt_gs_dm_relaxed_diff(relaxed_diff_dm_real, dm_gs, false, this->pot_hxc_gs.get());
        if (PARAM.inp.test_force)
            ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "H_GS-(T+Z) FORCE (without EXX) (eV/Angstrom)", force_hamiltgs_relaxed_diff, false);

        ModuleBase::matrix force_overlap_edm = lr_force.cal_force_overlap_edm(edm_real);    // "-" sign has been included in the force factor
        if (PARAM.inp.test_force)
            ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "OVERLAP-EDM FORCE (eV/Angstrom)", force_overlap_edm, false);

        if (PARAM.inp.test_force)
        {
            // test H[T] force (Z=0), non-EXX part
            elecstate::DensityMatrix<T, double> diff_dm_real(&this->paraMat_, 1, this->kv.kvec_d, this->nk);
            LR_Util::initialize_DMR(diff_dm_real, this->paraMat_, (*this->ucell_), this->gd(), this->orb_cutoff_);
            LR_Util::get_DMR_real_imag_part(diff_dm, diff_dm_real, 'R');

            GlobalV::ofs_running << "========== [TEST H_GS-(T) force (Z=0), non-EXX part] ===========" << std::endl;
            ModuleBase::matrix force_hamiltgs_diff = lr_force.cal_force_hamilt_gs_dm_relaxed_diff(diff_dm_real, dm_gs);
            ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "H_GS-T FORCE (without EXX) (eV/Angstrom)", force_hamiltgs_diff, false);
            GlobalV::ofs_running << "========== [\\TEST H_GS-(T) force (Z=0), non-EXX part] ===========" << std::endl;
        }


#ifdef __EXX
        const double& alpha = this->exx_info.info_global.hybrid_alpha;

        if (LR::exx_kernel_list().count(xc_kernel))
        {
            const auto& Ds_trans = LR_Util::get_exx_Ds_spin1(dm_trans, (*this->ucell_), this->kv, this->paraMat_);
            ModuleBase::matrix force_exx_dmtrans = lr_force.cal_force_exx_dm_trans(Ds_trans, alpha * 4.0);  // cancel the two 0.5s in Ds
            if (PARAM.inp.test_force)
                ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "EXX DMTRANS FORCE (eV/Angstrom)", force_exx_dmtrans, false);
            force_hxc_dmtrans += force_exx_dmtrans;

        }

        if (LR::gs_is_hybrid())
        {
            const auto& Ds_gs = LR_Util::get_exx_Ds_spin1(dm_gs, (*this->ucell_), this->kv, this->paraMat_);    // returns 0.5*D[0]
            const auto& Ds_relaxed_diff = LR_Util::get_exx_Ds_spin1(relaxed_diff_dm, (*this->ucell_), this->kv, this->paraMat_);   // returns 0.5*D[0]
            // LR_Util::print_CV(Ds_relaxed_diff, "Ds_relaxed_diff for EXX force");
            // `get_exx_Ds_spin1` feeds `split_m2D_ktoR(..., nspin=1)`, which reads only channel 0
            // with a 0.5 prefactor. For `dm_gs` that channel is $D^\text{gs}_\uparrow$ at nspin=2
            // but the spin-summed $D^\text{gs}$ at nspin=1, i.e. twice as large.
            ModuleBase::matrix force_exx_gs_relaxed_diff = lr_force.cal_force_exx_gs_dm_relaxed_diff(Ds_gs, Ds_relaxed_diff, alpha * 4.0) * gs_dm_channel_factor();  // cancel the two 0.5s in Ds
            if (PARAM.inp.test_force)
                ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "EXX GS-(T+Z) FORCE (eV/Angstrom)", force_exx_gs_relaxed_diff, false);
            force_hamiltgs_relaxed_diff += force_exx_gs_relaxed_diff;

            if (PARAM.inp.test_force)
            {
                // test H[T] force (Z=0), EXX part
                const auto& Ds_diff = LR_Util::get_exx_Ds_spin1(diff_dm, (*this->ucell_), this->kv, this->paraMat_);   // returns 0.5*D[0]
                GlobalV::ofs_running << "========== [TEST H_GS-(T) force (Z=0), EXX part] ===========" << std::endl;
                ModuleBase::matrix force_exx_gs_diff = lr_force.cal_force_exx_gs_dm_relaxed_diff(Ds_gs, Ds_diff, alpha * 4.0) * gs_dm_channel_factor();  // cancel the two 0.5s in Ds
                ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "H_GS-T EXX FORCE (Z=0) (eV/Angstrom)", force_exx_gs_diff, false);
                GlobalV::ofs_running << "========== [\\TEST H_GS-(T) force (Z=0), EXX part] ===========" << std::endl;
            }
        }
#endif
        forces[istate - ist_begin] = force_hxc_dmtrans + force_hamiltgs_relaxed_diff + force_overlap_edm;
    }
    ModuleBase::timer::end("ESolver_LR", "cal_force");
    // total force
    print_force(forces, std::cout, ist_begin);
    print_force(forces, GlobalV::ofs_running, ist_begin);
    return forces;
}

template<typename T, typename TR>
ct::Tensor ModuleESolver::ESolver_LR<T, TR>::pad_group_to_z_(const int ispin,
    const std::vector<int>& group) const
{
    const int d = static_cast<int>(group.size());
    const int nloc_g = this->nloc_per_state_z_;
    ct::Tensor Xz = LR_Util::newTensor<T>({ d, nloc_g });
    Xz.zero();
    // One at a time rather than as a range: `group` is sorted, but nothing guarantees its members
    // are contiguous in the state list.
    for (int k = 0; k < d; ++k)
    {
        const ct::Tensor one = this->pad_X_to_z_(ispin, group[k], 1);
        std::copy(one.template data<T>(), one.template data<T>() + nloc_g,
            Xz.template data<T>() + static_cast<size_t>(k) * nloc_g);
    }
    return Xz;
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::resolve_target_multiplet_()
{
    this->target_group_.clear();
    if (LR_Util::tolower(this->inp_->lr_relax_degen_mode) == "state") { return; }
    const int ekb_off = this->openshell ? 0 : this->target_is_ * this->nstates;
    std::vector<double> omega(this->nstates);
    for (int ist = 0; ist < this->nstates; ++ist) { omega[ist] = this->pelec->ekb.c[ekb_off + ist]; }
    const std::vector<std::vector<int>> groups
        = LR::group_degenerate_states(omega, this->inp_->lr_grad_degen_thr);
    for (const std::vector<int>& g : groups)
    {
        if (std::find(g.begin(), g.end(), this->inp_->lr_target_state) == g.end()) { continue; }
        // a one-member group means the target is not degenerate here, and averaging over it would
        // be the single-state path with extra steps
        if (g.size() > 1) { this->target_group_ = g; }
        break;
    }
}

template<typename T, typename TR>
double ModuleESolver::ESolver_LR<T, TR>::target_omega_() const
{
    if (this->target_group_.empty()) { return this->pelec->ekb.c[this->target_ekb_offset_()]; }
    const int ekb_off = this->openshell ? 0 : this->target_is_ * this->nstates;
    double sum = 0.0;
    for (const int ist : this->target_group_) { sum += this->pelec->ekb.c[ekb_off + ist]; }
    return sum / static_cast<double>(this->target_group_.size());
}

template<typename T, typename TR>
ModuleBase::matrix ModuleESolver::ESolver_LR<T, TR>::cal_lr_force_relax_(std::ofstream& ofs)
{
    this->resolve_target_multiplet_();
    if (this->target_group_.empty())
    {
        return this->cal_force(this->target_is_, this->inp_->lr_target_state)[0];
    }
    // Regime (b): follow the multiplet average, whose gradient is Tr(G)/d. Only the DIAGONAL of the
    // gradient matrix is needed for this -- the average is basis-independent by construction, so
    // the off-diagonal part (and the extra d(d-1)/2 solves it costs) is not involved. The
    // Jahn-Teller distortion is orthogonal to this direction and does need them.
    const int d = static_cast<int>(this->target_group_.size());
    const int ekb_off = this->openshell ? 0 : this->target_is_ * this->nstates;
    const ct::Tensor Xz = this->pad_group_to_z_(this->target_is_, this->target_group_);
    std::vector<double> omega(d);
    for (int k = 0; k < d; ++k) { omega[k] = this->pelec->ekb.c[ekb_off + this->target_group_[k]]; }
    const std::vector<ModuleBase::matrix> forces = this->openshell
        ? this->cal_force_openshell_Xz(Xz, omega, this->target_group_.front())
        : this->cal_force_Xz(this->target_is_, Xz, omega, this->target_group_.front());
    ofs << " Target state " << this->inp_->lr_target_state << " is degenerate with";
    for (const int ist : this->target_group_)
    {
        if (ist != this->inp_->lr_target_state) { ofs << " " << ist; }
    }
    ofs << "; lr_relax_degen_mode=average, so the relaxation follows the multiplet average"
        " (Omega_bar = " << this->target_omega_() << " Ry) and stays on the symmetric"
        " configuration." << std::endl;
    return average_forces(forces);
}

template<typename T, typename TR>
ModuleBase::matrix ModuleESolver::ESolver_LR<T, TR>::MultipletLVC::average_force() const
{
    assert(!this->g.empty());
    std::vector<ModuleBase::matrix> diag;
    for (int k = 0; k < this->dim(); ++k) { diag.push_back(this->g[k][k]); }
    return average_forces(diag);
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::cal_force_and_grad_matrix_(const int ispin, std::ofstream& ofs)
{
    const std::vector<ModuleBase::matrix> forces = this->cal_force(ispin);
    if (this->inp_->lr_grad_degen_thr <= 0.0) { return; }
    // Rebuilt from scratch for this channel: a relaxation calls this once per ionic step, and a
    // stale multiplet from the previous geometry must not survive into the next one.
    this->multiplet_lvc_.erase(
        std::remove_if(this->multiplet_lvc_.begin(), this->multiplet_lvc_.end(),
            [ispin](const MultipletLVC& m) { return m.ispin == ispin; }),
        this->multiplet_lvc_.end());
    // The per-state gradients above are the diagonal of the degenerate-subspace gradient matrix in
    // whatever basis the eigensolver returned, so inside a multiplet only their trace means
    // anything. `lr_grad_degen_thr` asks for the off-diagonal part too, which is the rest of the
    // first-order information.
    const int ekb_off = this->openshell ? 0 : ispin * this->nstates;
    std::vector<double> omega(this->nstates);
    for (int ist = 0; ist < this->nstates; ++ist) { omega[ist] = this->pelec->ekb.c[ekb_off + ist]; }
    const std::vector<std::vector<int>> groups
        = LR::group_degenerate_states(omega, this->inp_->lr_grad_degen_thr);
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
            for (int i = 0; i < nloc_g; ++i) { loc += xk[i] * xl[i]; }
            Parallel_Reduce::reduce_all(loc);
            const double ref = (k == l) ? 1.0 : 0.0;
            max_ovlp_err = std::max(max_ovlp_err, std::abs(std::real(loc) - ref));
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
    if (max_ovlp_err > 1e-6)
    {
        ofs << " WARNING: this multiplet's eigenvectors are not orthonormal to"
            " 1e-6, so (X_k+X_l)/sqrt(2) is not normalized and the assembled gradient matrix is"
            " wrong by that much." << std::endl;
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
    ModuleBase::timer::start("ESolver_LR", "cal_force");

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
    const int nst = static_cast<int>(omega.size());
    assert(static_cast<int>(Xz.shape().dim_size(0)) == nst);
    const int ist_begin_ = label_begin;
    const int ist_end_ = label_begin + nst;
    const std::vector<int>& nvirt_g = this->nvirt_z_;
    const std::vector<Parallel_2D>& paraX_g = this->paraX_z_;
    const int nloc_g = this->nloc_per_state_z_;

    const ct::Tensor& Z = this->solve_zvector_eqation(0, nst, Xz);

    const std::vector<int> ld_x = { static_cast<int>(this->nk * paraX_g[0].get_local_size()),
                                    static_cast<int>(this->nk * paraX_g[1].get_local_size()) };
    const std::vector<int> off_x = { 0, ld_x[0] };
    std::vector<psi::Psi<T>> c_spin;
    for (int is : {0, 1}) { c_spin.push_back(LR_Util::get_psi_spin(*this->psi_ks_z_, is, this->nk)); }

    LR_Force<T> lr_force((*this->ucell_), this->kv.kvec_d, this->paraMat_,
        *this->pw_rhod, *this->pw_rho, this->vloc(), this->sfac(), this->gd(), this->tcb()
#ifdef __EXX
        , std::weak_ptr<Exx_LRI<T>>(this->exx_lri), this->exx_info.info_global.hybrid_alpha
#endif
    );
    GlobalV::ofs_running << "Start to calculate excited-state force of updown (open shell)" << std::endl;

    const int ist_begin = ist_begin_;
    const int ist_end = ist_end_;
    std::vector<ModuleBase::matrix> forces(ist_end - ist_begin);
    for (int istate = ist_begin;istate < ist_end;++istate)
    {
        const int offset = (istate - ist_begin) * nloc_g;   // X widened into the Z window
        const T* const X_istate = Xz.data<T>() + offset;
        const T* const Z_istate = Z.template data<T>() + offset;

        // 1. the k-space blocks of each spin channel
        std::vector<std::vector<ct::Tensor>> dmx_k(2), dmdiff_k(2), relaxed_k(2);
        for (int is : {0, 1})
        {
            dmx_k[is] = cal_dm_trans_pblas(X_istate + off_x[is], paraX_g[is], c_spin[is], this->paraC_z_,
                this->nbasis, this->nocc[is], nvirt_g[is], this->paraMat_);
            dmdiff_k[is] = cal_dm_diff_pblas(X_istate + off_x[is], paraX_g[is], c_spin[is], this->paraC_z_,
                this->nbasis, this->nocc[is], nvirt_g[is], this->paraMat_);
            std::vector<ct::Tensor> dmz_k = cal_dm_trans_pblas(Z_istate + off_x[is], paraX_g[is], c_spin[is],
                this->paraC_z_, this->nbasis, this->nocc[is], nvirt_g[is], this->paraMat_);
            for (auto& d : dmz_k) { LR_Util::matsym(d.template data<T>(), this->nbasis, this->paraMat_); }
            relaxed_k[is] = dmdiff_k[is] + dmz_k;
        }

        // 2. $D^X$. Complex and UN-symmetrized first (the EXX kernel needs the full
        //    non-symmetric $D^X$), then the real symmetrized copy for the grid Hxc force --
        //    `build_dm_from_dmk_spin` symmetrizes IN PLACE, hence the ordering.
        auto dm_trans = LR_Util::build_dm_from_dmk_spin<T, T>(dmx_k,
            this->paraMat_, this->nk, this->kv.kvec_d, (*this->ucell_), this->gd(), this->orb_cutoff_);
        LR_Util::transpose_DMR(dm_trans, (*this->ucell_).nat);
        auto dm_trans_real = LR_Util::build_dm_from_dmk_spin<T, double>(dmx_k,
            this->paraMat_, this->nk, this->kv.kvec_d, (*this->ucell_), this->gd(), this->orb_cutoff_,
            /*symmetrize=*/true);
        LR_Util::transpose_DMR(dm_trans_real, (*this->ucell_).nat);

        // 3. the relaxed difference density matrix $T+D^Z$
        const elecstate::DensityMatrix<T, T>& relaxed_diff_dm =
            LR_Util::build_dm_from_dmk_spin<T, T>(relaxed_k,
                this->paraMat_, this->nk, this->kv.kvec_d, (*this->ucell_), this->gd(), this->orb_cutoff_);
        elecstate::DensityMatrix<T, double> relaxed_diff_dm_real(&this->paraMat_, 2, this->kv.kvec_d, this->nk);
        LR_Util::initialize_DMR(relaxed_diff_dm_real, this->paraMat_, (*this->ucell_), this->gd(), this->orb_cutoff_);
        LR_Util::get_DMR_real_imag_part(relaxed_diff_dm, relaxed_diff_dm_real, 'R');

        // 4. the energy-weighted density matrix
        std::weak_ptr<PotHxcLR> pot_weak = this->pot[0];
        std::weak_ptr<PotHxcLR> pot_hxc_gs_weak = this->pot_hxc_gs;
#ifdef __EXX
        std::weak_ptr<Exx_LRI<T>> exx_lri_weak = this->exx_lri;
#endif
        const std::vector<std::vector<ct::Tensor>>& edm_k =
            cal_edm_from_XZ_istate_openshell(X_istate, Z_istate,
                omega[istate - ist_begin], this->eig_ks_z_.c, dm_trans,
                *this->psi_ks_z_, this->nspin, this->nbasis, this->nocc, nvirt_g,
                (*this->ucell_), this->orb_cutoff_,
#ifdef __EXX
                exx_lri_weak, this->exx_info.info_global.hybrid_alpha,
#endif
                pot_weak, pot_hxc_gs_weak,
                this->kv, this->gd(), paraX_g, this->paraC_z_, this->paraMat_, this->xc_kernel);
        elecstate::DensityMatrix<T, double> edm_real = LR_Util::build_dm_from_dmk_spin<T, double>(edm_k,
            this->paraMat_, this->nk, this->kv.kvec_d, (*this->ucell_), this->gd(), this->orb_cutoff_,
            /*symmetrize=*/true);

        // 5. the force terms
        ModuleBase::matrix force_hxc_dmtrans = lr_force.cal_force_hxc_dmtrans(dm_trans_real, *this->pot[0]);
        if (PARAM.inp.test_force)
            ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "HXC DMTRANS FORCE (eV/Angstrom)", force_hxc_dmtrans, false);

        const elecstate::DensityMatrix<T, double>& dm_gs = this->cal_dm_gs();

        // the $g^{xc}$ half of $\partial_x K[D^X]D^X$, i.e. the derivative of the xc kernel
        // through the ground-state density. Only for local kernels.
        if (LR_Util::has_local_xc(this->xc_kernel))
        {
            PotGradXCLR pot_grad(this->pot_hxc_gs->xc_kernel_components(), this->pot_hxc_gs->get_rho_basis(),
                (*this->ucell_), this->pot_hxc_gs->nrxx, /*triplet=*/false);
            ModuleBase::matrix force_gxc_dmtrans =
                lr_force.cal_force_gxc_dmtrans_openshell(dm_trans_real, dm_gs, pot_grad);
            if (PARAM.inp.test_force)
                ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "GXC DMTRANS FORCE (eV/Angstrom)", force_gxc_dmtrans, false);
            force_hxc_dmtrans += force_gxc_dmtrans;
        }

        ModuleBase::matrix force_hamiltgs_relaxed_diff = lr_force.cal_force_hamilt_gs_dm_relaxed_diff(
            relaxed_diff_dm_real, dm_gs, false, this->pot_hxc_gs.get());
        if (PARAM.inp.test_force)
            ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "H_GS-(T+Z) FORCE (without EXX) (eV/Angstrom)", force_hamiltgs_relaxed_diff, false);

        ModuleBase::matrix force_overlap_edm = lr_force.cal_force_overlap_edm(edm_real);
        if (PARAM.inp.test_force)
            ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "OVERLAP-EDM FORCE (eV/Angstrom)", force_overlap_edm, false);

#ifdef __EXX
        const double& alpha = this->exx_info.info_global.hybrid_alpha;
        // Exchange is spin-diagonal, so each channel is done independently and summed.
        // `get_exx_Ds_gs` returns the channels unscaled (SPIN_multiple = 1 at nspin=2), unlike
        // the closed-shell `get_exx_Ds_spin1` which returns 0.5*D. Counting the closed-shell
        // `alpha*4.0` back down for both terms:
        //   $D^XD^X$: each slot goes 0.5*D^X_tot -> D^X_is, i.e. x2 each, and the explicit
        //     sum over is adds another x2 -- but $D^X_\text{tot}=\sqrt2 D^X_\sigma$ eats one,
        //     so 4/(2*2) * ... = `alpha`.
        //   $D^\text{gs}(T{+}D^Z)$: the left slot goes 0.5*D_up -> D_is (x2); the right slot
        //     goes 0.5*(T+Z)_tot = (T+Z)_up -> (T+Z)_is (x1, no sqrt2 here); the explicit sum
        //     over is adds x2. So 4/(2*1*2) = `alpha` as well -- NOT `2*alpha`: the earlier
        //     comment forgot that the closed-shell right slot is already the spin SUM, which is
        //     exactly what the `for (is)` loop below now supplies.
        if (LR::exx_kernel_list().count(this->xc_kernel))
        {
            const auto& Ds_trans = LR_Util::get_exx_Ds_gs(dm_trans, (*this->ucell_), this->kv, this->paraMat_);
            ModuleBase::matrix force_exx_dmtrans(this->ucell_->nat, 3);
            for (int is : {0, 1})
            {
                force_exx_dmtrans += lr_force.cal_force_exx_dm_trans(Ds_trans[is], alpha, std::to_string(is));
            }
            if (PARAM.inp.test_force)
                ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "EXX DMTRANS FORCE (eV/Angstrom)", force_exx_dmtrans, false);
            force_hxc_dmtrans += force_exx_dmtrans;
        }
        if (LR::gs_is_hybrid())
        {
            const auto& Ds_gs = LR_Util::get_exx_Ds_gs(dm_gs, (*this->ucell_), this->kv, this->paraMat_);
            const auto& Ds_relaxed_diff = LR_Util::get_exx_Ds_gs(relaxed_diff_dm, (*this->ucell_), this->kv, this->paraMat_);
            ModuleBase::matrix force_exx_gs_relaxed_diff(this->ucell_->nat, 3);
            for (int is : {0, 1})
            {
                force_exx_gs_relaxed_diff += lr_force.cal_force_exx_gs_dm_relaxed_diff(
                    Ds_gs[is], Ds_relaxed_diff[is], alpha, std::to_string(is));
            }
            if (PARAM.inp.test_force)
                ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "EXX GS-(T+Z) FORCE (eV/Angstrom)", force_exx_gs_relaxed_diff, false);
            force_hamiltgs_relaxed_diff += force_exx_gs_relaxed_diff;
        }
#endif
        forces[istate - ist_begin] = force_hxc_dmtrans + force_hamiltgs_relaxed_diff + force_overlap_edm;
    }
    ModuleBase::timer::end("ESolver_LR", "cal_force");
    print_force(forces, std::cout, ist_begin);
    print_force(forces, GlobalV::ofs_running, ist_begin);
    return forces;
}

template<typename T, typename TR>
elecstate::DensityMatrix<T, double> ModuleESolver::ESolver_LR<T, TR>::cal_dm_gs()
{
    elecstate::DensityMatrix<T, double> dm_gs(&this->paraMat_, this->nspin, this->kv.kvec_d, this->nk);
    elecstate::cal_dm_psi(&this->paraMat_all_, this->wg_ks_all, *this->psi_ks_all_, dm_gs);   // nbands is important here
    LR_Util::initialize_DMR(dm_gs, this->paraMat_, (*this->ucell_), this->gd(), this->orb_cutoff_);   // nbands is not important here
    dm_gs.cal_DMR();
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

    const elecstate::DensityMatrix<T, double>& dm_gs = this->cal_dm_gs();
    // LR_Util::print_DMR(dm_gs, "DM(R) of ground state");
    ///========================== test 1: reproduce the force of ground state =========================
    // energy density matrix of the ground state
    elecstate::DensityMatrix<T, double> edm_gs(&this->paraMat_, this->nspin, this->kv.kvec_d, this->nk);   //DX
    ModuleBase::matrix wg_ekb_ks_all(nspin, PARAM.inp.nbands);
    std::transform(this->wg_ks_all.c, this->wg_ks_all.c + nspin * PARAM.inp.nbands,
        this->eig_ks_all.c, wg_ekb_ks_all.c, std::multiplies<double>());
    elecstate::cal_dm_psi(&this->paraMat_all_, wg_ekb_ks_all, *this->psi_ks_all_, edm_gs);
    LR_Util::initialize_DMR(edm_gs, this->paraMat_, (*this->ucell_), this->gd(), this->orb_cutoff_);
    edm_gs.cal_DMR();
    // ground-state force
    ModuleBase::matrix force_gs = lr_force.reproduce_force_gs(kv, dm_gs, edm_gs);
    ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "Ground State FORCE (eV/Angstrom)", force_gs, false);
    /// ======================================= END test 1 =========================================
    ///========================== test 2: reproduce the DX Hartree term =========================
    ModuleBase::matrix f_hxc_potgs = lr_force.reproduce_force_gs_loc(dm_gs, *this->pot_gs_hartree);
    ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "GS Hartree force calculated by 'cal_pulay_fs' from potential (eV/Angstrom)", f_hxc_potgs, false);
    ModuleBase::matrix f_hxc_potlr = lr_force.cal_force_hxc_dmtrans(dm_gs, *this->pot[0]);
    ModuleIO::print_force(GlobalV::ofs_running, (*this->ucell_), "GS Hxc force calculated by 'LR_Force' from kernel (eV/Angstrom)", f_hxc_potlr, false);
    // `cal_force_hxc_dmtrans` now includes the Pulay -> Pulay+Hellmann-Feynman factor 2 itself,
    // so this must match the ground-state Hartree force directly (dm_gs is already symmetric).
    /// ======================================= END test 2 =========================================
    ///========================== test 3: H2 SZ 4-center gradients =========================
    if (this->nbasis == 2 && (*this->ucell_).nat == 2)
    {
        // lr_force.cal_H2_sz_center2_deriv(orb_cutoff_, kv);  // for gradient
        // lr_force.cal_H2_sz_center4(orb_cutoff_, kv, /*is_grad=*/false);  // for 4-center integrals
        // lr_force.cal_H2_sz_center4(orb_cutoff_, kv, /*is_grad=*/true);   // for gradient
        // exit(0);
    }
}

template class ModuleESolver::ESolver_LR<double, double>;
template class ModuleESolver::ESolver_LR<std::complex<double>, double>;
