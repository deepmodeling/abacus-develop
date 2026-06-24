#include "source_hsolver/diago_lobpcg.h"
#include "source_hsolver/diago_lobpcg_detail.h"

namespace hsolver {

// ============================================================================
// State and fallback helpers
// ============================================================================

template <typename T, typename Device>
bool DiagoLobpcg<T, Device>::update_best_state(State& state,
                                                StateQuality& quality,
                                                const int candidate_notconv,
                                                const Real candidate_residual)
{
    if (!lobpcg_detail::state_is_better(candidate_notconv,
                                        candidate_residual,
                                        quality.valid,
                                        quality.notconv,
                                        quality.residual)) {
        return false;
    }
    quality.residual = candidate_residual;
    quality.notconv = candidate_notconv;
    quality.valid = true;
    this->save_state(state);
    return true;
}

template <typename T, typename Device>
template <typename UpdateFunc, typename RepairFunc>
void DiagoLobpcg<T, Device>::update_subspace_with_fallback(const char* first_failure,
                                                            const char* retry_failure,
                                                            State& rollback_state,
                                                            const int used_iter,
                                                            const UpdateFunc& update_func,
                                                            const RepairFunc& repair_func)
{
    try {
        update_func();
    } catch (const std::exception& e1) {
        this->diag_log(std::string(first_failure) + ": " + e1.what(),
                       "retry without previous search direction",
                       "iteration=" + std::to_string(used_iter));
        this->restore_state(rollback_state, true);

        try {
            update_func();
        } catch (const std::exception& e2) {
            this->diag_log(std::string(retry_failure) + ": " + e2.what(),
                           "fallback to Rayleigh-Ritz repair",
                           "iteration=" + std::to_string(used_iter));
            this->restore_state(rollback_state, false);
            repair_func();
        }
    }
}

template <typename T, typename Device>
template <typename RecomputeFunc>
void DiagoLobpcg<T, Device>::restore_best_state_if_needed(
    const char* diag_name,
    State& best_state,
    const StateQuality& best_quality,
    const int used_iter,
    Real& final_residual,
    int& final_notconv,
    const RecomputeFunc& recompute_final_quality)
{
    if (!lobpcg_detail::should_restore_best_state(best_quality.valid,
                                                  final_notconv,
                                                  final_residual,
                                                  best_quality.notconv,
                                                  best_quality.residual)) {
        return;
    }

    this->diag_log(std::string(diag_name) + " restored best residual state",
                   "final=" + std::to_string(final_residual)
                       + ", final_notconv=" + std::to_string(final_notconv),
                   "best=" + std::to_string(best_quality.residual)
                       + ", best_notconv=" + std::to_string(best_quality.notconv),
                   "context={" + this->diag_context + "}");
    this->restore_state(best_state, true);
    recompute_final_quality(final_residual, final_notconv);
}

template <typename T, typename Device>
int DiagoLobpcg<T, Device>::run_lobpcg_loop(
    const char* problem_type,
    const char* diag_name,
    const HPsiFunc& hpsi_func,
    const SPsiFunc& spsi_func,
    const std::vector<double>& effective_ethr_band,
    const int scf_iter)
{
    const int default_max_iter = (scf_iter > 1) ? this->nline : (this->nline * 20);
    const int max_iter = (this->max_iter > 0) ? this->max_iter : default_max_iter;
    int used_iter = 0;
    State rollback_state;
    State best_state;
    StateQuality best_quality;
    std::string stop_reason;

    auto compute_residual = [&]() {
        this->compute_residual_s(this->psi, this->hpsi, this->spsi, this->eigen,
                                 this->prec, this->grad, this->err_st);
    };

    for (int ntry = 0; ntry < max_iter; ++ntry) {
        used_iter = ntry + 1;
        compute_residual();
        const Real residual_before_update = this->max_error(this->err_st);
        const int notconv_before_update = this->count_not_converged(this->err_st, effective_ethr_band);
        const bool continues_for_notconv =
            should_continue_for_notconv(notconv_before_update, this->notconv_max);
        this->update_best_state(best_state,
                                best_quality,
                                notconv_before_update,
                                residual_before_update);
        const std::vector<char> soft_lock_mask =
            this->make_soft_lock_mask(this->err_st,
                                      effective_ethr_band,
                                      notconv_before_update);
        if (!continues_for_notconv) {
            break;
        }

        this->calc_spsi_with_block(spsi_func, this->grad.data<T>(), this->sgrad);
        this->orth_projection_s(this->psi, this->spsi, this->tmp_hsub,
                                this->sgrad, this->grad);
        this->calc_hpsi_with_block(hpsi_func, this->grad.data<T>(), this->hgrad);
        this->save_state(rollback_state);

        this->update_subspace_with_fallback(
            "lobpcg_update_s failed",
            "lobpcg_update_s retry failed",
            rollback_state,
            used_iter,
            [&]() {
                this->update_generalized_subspace(false);
            },
            [&]() {
                this->repair_generalized_subspace(hpsi_func, spsi_func);
            });

        bool update_rejected = false;
        compute_residual();
        Real residual_after_update = this->max_error(this->err_st);
        int notconv_after_update = this->count_not_converged(this->err_st, effective_ethr_band);
        if (!soft_lock_mask.empty()
            && this->restore_soft_locked_bands(rollback_state,
                                               soft_lock_mask,
                                               this->err_st,
                                               effective_ethr_band)) {
            compute_residual();
            residual_after_update = this->max_error(this->err_st);
            notconv_after_update = this->count_not_converged(this->err_st, effective_ethr_band);
        }
        const Real residual_growth_limit = lobpcg_detail::generalized_residual_growth_limit<Real>();
        const Real residual_limit = lobpcg_detail::residual_guard_limit(
            residual_before_update,
            residual_growth_limit);
        update_rejected = lobpcg_detail::should_reject_residual_update(
            notconv_before_update,
            residual_before_update,
            notconv_after_update,
            residual_after_update,
            residual_growth_limit);
        lobpcg_reduce_bool_or(update_rejected);
        bool stop_after_rejected_update = false;
        if (update_rejected) {
            stop_after_rejected_update = this->handle_generalized_rejected_update(
                rollback_state,
                best_state,
                best_quality,
                update_rejected,
                used_iter,
                notconv_before_update,
                residual_before_update,
                residual_after_update,
                residual_limit,
                residual_growth_limit,
                effective_ethr_band);
        } else {
            this->update_best_state(best_state, best_quality,
                                    notconv_after_update, residual_after_update);
        }
        if (stop_after_rejected_update) {
            this->diag_log("lobpcg_update_s stopped after residual-guard rollback",
                           "no accepted updated state; returning control to outer SCF",
                           "iteration=" + std::to_string(used_iter));
            stop_reason = "stopped after residual-guard rollback";
            break;
        }

        const bool has_next_iteration = (ntry + 1) < max_iter;
        const bool restart_next = has_next_iteration && scf_iter == 1 && ((ntry + 1) % this->nline == 0);
        if (has_next_iteration && !restart_next && !update_rejected) {
            try {
                this->project_generalized_search_direction(hpsi_func, spsi_func);
            } catch (const std::exception&) {
                this->restore_state(rollback_state, true);
            }
        }
        if (restart_next) {
            this->clear_search_directions();
        }
    }

    compute_residual();
    Real final_residual = this->max_error(this->err_st);
    int final_notconv = this->count_not_converged(this->err_st, effective_ethr_band);
    this->restore_best_state_if_needed(
        diag_name,
        best_state,
        best_quality,
        used_iter,
        final_residual,
        final_notconv,
        [&](Real& residual, int& notconv) {
            compute_residual();
            residual = this->max_error(this->err_st);
            notconv = this->count_not_converged(this->err_st, effective_ethr_band);
        });
    if (stop_reason.empty()
        && !should_continue_for_notconv(final_notconv, this->notconv_max)) {
        stop_reason = "stopped with allowed unconverged bands";
    }
    this->report_not_converged(problem_type,
                               used_iter,
                               max_iter,
                               effective_ethr_band,
                               stop_reason);
    return used_iter;
}

template <typename T, typename Device>
bool DiagoLobpcg<T, Device>::handle_generalized_rejected_update(
    const State& rollback_state,
    State& best_state,
    StateQuality& best_quality,
    bool& update_rejected,
    const int used_iter,
    const int notconv_before_update,
    const Real residual_before_update,
    const Real residual_after_update,
    const Real residual_limit,
    const Real residual_growth_limit,
    const std::vector<double>& effective_ethr_band)
{
    auto restore_backup_state = [&]() {
        this->restore_state(rollback_state, true);
    };
    this->diag_log("lobpcg_update_s rejected residual-increasing step",
                   "before=" + std::to_string(residual_before_update)
                       + ", after=" + std::to_string(residual_after_update)
                       + ", limit=" + std::to_string(residual_limit),
                   this->n_band_l != this->n_band ? "retry with compressed S_sub path"
                                                  : "fallback to backed-up Rayleigh-Ritz state",
                   "iteration=" + std::to_string(used_iter));
    restore_backup_state();

    Real guarded_residual = std::numeric_limits<Real>::max();
    bool stop_after_rejected_update = false;
    if (this->n_band_l != this->n_band) {
        bool compressed_ok = false;
        try {
            this->update_generalized_subspace(true);
            this->compute_residual_s(this->psi, this->hpsi, this->spsi, this->eigen,
                                     this->prec, this->grad, this->err_st);
            guarded_residual = this->max_error(this->err_st);
            const int guarded_notconv = this->count_not_converged(this->err_st, effective_ethr_band);
            compressed_ok = lobpcg_detail::compressed_guard_is_acceptable(
                notconv_before_update,
                residual_before_update,
                guarded_notconv,
                guarded_residual,
                residual_growth_limit);
            lobpcg_reduce_bool_and(compressed_ok);
        } catch (const std::exception& e) {
            this->diag_log("lobpcg_update_s compressed guard failed: " + std::string(e.what()),
                           "fallback to backed-up Rayleigh-Ritz state",
                           "iteration=" + std::to_string(used_iter));
        }
        if (!compressed_ok) {
            restore_backup_state();
            this->compute_residual_s(this->psi, this->hpsi, this->spsi, this->eigen,
                                     this->prec, this->grad, this->err_st);
            guarded_residual = this->max_error(this->err_st);
            stop_after_rejected_update = true;
        } else {
            update_rejected = false;
        }
    } else {
        this->compute_residual_s(this->psi, this->hpsi, this->spsi, this->eigen,
                                 this->prec, this->grad, this->err_st);
        guarded_residual = this->max_error(this->err_st);
        stop_after_rejected_update = true;
    }

    const int guarded_notconv = this->count_not_converged(this->err_st, effective_ethr_band);
    this->update_best_state(best_state, best_quality,
                            guarded_notconv, guarded_residual);
    return stop_after_rejected_update;
}


template <typename T, typename Device>
void DiagoLobpcg<T, Device>::diag_log(const std::string& context,
                                      const std::string& line1,
                                      const std::string& line2,
                                      const std::string& line3) const
{
    const std::string full_context = this->diag_context.empty()
        ? context
        : context + " [" + this->diag_context + "]";
    std::ostringstream oss;
    oss << " LOBPCG_DIAG " << full_context << '\n'
        << "   " << line1 << '\n'
        << "   " << line2 << '\n';
    if (!line3.empty()) {
        oss << "   " << line3 << '\n';
    }
    if (GlobalV::ofs_running.good()) {
        GlobalV::ofs_running << oss.str();
        GlobalV::ofs_running.flush();
    }
}

template <typename T, typename Device>
std::string DiagoLobpcg<T, Device>::error_with_context(const std::string& message) const
{
    if (this->diag_context.empty()) {
        return message;
    }
    return message + "; context={" + this->diag_context + "}";
}

// ============================================================================
// Constructor / Destructor
// ============================================================================

template <typename T, typename Device>
DiagoLobpcg<T, Device>::DiagoLobpcg(const Real* precondition)
{
    this->r_type   = ct::DataTypeToEnum<Real>::value;
    this->t_type   = ct::DataTypeToEnum<T>::value;
    this->dev_type = ct::DeviceTypeToEnum<Device>::value;
    this->h_prec_ptr = precondition;
    this->one     = &one_;
    this->zero    = &zero_;
    this->neg_one = &neg_one_;
}

template <typename T, typename Device>
DiagoLobpcg<T, Device>::~DiagoLobpcg() {}

// ============================================================================
// init_iter
// ============================================================================

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::init_iter(int nband, int nband_l,
                                        int nbasis, int ndim)
{
    this->n_band   = nband;
    this->n_band_l = nband_l;
    this->n_basis  = nbasis;
    this->n_dim    = ndim;
    this->nsub     = 3 * n_band;
    this->has_pdir = false;

    this->eigen     = ct::Tensor(r_type, dev_type, {this->n_band});
    this->sub_eigen = ct::Tensor(r_type, dev_type, {this->nsub});
    this->err_st    = ct::Tensor(r_type, dev_type, {this->n_band_l});

    this->hsub = ct::Tensor(t_type, dev_type, {this->nsub, this->nsub});
    this->ssub = ct::Tensor(t_type, dev_type, {this->nsub, this->nsub});
    this->tmp_hsub = ct::Tensor(t_type, dev_type, {this->n_band, this->n_band});
    this->tmp_ssub = ct::Tensor(t_type, dev_type, {this->n_band, this->n_band});

    this->hpsi  = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});
    this->spsi  = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});
    this->grad  = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});
    this->hgrad = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});
    this->sgrad = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});
    this->pdir  = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});
    this->hpdir = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});
    this->spdir = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});

    this->work   = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});
    this->hwork  = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});
    this->swork  = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});
    this->pwork  = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});
    this->hpwork = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});
    this->spwork = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});

    this->prec = ct::Tensor(r_type, dev_type, {this->n_basis});
    this->h_prec = ct::TensorMap(
        (void*)this->h_prec_ptr, this->r_type,
        ct::DeviceType::CpuDevice, {this->n_basis});

#ifdef __MPI
    this->pmmcn.set_dimension(BP_WORLD, POOL_WORLD,
                              n_band_l, n_basis,
                              n_band_l, n_basis,
                              n_dim, n_band);
    this->plintrans.set_dimension(n_dim, nband_l, n_band_l, n_basis,
                                  BP_WORLD, false);
    this->plintrans_batch2.set_dimension(2 * n_dim, nband_l, n_band_l, 2 * n_dim,
                                         BP_WORLD, false);
    this->plintrans_batch3.set_dimension(3 * n_dim, nband_l, n_band_l, 3 * n_dim,
                                         BP_WORLD, false);
#else
    this->pmmcn.set_dimension(n_band_l, n_basis,
                              n_band_l, n_basis,
                              n_dim, n_band);
    this->plintrans.set_dimension(n_dim, nband_l, n_band_l, n_basis,
                                  false);
    this->plintrans_batch2.set_dimension(2 * n_dim, nband_l, n_band_l, 2 * n_dim,
                                         false);
    this->plintrans_batch3.set_dimension(3 * n_dim, nband_l, n_band_l, 3 * n_dim,
                                         false);
#endif
}

template <typename T, typename Device>
int DiagoLobpcg<T, Device>::local_band_start() const
{
#ifdef __MPI
    if (this->plintrans.nproc_col > 1) {
        return this->plintrans.start_colB[this->plintrans.rank_col];
    }
#endif
    return 0;
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::clear_search_directions()
{
    const int psi_sz = this->n_basis * this->n_band_l;
    setmem_complex_op()(this->pdir.data<T>(), static_cast<T>(0.0), psi_sz);
    setmem_complex_op()(this->hpdir.data<T>(), static_cast<T>(0.0), psi_sz);
    setmem_complex_op()(this->spdir.data<T>(), static_cast<T>(0.0), psi_sz);
    this->has_pdir = false;
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::save_state(std::vector<T>& psi_out,
                                         std::vector<T>& hpsi_out,
                                         std::vector<Real>& eigen_out)
{
    const int psi_sz = this->n_basis * this->n_band_l;
    const int eig_sz = this->n_band;
    psi_out.resize(psi_sz);
    hpsi_out.resize(psi_sz);
    eigen_out.resize(eig_sz);

    syncmem_complex_d2h_op()(psi_out.data(), this->psi.data<T>(), psi_sz);
    syncmem_complex_d2h_op()(hpsi_out.data(), this->hpsi.data<T>(), psi_sz);
    syncmem_var_d2h_op()(eigen_out.data(), this->eigen.data<Real>(), eig_sz);
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::restore_state(const std::vector<T>& psi_in,
                                            const std::vector<T>& hpsi_in,
                                            const std::vector<Real>& eigen_in,
                                            const bool reset_search)
{
    const int psi_sz = this->n_basis * this->n_band_l;
    const int eig_sz = this->n_band;
    syncmem_complex_h2d_op()(this->psi.data<T>(), psi_in.data(), psi_sz);
    syncmem_complex_h2d_op()(this->hpsi.data<T>(), hpsi_in.data(), psi_sz);
    syncmem_var_h2d_op()(this->eigen.data<Real>(), eigen_in.data(), eig_sz);
    if (reset_search) {
        this->clear_search_directions();
    }
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::save_state(State& state)
{
    this->save_state(state.psi, state.hpsi, state.eigen);
    const int psi_sz = this->n_basis * this->n_band_l;
    state.spsi.resize(psi_sz);
    syncmem_complex_d2h_op()(state.spsi.data(), this->spsi.data<T>(), psi_sz);
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::restore_state(const State& state,
                                            const bool reset_search)
{
    this->restore_state(state.psi, state.hpsi, state.eigen, false);
    const int psi_sz = this->n_basis * this->n_band_l;
    syncmem_complex_h2d_op()(this->spsi.data<T>(), state.spsi.data(), psi_sz);
    if (reset_search) {
        this->clear_search_directions();
    }
}

template <typename T, typename Device>
std::vector<char> DiagoLobpcg<T, Device>::make_soft_lock_mask(
    const ct::Tensor& err_in,
    const std::vector<double>& ethr_band,
    const int notconv) const
{
    // Soft locking is only useful when at least one global band has converged.
    if (this->n_band_l == this->n_band || notconv <= 0 || notconv >= this->n_band) {
        return {};
    }

    std::vector<char> soft_lock_mask(this->n_band_l, 0);
    const Real* err_d = err_in.data<Real>();
    for (int ib = 0; ib < this->n_band_l; ++ib) {
        soft_lock_mask[ib] = err_d[ib] <= static_cast<Real>(ethr_band[ib]) ? 1 : 0;
    }
    return soft_lock_mask;
}

template <typename T, typename Device>
bool DiagoLobpcg<T, Device>::restore_soft_locked_bands(
    const State& state,
    const std::vector<char>& soft_lock_mask,
    const ct::Tensor& err_in,
    const std::vector<double>& ethr_band)
{
    if (soft_lock_mask.empty()) {
        return false;
    }

    const Real* err_d = err_in.data<Real>();
    Real* eigen_d = this->eigen.data<Real>();
    T* psi_d = this->psi.data<T>();
    T* hpsi_d = this->hpsi.data<T>();
    T* spsi_d = this->spsi.data<T>();
    int restored_local = 0;
    const int local_start = this->local_band_start();

    for (int ib = 0; ib < this->n_band_l; ++ib) {
        if (soft_lock_mask[ib] == 0) {
            continue;
        }
        const bool damaged = !std::isfinite(err_d[ib])
            || err_d[ib] > static_cast<Real>(ethr_band[ib]);
        if (!damaged) {
            continue;
        }

        const int offset = ib * this->n_basis;
        std::copy(state.psi.begin() + offset,
                  state.psi.begin() + offset + this->n_basis,
                  psi_d + offset);
        std::copy(state.hpsi.begin() + offset,
                  state.hpsi.begin() + offset + this->n_basis,
                  hpsi_d + offset);
        std::copy(state.spsi.begin() + offset,
                  state.spsi.begin() + offset + this->n_basis,
                  spsi_d + offset);
        eigen_d[local_start + ib] = state.eigen[local_start + ib];
        ++restored_local;
    }

    int restored_global = restored_local;
    lobpcg_reduce_sum(restored_global);
    if (restored_global > 0) {
        this->clear_search_directions();
    }
    return restored_global > 0;
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::ensure_generalized_update_finite(const ct::Tensor& psi,
                                                               const ct::Tensor& hpsi,
                                                               const ct::Tensor& spsi,
                                                               const ct::Tensor& eigen,
                                                               const int nbs,
                                                               const int n,
                                                               const char* log_context,
                                                               const char* error_message) const
{
    const int invalid_mask = update_invalid_mask(psi.data<T>(), hpsi.data<T>(), spsi.data<T>(),
                                                 eigen.data<Real>(), this->n_band_l, nbs, this->n_dim, n);
    if (invalid_mask == 0) {
        return;
    }
    this->diag_log(log_context,
                   vector_block_diagnostics("psi", psi.data<T>(), this->n_band_l, nbs, this->n_dim),
                   vector_block_diagnostics("hpsi", hpsi.data<T>(), this->n_band_l, nbs, this->n_dim),
                   vector_block_diagnostics("spsi", spsi.data<T>(), this->n_band_l, nbs, this->n_dim)
                       + ", eigen_invalid=" + std::to_string((invalid_mask & 8) != 0)
                       + ", invalid_mask=" + std::to_string(invalid_mask));
    throw std::runtime_error(error_message);
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::update_generalized_subspace(const bool force_compressed)
{
    if (this->n_band_l != this->n_band) {
        this->lobpcg_update_s_parallel(this->psi, this->hpsi, this->spsi,
                                       this->grad, this->hgrad, this->sgrad,
                                       this->pdir, this->hpdir, this->spdir,
                                       this->eigen, force_compressed);
    } else {
        this->lobpcg_update_s(this->psi, this->hpsi, this->spsi,
                              this->grad, this->hgrad, this->sgrad,
                              this->pdir, this->hpdir, this->spdir,
                              this->eigen);
    }
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::repair_generalized_subspace(const HPsiFunc& hpsi_func,
                                                          const SPsiFunc& spsi_func)
{
    this->calc_hpsi_with_block(hpsi_func, this->psi.data<T>(), this->hpsi);
    this->calc_spsi_with_block(spsi_func, this->psi.data<T>(), this->spsi);
    if (this->n_band_l != this->n_band) {
        this->generalized_rayleigh_ritz_parallel(this->psi, this->hpsi, this->spsi, this->eigen);
    } else {
        this->generalized_rayleigh_ritz(this->psi, this->hpsi, this->spsi, this->eigen);
    }
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::project_generalized_search_direction(const HPsiFunc& hpsi_func,
                                                                   const SPsiFunc& spsi_func)
{
    if (this->n_band_l != this->n_band) {
        this->orth_projection_s_with_h(this->psi, this->hpsi, this->spsi,
                                       this->tmp_hsub, this->hpdir,
                                       this->spdir, this->pdir);
    } else {
        this->calc_spsi_with_block(spsi_func, this->pdir.data<T>(), this->spdir);
        this->orth_projection_s(this->psi, this->spsi, this->tmp_hsub,
                                this->spdir, this->pdir);
        this->calc_hpsi_with_block(hpsi_func, this->pdir.data<T>(), this->hpdir);
    }
}

// ============================================================================
// calc_prec
// ============================================================================

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::calc_prec()
{
    syncmem_var_h2d_op()(this->prec.template data<Real>(),
                         this->h_prec.template data<Real>(),
                         this->n_basis);
}

// ============================================================================
// calc_hpsi_with_block / calc_spsi_with_block
// ============================================================================

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::calc_hpsi_with_block(
    const HPsiFunc& hpsi_func, T* psi_in, ct::Tensor& hpsi_out)
{
    hpsi_func(psi_in, hpsi_out.data<T>(), this->n_basis, this->n_band_l);
    int invalid_mask =
        (!finite_vector_block(psi_in, this->n_band_l, this->n_basis, this->n_dim) ? 1 : 0)
      | (!finite_vector_block(hpsi_out.data<T>(), this->n_band_l, this->n_basis, this->n_dim) ? 2 : 0);
    lobpcg_reduce_bit_or(invalid_mask);
    if (invalid_mask != 0) {
        this->diag_log("calc_hpsi_with_block non-finite",
                        vector_block_diagnostics("psi_in",
                                                 psi_in,
                                                 this->n_band_l,
                                                 this->n_basis,
                                                 this->n_dim),
                        vector_block_diagnostics("hpsi_out",
                                                 hpsi_out.data<T>(),
                                                 this->n_band_l,
                                                 this->n_basis,
                                                 this->n_dim));
        throw std::runtime_error(this->error_with_context("LOBPCG hPsi produced non-finite values"));
    }
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::calc_spsi_with_block(
    const SPsiFunc& spsi_func, const T* psi_in, ct::Tensor& spsi_out)
{
    spsi_func(psi_in, spsi_out.data<T>(), this->n_basis, this->n_band_l);
    int invalid_mask =
        (!finite_vector_block(psi_in, this->n_band_l, this->n_basis, this->n_dim) ? 1 : 0)
      | (!finite_vector_block(spsi_out.data<T>(), this->n_band_l, this->n_basis, this->n_dim) ? 2 : 0);
    lobpcg_reduce_bit_or(invalid_mask);
    if (invalid_mask != 0) {
        this->diag_log("calc_spsi_with_block non-finite",
                        vector_block_diagnostics("psi_in",
                                                 psi_in,
                                                 this->n_band_l,
                                                 this->n_basis,
                                                 this->n_dim),
                        vector_block_diagnostics("spsi_out",
                                                 spsi_out.data<T>(),
                                                 this->n_band_l,
                                                 this->n_basis,
                                                 this->n_dim));
        throw std::runtime_error(this->error_with_context("LOBPCG sPsi produced non-finite values"));
    }
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::repair_initial_subspace_s(
    const HPsiFunc& hpsi_func, const SPsiFunc& spsi_func)
{
    const int nb = this->n_band;
    const int lda = this->n_basis;
    const int nvalid = this->n_dim;
    const Real eps = static_cast<Real>(100)
                   * std::numeric_limits<Real>::epsilon();

    T* psi_d = this->psi.data<T>();
    T* spsi_d = this->spsi.data<T>();
    bool repaired = false;

    for (int ib = 0; ib < nb; ++ib) {
        bool ready = false;
        for (int attempt = -1; attempt < nvalid && !ready; ++attempt) {
            if (attempt >= 0) {
                repaired = true;
                for (int ig = 0; ig < lda; ++ig) {
                    psi_d[ib * lda + ig] = static_cast<T>(0.0);
                    spsi_d[ib * lda + ig] = static_cast<T>(0.0);
                }
                psi_d[ib * lda + ((ib + attempt) % nvalid)] = static_cast<T>(1.0);
                spsi_func(psi_d + ib * lda, spsi_d + ib * lda, lda, 1);
            }

            bool finite_vec = true;
            for (int ig = 0; ig < nvalid; ++ig) {
                const T qv = psi_d[ib * lda + ig];
                const T sqv = spsi_d[ib * lda + ig];
                if (!std::isfinite(std::real(qv)) || !std::isfinite(std::imag(qv))
                    || !std::isfinite(std::real(sqv)) || !std::isfinite(std::imag(sqv))) {
                    finite_vec = false;
                    break;
                }
            }
            if (!finite_vec) {
                repaired = true;
                continue;
            }

            bool finite_projection = true;
            for (int pass = 0; pass < 2; ++pass) {
                for (int jb = 0; jb < ib; ++jb) {
                    T dot = static_cast<T>(0.0);
                    for (int ig = 0; ig < nvalid; ++ig) {
                        dot += std::conj(psi_d[jb * lda + ig])
                             * spsi_d[ib * lda + ig];
                    }
#ifdef __MPI
                    Parallel_Reduce::reduce_pool(&dot, 1);
#endif
                    if (!std::isfinite(std::real(dot)) || !std::isfinite(std::imag(dot))) {
                        finite_projection = false;
                        break;
                    }
                    for (int ig = 0; ig < nvalid; ++ig) {
                        psi_d[ib * lda + ig] -= dot * psi_d[jb * lda + ig];
                        spsi_d[ib * lda + ig] -= dot * spsi_d[jb * lda + ig];
                    }
                }
                if (!finite_projection) {
                    break;
                }
            }
            if (!finite_projection) {
                repaired = true;
                continue;
            }

            Real norm2 = static_cast<Real>(0.0);
            for (int ig = 0; ig < nvalid; ++ig) {
                norm2 += std::real(std::conj(psi_d[ib * lda + ig])
                                  * spsi_d[ib * lda + ig]);
            }
#ifdef __MPI
            Parallel_Reduce::reduce_pool(&norm2, 1);
#endif
            if (!std::isfinite(norm2) || norm2 <= eps) {
                repaired = true;
                continue;
            }

            const Real inv_norm = static_cast<Real>(1.0) / std::sqrt(norm2);
            for (int ig = 0; ig < nvalid; ++ig) {
                psi_d[ib * lda + ig] *= inv_norm;
                spsi_d[ib * lda + ig] *= inv_norm;
            }
            for (int ig = nvalid; ig < lda; ++ig) {
                psi_d[ib * lda + ig] = static_cast<T>(0.0);
                spsi_d[ib * lda + ig] = static_cast<T>(0.0);
            }
            ready = true;
        }

        if (!ready) {
            throw std::runtime_error("LOBPCG failed to repair rank-deficient initial subspace");
        }
    }

    if (repaired) {
        this->calc_hpsi_with_block(hpsi_func, this->psi.data<T>(), this->hpsi);
        this->calc_spsi_with_block(spsi_func, this->psi.data<T>(), this->spsi);
    } else {
        this->calc_hpsi_with_block(hpsi_func, this->psi.data<T>(), this->hpsi);
    }

    for (int ib = 0; ib < nb; ++ib) {
        Real norm2 = static_cast<Real>(0.0);
        for (int ig = 0; ig < nvalid; ++ig) {
            norm2 += std::real(std::conj(psi_d[ib * lda + ig])
                              * spsi_d[ib * lda + ig]);
        }
#ifdef __MPI
        Parallel_Reduce::reduce_pool(&norm2, 1);
#endif
        if (!std::isfinite(norm2) || !(norm2 > eps)) {
            throw std::runtime_error("LOBPCG repaired initial subspace has invalid S-norm at band "
                                     + std::to_string(ib)
                                     + ", norm2=" + std::to_string(norm2));
        }
    }
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::generalized_rayleigh_ritz(
    ct::Tensor& psi_inout, ct::Tensor& hpsi_inout,
    ct::Tensor& spsi_inout, ct::Tensor& eigen_out)
{
    const int nb = this->n_band;
    const int nbs = this->n_basis;
    const int nvalid = this->n_dim;
    const int local_sz = this->n_band_l * nbs;

    this->pgemm_multiply(this->one_,
                         psi_inout.data<T>(),
                         hpsi_inout.data<T>(),
                         this->zero_,
                         this->tmp_hsub.data<T>());
    this->pgemm_multiply(this->one_,
                         psi_inout.data<T>(),
                         spsi_inout.data<T>(),
                         this->zero_,
                         this->tmp_ssub.data<T>());
    hermitize(this->tmp_hsub.data<T>(), nb, nb);
    hermitize(this->tmp_ssub.data<T>(), nb, nb);

    bool rr_ok = false;
    std::string rr_error;
    try {
        ct::kernels::lapack_hegvd<T, ct_Device>()(
            nb, nb,
            this->tmp_hsub.data<T>(),
            this->tmp_ssub.data<T>(),
            eigen_out.data<Real>(),
            this->hsub.data<T>());
        rr_ok = finite_real_block(eigen_out.data<Real>(), nb)
             && finite_scalar_block(this->hsub.data<T>(), nb * nb);
        if (!rr_ok) {
            this->diag_log("generalized_rayleigh_ritz hegvd returned non-finite eigen data",
                            hermitian_matrix_diagnostics("H_sub", this->tmp_hsub.data<T>(), nb, nb),
                            hermitian_matrix_diagnostics("S_sub", this->tmp_ssub.data<T>(), nb, nb),
                            vector_block_diagnostics("eigvec",
                                                     this->hsub.data<T>(),
                                                     nb,
                                                     nb,
                                                     nb));
        }
    } catch (const std::exception& e) {
        rr_error = e.what();
        rr_ok = false;
    }

    if (!rr_ok) {
        this->diag_log("generalized_rayleigh_ritz hegvd failed"
                        + (rr_error.empty() ? std::string() : (": " + rr_error)),
                        hermitian_matrix_diagnostics("H_sub", this->tmp_hsub.data<T>(), nb, nb),
                        hermitian_matrix_diagnostics("S_sub", this->tmp_ssub.data<T>(), nb, nb),
                        s_overlap_diagnostics(this->tmp_ssub.data<T>(), nb, nb));
        throw std::runtime_error("LOBPCG generalized Rayleigh-Ritz failed in hegvd");
    }

    const T* rotate_in[3] = {psi_inout.data<T>(), hpsi_inout.data<T>(), spsi_inout.data<T>()};
    T* rotate_out[3] = {this->work.data<T>(), this->hwork.data<T>(), this->swork.data<T>()};
    this->plintrans_batched_act(this->one_, rotate_in, 3,
                                this->hsub.data<T>(), this->zero_, rotate_out);
    syncmem_complex_op()(psi_inout.data<T>(), this->work.data<T>(), local_sz);
    syncmem_complex_op()(hpsi_inout.data<T>(), this->hwork.data<T>(), local_sz);
    syncmem_complex_op()(spsi_inout.data<T>(), this->swork.data<T>(), local_sz);

    if (!finite_vector_block(psi_inout.data<T>(), nb, nbs, nvalid)
        || !finite_vector_block(hpsi_inout.data<T>(), nb, nbs, nvalid)
        || !finite_vector_block(spsi_inout.data<T>(), nb, nbs, nvalid)) {
        this->diag_log("generalized_rayleigh_ritz rotation produced non-finite vectors",
                        vector_block_diagnostics("psi", psi_inout.data<T>(), nb, nbs, nvalid),
                        vector_block_diagnostics("hpsi", hpsi_inout.data<T>(), nb, nbs, nvalid),
                        vector_block_diagnostics("spsi", spsi_inout.data<T>(), nb, nbs, nvalid));
        throw std::runtime_error("LOBPCG generalized Rayleigh-Ritz rotation produced non-finite vectors");
    }
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::generalized_rayleigh_ritz_parallel(
    ct::Tensor& psi_inout, ct::Tensor& hpsi_inout,
    ct::Tensor& spsi_inout, ct::Tensor& eigen_out)
{
    const int nb = this->n_band;
    const int nbs = this->n_basis;
    const int local_sz = this->n_band_l * nbs;

    T* hsub_d = this->hsub.data<T>();
    T* ssub_d = this->ssub.data<T>();
    setmem_complex_op()(hsub_d, static_cast<T>(0.0), this->nsub * this->nsub);
    setmem_complex_op()(ssub_d, static_cast<T>(0.0), this->nsub * this->nsub);

    this->pgemm_multiply(this->one_,
                         psi_inout.data<T>(),
                         hpsi_inout.data<T>(),
                         this->zero_,
                         this->tmp_hsub.data<T>());
    this->pgemm_multiply(this->one_,
                         psi_inout.data<T>(),
                         spsi_inout.data<T>(),
                         this->zero_,
                         this->tmp_ssub.data<T>());
    ModuleBase::matrixCopy<T, Device>()(nb, nb, this->tmp_hsub.data<T>(), nb, hsub_d, this->nsub);
    ModuleBase::matrixCopy<T, Device>()(nb, nb, this->tmp_ssub.data<T>(), nb, ssub_d, this->nsub);
    hermitize(hsub_d, this->nsub, nb);
    hermitize(ssub_d, this->nsub, nb);

    std::vector<Real> inv_subspace_norm;
    std::string scale_error;
    if (!scale_subspace_by_overlap_diag(hsub_d, ssub_d, this->nsub, nb,
                                        inv_subspace_norm, scale_error)) {
        this->diag_log("generalized_rayleigh_ritz_parallel failed to normalize S_sub before hegvd",
                       hermitian_matrix_diagnostics("S_sub", ssub_d, this->nsub, nb),
                       s_overlap_diagnostics(ssub_d, this->nsub, nb),
                       scale_error);
        throw std::runtime_error("LOBPCG generalized parallel initial overlap has invalid diagonal");
    }

    const auto spd_check = check_subspace_spd<T, ct_Device>(ssub_d, this->nsub, nb);
    if (!spd_check.ok) {
        this->diag_log("generalized_rayleigh_ritz_parallel rejected ill-conditioned S_sub before hegvd",
                       hermitian_matrix_diagnostics("S_sub", ssub_d, this->nsub, nb),
                       s_overlap_diagnostics(ssub_d, this->nsub, nb),
                       subspace_spd_diagnostics(spd_check));
        throw std::runtime_error("LOBPCG generalized parallel initial overlap is ill-conditioned before hegvd");
    }

    try {
        ct::kernels::lapack_hegvd<T, ct_Device>()(
            nb, this->nsub, hsub_d, ssub_d, eigen_out.data<Real>(), hsub_d);
    } catch (const std::exception& e) {
        this->diag_log("generalized_rayleigh_ritz_parallel hegvd failed: " + std::string(e.what()),
                       hermitian_matrix_diagnostics("H_sub", hsub_d, this->nsub, nb),
                       hermitian_matrix_diagnostics("S_sub", ssub_d, this->nsub, nb),
                       subspace_spd_diagnostics(spd_check));
        throw;
    }

    if (!finite_real_block(eigen_out.data<Real>(), nb)
        || !finite_scalar_block(hsub_d, nb * this->nsub)) {
        this->diag_log("generalized_rayleigh_ritz_parallel hegvd produced non-finite values",
                       hermitian_matrix_diagnostics("H_sub", hsub_d, this->nsub, nb),
                       hermitian_matrix_diagnostics("S_sub", ssub_d, this->nsub, nb),
                       subspace_spd_diagnostics(spd_check));
        throw std::runtime_error("LOBPCG generalized parallel initial diagonalization produced non-finite values");
    }

    for (int jc = 0; jc < nb; ++jc) {
        for (int ir = 0; ir < nb; ++ir) {
            hsub_d[jc * this->nsub + ir] *= inv_subspace_norm[ir];
        }
    }

    ModuleBase::matrixCopy<T, Device>()(nb, nb, hsub_d, this->nsub, this->tmp_hsub.data<T>(), nb);

    const T* rotate_in[3] = {psi_inout.data<T>(), hpsi_inout.data<T>(), spsi_inout.data<T>()};
    T* rotate_out[3] = {this->work.data<T>(), this->hwork.data<T>(), this->swork.data<T>()};
    this->plintrans_batched_act(this->one_, rotate_in, 3,
                                this->tmp_hsub.data<T>(), this->zero_, rotate_out);
    syncmem_complex_op()(psi_inout.data<T>(), this->work.data<T>(), local_sz);
    syncmem_complex_op()(hpsi_inout.data<T>(), this->hwork.data<T>(), local_sz);
    syncmem_complex_op()(spsi_inout.data<T>(), this->swork.data<T>(), local_sz);

    if (!finite_vector_block(psi_inout.data<T>(), this->n_band_l, nbs, this->n_dim)
        || !finite_vector_block(hpsi_inout.data<T>(), this->n_band_l, nbs, this->n_dim)
        || !finite_vector_block(spsi_inout.data<T>(), this->n_band_l, nbs, this->n_dim)) {
        this->diag_log("generalized_rayleigh_ritz_parallel rotation produced non-finite vectors",
                       vector_block_diagnostics("psi", psi_inout.data<T>(), this->n_band_l, nbs, this->n_dim),
                       vector_block_diagnostics("hpsi", hpsi_inout.data<T>(), this->n_band_l, nbs, this->n_dim),
                       vector_block_diagnostics("spsi", spsi_inout.data<T>(), this->n_band_l, nbs, this->n_dim));
        throw std::runtime_error("LOBPCG generalized parallel initial rotation produced non-finite vectors");
    }
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::compute_residual_s(
    const ct::Tensor& psi_in, const ct::Tensor& hpsi_in,
    const ct::Tensor& spsi_in, const ct::Tensor& eigen_in,
    const ct::Tensor& prec_in, ct::Tensor& grad_out, ct::Tensor& err_out)
{
    const Real* _prec  = prec_in.data<Real>();
    const Real* _eigen = eigen_in.data<Real>();
    const T*    _hpsi  = hpsi_in.data<T>();
    const T*    _spsi  = spsi_in.data<T>();
    T*          _grad  = grad_out.data<T>();
    Real*       _err   = err_out.data<Real>();
    const int band_start = this->local_band_start();
    int invalid_grad = 0;
    std::string invalid_grad_context;

    for (int ib = 0; ib < this->n_band_l; ib++) {
        const int  ioff   = ib * this->n_basis;
        const Real lambda = _eigen[band_start + ib];
        Real       err_j  = 0.0;
        for (int ig = 0; ig < this->n_dim; ig++) {
            const int idx = ioff + ig;
            const T   r   = _hpsi[idx] - lambda * _spsi[idx];
            const Real denom = std::max(_prec[ig], static_cast<Real>(1e-8));
            _grad[idx] = r / denom;
            if (!std::isfinite(std::real(_grad[idx])) || !std::isfinite(std::imag(_grad[idx]))) {
                invalid_grad = 1;
                if (invalid_grad_context.empty()) {
                    std::ostringstream oss;
                    oss << "ib=" << ib
                        << " ig=" << ig
                        << " lambda=" << lambda
                        << " prec=" << _prec[ig]
                        << " denom=" << denom
                        << " hpsi=(" << std::real(_hpsi[idx]) << "," << std::imag(_hpsi[idx]) << ")"
                        << " spsi=(" << std::real(_spsi[idx]) << "," << std::imag(_spsi[idx]) << ")"
                        << " residual=(" << std::real(r) << "," << std::imag(r) << ")"
                        << " grad=(" << std::real(_grad[idx]) << "," << std::imag(_grad[idx]) << ")";
                    invalid_grad_context = oss.str();
                }
                _grad[idx] = static_cast<T>(0.0);
            } else {
                err_j += std::norm(r);
            }
        }
        for (int ig = this->n_dim; ig < this->n_basis; ig++)
            _grad[ioff + ig] = static_cast<T>(0.0);
        _err[ib] = invalid_grad == 0 ? err_j : static_cast<Real>(0.0);
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(_err, this->n_band_l);
#endif
    lobpcg_reduce_bit_or(invalid_grad);
    if (invalid_grad != 0) {
        this->diag_log("compute_residual_s non-finite grad",
                       invalid_grad_context.empty() ? "non-finite gradient on another rank" : invalid_grad_context,
                       vector_block_diagnostics("hpsi", _hpsi, this->n_band_l, this->n_basis, this->n_dim),
                       vector_block_diagnostics("spsi", _spsi, this->n_band_l, this->n_basis, this->n_dim));
        throw std::runtime_error("LOBPCG generalized residual produced non-finite gradient");
    }
    for (int ib = 0; ib < this->n_band_l; ib++)
        _err[ib] = std::sqrt(_err[ib]);
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::orth_projection_s(
    const ct::Tensor& psi_in, const ct::Tensor& spsi_in,
    ct::Tensor& hsub_work, ct::Tensor& sgrad_out, ct::Tensor& grad_out)
{
    const int nbs = this->n_basis;
    const int nvalid = this->n_dim;

    this->pgemm_multiply(this->one_,
                         psi_in.data<T>(),
                         sgrad_out.data<T>(),
                         this->zero_,
                         hsub_work.data<T>());
    const T* proj_in[2] = {psi_in.data<T>(), spsi_in.data<T>()};
    T* proj_out[2] = {grad_out.data<T>(), sgrad_out.data<T>()};
    this->plintrans_batched_act(this->neg_one_, proj_in, 2,
                                hsub_work.data<T>(), this->one_, proj_out);

    T* grad = grad_out.data<T>();
    T* sgrad = sgrad_out.data<T>();
    for (int jb = 0; jb < this->n_band_l; jb++) {
        for (int ig = nvalid; ig < nbs; ig++) {
            grad[jb * nbs + ig] = static_cast<T>(0.0);
            sgrad[jb * nbs + ig] = static_cast<T>(0.0);
        }
    }
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::orth_projection_s_with_h(
    const ct::Tensor& psi_in, const ct::Tensor& hpsi_in,
    const ct::Tensor& spsi_in, ct::Tensor& hsub_work,
    ct::Tensor& hpdir_out, ct::Tensor& spdir_out, ct::Tensor& pdir_out)
{
    const int nbs = this->n_basis;
    const int nvalid = this->n_dim;

    this->pgemm_multiply(this->one_,
                         psi_in.data<T>(),
                         spdir_out.data<T>(),
                         this->zero_,
                         hsub_work.data<T>());
    const T* proj_in[3] = {psi_in.data<T>(), spsi_in.data<T>(), hpsi_in.data<T>()};
    T* proj_out[3] = {pdir_out.data<T>(), spdir_out.data<T>(), hpdir_out.data<T>()};
    this->plintrans_batched_act(this->neg_one_, proj_in, 3,
                                hsub_work.data<T>(), this->one_, proj_out);

    T* pdir = pdir_out.data<T>();
    T* spdir = spdir_out.data<T>();
    T* hpdir = hpdir_out.data<T>();
    for (int jb = 0; jb < this->n_band_l; jb++) {
        for (int ig = nvalid; ig < nbs; ig++) {
            pdir[jb * nbs + ig] = static_cast<T>(0.0);
            spdir[jb * nbs + ig] = static_cast<T>(0.0);
            hpdir[jb * nbs + ig] = static_cast<T>(0.0);
        }
    }
}

// ============================================================================
// rotate_wf
// ============================================================================

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::rotate_wf(
    const ct::Tensor& hsub_in, ct::Tensor& psi_out, ct::Tensor& workspace_in)
{
    const int nbs = this->n_basis;
    const int nvalid = this->n_dim;
    this->plintrans_act(1.0, psi_out.data<T>(), hsub_in.data<T>(),
                        0.0, workspace_in.data<T>());
    T* workspace = workspace_in.data<T>();
    for (int ib = 0; ib < this->n_band_l; ++ib) {
        for (int ig = nvalid; ig < nbs; ++ig) {
            workspace[ib * nbs + ig] = static_cast<T>(0.0);
        }
    }
    syncmem_complex_op()(psi_out.data<T>(), workspace_in.data<T>(),
                         this->n_band_l * nbs);
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::orth_cholesky_s(
    ct::Tensor&, ct::Tensor& psi_out, ct::Tensor& hpsi_out,
    ct::Tensor& spsi_out, ct::Tensor&)
{
    const int nb  = this->n_band;
    const int nbs = this->n_basis;
    const int nvalid = this->n_dim;

    T* psi_d  = psi_out.data<T>();
    T* hpsi_d = hpsi_out.data<T>();
    T* spsi_d = spsi_out.data<T>();
    const Real eps = static_cast<Real>(100)
                   * std::numeric_limits<Real>::epsilon();

    for (int ib = 0; ib < nb; ++ib) {
        for (int jb = 0; jb < ib; ++jb) {
            T dot = static_cast<T>(0.0);
            for (int ig = 0; ig < nvalid; ++ig) {
                dot += std::conj(psi_d[jb * nbs + ig])
                     * spsi_d[ib * nbs + ig];
            }
#ifdef __MPI
            Parallel_Reduce::reduce_pool(&dot, 1);
#endif
            for (int ig = 0; ig < nvalid; ++ig) {
                psi_d[ib * nbs + ig]  -= dot * psi_d[jb * nbs + ig];
                hpsi_d[ib * nbs + ig] -= dot * hpsi_d[jb * nbs + ig];
                spsi_d[ib * nbs + ig] -= dot * spsi_d[jb * nbs + ig];
            }
        }

        Real norm2 = static_cast<Real>(0.0);
        for (int ig = 0; ig < nvalid; ++ig) {
            norm2 += std::real(std::conj(psi_d[ib * nbs + ig])
                              * spsi_d[ib * nbs + ig]);
        }
#ifdef __MPI
        Parallel_Reduce::reduce_pool(&norm2, 1);
#endif
        if (!(norm2 > eps)) {
            throw std::runtime_error("orth_cholesky_s failed: dependent vectors at band "
                                     + std::to_string(ib)
                                     + ", norm2=" + std::to_string(norm2)
                                     + ", nvalid=" + std::to_string(nvalid)
                                     + ", lda=" + std::to_string(nbs));
        }
        const Real inv_norm = static_cast<Real>(1.0) / std::sqrt(norm2);
        for (int ig = 0; ig < nvalid; ++ig) {
            psi_d[ib * nbs + ig]  *= inv_norm;
            hpsi_d[ib * nbs + ig] *= inv_norm;
            spsi_d[ib * nbs + ig] *= inv_norm;
        }
        for (int ig = nvalid; ig < nbs; ++ig) {
            psi_d[ib * nbs + ig] = static_cast<T>(0.0);
            hpsi_d[ib * nbs + ig] = static_cast<T>(0.0);
            spsi_d[ib * nbs + ig] = static_cast<T>(0.0);
        }
    }
}

// ============================================================================
// test_error
// ============================================================================

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::validate_ethr_band(const std::vector<double>& ethr_band) const
{
    if (ethr_band.size() != static_cast<size_t>(this->n_band_l)) {
        std::ostringstream oss;
        oss << "LOBPCG local ethr_band size mismatch: size=" << ethr_band.size()
            << ", required local bands=" << this->n_band_l
            << ", global bands=" << this->n_band;
        if (!this->diag_context.empty()) {
            oss << ", context={" << this->diag_context << "}";
        }
        throw std::invalid_argument(oss.str());
    }
}

template <typename T, typename Device>
typename DiagoLobpcg<T, Device>::Real DiagoLobpcg<T, Device>::max_error(
    const ct::Tensor& err_in) const
{
    Real* err = err_in.data<Real>();
    std::vector<Real> tmp_cpu;
    if (err_in.device_type() == ct::DeviceType::GpuDevice) {
        tmp_cpu.resize(this->n_band_l);
        err = tmp_cpu.data();
        syncmem_var_d2h_op()(err, err_in.data<Real>(), this->n_band_l);
    }

    Real max_residual = max_residual_local(err, this->n_band_l);
    lobpcg_reduce_max(max_residual);
    return max_residual;
}

template <typename T, typename Device>
int DiagoLobpcg<T, Device>::count_not_converged(
    const ct::Tensor& err_in,
    const std::vector<double>& ethr_band) const
{
    Real* err = err_in.data<Real>();
    std::vector<Real> tmp_cpu;
    if (err_in.device_type() == ct::DeviceType::GpuDevice) {
        tmp_cpu.resize(this->n_band_l);
        err = tmp_cpu.data();
        syncmem_var_d2h_op()(err, err_in.data<Real>(), this->n_band_l);
    }

    int notconv = count_not_converged_local(err,
                                            ethr_band.data(),
                                            this->n_band_l);
    lobpcg_reduce_sum(notconv);
    return notconv;
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::report_not_converged(
    const char* problem_type,
    const int used_iter,
    const int max_iter,
    const std::vector<double>& ethr_band,
    const std::string& stop_reason) const
{
    const Real* err = this->err_st.data<Real>();
    int notconv = 0;
    Real max_residual = static_cast<Real>(0.0);
    for (int ib = 0; ib < this->n_band_l; ++ib) {
        max_residual = std::max(max_residual, err[ib]);
        if (err[ib] > ethr_band[ib]) {
            ++notconv;
        }
    }
    lobpcg_reduce_sum(notconv);
    lobpcg_reduce_max(max_residual);
    if (notconv <= 0) {
        return;
    }

    std::ostringstream msg;
    msg << "DiagoLobpcg::diag(" << problem_type << ") "
        << (stop_reason.empty() ? std::string("reached max_iter=") + std::to_string(max_iter)
                                : stop_reason)
        << " after " << used_iter
        << " iterations; notconv=" << notconv
        << ", max_residual=" << max_residual;
    if (!this->diag_context.empty()) {
        msg << ", context={" << this->diag_context << "}";
    }
    print_lobpcg_diag_message(msg.str());
    if (this->throw_on_notconv_exceed
        && notconv_exceeds_limit(notconv, this->notconv_max)) {
        throw std::runtime_error(msg.str());
    }
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::pgemm_multiply(const T alpha,
                                            const T* a,
                                            const T* b,
                                            const T beta,
                                            T* c)
{
    this->pmmcn.multiply(alpha, a, b, beta, c);
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::plintrans_act(const T alpha,
                                           const T* a,
                                           const T* u,
                                           const T beta,
                                           T* b)
{
    this->plintrans.act(alpha, a, u, beta, b);
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::plintrans_batched_act(const T alpha,
                                                   const T* const* a_blocks,
                                                   const int block_count,
                                                   const T* u,
                                                   const T beta,
                                                   T* const* b_blocks)
{
    if (block_count <= 1) {
        this->plintrans_act(alpha, a_blocks[0], u, beta, b_blocks[0]);
        return;
    }
    const int batch_rows = block_count * this->n_dim;
    const int batch_size = batch_rows * this->n_band_l;
    this->plintrans_batch_in.resize(batch_size);
    this->plintrans_batch_out.resize(batch_size);

    const bool use_beta = beta != static_cast<T>(0.0);
    for (int block = 0; block < block_count; ++block) {
        ModuleBase::matrixCopy<T, Device>()(this->n_band_l,
                                            this->n_dim,
                                            a_blocks[block],
                                            this->n_basis,
                                            this->plintrans_batch_in.data() + block * this->n_dim,
                                            batch_rows);
        if (use_beta) {
            ModuleBase::matrixCopy<T, Device>()(this->n_band_l,
                                                this->n_dim,
                                                b_blocks[block],
                                                this->n_basis,
                                                this->plintrans_batch_out.data() + block * this->n_dim,
                                                batch_rows);
        }
    }
    if (!use_beta) {
        std::fill(this->plintrans_batch_out.begin(), this->plintrans_batch_out.end(), static_cast<T>(0.0));
    }
    if (block_count == 2) {
        this->plintrans_batch2.act(alpha, this->plintrans_batch_in.data(), u, beta, this->plintrans_batch_out.data());
    } else if (block_count == 3) {
        this->plintrans_batch3.act(alpha, this->plintrans_batch_in.data(), u, beta, this->plintrans_batch_out.data());
    } else {
        throw std::runtime_error("LOBPCG batched linear transform supports only 2 or 3 blocks");
    }

    for (int block = 0; block < block_count; ++block) {
        ModuleBase::matrixCopy<T, Device>()(this->n_band_l,
                                            this->n_dim,
                                            this->plintrans_batch_out.data() + block * this->n_dim,
                                            batch_rows,
                                            b_blocks[block],
                                            this->n_basis);
        for (int ib = 0; ib < this->n_band_l; ++ib) {
            T* dst = b_blocks[block] + ib * this->n_basis;
            std::fill(dst + this->n_dim, dst + this->n_basis, static_cast<T>(0.0));
        }
    }
}

// ============================================================================
// lobpcg_update_s — generalized R-R on S-orthonormalized W = [X, Z, P]
// ============================================================================

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::lobpcg_update_s_parallel(
    ct::Tensor& psi, ct::Tensor& hpsi, ct::Tensor& spsi,
    ct::Tensor& grad, ct::Tensor& hgrad, ct::Tensor& sgrad,
    ct::Tensor& pdir, ct::Tensor& hpdir, ct::Tensor& spdir,
    ct::Tensor& eigen,
    const bool force_compressed)
{
    const int n = this->n_band;
    const int nbs = this->n_basis;
    const int local_sz = this->n_band_l * nbs;
    int block_count = this->has_pdir ? 3 : 2;
    int m = block_count * n;

    const T* basis_blocks[3] = {psi.data<T>(), grad.data<T>(), pdir.data<T>()};
    const T* hbasis_blocks[3] = {hpsi.data<T>(), hgrad.data<T>(), hpdir.data<T>()};
    const T* sbasis_blocks[3] = {spsi.data<T>(), sgrad.data<T>(), spdir.data<T>()};

    T* hsub_d = this->hsub.data<T>();
    T* ssub_d = this->ssub.data<T>();

    auto store_block = [=](const T* src, T* dst, const int iblock, const int jblock) {
        ModuleBase::matrixCopy<T, Device>()(n,
                                            n,
                                            src,
                                            n,
                                            dst + jblock * n * this->nsub + iblock * n,
                                            this->nsub);
    };
    auto store_hermitian_block = [=](const T* src, T* dst, const int iblock, const int jblock) {
        store_block(src, dst, iblock, jblock);
        if (iblock == jblock) {
            return;
        }
        for (int jc = 0; jc < n; ++jc) {
            for (int ir = 0; ir < n; ++ir) {
                dst[(iblock * n + jc) * this->nsub + jblock * n + ir]
                    = std::conj(src[ir * n + jc]);
            }
        }
    };

    auto build_subspace = [&](const int active_blocks) {
        setmem_complex_op()(hsub_d, static_cast<T>(0.0), this->nsub * this->nsub);
        setmem_complex_op()(ssub_d, static_cast<T>(0.0), this->nsub * this->nsub);
        for (int jb = 0; jb < active_blocks; ++jb) {
            for (int ib = 0; ib <= jb; ++ib) {
                this->pgemm_multiply(this->one_,
                                     basis_blocks[ib],
                                     hbasis_blocks[jb],
                                     this->zero_,
                                     this->tmp_hsub.data<T>());
                store_hermitian_block(this->tmp_hsub.data<T>(), hsub_d, ib, jb);

                this->pgemm_multiply(this->one_,
                                     basis_blocks[ib],
                                     sbasis_blocks[jb],
                                     this->zero_,
                                     this->tmp_ssub.data<T>());
                store_hermitian_block(this->tmp_ssub.data<T>(), ssub_d, ib, jb);
            }
        }
        hermitize(hsub_d, this->nsub, active_blocks * n);
        hermitize(ssub_d, this->nsub, active_blocks * n);
    };

    build_subspace(block_count);
    std::vector<Real> inv_subspace_norm;
    std::string scale_error;

    auto solve_scaled_hegvd = [&]() -> bool {
        if (!scale_subspace_by_overlap_diag(hsub_d, ssub_d, this->nsub, m,
                                            inv_subspace_norm, scale_error)) {
            return false;
        }
        auto spd_check = check_subspace_spd<T, ct_Device>(ssub_d, this->nsub, m);
        if (!spd_check.ok) {
            scale_error = subspace_spd_diagnostics(spd_check);
            return false;
        }
        try {
            ct::kernels::lapack_hegvd<T, ct_Device>()(
                m, this->nsub, hsub_d, ssub_d,
                this->sub_eigen.data<Real>(), hsub_d);
        } catch (const std::exception& e) {
            scale_error = e.what();
            return false;
        }
        if (!finite_real_block(this->sub_eigen.data<Real>(), m)
            || !finite_scalar_block(hsub_d, m * this->nsub)) {
            scale_error = "hegvd produced non-finite values";
            return false;
        }
        for (int jc = 0; jc < n; ++jc) {
            for (int ir = 0; ir < m; ++ir) {
                hsub_d[jc * this->nsub + ir] *= inv_subspace_norm[ir];
            }
        }
        return true;
    };

    auto solve_compressed_heevd = [&]() {
        std::vector<T> h_scaled(m * m, static_cast<T>(0.0));
        std::vector<T> s_scaled(m * m, static_cast<T>(0.0));
        for (int jc = 0; jc < m; ++jc) {
            for (int ir = 0; ir < m; ++ir) {
                h_scaled[jc * m + ir] = hsub_d[jc * this->nsub + ir];
                s_scaled[jc * m + ir] = ssub_d[jc * this->nsub + ir];
            }
        }

        std::vector<T> s_evec = s_scaled;
        std::vector<Real> s_eval(m, static_cast<Real>(0.0));
        ct::kernels::lapack_heevd<T, ct_Device>()(m, s_evec.data(), m, s_eval.data());
        if (!finite_real_block(s_eval.data(), m) || !finite_scalar_block(s_evec.data(), m * m)) {
            throw std::runtime_error("LOBPCG generalized parallel overlap diagonalization produced non-finite values");
        }

        const Real s_max = s_eval.empty() ? static_cast<Real>(0.0)
                                          : *std::max_element(s_eval.begin(), s_eval.end());
        const Real eps = std::numeric_limits<Real>::epsilon();
        const Real rank_floor = std::max(std::abs(s_max) * static_cast<Real>(1.0e-10),
                                         static_cast<Real>(100.0) * eps * std::max(static_cast<Real>(1.0),
                                                                                   std::abs(s_max)));
        int first_kept = 0;
        while (first_kept < m && s_eval[first_kept] <= rank_floor) {
            ++first_kept;
        }
        const int rank = m - first_kept;
        if (rank < n) {
            throw std::runtime_error("LOBPCG generalized parallel compressed subspace lost rank");
        }

        std::vector<T> q(m * rank, static_cast<T>(0.0));
        for (int jc = 0; jc < rank; ++jc) {
            const int src_col = first_kept + jc;
            const Real inv_sqrt = static_cast<Real>(1.0) / std::sqrt(s_eval[src_col]);
            for (int ir = 0; ir < m; ++ir) {
                q[jc * m + ir] = s_evec[src_col * m + ir] * inv_sqrt;
            }
        }

        std::vector<T> hq(m * rank, static_cast<T>(0.0));
        ModuleBase::gemm_op<T, Device>()('N', 'N',
                                         m, rank, m,
                                         this->one,
                                         h_scaled.data(), m,
                                         q.data(), m,
                                         this->zero,
                                         hq.data(), m);
        std::vector<T> h_comp(rank * rank, static_cast<T>(0.0));
        ModuleBase::gemm_op<T, Device>()('C', 'N',
                                         rank, rank, m,
                                         this->one,
                                         q.data(), m,
                                         hq.data(), m,
                                         this->zero,
                                         h_comp.data(), rank);
        hermitize(h_comp.data(), rank, rank);

        ct::kernels::lapack_heevd<T, ct_Device>()(
            rank, h_comp.data(), rank, this->sub_eigen.data<Real>());
        if (!finite_real_block(this->sub_eigen.data<Real>(), rank)
            || !finite_scalar_block(h_comp.data(), rank * rank)) {
            throw std::runtime_error("LOBPCG generalized parallel compressed diagonalization produced non-finite values");
        }

        std::vector<T> coeff_scaled(m * n, static_cast<T>(0.0));
        ModuleBase::gemm_op<T, Device>()('N', 'N',
                                         m, n, rank,
                                         this->one,
                                         q.data(), m,
                                         h_comp.data(), rank,
                                         this->zero,
                                         coeff_scaled.data(), m);
        setmem_complex_op()(hsub_d, static_cast<T>(0.0), this->nsub * this->nsub);
        ModuleBase::matrixCopy<T, Device>()(n, m, coeff_scaled.data(), m, hsub_d, this->nsub);
    };

    if (force_compressed || !solve_scaled_hegvd()) {
        bool solved = false;
        if (!force_compressed && block_count == 3) {
            block_count = 2;
            m = block_count * n;
            build_subspace(block_count);
            scale_error.clear();
            solved = solve_scaled_hegvd();
        }
        if (!solved) {
            build_subspace(block_count);
            try {
                solve_compressed_heevd();
            } catch (const std::exception& e) {
                this->diag_log("lobpcg_update_s_parallel compressed fallback failed: " + std::string(e.what()),
                               hermitian_matrix_diagnostics("H_sub", hsub_d, this->nsub, m),
                               hermitian_matrix_diagnostics("S_sub", ssub_d, this->nsub, m),
                               s_overlap_diagnostics(ssub_d, this->nsub, m)
                                   + ", scaled-path reason=" + scale_error);
                throw;
            }
        }
    }

    copy_lowest_subspace_eigenvalues(this->sub_eigen.data<Real>(), eigen.data<Real>(), n);

    auto buffers = make_generalized_update_buffers(this->work.data<T>(),
                                                   this->hwork.data<T>(),
                                                   this->swork.data<T>(),
                                                   this->pwork.data<T>(),
                                                   this->hpwork.data<T>(),
                                                   this->spwork.data<T>());
    zero_update_buffers<T, Device>(buffers, local_sz);

    auto copy_coeff_block = [=](const int block, T* coeff) {
        setmem_complex_op()(coeff, static_cast<T>(0.0), n * n);
        ModuleBase::matrixCopy<T, Device>()(n,
                                            n,
                                            hsub_d + block * n,
                                            this->nsub,
                                            coeff,
                                            n);
    };

    for (int ib = 0; ib < block_count; ++ib) {
        copy_coeff_block(ib, this->tmp_hsub.data<T>());
        const T* update_in[3] = {basis_blocks[ib], hbasis_blocks[ib], sbasis_blocks[ib]};
        if (ib == 0) {
            T* update_out[3] = {buffers.x, buffers.hx, buffers.sx};
            this->plintrans_batched_act(this->one_, update_in, 3,
                                        this->tmp_hsub.data<T>(),
                                        this->zero_,
                                        update_out);
        } else {
            T* update_out[3] = {buffers.p, buffers.hp, buffers.sp};
            this->plintrans_batched_act(this->one_, update_in, 3,
                                        this->tmp_hsub.data<T>(),
                                        ib == 1 ? this->zero_ : this->one_,
                                        update_out);
        }
    }
    add_block_to<T, Device>(buffers.p, buffers.x, local_sz, this->one);
    add_block_to<T, Device>(buffers.hp, buffers.hx, local_sz, this->one);
    add_block_to<T, Device>(buffers.sp, buffers.sx, local_sz, this->one);

    sync_generalized_update<T, Device>(psi, hpsi, spsi, pdir, hpdir, spdir, buffers, local_sz);
    this->ensure_generalized_update_finite(
        psi, hpsi, spsi, eigen, nbs, n,
        "lobpcg_update_s_parallel produced non-finite values",
        "LOBPCG generalized parallel update produced non-finite values");

    this->has_pdir = true;
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::lobpcg_update_s(
    ct::Tensor& psi, ct::Tensor& hpsi, ct::Tensor& spsi,
    ct::Tensor& grad, ct::Tensor& hgrad, ct::Tensor& sgrad,
    ct::Tensor& pdir, ct::Tensor& hpdir, ct::Tensor& spdir,
    ct::Tensor& eigen)
{
    const int n    = this->n_band;
    const int nbs  = this->n_basis;
    const int nvalid = this->n_dim;
    const int local_sz = this->n_band_l * nbs;
    const Real eps = static_cast<Real>(100)
                   * std::numeric_limits<Real>::epsilon();

    std::vector<T> basis;
    std::vector<T> hbasis;
    std::vector<T> sbasis;
    const int block_capacity = this->has_pdir ? 3 : 2;
    basis.reserve(block_capacity * local_sz);
    hbasis.reserve(block_capacity * local_sz);
    sbasis.reserve(block_capacity * local_sz);

    append_s_orthonormal_block<T>(n, nbs, nvalid, eps,
                                  psi.data<T>(), hpsi.data<T>(), spsi.data<T>(),
                                  basis, hbasis, sbasis);
    append_s_orthonormal_block<T>(n, nbs, nvalid, eps,
                                  grad.data<T>(), hgrad.data<T>(), sgrad.data<T>(),
                                  basis, hbasis, sbasis);
    if (this->has_pdir) {
        append_s_orthonormal_block<T>(n, nbs, nvalid, eps,
                                      pdir.data<T>(), hpdir.data<T>(), spdir.data<T>(),
                                      basis, hbasis, sbasis);
    }

    const int m = static_cast<int>(basis.size() / nbs);
    if (m < n) {
        throw std::runtime_error("LOBPCG generalized subspace lost rank");
    }

    T* hsub_d = this->hsub.data<T>();
    T* ssub_d = this->ssub.data<T>();
    setmem_complex_op()(hsub_d, static_cast<T>(0.0), this->nsub * this->nsub);
    setmem_complex_op()(ssub_d, static_cast<T>(0.0), this->nsub * this->nsub);

    ModuleBase::gemm_op<T, Device>()('C', 'N',
                                     m, m, nvalid,
                                     this->one,
                                     basis.data(), nbs,
                                     hbasis.data(), nbs,
                                     this->zero,
                                     hsub_d, this->nsub);
    ModuleBase::gemm_op<T, Device>()('C', 'N',
                                     m, m, nvalid,
                                     this->one,
                                     basis.data(), nbs,
                                     sbasis.data(), nbs,
                                     this->zero,
                                     ssub_d, this->nsub);
#ifdef __MPI
    for (int jc = 0; jc < m; ++jc) {
        Parallel_Reduce::reduce_pool(hsub_d + jc * this->nsub, m);
        Parallel_Reduce::reduce_pool(ssub_d + jc * this->nsub, m);
    }
#endif
    hermitize(hsub_d, this->nsub, m);
    hermitize(ssub_d, this->nsub, m);

    try {
        ct::kernels::lapack_hegvd<T, ct_Device>()(
            m, this->nsub, hsub_d, ssub_d,
            this->sub_eigen.data<Real>(), hsub_d);
    } catch (const std::exception& e) {
        this->diag_log("lobpcg_update_s hegvd failed: " + std::string(e.what()),
                        hermitian_matrix_diagnostics("H_sub", hsub_d, this->nsub, m),
                        hermitian_matrix_diagnostics("S_sub", ssub_d, this->nsub, m),
                        s_overlap_diagnostics(ssub_d, this->nsub, m));
        throw;
    }

    if (!finite_real_block(this->sub_eigen.data<Real>(), m)
        || !finite_scalar_block(hsub_d, m * this->nsub)) {
        this->diag_log("lobpcg_update_s hegvd produced non-finite values",
                        hermitian_matrix_diagnostics("H_sub", hsub_d, this->nsub, m),
                        hermitian_matrix_diagnostics("S_sub", ssub_d, this->nsub, m),
                        s_overlap_diagnostics(ssub_d, this->nsub, m));
        throw std::runtime_error("LOBPCG generalized subspace diagonalization produced non-finite values");
    }

    copy_lowest_subspace_eigenvalues(this->sub_eigen.data<Real>(), eigen.data<Real>(), n);

    auto buffers = make_generalized_update_buffers(this->work.data<T>(),
                                                   this->hwork.data<T>(),
                                                   this->swork.data<T>(),
                                                   this->pwork.data<T>(),
                                                   this->hpwork.data<T>(),
                                                   this->spwork.data<T>());
    zero_update_buffers<T, Device>(buffers, local_sz);

    ModuleBase::gemm_op<T, Device>()('N', 'N',
                                     nvalid, n, n,
                                     this->one,
                                     basis.data(), nbs,
                                     hsub_d, this->nsub,
                                     this->zero,
                                     buffers.x, nbs);
    ModuleBase::gemm_op<T, Device>()('N', 'N',
                                     nvalid, n, n,
                                     this->one,
                                     hbasis.data(), nbs,
                                     hsub_d, this->nsub,
                                     this->zero,
                                     buffers.hx, nbs);
    ModuleBase::gemm_op<T, Device>()('N', 'N',
                                     nvalid, n, n,
                                     this->one,
                                     sbasis.data(), nbs,
                                     hsub_d, this->nsub,
                                     this->zero,
                                     buffers.sx, nbs);

    const int tail_cols = m - n;
    if (tail_cols > 0) {
        ModuleBase::gemm_op<T, Device>()('N', 'N',
                                         nvalid, n, tail_cols,
                                         this->one,
                                         basis.data() + n * nbs, nbs,
                                         hsub_d + n, this->nsub,
                                         this->zero,
                                         buffers.p, nbs);
        ModuleBase::gemm_op<T, Device>()('N', 'N',
                                         nvalid, n, tail_cols,
                                         this->one,
                                         hbasis.data() + n * nbs, nbs,
                                         hsub_d + n, this->nsub,
                                         this->zero,
                                         buffers.hp, nbs);
        ModuleBase::gemm_op<T, Device>()('N', 'N',
                                         nvalid, n, tail_cols,
                                         this->one,
                                         sbasis.data() + n * nbs, nbs,
                                         hsub_d + n, this->nsub,
                                         this->zero,
                                         buffers.sp, nbs);
        add_block_to<T, Device>(buffers.p, buffers.x, local_sz, this->one);
        add_block_to<T, Device>(buffers.hp, buffers.hx, local_sz, this->one);
        add_block_to<T, Device>(buffers.sp, buffers.sx, local_sz, this->one);
    }

    sync_generalized_update<T, Device>(psi, hpsi, spsi, pdir, hpdir, spdir, buffers, local_sz);
    this->ensure_generalized_update_finite(
        psi, hpsi, spsi, eigen, nbs, n,
        "lobpcg_update_s produced non-finite values",
        "LOBPCG generalized update produced non-finite values");

    this->has_pdir = true;
}


template <typename T, typename Device>
int DiagoLobpcg<T, Device>::diag(
    const HPsiFunc& hpsi_func, const SPsiFunc& spsi_func, T* psi_in,
    Real* eigenvalue_in, const std::vector<double>& ethr_band)
{
    this->validate_ethr_band(ethr_band);

    this->has_pdir = false;
    this->psi = ct::TensorMap(psi_in, t_type, dev_type,
                               {this->n_band_l, this->n_basis});

    this->calc_spsi_with_block(spsi_func, psi_in, this->spsi);

    const int scf_iter = DiagoIterAssist<T, Device>::SCF_ITER;

    this->calc_prec();
    if (this->n_band_l != this->n_band) {
        this->calc_hpsi_with_block(hpsi_func, this->psi.data<T>(), this->hpsi);
        this->generalized_rayleigh_ritz_parallel(this->psi, this->hpsi, this->spsi, this->eigen);
    } else {
        this->repair_initial_subspace_s(hpsi_func, spsi_func);
        this->generalized_rayleigh_ritz(this->psi, this->hpsi, this->spsi, this->eigen);
    }

    this->clear_search_directions();

    std::vector<double> effective_ethr_band = ethr_band;
    if (this->notconv_max < 0) {
        // SCF can refine the density across outer iterations; avoid chasing a tiny diagonalization threshold.
        constexpr double scf_generalized_residual_floor = 1.0e-5;
        for (double& ethr : effective_ethr_band) {
            ethr = std::max(ethr, scf_generalized_residual_floor);
        }
    }
    const int used_iter = this->run_lobpcg_loop(
        LOBPCG_PROBLEM_GENERALIZED,
        "DiagoLobpcg::diag(S!=I)",
        hpsi_func,
        spsi_func,
        effective_ethr_band,
        scf_iter);

    syncmem_var_d2h_op()(eigenvalue_in,
                          this->eigen.data<Real>() + this->local_band_start(),
                          this->n_band_l);
    return used_iter;
}

template class DiagoLobpcg<std::complex<double>, base_device::DEVICE_CPU>;

} // namespace hsolver
