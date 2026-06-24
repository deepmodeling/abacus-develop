#include "source_hsolver/diago_lobpcg.h"

#include "source_base/kernels/math_kernel_op.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_reduce.h"
#include "source_hsolver/diago_iter_assist.h"
#include <ATen/core/tensor_map.h>
#include <ATen/kernels/lapack.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace hsolver {

namespace detail {

static void lobpcg_reduce_sum(int& value)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_SUM, BP_WORLD);
#endif
}

static void lobpcg_reduce_max(double& value)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_DOUBLE, MPI_MAX, BP_WORLD);
#endif
}

static void lobpcg_reduce_bool_or(bool& value)
{
#ifdef __MPI
    int reduced = value ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &reduced, 1, MPI_INT, MPI_LOR, BP_WORLD);
    value = (reduced != 0);
#endif
}

static void lobpcg_reduce_bool_and(bool& value)
{
#ifdef __MPI
    int reduced = value ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &reduced, 1, MPI_INT, MPI_LAND, BP_WORLD);
    value = (reduced != 0);
#endif
}

template <typename T>
static void hermitize(T* mat, int ld, int active_sub)
{
    using Real = typename GetTypeReal<T>::type;
    for (int i = 0; i < active_sub; i++) {
        const int idx = i * ld + i;
        mat[idx] = T(std::real(mat[idx]), static_cast<Real>(0));
    }
    for (int jc = 0; jc < active_sub; ++jc) {
        for (int ir = jc + 1; ir < active_sub; ++ir) {
            const T avg = static_cast<T>(0.5)
                        * (mat[jc * ld + ir]
                        +  std::conj(mat[ir * ld + jc]));
            mat[jc * ld + ir] = avg;
            mat[ir * ld + jc] = std::conj(avg);
        }
    }
}

template <typename Real>
static bool finite_real_block(const Real* data, int n)
{
    for (int i = 0; i < n; ++i) {
        if (!std::isfinite(data[i])) {
            return false;
        }
    }
    return true;
}

template <typename T>
static bool finite_scalar_block(const T* data, int n)
{
    for (int i = 0; i < n; ++i) {
        if (!std::isfinite(std::real(data[i])) || !std::isfinite(std::imag(data[i]))) {
            return false;
        }
    }
    return true;
}

template <typename T, typename CtDevice>
static bool check_subspace_spd(const T* ssub, int ld, int n)
{
    using Real = typename GetTypeReal<T>::type;
    if (n <= 0) {
        return false;
    }

    std::vector<T> smat(n * n, static_cast<T>(0.0));
    std::vector<Real> eval(n, static_cast<Real>(0.0));
    for (int jc = 0; jc < n; ++jc) {
        for (int ir = 0; ir < n; ++ir) {
            smat[jc * n + ir] = ssub[jc * ld + ir];
        }
    }

    try {
        ct::kernels::lapack_heevd<T, CtDevice>()(n, smat.data(), n, eval.data());
    } catch (const std::exception&) {
        return false;
    }

    if (!finite_real_block(eval.data(), n)) {
        return false;
    }

    const Real min_eval = eval.front();
    const Real max_eval = eval.back();
    const Real eps = std::numeric_limits<Real>::epsilon();
    const Real cond_limit = static_cast<Real>(1.0e10);
    const Real scale = std::max(static_cast<Real>(1.0), std::abs(max_eval));
    const Real floor = std::max(scale / cond_limit,
                                static_cast<Real>(100.0) * eps * scale);

    if (!std::isfinite(min_eval) || !std::isfinite(max_eval)
        || max_eval <= static_cast<Real>(0.0)
        || min_eval <= floor) {
        return false;
    }

    const Real cond = max_eval / min_eval;
    return std::isfinite(cond) && cond <= cond_limit;
}

template <typename T>
static bool scale_subspace_by_overlap_diag(T* hsub,
                                           T* ssub,
                                           int ld,
                                           int n,
                                           std::vector<typename GetTypeReal<T>::type>& inv_norm)
{
    using Real = typename GetTypeReal<T>::type;
    inv_norm.assign(n, static_cast<Real>(0.0));
    const Real diag_floor = std::numeric_limits<Real>::min();

    for (int i = 0; i < n; ++i) {
        const T diag = ssub[i * ld + i];
        const Real diag_real = static_cast<Real>(std::real(diag));
        if (!std::isfinite(diag_real) || diag_real <= diag_floor) {
            return false;
        }
        inv_norm[i] = static_cast<Real>(1.0) / std::sqrt(diag_real);
    }

    for (int jc = 0; jc < n; ++jc) {
        for (int ir = 0; ir < n; ++ir) {
            const Real scale = inv_norm[ir] * inv_norm[jc];
            hsub[jc * ld + ir] *= scale;
            ssub[jc * ld + ir] *= scale;
        }
    }
    hermitize(hsub, ld, n);
    hermitize(ssub, ld, n);
    return true;
}

template <typename T>
static bool finite_vector_block(const T* data, int nvec, int lda, int nvalid)
{
    for (int ib = 0; ib < nvec; ++ib) {
        for (int ig = 0; ig < nvalid; ++ig) {
            const T value = data[ib * lda + ig];
            if (!std::isfinite(std::real(value)) || !std::isfinite(std::imag(value))) {
                return false;
            }
        }
    }
    return true;
}

template <typename T>
static void append_s_orthonormal_block(
    const int nvec, const int lda, const int nvalid,
    const typename GetTypeReal<T>::type thresh,
    const T* block, const T* hblock, const T* sblock,
    std::vector<T>& basis, std::vector<T>& hbasis, std::vector<T>& sbasis)
{
    using Real = typename GetTypeReal<T>::type;

    for (int ib = 0; ib < nvec; ++ib) {
        std::vector<T> q(lda, static_cast<T>(0.0));
        std::vector<T> hq(lda, static_cast<T>(0.0));
        std::vector<T> sq(lda, static_cast<T>(0.0));
        for (int ig = 0; ig < nvalid; ++ig) {
            q[ig] = block[ib * lda + ig];
            hq[ig] = hblock[ib * lda + ig];
            sq[ig] = sblock[ib * lda + ig];
        }

        const int nold = static_cast<int>(basis.size() / lda);
        for (int pass = 0; pass < 2; ++pass) {
            for (int jq = 0; jq < nold; ++jq) {
                T dot = static_cast<T>(0.0);
                for (int ig = 0; ig < nvalid; ++ig) {
                    dot += std::conj(basis[jq * lda + ig]) * sq[ig];
                }
#ifdef __MPI
                Parallel_Reduce::reduce_pool(&dot, 1);
#endif
                for (int ig = 0; ig < nvalid; ++ig) {
                    q[ig] -= dot * basis[jq * lda + ig];
                    hq[ig] -= dot * hbasis[jq * lda + ig];
                    sq[ig] -= dot * sbasis[jq * lda + ig];
                }
            }
        }

        Real norm2 = static_cast<Real>(0.0);
        for (int ig = 0; ig < nvalid; ++ig) {
            norm2 += std::real(std::conj(q[ig]) * sq[ig]);
        }
#ifdef __MPI
        Parallel_Reduce::reduce_pool(&norm2, 1);
#endif
        if (!std::isfinite(norm2) || norm2 <= thresh) {
            continue;
        }

        const Real inv_norm = static_cast<Real>(1.0) / std::sqrt(norm2);
        for (int ig = 0; ig < lda; ++ig) {
            basis.push_back(q[ig] * inv_norm);
            hbasis.push_back(hq[ig] * inv_norm);
            sbasis.push_back(sq[ig] * inv_norm);
        }
    }
}

} // namespace detail

template <typename T, typename Device>
bool DiagoLobpcg<T, Device>::update_best_state(State& state,
                                                StateQuality& quality,
                                                const int candidate_notconv,
                                                const Real candidate_residual)
{
    if (!std::isfinite(candidate_residual)) {
        return false;
    }
    bool better = !quality.valid || candidate_notconv < quality.notconv;
    if (quality.valid && candidate_notconv == quality.notconv) {
        const Real best_tol = std::max(static_cast<Real>(1.0e-12),
                                       std::abs(quality.residual) * static_cast<Real>(1.0e-8));
        better = candidate_residual < quality.residual - best_tol;
    }
    if (!better) {
        return false;
    }

    quality.residual = candidate_residual;
    quality.notconv = candidate_notconv;
    quality.valid = true;
    this->save_state(state);
    return true;
}

template <typename T, typename Device>
int DiagoLobpcg<T, Device>::run_lobpcg_loop(
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
            this->notconv_max >= 0 ? notconv_before_update > this->notconv_max
                                   : notconv_before_update > 0;
        this->update_best_state(best_state,
                                best_quality,
                                notconv_before_update,
                                residual_before_update);
        std::vector<char> soft_lock_mask;
        if (this->n_band_l != this->n_band
            && notconv_before_update > 0
            && notconv_before_update < this->n_band) {
            const Real* err_d = this->err_st.data<Real>();
            soft_lock_mask.resize(this->n_band_l, 0);
            for (int ib = 0; ib < this->n_band_l; ++ib) {
                soft_lock_mask[ib] = err_d[ib] <= static_cast<Real>(effective_ethr_band[ib]) ? 1 : 0;
            }
        }
        if (!continues_for_notconv) {
            break;
        }

        this->calc_spsi_with_block(spsi_func, this->grad.data<T>(), this->sgrad);
        this->orth_projection_s(this->psi, this->spsi, this->tmp_hsub,
                                this->sgrad, this->grad);
        this->calc_hpsi_with_block(hpsi_func, this->grad.data<T>(), this->hgrad);
        this->save_state(rollback_state);

        try {
            this->update_generalized_subspace(false);
        } catch (const std::exception&) {
            this->restore_state(rollback_state, true);
            try {
                this->update_generalized_subspace(false);
            } catch (const std::exception&) {
                this->restore_state(rollback_state, false);
                this->calc_hpsi_with_block(hpsi_func, this->psi.data<T>(), this->hpsi);
                this->calc_spsi_with_block(spsi_func, this->psi.data<T>(), this->spsi);
                if (this->n_band_l != this->n_band) {
                    this->generalized_rayleigh_ritz_parallel(this->psi, this->hpsi, this->spsi, this->eigen);
                } else {
                    this->generalized_rayleigh_ritz(this->psi, this->hpsi, this->spsi, this->eigen);
                }
            }
        }

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
        const Real residual_growth_limit = static_cast<Real>(10.0);
        const Real residual_limit = std::max(static_cast<Real>(1.0e-8),
                                             residual_before_update * residual_growth_limit);
        update_rejected = !std::isfinite(residual_after_update)
            || (residual_after_update > residual_limit
                && notconv_after_update >= notconv_before_update);
        detail::lobpcg_reduce_bool_or(update_rejected);
        bool stop_after_rejected_update = false;
        if (update_rejected) {
            stop_after_rejected_update = this->handle_generalized_rejected_update(
                rollback_state,
                best_state,
                best_quality,
                update_rejected,
                notconv_before_update,
                residual_before_update,
                residual_growth_limit,
                effective_ethr_band);
        } else {
            this->update_best_state(best_state, best_quality,
                                    notconv_after_update, residual_after_update);
        }
        if (stop_after_rejected_update) {
            stop_reason = "stopped after residual-guard rollback";
            break;
        }

        const bool has_next_iteration = (ntry + 1) < max_iter;
        const bool restart_next = has_next_iteration && scf_iter == 1 && ((ntry + 1) % this->nline == 0);
        if (has_next_iteration && !restart_next && !update_rejected) {
            try {
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
    const Real best_restore_tol = std::max(static_cast<Real>(1.0e-12),
                                           std::abs(best_quality.residual) * static_cast<Real>(1.0e-8));
    const bool should_restore_best =
        best_quality.valid
        && (final_notconv > best_quality.notconv
            || !std::isfinite(final_residual)
            || (final_notconv == best_quality.notconv
                && final_residual > best_quality.residual + best_restore_tol));
    if (should_restore_best) {
        this->restore_state(best_state, true);
        compute_residual();
        final_residual = this->max_error(this->err_st);
        final_notconv = this->count_not_converged(this->err_st, effective_ethr_band);
    }
    if (stop_reason.empty()
        && !(this->notconv_max >= 0 ? final_notconv > this->notconv_max
                                    : final_notconv > 0)) {
        stop_reason = "stopped with allowed unconverged bands";
    }
    this->report_not_converged(used_iter,
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
    const int notconv_before_update,
    const Real residual_before_update,
    const Real residual_growth_limit,
    const std::vector<double>& effective_ethr_band)
{
    auto restore_backup_state = [&]() {
        this->restore_state(rollback_state, true);
    };
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
            const Real guarded_limit = std::max(static_cast<Real>(1.0e-8),
                                                residual_before_update * residual_growth_limit);
            compressed_ok = std::isfinite(guarded_residual)
                && (guarded_residual <= guarded_limit || guarded_notconv < notconv_before_update);
            detail::lobpcg_reduce_bool_and(compressed_ok);
        } catch (const std::exception&) {
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
    T* blocks[3] = {this->pdir.data<T>(), this->hpdir.data<T>(), this->spdir.data<T>()};
    for (T* block : blocks)
        std::fill(block, block + psi_sz, static_cast<T>(0.0));
    this->has_pdir = false;
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::save_state(State& state)
{
    const int psi_sz = this->n_basis * this->n_band_l;
    state.psi.assign(this->psi.data<T>(), this->psi.data<T>() + psi_sz);
    state.hpsi.assign(this->hpsi.data<T>(), this->hpsi.data<T>() + psi_sz);
    state.spsi.assign(this->spsi.data<T>(), this->spsi.data<T>() + psi_sz);
    state.eigen.assign(this->eigen.data<Real>(), this->eigen.data<Real>() + this->n_band);
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::restore_state(const State& state,
                                            const bool reset_search)
{
    const int psi_sz = this->n_basis * this->n_band_l;
    std::copy(state.psi.data(), state.psi.data() + psi_sz, this->psi.data<T>());
    std::copy(state.hpsi.data(), state.hpsi.data() + psi_sz, this->hpsi.data<T>());
    std::copy(state.spsi.data(), state.spsi.data() + psi_sz, this->spsi.data<T>());
    std::copy(state.eigen.data(), state.eigen.data() + this->n_band, this->eigen.data<Real>());
    if (reset_search) {
        this->clear_search_directions();
    }
}

template <typename T, typename Device>
bool DiagoLobpcg<T, Device>::restore_soft_locked_bands(
    const State& state,
    const std::vector<char>& soft_lock_mask,
    const ct::Tensor& err_in,
    const std::vector<double>& ethr_band)
{
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
    detail::lobpcg_reduce_sum(restored_global);
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
                                                               const char* error_message) const
{
    bool invalid_update =
        !detail::finite_vector_block(psi.data<T>(), this->n_band_l, this->n_basis, this->n_dim)
        || !detail::finite_vector_block(hpsi.data<T>(), this->n_band_l, this->n_basis, this->n_dim)
        || !detail::finite_vector_block(spsi.data<T>(), this->n_band_l, this->n_basis, this->n_dim)
        || !detail::finite_real_block(eigen.data<Real>(), this->n_band);
    detail::lobpcg_reduce_bool_or(invalid_update);
    if (!invalid_update) {
        return;
    }
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
void DiagoLobpcg<T, Device>::calc_hpsi_with_block(
    const HPsiFunc& hpsi_func, T* psi_in, ct::Tensor& hpsi_out)
{
    hpsi_func(psi_in, hpsi_out.data<T>(), this->n_basis, this->n_band_l);
    bool invalid_hpsi =
        !detail::finite_vector_block(psi_in, this->n_band_l, this->n_basis, this->n_dim)
        || !detail::finite_vector_block(hpsi_out.data<T>(), this->n_band_l, this->n_basis, this->n_dim);
    detail::lobpcg_reduce_bool_or(invalid_hpsi);
    if (invalid_hpsi) {
        std::string message = "LOBPCG hPsi produced non-finite values";
        if (!this->diag_context.empty()) {
            message += "; context={" + this->diag_context + "}";
        }
        throw std::runtime_error(message);
    }
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::calc_spsi_with_block(
    const SPsiFunc& spsi_func, const T* psi_in, ct::Tensor& spsi_out)
{
    spsi_func(psi_in, spsi_out.data<T>(), this->n_basis, this->n_band_l);
    bool invalid_spsi =
        !detail::finite_vector_block(psi_in, this->n_band_l, this->n_basis, this->n_dim)
        || !detail::finite_vector_block(spsi_out.data<T>(), this->n_band_l, this->n_basis, this->n_dim);
    detail::lobpcg_reduce_bool_or(invalid_spsi);
    if (invalid_spsi) {
        std::string message = "LOBPCG sPsi produced non-finite values";
        if (!this->diag_context.empty()) {
            message += "; context={" + this->diag_context + "}";
        }
        throw std::runtime_error(message);
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
                std::fill(psi_d + ib * lda, psi_d + (ib + 1) * lda, static_cast<T>(0.0));
                std::fill(spsi_d + ib * lda, spsi_d + (ib + 1) * lda, static_cast<T>(0.0));
                psi_d[ib * lda + ((ib + attempt) % nvalid)] = static_cast<T>(1.0);
                spsi_func(psi_d + ib * lda, spsi_d + ib * lda, lda, 1);
            }

            const bool finite_vec =
                detail::finite_vector_block(psi_d + ib * lda, 1, lda, nvalid)
                && detail::finite_vector_block(spsi_d + ib * lda, 1, lda, nvalid);
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
            std::fill(psi_d + ib * lda + nvalid, psi_d + (ib + 1) * lda, static_cast<T>(0.0));
            std::fill(spsi_d + ib * lda + nvalid, spsi_d + (ib + 1) * lda, static_cast<T>(0.0));
            ready = true;
        }

        if (!ready) {
            throw std::runtime_error("LOBPCG failed to repair rank-deficient initial subspace");
        }
    }

    this->calc_hpsi_with_block(hpsi_func, this->psi.data<T>(), this->hpsi);
    if (repaired) {
        this->calc_spsi_with_block(spsi_func, this->psi.data<T>(), this->spsi);
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

    this->pmmcn.multiply(this->one_, psi_inout.data<T>(), hpsi_inout.data<T>(),
                         this->zero_, this->tmp_hsub.data<T>());
    this->pmmcn.multiply(this->one_, psi_inout.data<T>(), spsi_inout.data<T>(),
                         this->zero_, this->tmp_ssub.data<T>());
    detail::hermitize(this->tmp_hsub.data<T>(), nb, nb);
    detail::hermitize(this->tmp_ssub.data<T>(), nb, nb);

    try {
        ct::kernels::lapack_hegvd<T, ct_Device>()(
            nb, nb,
            this->tmp_hsub.data<T>(),
            this->tmp_ssub.data<T>(),
            eigen_out.data<Real>(),
            this->hsub.data<T>());
    } catch (const std::exception&) {
        throw std::runtime_error("LOBPCG generalized Rayleigh-Ritz failed in hegvd");
    }

    if (!detail::finite_real_block(eigen_out.data<Real>(), nb)
        || !detail::finite_scalar_block(this->hsub.data<T>(), nb * nb)) {
        throw std::runtime_error("LOBPCG generalized Rayleigh-Ritz failed in hegvd");
    }

    const T* rotate_in[3] = {psi_inout.data<T>(), hpsi_inout.data<T>(), spsi_inout.data<T>()};
    T* rotate_out[3] = {this->work.data<T>(), this->hwork.data<T>(), this->swork.data<T>()};
    this->plintrans_batched_act(this->one_, rotate_in, 3,
                                this->hsub.data<T>(), this->zero_, rotate_out);
    std::copy(this->work.data<T>(), this->work.data<T>() + local_sz, psi_inout.data<T>());
    std::copy(this->hwork.data<T>(), this->hwork.data<T>() + local_sz, hpsi_inout.data<T>());
    std::copy(this->swork.data<T>(), this->swork.data<T>() + local_sz, spsi_inout.data<T>());

    if (!detail::finite_vector_block(psi_inout.data<T>(), nb, nbs, nvalid)
        || !detail::finite_vector_block(hpsi_inout.data<T>(), nb, nbs, nvalid)
        || !detail::finite_vector_block(spsi_inout.data<T>(), nb, nbs, nvalid)) {
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
    std::fill(hsub_d, hsub_d + this->nsub * this->nsub, static_cast<T>(0.0));
    std::fill(ssub_d, ssub_d + this->nsub * this->nsub, static_cast<T>(0.0));

    this->pmmcn.multiply(this->one_, psi_inout.data<T>(), hpsi_inout.data<T>(),
                         this->zero_, this->tmp_hsub.data<T>());
    this->pmmcn.multiply(this->one_, psi_inout.data<T>(), spsi_inout.data<T>(),
                         this->zero_, this->tmp_ssub.data<T>());
    ModuleBase::matrixCopy<T, Device>()(nb, nb, this->tmp_hsub.data<T>(), nb, hsub_d, this->nsub);
    ModuleBase::matrixCopy<T, Device>()(nb, nb, this->tmp_ssub.data<T>(), nb, ssub_d, this->nsub);
    detail::hermitize(hsub_d, this->nsub, nb);
    detail::hermitize(ssub_d, this->nsub, nb);

    std::vector<Real> inv_subspace_norm;
    if (!detail::scale_subspace_by_overlap_diag(hsub_d, ssub_d, this->nsub, nb,
                                        inv_subspace_norm)) {
        throw std::runtime_error("LOBPCG generalized parallel initial overlap has invalid diagonal");
    }

    if (!detail::check_subspace_spd<T, ct_Device>(ssub_d, this->nsub, nb)) {
        throw std::runtime_error("LOBPCG generalized parallel initial overlap is ill-conditioned before hegvd");
    }

    ct::kernels::lapack_hegvd<T, ct_Device>()(
        nb, this->nsub, hsub_d, ssub_d, eigen_out.data<Real>(), hsub_d);

    if (!detail::finite_real_block(eigen_out.data<Real>(), nb)
        || !detail::finite_scalar_block(hsub_d, nb * this->nsub)) {
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
    std::copy(this->work.data<T>(), this->work.data<T>() + local_sz, psi_inout.data<T>());
    std::copy(this->hwork.data<T>(), this->hwork.data<T>() + local_sz, hpsi_inout.data<T>());
    std::copy(this->swork.data<T>(), this->swork.data<T>() + local_sz, spsi_inout.data<T>());

    if (!detail::finite_vector_block(psi_inout.data<T>(), this->n_band_l, nbs, this->n_dim)
        || !detail::finite_vector_block(hpsi_inout.data<T>(), this->n_band_l, nbs, this->n_dim)
        || !detail::finite_vector_block(spsi_inout.data<T>(), this->n_band_l, nbs, this->n_dim)) {
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
    bool invalid_grad = false;

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
                invalid_grad = true;
                _grad[idx] = static_cast<T>(0.0);
            } else {
                err_j += std::norm(r);
            }
        }
        for (int ig = this->n_dim; ig < this->n_basis; ig++)
            _grad[ioff + ig] = static_cast<T>(0.0);
        _err[ib] = invalid_grad ? static_cast<Real>(0.0) : err_j;
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(_err, this->n_band_l);
#endif
    detail::lobpcg_reduce_bool_or(invalid_grad);
    if (invalid_grad) {
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

    this->pmmcn.multiply(this->one_, psi_in.data<T>(), sgrad_out.data<T>(),
                         this->zero_, hsub_work.data<T>());
    const T* proj_in[2] = {psi_in.data<T>(), spsi_in.data<T>()};
    T* proj_out[2] = {grad_out.data<T>(), sgrad_out.data<T>()};
    this->plintrans_batched_act(this->neg_one_, proj_in, 2,
                                hsub_work.data<T>(), this->one_, proj_out);

    T* grad = grad_out.data<T>();
    T* sgrad = sgrad_out.data<T>();
    for (int jb = 0; jb < this->n_band_l; jb++) {
        std::fill(grad + jb * nbs + nvalid, grad + (jb + 1) * nbs, static_cast<T>(0.0));
        std::fill(sgrad + jb * nbs + nvalid, sgrad + (jb + 1) * nbs, static_cast<T>(0.0));
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

    this->pmmcn.multiply(this->one_, psi_in.data<T>(), spdir_out.data<T>(),
                         this->zero_, hsub_work.data<T>());
    const T* proj_in[3] = {psi_in.data<T>(), spsi_in.data<T>(), hpsi_in.data<T>()};
    T* proj_out[3] = {pdir_out.data<T>(), spdir_out.data<T>(), hpdir_out.data<T>()};
    this->plintrans_batched_act(this->neg_one_, proj_in, 3,
                                hsub_work.data<T>(), this->one_, proj_out);

    T* pdir = pdir_out.data<T>();
    T* spdir = spdir_out.data<T>();
    T* hpdir = hpdir_out.data<T>();
    for (int jb = 0; jb < this->n_band_l; jb++) {
        std::fill(pdir + jb * nbs + nvalid, pdir + (jb + 1) * nbs, static_cast<T>(0.0));
        std::fill(spdir + jb * nbs + nvalid, spdir + (jb + 1) * nbs, static_cast<T>(0.0));
        std::fill(hpdir + jb * nbs + nvalid, hpdir + (jb + 1) * nbs, static_cast<T>(0.0));
    }
}

template <typename T, typename Device>
typename DiagoLobpcg<T, Device>::Real DiagoLobpcg<T, Device>::max_error(
    const ct::Tensor& err_in) const
{
    Real* err = err_in.data<Real>();
    Real max_residual = static_cast<Real>(0.0);
    for (int ib = 0; ib < this->n_band_l; ++ib) {
        max_residual = std::max(max_residual, err[ib]);
    }
    detail::lobpcg_reduce_max(max_residual);
    return max_residual;
}

template <typename T, typename Device>
int DiagoLobpcg<T, Device>::count_not_converged(
    const ct::Tensor& err_in,
    const std::vector<double>& ethr_band) const
{
    Real* err = err_in.data<Real>();
    int notconv = 0;
    for (int ib = 0; ib < this->n_band_l; ++ib) {
        if (err[ib] > static_cast<Real>(ethr_band[ib])) {
            ++notconv;
        }
    }
    detail::lobpcg_reduce_sum(notconv);
    return notconv;
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::report_not_converged(
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
    detail::lobpcg_reduce_sum(notconv);
    detail::lobpcg_reduce_max(max_residual);
    if (notconv <= 0) {
        return;
    }

    std::ostringstream msg;
    msg << "DiagoLobpcg::diag(S!=I) "
        << (stop_reason.empty() ? std::string("reached max_iter=") + std::to_string(max_iter)
                                : stop_reason)
        << " after " << used_iter
        << " iterations; notconv=" << notconv
        << ", max_residual=" << max_residual;
    if (!this->diag_context.empty()) {
        msg << ", context={" << this->diag_context << "}";
    }
    bool print_message = true;
#ifdef __MPI
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    print_message = (rank == 0);
#endif
    if (print_message) {
        std::cout << "\n " << msg.str() << std::endl;
        if (GlobalV::ofs_running.good()) {
            GlobalV::ofs_running << " " << msg.str() << std::endl;
            GlobalV::ofs_running.flush();
        }
    }
    if (this->throw_on_notconv_exceed && this->notconv_max >= 0 && notconv > this->notconv_max) {
        throw std::runtime_error(msg.str());
    }
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::plintrans_batched_act(const T alpha,
                                                   const T* const* a_blocks,
                                                   const int block_count,
                                                   const T* u,
                                                   const T beta,
                                                   T* const* b_blocks)
{
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
        std::fill(hsub_d, hsub_d + this->nsub * this->nsub, static_cast<T>(0.0));
        std::fill(ssub_d, ssub_d + this->nsub * this->nsub, static_cast<T>(0.0));
        for (int jb = 0; jb < active_blocks; ++jb) {
            for (int ib = 0; ib <= jb; ++ib) {
                this->pmmcn.multiply(this->one_, basis_blocks[ib], hbasis_blocks[jb],
                                     this->zero_, this->tmp_hsub.data<T>());
                store_hermitian_block(this->tmp_hsub.data<T>(), hsub_d, ib, jb);

                this->pmmcn.multiply(this->one_, basis_blocks[ib], sbasis_blocks[jb],
                                     this->zero_, this->tmp_ssub.data<T>());
                store_hermitian_block(this->tmp_ssub.data<T>(), ssub_d, ib, jb);
            }
        }
        detail::hermitize(hsub_d, this->nsub, active_blocks * n);
        detail::hermitize(ssub_d, this->nsub, active_blocks * n);
    };

    build_subspace(block_count);
    std::vector<Real> inv_subspace_norm;

    auto solve_scaled_hegvd = [&]() -> bool {
        if (!detail::scale_subspace_by_overlap_diag(hsub_d, ssub_d, this->nsub, m,
                                            inv_subspace_norm)) {
            return false;
        }
        if (!detail::check_subspace_spd<T, ct_Device>(ssub_d, this->nsub, m)) {
            return false;
        }
        try {
            ct::kernels::lapack_hegvd<T, ct_Device>()(
                m, this->nsub, hsub_d, ssub_d,
                this->sub_eigen.data<Real>(), hsub_d);
        } catch (const std::exception&) {
            return false;
        }
        if (!detail::finite_real_block(this->sub_eigen.data<Real>(), m)
            || !detail::finite_scalar_block(hsub_d, m * this->nsub)) {
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
        if (!detail::finite_real_block(s_eval.data(), m) || !detail::finite_scalar_block(s_evec.data(), m * m)) {
            throw std::runtime_error("LOBPCG generalized parallel overlap diagonalization produced non-finite values");
        }

        const Real s_max = s_eval.back();
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
        detail::hermitize(h_comp.data(), rank, rank);

        ct::kernels::lapack_heevd<T, ct_Device>()(
            rank, h_comp.data(), rank, this->sub_eigen.data<Real>());
        if (!detail::finite_real_block(this->sub_eigen.data<Real>(), rank)
            || !detail::finite_scalar_block(h_comp.data(), rank * rank)) {
            throw std::runtime_error(
                "LOBPCG generalized parallel compressed diagonalization produced non-finite values");
        }

        std::vector<T> coeff_scaled(m * n, static_cast<T>(0.0));
        ModuleBase::gemm_op<T, Device>()('N', 'N',
                                         m, n, rank,
                                         this->one,
                                         q.data(), m,
                                         h_comp.data(), rank,
                                         this->zero,
                                         coeff_scaled.data(), m);
        std::fill(hsub_d, hsub_d + this->nsub * this->nsub, static_cast<T>(0.0));
        ModuleBase::matrixCopy<T, Device>()(n, m, coeff_scaled.data(), m, hsub_d, this->nsub);
    };

    if (force_compressed || !solve_scaled_hegvd()) {
        bool solved = false;
        if (!force_compressed && block_count == 3) {
            block_count = 2;
            m = block_count * n;
            build_subspace(block_count);
            solved = solve_scaled_hegvd();
        }
        if (!solved) {
            build_subspace(block_count);
            solve_compressed_heevd();
        }
    }

    std::copy(this->sub_eigen.data<Real>(), this->sub_eigen.data<Real>() + n, eigen.data<Real>());

    T* x = this->work.data<T>();
    T* hx = this->hwork.data<T>();
    T* sx = this->swork.data<T>();
    T* p = this->pwork.data<T>();
    T* hp = this->hpwork.data<T>();
    T* sp = this->spwork.data<T>();
    T* update_blocks[6] = {x, hx, sx, p, hp, sp};
    for (T* block : update_blocks)
        std::fill(block, block + local_sz, static_cast<T>(0.0));

    auto copy_coeff_block = [=](const int block, T* coeff) {
        std::fill(coeff, coeff + n * n, static_cast<T>(0.0));
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
            T* update_out[3] = {x, hx, sx};
            this->plintrans_batched_act(this->one_, update_in, 3,
                                        this->tmp_hsub.data<T>(),
                                        this->zero_,
                                        update_out);
        } else {
            T* update_out[3] = {p, hp, sp};
            this->plintrans_batched_act(this->one_, update_in, 3,
                                        this->tmp_hsub.data<T>(),
                                        ib == 1 ? this->zero_ : this->one_,
                                        update_out);
        }
    }
    ModuleBase::axpy_op<T, Device>()(local_sz, this->one, p, 1, x, 1);
    ModuleBase::axpy_op<T, Device>()(local_sz, this->one, hp, 1, hx, 1);
    ModuleBase::axpy_op<T, Device>()(local_sz, this->one, sp, 1, sx, 1);

    T* dst_blocks[6] = {psi.data<T>(), hpsi.data<T>(), spsi.data<T>(),
                        pdir.data<T>(), hpdir.data<T>(), spdir.data<T>()};
    for (int i = 0; i < 6; ++i)
        std::copy(update_blocks[i], update_blocks[i] + local_sz, dst_blocks[i]);
    this->ensure_generalized_update_finite(
        psi, hpsi, spsi, eigen, "LOBPCG generalized parallel update produced non-finite values");

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

    detail::append_s_orthonormal_block<T>(n, nbs, nvalid, eps,
                                          psi.data<T>(), hpsi.data<T>(), spsi.data<T>(),
                                          basis, hbasis, sbasis);
    detail::append_s_orthonormal_block<T>(n, nbs, nvalid, eps,
                                          grad.data<T>(), hgrad.data<T>(), sgrad.data<T>(),
                                          basis, hbasis, sbasis);
    if (this->has_pdir) {
        detail::append_s_orthonormal_block<T>(n, nbs, nvalid, eps,
                                              pdir.data<T>(), hpdir.data<T>(), spdir.data<T>(),
                                              basis, hbasis, sbasis);
    }

    const int m = static_cast<int>(basis.size() / nbs);
    if (m < n) {
        throw std::runtime_error("LOBPCG generalized subspace lost rank");
    }

    T* hsub_d = this->hsub.data<T>();
    T* ssub_d = this->ssub.data<T>();
    std::fill(hsub_d, hsub_d + this->nsub * this->nsub, static_cast<T>(0.0));
    std::fill(ssub_d, ssub_d + this->nsub * this->nsub, static_cast<T>(0.0));

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
    detail::hermitize(hsub_d, this->nsub, m);
    detail::hermitize(ssub_d, this->nsub, m);

    ct::kernels::lapack_hegvd<T, ct_Device>()(
        m, this->nsub, hsub_d, ssub_d,
        this->sub_eigen.data<Real>(), hsub_d);

    if (!detail::finite_real_block(this->sub_eigen.data<Real>(), m)
        || !detail::finite_scalar_block(hsub_d, m * this->nsub)) {
        throw std::runtime_error("LOBPCG generalized subspace diagonalization produced non-finite values");
    }

    std::copy(this->sub_eigen.data<Real>(), this->sub_eigen.data<Real>() + n, eigen.data<Real>());

    T* x = this->work.data<T>();
    T* hx = this->hwork.data<T>();
    T* sx = this->swork.data<T>();
    T* p = this->pwork.data<T>();
    T* hp = this->hpwork.data<T>();
    T* sp = this->spwork.data<T>();
    T* update_blocks[6] = {x, hx, sx, p, hp, sp};
    for (T* block : update_blocks)
        std::fill(block, block + local_sz, static_cast<T>(0.0));

    ModuleBase::gemm_op<T, Device>()('N', 'N',
                                     nvalid, n, n,
                                     this->one,
                                     basis.data(), nbs,
                                     hsub_d, this->nsub,
                                     this->zero,
                                     x, nbs);
    ModuleBase::gemm_op<T, Device>()('N', 'N',
                                     nvalid, n, n,
                                     this->one,
                                     hbasis.data(), nbs,
                                     hsub_d, this->nsub,
                                     this->zero,
                                     hx, nbs);
    ModuleBase::gemm_op<T, Device>()('N', 'N',
                                     nvalid, n, n,
                                     this->one,
                                     sbasis.data(), nbs,
                                     hsub_d, this->nsub,
                                     this->zero,
                                     sx, nbs);

    const int tail_cols = m - n;
    if (tail_cols > 0) {
        ModuleBase::gemm_op<T, Device>()('N', 'N',
                                         nvalid, n, tail_cols,
                                         this->one,
                                         basis.data() + n * nbs, nbs,
                                         hsub_d + n, this->nsub,
                                         this->zero,
                                         p, nbs);
        ModuleBase::gemm_op<T, Device>()('N', 'N',
                                         nvalid, n, tail_cols,
                                         this->one,
                                         hbasis.data() + n * nbs, nbs,
                                         hsub_d + n, this->nsub,
                                         this->zero,
                                         hp, nbs);
        ModuleBase::gemm_op<T, Device>()('N', 'N',
                                         nvalid, n, tail_cols,
                                         this->one,
                                         sbasis.data() + n * nbs, nbs,
                                         hsub_d + n, this->nsub,
                                         this->zero,
                                         sp, nbs);
        ModuleBase::axpy_op<T, Device>()(local_sz, this->one, p, 1, x, 1);
        ModuleBase::axpy_op<T, Device>()(local_sz, this->one, hp, 1, hx, 1);
        ModuleBase::axpy_op<T, Device>()(local_sz, this->one, sp, 1, sx, 1);
    }

    T* dst_blocks[6] = {psi.data<T>(), hpsi.data<T>(), spsi.data<T>(),
                        pdir.data<T>(), hpdir.data<T>(), spdir.data<T>()};
    for (int i = 0; i < 6; ++i)
        std::copy(update_blocks[i], update_blocks[i] + local_sz, dst_blocks[i]);
    this->ensure_generalized_update_finite(
        psi, hpsi, spsi, eigen, "LOBPCG generalized update produced non-finite values");

    this->has_pdir = true;
}


template <typename T, typename Device>
int DiagoLobpcg<T, Device>::diag(
    const HPsiFunc& hpsi_func, const SPsiFunc& spsi_func, T* psi_in,
    Real* eigenvalue_in, const std::vector<double>& ethr_band)
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

    this->has_pdir = false;
    this->psi = ct::TensorMap(psi_in, t_type, dev_type,
                               {this->n_band_l, this->n_basis});

    this->calc_spsi_with_block(spsi_func, psi_in, this->spsi);

    const int scf_iter = DiagoIterAssist<T, Device>::SCF_ITER;

    std::copy(this->h_prec_ptr,
              this->h_prec_ptr + this->n_basis,
              this->prec.data<Real>());
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
        hpsi_func,
        spsi_func,
        effective_ethr_band,
        scf_iter);

    std::copy(this->eigen.data<Real>() + this->local_band_start(),
              this->eigen.data<Real>() + this->local_band_start() + this->n_band_l,
              eigenvalue_in);
    return used_iter;
}

template class DiagoLobpcg<std::complex<double>, base_device::DEVICE_CPU>;

} // namespace hsolver
