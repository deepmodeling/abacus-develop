#ifndef DIAGO_LOBPCG_DETAIL_H_
#define DIAGO_LOBPCG_DETAIL_H_

#include "source_hsolver/diago_lobpcg.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_comm.h"
#include <ATen/kernels/lapack.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace hsolver {
namespace lobpcg_detail {

// ============================================================================
// Internal helpers
// ============================================================================

template <typename Real>
inline Real generalized_residual_growth_limit()
{
    return static_cast<Real>(10.0);
}

template <typename Real>
inline Real residual_guard_limit(const Real residual_before,
                                 const Real residual_growth_limit)
{
    return std::max(static_cast<Real>(1.0e-8),
                    residual_before * residual_growth_limit);
}

template <typename Real>
inline bool state_is_better(const int candidate_notconv,
                            const Real candidate_residual,
                            const bool has_best_state,
                            const int best_notconv,
                            const Real best_residual)
{
    if (!std::isfinite(candidate_residual)) {
        return false;
    }
    if (!has_best_state || candidate_notconv < best_notconv) {
        return true;
    }
    if (candidate_notconv > best_notconv) {
        return false;
    }
    const Real best_tol = std::max(static_cast<Real>(1.0e-12),
                                   std::abs(best_residual) * static_cast<Real>(1.0e-8));
    return candidate_residual < best_residual - best_tol;
}

template <typename Real>
inline bool should_reject_residual_update(const int notconv_before,
                                          const Real residual_before,
                                          const int notconv_after,
                                          const Real residual_after,
                                          const Real residual_growth_limit)
{
    const Real residual_limit = residual_guard_limit(residual_before, residual_growth_limit);
    return !std::isfinite(residual_after)
        || (residual_after > residual_limit && notconv_after >= notconv_before);
}

template <typename Real>
inline bool compressed_guard_is_acceptable(const int notconv_before,
                                           const Real residual_before,
                                           const int guarded_notconv,
                                           const Real guarded_residual,
                                           const Real residual_growth_limit)
{
    const Real guarded_limit = residual_guard_limit(residual_before, residual_growth_limit);
    return std::isfinite(guarded_residual)
        && (guarded_residual <= guarded_limit || guarded_notconv < notconv_before);
}

template <typename Real>
inline bool should_restore_best_state(const bool has_best_state,
                                      const int final_notconv,
                                      const Real final_residual,
                                      const int best_notconv,
                                      const Real best_residual)
{
    const Real best_restore_tol = std::max(static_cast<Real>(1.0e-12),
                                           std::abs(best_residual) * static_cast<Real>(1.0e-8));
    return has_best_state
        && (final_notconv > best_notconv
            || !std::isfinite(final_residual)
            || (final_notconv == best_notconv
                && final_residual > best_residual + best_restore_tol));
}

static constexpr const char* LOBPCG_PROBLEM_GENERALIZED = "S!=I";

template <typename T>
struct LobpcgGeneralizedUpdateBuffers
{
    T* x = nullptr;
    T* hx = nullptr;
    T* sx = nullptr;
    T* p = nullptr;
    T* hp = nullptr;
    T* sp = nullptr;
};

template <typename Real>
static void copy_lowest_subspace_eigenvalues(const Real* sub, Real* eig, const int n)
{
    for (int ib = 0; ib < n; ++ib) {
        eig[ib] = sub[ib];
    }
}

static bool should_continue_for_notconv(const int notconv, const int notconv_max)
{
    return notconv_max >= 0 ? notconv > notconv_max : notconv > 0;
}

static bool notconv_exceeds_limit(const int notconv, const int notconv_max)
{
    return notconv_max >= 0 && notconv > notconv_max;
}

template <typename Real>
static Real max_residual_local(const Real* residual, const int local_size)
{
    Real max_residual = static_cast<Real>(0.0);
    for (int ib = 0; ib < local_size; ++ib) {
        max_residual = std::max(max_residual, residual[ib]);
    }
    return max_residual;
}

template <typename Real>
static int count_not_converged_local(
    const Real* residual,
    const double* threshold,
    const int local_size)
{
    int notconv = 0;
    for (int ib = 0; ib < local_size; ++ib) {
        if (residual[ib] > static_cast<Real>(threshold[ib])) {
            ++notconv;
        }
    }
    return notconv;
}

static void print_lobpcg_diag_message(const std::string& message)
{
    bool print_message = true;
#ifdef __MPI
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    print_message = (rank == 0);
#endif
    if (!print_message) {
        return;
    }
    std::cout << "\n " << message << std::endl;
    if (GlobalV::ofs_running.good()) {
        GlobalV::ofs_running << " " << message << std::endl;
        GlobalV::ofs_running.flush();
    }
}

static void lobpcg_reduce_sum(int& value)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_SUM, BP_WORLD);
#endif
}

static void lobpcg_reduce_max(float& value)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_FLOAT, MPI_MAX, BP_WORLD);
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

static void lobpcg_reduce_bit_or(int& value)
{
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_BOR, BP_WORLD);
#endif
}

template <typename T>
static LobpcgGeneralizedUpdateBuffers<T> make_generalized_update_buffers(T* x, T* hx, T* sx,
                                                                         T* p, T* hp, T* sp)
{
    LobpcgGeneralizedUpdateBuffers<T> buffers;
    buffers.x = x;
    buffers.hx = hx;
    buffers.sx = sx;
    buffers.p = p;
    buffers.hp = hp;
    buffers.sp = sp;
    return buffers;
}

template <typename T, typename Device>
static void zero_update_buffers(const LobpcgGeneralizedUpdateBuffers<T>& buffers, const int local_sz)
{
    using ct_Device = typename ct::PsiToContainer<Device>::type;
    using setmem_complex_op = ct::kernels::set_memory<T, ct_Device>;
    setmem_complex_op()(buffers.x, static_cast<T>(0.0), local_sz);
    setmem_complex_op()(buffers.hx, static_cast<T>(0.0), local_sz);
    setmem_complex_op()(buffers.sx, static_cast<T>(0.0), local_sz);
    setmem_complex_op()(buffers.p, static_cast<T>(0.0), local_sz);
    setmem_complex_op()(buffers.hp, static_cast<T>(0.0), local_sz);
    setmem_complex_op()(buffers.sp, static_cast<T>(0.0), local_sz);
}

template <typename T, typename Device>
static void add_block_to(const T* src, T* dst, const int local_sz, const T* one)
{
    ModuleBase::axpy_op<T, Device>()(local_sz, one, src, 1, dst, 1);
}

template <typename T, typename Device>
static void sync_generalized_update(ct::Tensor& psi,
                                    ct::Tensor& hpsi,
                                    ct::Tensor& spsi,
                                    ct::Tensor& pdir,
                                    ct::Tensor& hpdir,
                                    ct::Tensor& spdir,
                                    const LobpcgGeneralizedUpdateBuffers<T>& buffers,
                                    const int local_sz)
{
    using ct_Device = typename ct::PsiToContainer<Device>::type;
    using syncmem_complex_op = ct::kernels::synchronize_memory<T, ct_Device, ct_Device>;
    syncmem_complex_op()(psi.data<T>(),   buffers.x,  local_sz);
    syncmem_complex_op()(hpsi.data<T>(),  buffers.hx, local_sz);
    syncmem_complex_op()(spsi.data<T>(),  buffers.sx, local_sz);
    syncmem_complex_op()(pdir.data<T>(),  buffers.p,  local_sz);
    syncmem_complex_op()(hpdir.data<T>(), buffers.hp, local_sz);
    syncmem_complex_op()(spdir.data<T>(), buffers.sp, local_sz);
}

template <typename T>
static void clean_hermitian_diag(T* mat, int ld, int active_sub)
{
    using Real = typename GetTypeReal<T>::type;
    for (int i = 0; i < active_sub; i++) {
        int idx = i * ld + i;
        mat[idx] = T(std::real(mat[idx]), static_cast<Real>(0));
    }
}

template <typename T>
static void hermitize(T* mat, int ld, int active_sub)
{
    clean_hermitian_diag(mat, ld, active_sub);
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

template <typename T>
static std::string hermitian_matrix_diagnostics(const char* name, const T* mat, int ld, int n)
{
    using Real = typename GetTypeReal<T>::type;
    int nonfinite = 0;
    Real max_abs = static_cast<Real>(0.0);
    Real max_antiherm = static_cast<Real>(0.0);
    Real min_diag = std::numeric_limits<Real>::max();
    Real max_diag = -std::numeric_limits<Real>::max();

    for (int j = 0; j < n; ++j) {
        const T diag = mat[j * ld + j];
        if (!std::isfinite(std::real(diag)) || !std::isfinite(std::imag(diag))) {
            ++nonfinite;
        } else {
            min_diag = std::min(min_diag, static_cast<Real>(std::real(diag)));
            max_diag = std::max(max_diag, static_cast<Real>(std::real(diag)));
        }
        for (int i = 0; i < n; ++i) {
            const T a = mat[j * ld + i];
            if (!std::isfinite(std::real(a)) || !std::isfinite(std::imag(a))) {
                ++nonfinite;
                continue;
            }
            max_abs = std::max(max_abs, static_cast<Real>(std::abs(a)));
            const T b = mat[i * ld + j];
            if (std::isfinite(std::real(b)) && std::isfinite(std::imag(b))) {
                max_antiherm = std::max(max_antiherm, static_cast<Real>(std::abs(a - std::conj(b))));
            }
        }
    }

    std::ostringstream oss;
    oss << name
        << " nonfinite=" << nonfinite
        << " max_abs=" << std::setprecision(12) << max_abs
        << " max_antiherm=" << max_antiherm
        << " diag_min=" << (min_diag == std::numeric_limits<Real>::max() ? static_cast<Real>(0.0) : min_diag)
        << " diag_max=" << (max_diag == -std::numeric_limits<Real>::max() ? static_cast<Real>(0.0) : max_diag);
    return oss.str();
}

template <typename T>
static std::string s_overlap_diagnostics(const T* ssub, int ld, int n)
{
    using Real = typename GetTypeReal<T>::type;
    Real max_dev = static_cast<Real>(0.0);
    Real max_offdiag = static_cast<Real>(0.0);
    Real min_diag = std::numeric_limits<Real>::max();
    Real max_diag = -std::numeric_limits<Real>::max();
    int nonfinite = 0;

    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            const T value = ssub[j * ld + i];
            if (!std::isfinite(std::real(value)) || !std::isfinite(std::imag(value))) {
                ++nonfinite;
                continue;
            }
            const Real abs_value = static_cast<Real>(std::abs(value));
            if (i == j) {
                const Real diag = static_cast<Real>(std::real(value));
                min_diag = std::min(min_diag, diag);
                max_diag = std::max(max_diag, diag);
                max_dev = std::max(max_dev, static_cast<Real>(std::abs(diag - static_cast<Real>(1.0))));
            } else {
                max_offdiag = std::max(max_offdiag, abs_value);
                max_dev = std::max(max_dev, abs_value);
            }
        }
    }

    std::ostringstream oss;
    oss << "S_orth nonfinite=" << nonfinite
        << " max_abs(S-I)=" << std::setprecision(12) << max_dev
        << " max_offdiag=" << max_offdiag
        << " diag_min=" << (min_diag == std::numeric_limits<Real>::max() ? static_cast<Real>(0.0) : min_diag)
        << " diag_max=" << (max_diag == -std::numeric_limits<Real>::max() ? static_cast<Real>(0.0) : max_diag);
    return oss.str();
}

template <typename Real>
struct SubspaceSpdCheck
{
    bool ok = false;
    Real min_eval = static_cast<Real>(0.0);
    Real max_eval = static_cast<Real>(0.0);
    Real cond = std::numeric_limits<Real>::infinity();
    Real floor = static_cast<Real>(0.0);
    std::string error;
};

template <typename T, typename CtDevice>
static SubspaceSpdCheck<typename GetTypeReal<T>::type>
check_subspace_spd(const T* ssub, int ld, int n)
{
    using Real = typename GetTypeReal<T>::type;
    SubspaceSpdCheck<Real> result;
    if (n <= 0) {
        result.error = "empty matrix";
        return result;
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
    } catch (const std::exception& e) {
        result.error = e.what();
        return result;
    }

    if (!finite_real_block(eval.data(), n)) {
        result.error = "non-finite overlap eigenvalues";
        return result;
    }

    result.min_eval = eval.front();
    result.max_eval = eval.back();
    for (int i = 0; i < n; ++i) {
        result.min_eval = std::min(result.min_eval, eval[i]);
        result.max_eval = std::max(result.max_eval, eval[i]);
    }

    const Real eps = std::numeric_limits<Real>::epsilon();
    const Real cond_limit = std::is_same<Real, double>::value
        ? static_cast<Real>(1.0e10)
        : static_cast<Real>(1.0e5);
    const Real scale = std::max(static_cast<Real>(1.0), std::abs(result.max_eval));
    result.floor = std::max(scale / cond_limit,
                            static_cast<Real>(100.0) * eps * scale);

    if (!std::isfinite(result.min_eval) || !std::isfinite(result.max_eval)
        || result.max_eval <= static_cast<Real>(0.0)
        || result.min_eval <= result.floor) {
        result.error = "ill-conditioned overlap";
        return result;
    }

    result.cond = result.max_eval / result.min_eval;
    if (!std::isfinite(result.cond) || result.cond > cond_limit) {
        result.error = "overlap condition exceeds limit";
        result.ok = false;
        return result;
    }

    result.ok = true;
    return result;
}

template <typename Real>
static std::string subspace_spd_diagnostics(const SubspaceSpdCheck<Real>& check)
{
    std::ostringstream oss;
    oss << "S_sub spectrum min=" << std::setprecision(12) << check.min_eval
        << " max=" << check.max_eval
        << " cond=" << check.cond
        << " floor=" << check.floor;
    if (!check.error.empty()) {
        oss << " error=" << check.error;
    }
    return oss.str();
}

template <typename T>
static bool scale_subspace_by_overlap_diag(T* hsub,
                                           T* ssub,
                                           int ld,
                                           int n,
                                           std::vector<typename GetTypeReal<T>::type>& inv_norm,
                                           std::string& error)
{
    using Real = typename GetTypeReal<T>::type;
    inv_norm.assign(n, static_cast<Real>(0.0));
    const Real diag_floor = std::numeric_limits<Real>::min();

    for (int i = 0; i < n; ++i) {
        const T diag = ssub[i * ld + i];
        const Real diag_real = static_cast<Real>(std::real(diag));
        if (!std::isfinite(diag_real) || diag_real <= diag_floor) {
            std::ostringstream oss;
            oss << "invalid overlap diagonal at column " << i
                << ": " << std::setprecision(12) << diag_real;
            error = oss.str();
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

template <typename T, typename Real>
static int update_invalid_mask(const T* psi,
                               const T* hpsi,
                               const T* spsi,
                               const Real* eigen,
                               const int nvec,
                               const int lda,
                               const int nvalid,
                               const int n)
{
    int invalid_mask =
        (!finite_vector_block(psi, nvec, lda, nvalid) ? 1 : 0)
      | (!finite_vector_block(hpsi, nvec, lda, nvalid) ? 2 : 0);
    if (spsi != nullptr) {
        invalid_mask |= !finite_vector_block(spsi, nvec, lda, nvalid) ? 4 : 0;
    }
    invalid_mask |= !finite_real_block(eigen, n) ? (spsi != nullptr ? 8 : 4) : 0;
    lobpcg_reduce_bit_or(invalid_mask);
    return invalid_mask;
}

template <typename T>
static std::string vector_block_diagnostics(const char* name,
                                            const T* data,
                                            int nvec,
                                            int lda,
                                            int nvalid)
{
    using Real = typename GetTypeReal<T>::type;
    int nonfinite = 0;
    int nonfinite_bands = 0;
    int first_band = -1;
    int first_ig = -1;
    Real min_norm = std::numeric_limits<Real>::max();
    Real max_norm = static_cast<Real>(0.0);
    Real max_abs = static_cast<Real>(0.0);

    for (int ib = 0; ib < nvec; ++ib) {
        bool band_bad = false;
        Real norm2 = static_cast<Real>(0.0);
        for (int ig = 0; ig < nvalid; ++ig) {
            const T value = data[ib * lda + ig];
            if (!std::isfinite(std::real(value)) || !std::isfinite(std::imag(value))) {
                ++nonfinite;
                band_bad = true;
                if (first_band < 0) {
                    first_band = ib;
                    first_ig = ig;
                }
                continue;
            }
            const Real abs_value = static_cast<Real>(std::abs(value));
            max_abs = std::max(max_abs, abs_value);
            norm2 += std::norm(value);
        }
        if (band_bad) {
            ++nonfinite_bands;
        } else {
            const Real norm = std::sqrt(norm2);
            min_norm = std::min(min_norm, norm);
            max_norm = std::max(max_norm, norm);
        }
    }

    std::ostringstream oss;
    oss << name
        << " nonfinite=" << nonfinite
        << " nonfinite_bands=" << nonfinite_bands
        << " first=(" << first_band << "," << first_ig << ")"
        << " finite_norm_min="
        << (min_norm == std::numeric_limits<Real>::max() ? static_cast<Real>(0.0) : min_norm)
        << " finite_norm_max=" << max_norm
        << " max_abs=" << max_abs;
    return oss.str();
}

template <typename T>
static bool append_s_orthonormal_block(
    const int nvec, const int lda, const int nvalid,
    const typename GetTypeReal<T>::type thresh,
    const T* block, const T* hblock, const T* sblock,
    std::vector<T>& basis, std::vector<T>& hbasis, std::vector<T>& sbasis)
{
    using Real = typename GetTypeReal<T>::type;
    bool appended = false;

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
        appended = true;
    }

    return appended;
}

} // namespace lobpcg_detail
} // namespace hsolver

#endif // DIAGO_LOBPCG_DETAIL_H_
