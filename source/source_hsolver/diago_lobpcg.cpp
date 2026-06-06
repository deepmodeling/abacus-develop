#include "source_hsolver/diago_lobpcg.h"

#include "diago_iter_assist.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_comm.h"

#include <ATen/kernels/lapack.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace hsolver {

using LobpcgClock = std::chrono::steady_clock;

// ============================================================================
// Band-major explicit-loop helpers used by serial generalized fallback paths.
//
// Psi is stored band-major:  psi_data[ib * n_basis + ig],
// shape [n_band_l, n_basis].  ig >= n_dim must be zero-padded.
//
// Subspace matrices (C, V) are column-major for direct LAPACK use:
//   C[col * ld + row]  =  C[j * nb + i]   (nb = leading dimension).
// ============================================================================

/// C(i,j) = sum_ig  conj( A(i,ig) ) * B(j,ig)    standard inner-product <A|B>
template <typename T>
static void inner_product_loop(int nb, int lda, int nvalid,
                               T alpha, const T* A, const T* B,
                               T beta, T* C)
{
    for (int j = 0; j < nb; ++j) {
        for (int i = 0; i < nb; ++i) {
            T sum = static_cast<T>(0.0);
            for (int ig = 0; ig < nvalid; ++ig) {
                sum += std::conj(A[i * lda + ig]) * B[j * lda + ig];
            }
            C[j * nb + i] = alpha * sum + beta * C[j * nb + i];
        }
    }
}

/// newRow_i = sum_k  V(k,i) * oldRow_k
/// V col-major:  V(k,i) = V[i * nb + k]
template <typename T>
static void rotate_loop(int nb, int lda, int nvalid,
                        int ldv, T alpha, const T* V, const T* A,
                        T beta, T* C)
{
    for (int i = 0; i < nb; ++i) {
        for (int ig = 0; ig < nvalid; ++ig) {
            T sum = static_cast<T>(0.0);
            for (int k = 0; k < nb; ++k)
                sum += V[i * ldv + k] * A[k * lda + ig];
            C[i * lda + ig] = (beta == static_cast<T>(0.0))
                             ? alpha * sum
                             : alpha * sum + beta * C[i * lda + ig];
        }
        for (int ig = nvalid; ig < lda; ++ig) {
            C[i * lda + ig] = static_cast<T>(0.0);
        }
    }
}

// ============================================================================
// File-static helpers
// ============================================================================

template <typename T>
static void mirror_lower(T* mat, int ld, int active_sub)
{
    for (int c = 0; c < active_sub; c++)
        for (int r = c + 1; r < active_sub; r++)
            mat[c * ld + r] = std::conj(mat[r * ld + c]);
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

static void lobpcg_diag_log(const std::string& context,
                            const std::string& line1,
                            const std::string& line2,
                            const std::string& line3 = std::string())
{
    std::ostringstream oss;
    oss << " LOBPCG_DIAG " << context << '\n'
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
void DiagoLobpcg<T, Device>::diag_log(const std::string& context,
                                      const std::string& line1,
                                      const std::string& line2,
                                      const std::string& line3) const
{
    const std::string full_context = this->diag_context.empty()
        ? context
        : context + " [" + this->diag_context + "]";
    lobpcg_diag_log(full_context, line1, line2, line3);
}

template <typename T>
static bool append_orthonormal_block(
    const int nvec, const int lda, const int nvalid,
    const typename GetTypeReal<T>::type thresh,
    const T* block, const T* hblock,
    std::vector<T>& basis, std::vector<T>& hbasis)
{
    using Real = typename GetTypeReal<T>::type;
    bool appended = false;
    const Real thresh2 = thresh * thresh;

    for (int ib = 0; ib < nvec; ++ib) {
        std::vector<T> q(lda, static_cast<T>(0.0));
        std::vector<T> hq(lda, static_cast<T>(0.0));
        for (int ig = 0; ig < nvalid; ++ig) {
            q[ig] = block[ib * lda + ig];
            hq[ig] = hblock[ib * lda + ig];
        }

        const int nold = static_cast<int>(basis.size() / lda);
        for (int pass = 0; pass < 2; ++pass) {
            for (int jq = 0; jq < nold; ++jq) {
                T dot = static_cast<T>(0.0);
                for (int ig = 0; ig < nvalid; ++ig) {
                    dot += std::conj(basis[jq * lda + ig]) * q[ig];
                }
#ifdef __MPI
                Parallel_Reduce::reduce_pool(&dot, 1);
#endif
                for (int ig = 0; ig < nvalid; ++ig) {
                    q[ig] -= dot * basis[jq * lda + ig];
                    hq[ig] -= dot * hbasis[jq * lda + ig];
                }
            }
        }

        Real norm2 = static_cast<Real>(0.0);
        for (int ig = 0; ig < nvalid; ++ig) {
            norm2 += std::norm(q[ig]);
        }
#ifdef __MPI
        Parallel_Reduce::reduce_pool(&norm2, 1);
#endif
        if (!std::isfinite(norm2) || norm2 <= thresh2) {
            continue;
        }

        const Real inv_norm = static_cast<Real>(1.0) / std::sqrt(norm2);
        for (int ig = 0; ig < lda; ++ig) {
            basis.push_back(q[ig] * inv_norm);
            hbasis.push_back(hq[ig] * inv_norm);
        }
        appended = true;
    }

    return appended;
}

template <typename T>
static bool append_normalized_block(
    const int nvec, const int lda, const int nvalid,
    const typename GetTypeReal<T>::type thresh,
    const T* block, const T* hblock,
    std::vector<T>& basis, std::vector<T>& hbasis)
{
    using Real = typename GetTypeReal<T>::type;
    bool appended = false;
    const Real thresh2 = thresh * thresh;

    for (int ib = 0; ib < nvec; ++ib) {
        const T* src = block + ib * lda;
        const T* hsrc = hblock + ib * lda;
        Real norm2 = static_cast<Real>(0.0);
        for (int ig = 0; ig < nvalid; ++ig) {
            norm2 += std::norm(src[ig]);
        }
#ifdef __MPI
        Parallel_Reduce::reduce_pool(&norm2, 1);
#endif
        if (!std::isfinite(norm2) || norm2 <= thresh2) {
            continue;
        }

        const Real inv_norm = static_cast<Real>(1.0) / std::sqrt(norm2);
        for (int ig = 0; ig < nvalid; ++ig) {
            basis.push_back(src[ig] * inv_norm);
            hbasis.push_back(hsrc[ig] * inv_norm);
        }
        for (int ig = nvalid; ig < lda; ++ig) {
            basis.push_back(static_cast<T>(0.0));
            hbasis.push_back(static_cast<T>(0.0));
        }
        appended = true;
    }

    return appended;
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
#else
    this->pmmcn.set_dimension(n_band_l, n_basis,
                              n_band_l, n_basis,
                              n_dim, n_band);
    this->plintrans.set_dimension(n_dim, nband_l, n_band_l, n_basis,
                                  false);
#endif
}

template <typename T, typename Device>
int DiagoLobpcg<T, Device>::local_band_start() const
{
#ifdef __MPI
    if (this->plintrans.nproc_col > 1) {
        return this->plintrans.start_colB[GlobalV::MY_BNDGROUP];
    }
#endif
    return 0;
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
    if (!finite_vector_block(psi_in, this->n_band_l, this->n_basis, this->n_dim)
        || !finite_vector_block(hpsi_out.data<T>(), this->n_band_l, this->n_basis, this->n_dim)) {
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
        throw std::runtime_error("LOBPCG hPsi produced non-finite values");
    }
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::calc_spsi_with_block(
    const SPsiFunc& spsi_func, const T* psi_in, ct::Tensor& spsi_out)
{
    spsi_func(psi_in, spsi_out.data<T>(), this->n_basis, this->n_band_l);
    if (!finite_vector_block(psi_in, this->n_band_l, this->n_basis, this->n_dim)
        || !finite_vector_block(spsi_out.data<T>(), this->n_band_l, this->n_basis, this->n_dim)) {
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
        throw std::runtime_error("LOBPCG sPsi produced non-finite values");
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

// ============================================================================
// rayleigh_ritz (NC, S=I)
//
// psi_in is assumed orthonormal.  H_sub = <psi|H|psi>,  heevd,  rotate.
// ============================================================================

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::rayleigh_ritz(
    ct::Tensor& psi_inout, ct::Tensor& hpsi_inout, ct::Tensor& eigen_out)
{
    const int nb = this->n_band;
    const int nbs = this->n_basis;
    const int nvalid = this->n_dim;
    const int local_sz = this->n_band_l * nbs;

    this->pmmcn.multiply(1.0, psi_inout.data<T>(), hpsi_inout.data<T>(),
                         0.0, this->tmp_hsub.data<T>());
    mirror_lower(this->tmp_hsub.data<T>(), nb, nb);
    clean_hermitian_diag(this->tmp_hsub.data<T>(), nb, nb);

    // Force exact Hermitian symmetrization.
    {
        T* hsub = this->tmp_hsub.data<T>();
        for (int jj = 0; jj < nb; ++jj) {
            hsub[jj * nb + jj] = T(std::real(hsub[jj * nb + jj]), 0.0);
            for (int ii = jj + 1; ii < nb; ++ii) {
                T a = hsub[jj * nb + ii];
                T b = std::conj(hsub[ii * nb + jj]);
                T avg = static_cast<T>(0.5) * (a + b);
                hsub[jj * nb + ii] = avg;
                hsub[ii * nb + jj] = std::conj(avg);
            }
        }
    }

    try {
        ct::kernels::lapack_heevd<T, ct_Device>()(
            nb, this->tmp_hsub.data<T>(), nb, eigen_out.data<Real>());
    } catch (const std::exception& e) {
        this->diag_log("rayleigh_ritz heevd failed: " + std::string(e.what()),
                        hermitian_matrix_diagnostics("H_sub", this->tmp_hsub.data<T>(), nb, nb),
                        "S_sub unavailable for standard Rayleigh-Ritz");
        throw;
    }

    this->rotate_wf(this->tmp_hsub, psi_inout, this->work);
    syncmem_complex_op()(psi_inout.data<T>(), this->work.data<T>(), local_sz);

    this->rotate_wf(this->tmp_hsub, hpsi_inout, this->work);
    syncmem_complex_op()(hpsi_inout.data<T>(), this->work.data<T>(), local_sz);
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

    for (int ii = 0; ii < nb * nb; ++ii) {
        this->tmp_hsub.data<T>()[ii] = static_cast<T>(0.0);
        this->tmp_ssub.data<T>()[ii] = static_cast<T>(0.0);
    }

    inner_product_loop<T>(nb, nbs, nvalid, this->one_,
                          psi_inout.data<T>(), hpsi_inout.data<T>(),
                          this->zero_, this->tmp_hsub.data<T>());
    inner_product_loop<T>(nb, nbs, nvalid, this->one_,
                          psi_inout.data<T>(), spsi_inout.data<T>(),
                          this->zero_, this->tmp_ssub.data<T>());
#ifdef __MPI
    Parallel_Reduce::reduce_pool(this->tmp_hsub.data<T>(), nb * nb);
    Parallel_Reduce::reduce_pool(this->tmp_ssub.data<T>(), nb * nb);
#endif
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

    if (std::isfinite(std::real(this->hsub.data<T>()[0]))) {
        bool large_eigvec = false;
        const T* eigvec = this->hsub.data<T>();
        for (int ii = 0; ii < nb * nb; ++ii) {
            if (std::isfinite(std::real(eigvec[ii])) && std::isfinite(std::imag(eigvec[ii]))
                && std::abs(eigvec[ii]) > static_cast<Real>(1.0e100)) {
                large_eigvec = true;
                break;
            }
        }
        if (large_eigvec) {
            this->diag_log("generalized_rayleigh_ritz huge eigenvectors before rotation",
                            hermitian_matrix_diagnostics("H_sub", this->tmp_hsub.data<T>(), nb, nb),
                            hermitian_matrix_diagnostics("S_sub", this->tmp_ssub.data<T>(), nb, nb),
                            vector_block_diagnostics("eigvec", this->hsub.data<T>(), nb, nb, nb));
        }
    }

    rotate_loop<T>(nb, nbs, nvalid, nb, this->one_,
                   this->hsub.data<T>(), psi_inout.data<T>(),
                   this->zero_, this->work.data<T>());
    syncmem_complex_op()(psi_inout.data<T>(), this->work.data<T>(), local_sz);

    rotate_loop<T>(nb, nbs, nvalid, nb, this->one_,
                   this->hsub.data<T>(), hpsi_inout.data<T>(),
                   this->zero_, this->work.data<T>());
    syncmem_complex_op()(hpsi_inout.data<T>(), this->work.data<T>(), local_sz);

    rotate_loop<T>(nb, nbs, nvalid, nb, this->one_,
                   this->hsub.data<T>(), spsi_inout.data<T>(),
                   this->zero_, this->work.data<T>());
    syncmem_complex_op()(spsi_inout.data<T>(), this->work.data<T>(), local_sz);

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

    this->pmmcn.multiply(this->one_,
                         psi_inout.data<T>(),
                         hpsi_inout.data<T>(),
                         this->zero_,
                         this->tmp_hsub.data<T>());
    this->pmmcn.multiply(this->one_,
                         psi_inout.data<T>(),
                         spsi_inout.data<T>(),
                         this->zero_,
                         this->tmp_ssub.data<T>());
    for (int jc = 0; jc < nb; ++jc) {
        std::copy(this->tmp_hsub.data<T>() + jc * nb,
                  this->tmp_hsub.data<T>() + jc * nb + nb,
                  hsub_d + jc * this->nsub);
        std::copy(this->tmp_ssub.data<T>() + jc * nb,
                  this->tmp_ssub.data<T>() + jc * nb + nb,
                  ssub_d + jc * this->nsub);
    }
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

    setmem_complex_op()(this->tmp_hsub.data<T>(), static_cast<T>(0.0), nb * nb);
    for (int jc = 0; jc < nb; ++jc) {
        std::copy(hsub_d + jc * this->nsub,
                  hsub_d + jc * this->nsub + nb,
                  this->tmp_hsub.data<T>() + jc * nb);
    }

    this->plintrans.act(this->one_,
                        psi_inout.data<T>(),
                        this->tmp_hsub.data<T>(),
                        this->zero_,
                        this->work.data<T>());
    syncmem_complex_op()(psi_inout.data<T>(), this->work.data<T>(), local_sz);

    this->plintrans.act(this->one_,
                        hpsi_inout.data<T>(),
                        this->tmp_hsub.data<T>(),
                        this->zero_,
                        this->work.data<T>());
    syncmem_complex_op()(hpsi_inout.data<T>(), this->work.data<T>(), local_sz);

    this->plintrans.act(this->one_,
                        spsi_inout.data<T>(),
                        this->tmp_hsub.data<T>(),
                        this->zero_,
                        this->work.data<T>());
    syncmem_complex_op()(spsi_inout.data<T>(), this->work.data<T>(), local_sz);

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

// ============================================================================
// compute_residual — NC: R = HX - lambda*X,  grad = R ./ prec
// ============================================================================

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::compute_residual(
    const ct::Tensor& psi_in, const ct::Tensor& hpsi_in,
    const ct::Tensor& eigen_in, const ct::Tensor& prec_in,
    ct::Tensor& grad_out, ct::Tensor& err_out)
{
    const Real* _prec  = prec_in.data<Real>();
    const Real* _eigen = eigen_in.data<Real>();
    const T*    _psi   = psi_in.data<T>();
    const T*    _hpsi  = hpsi_in.data<T>();
    T*          _grad  = grad_out.data<T>();
    Real*       _err   = err_out.data<Real>();
    const int band_start = this->local_band_start();

    for (int ib = 0; ib < this->n_band_l; ib++) {
        const int  ioff   = ib * this->n_basis;
        const Real lambda = _eigen[band_start + ib];
        Real       err_j  = 0.0;
        for (int ig = 0; ig < this->n_dim; ig++) {
            const int idx = ioff + ig;
            const T   r   = _hpsi[idx] - lambda * _psi[idx];
            _grad[idx]    = r / std::max(_prec[ig],
                                             static_cast<Real>(1e-8));
            err_j        += std::norm(r);
        }
        for (int ig = this->n_dim; ig < this->n_basis; ig++)
            _grad[ioff + ig] = static_cast<T>(0.0);
        _err[ib] = err_j;
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(_err, this->n_band_l);
#endif
    for (int ib = 0; ib < this->n_band_l; ib++)
        _err[ib] = std::sqrt(_err[ib]);
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
                this->diag_log("compute_residual_s non-finite grad",
                                oss.str(),
                                vector_block_diagnostics("hpsi", _hpsi, this->n_band_l, this->n_basis, this->n_dim),
                                vector_block_diagnostics("spsi", _spsi, this->n_band_l, this->n_basis, this->n_dim));
                throw std::runtime_error("LOBPCG generalized residual produced non-finite gradient");
            }
            err_j        += std::norm(r);
        }
        for (int ig = this->n_dim; ig < this->n_basis; ig++)
            _grad[ioff + ig] = static_cast<T>(0.0);
        _err[ib] = err_j;
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(_err, this->n_band_l);
#endif
    for (int ib = 0; ib < this->n_band_l; ib++)
        _err[ib] = std::sqrt(_err[ib]);
}

// ============================================================================
// orth_projection — grad -= psi * <psi|grad>   [S=I]
//   inner(i,j) = <psi_i|grad_j>   (col-major)
//   grad_j    -= sum_i psi_i * inner(i,j)
// ============================================================================

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::orth_projection(
    const ct::Tensor& psi_in, ct::Tensor& hsub_work, ct::Tensor& grad_out)
{
    const int nbs = this->n_basis;
    const int nvalid = this->n_dim;

    this->pmmcn.multiply(1.0, psi_in.data<T>(), grad_out.data<T>(),
                         0.0, hsub_work.data<T>());
    this->plintrans.act(-1.0, psi_in.data<T>(), hsub_work.data<T>(),
                        1.0, grad_out.data<T>());
    T* grad = grad_out.data<T>();
    for (int jb = 0; jb < this->n_band_l; ++jb) {
        for (int ig = nvalid; ig < nbs; ++ig) {
            grad[jb * nbs + ig] = static_cast<T>(0.0);
        }
    }
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::orth_projection_s(
    const ct::Tensor& psi_in, const ct::Tensor& spsi_in,
    ct::Tensor& hsub_work, ct::Tensor& sgrad_out, ct::Tensor& grad_out)
{
    const int nbs = this->n_basis;
    const int nvalid = this->n_dim;

    this->pmmcn.multiply(this->one_,
                         psi_in.data<T>(),
                         sgrad_out.data<T>(),
                         this->zero_,
                         hsub_work.data<T>());
    this->plintrans.act(this->neg_one_,
                        psi_in.data<T>(),
                        hsub_work.data<T>(),
                        this->one_,
                        grad_out.data<T>());
    this->plintrans.act(this->neg_one_,
                        spsi_in.data<T>(),
                        hsub_work.data<T>(),
                        this->one_,
                        sgrad_out.data<T>());

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

    this->pmmcn.multiply(this->one_,
                         psi_in.data<T>(),
                         spdir_out.data<T>(),
                         this->zero_,
                         hsub_work.data<T>());
    this->plintrans.act(this->neg_one_,
                        psi_in.data<T>(),
                        hsub_work.data<T>(),
                        this->one_,
                        pdir_out.data<T>());
    this->plintrans.act(this->neg_one_,
                        spsi_in.data<T>(),
                        hsub_work.data<T>(),
                        this->one_,
                        spdir_out.data<T>());
    this->plintrans.act(this->neg_one_,
                        hpsi_in.data<T>(),
                        hsub_work.data<T>(),
                        this->one_,
                        hpdir_out.data<T>());

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
    this->plintrans.act(1.0, psi_out.data<T>(), hsub_in.data<T>(),
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

// ============================================================================
// orth_cholesky — S=I Cholesky orthonormalization
// ============================================================================

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::orth_cholesky(
    ct::Tensor& workspace_in, ct::Tensor& psi_out,
    ct::Tensor& hpsi_out, ct::Tensor& hsub_out)
{
    const int nb  = this->n_band;
    const int nbs = this->n_basis;
    const int nvalid = this->n_dim;

    try {
        this->pmmcn.multiply(1.0, psi_out.data<T>(), psi_out.data<T>(),
                             0.0, hsub_out.data<T>());
        ct::kernels::set_matrix<T, ct_Device>()('L', hsub_out.data<T>(), nb);
        ct::kernels::lapack_potrf<T, ct_Device>()('U', nb, hsub_out.data<T>(), nb);
        ct::kernels::lapack_trtri<T, ct_Device>()('U', 'N', nb, hsub_out.data<T>(), nb);
        this->rotate_wf(hsub_out, psi_out, workspace_in);
        this->rotate_wf(hsub_out, hpsi_out, workspace_in);
        return;
    } catch (const std::exception&) {
        // Fall back to modified Gram-Schmidt when Cholesky sees a near-dependent block.
    }

    T* psi_d  = psi_out.data<T>();
    T* hpsi_d = hpsi_out.data<T>();
    const Real eps = static_cast<Real>(100)
                   * std::numeric_limits<Real>::epsilon();

    for (int ib = 0; ib < nb; ++ib) {
        for (int jb = 0; jb < ib; ++jb) {
            T dot = static_cast<T>(0.0);
            for (int ig = 0; ig < nvalid; ++ig) {
                dot += std::conj(psi_d[jb * nbs + ig])
                     * psi_d[ib * nbs + ig];
            }
#ifdef __MPI
            Parallel_Reduce::reduce_pool(&dot, 1);
#endif
            for (int ig = 0; ig < nvalid; ++ig) {
                psi_d[ib * nbs + ig]  -= dot * psi_d[jb * nbs + ig];
                hpsi_d[ib * nbs + ig] -= dot * hpsi_d[jb * nbs + ig];
            }
        }

        Real norm2 = static_cast<Real>(0.0);
        for (int ig = 0; ig < nvalid; ++ig) {
            norm2 += std::norm(psi_d[ib * nbs + ig]);
        }
#ifdef __MPI
        Parallel_Reduce::reduce_pool(&norm2, 1);
#endif
        if (!(norm2 > eps)) {
            throw std::runtime_error("orth_cholesky failed: dependent vectors");
        }
        const Real inv_norm = static_cast<Real>(1.0) / std::sqrt(norm2);
        for (int ig = 0; ig < nvalid; ++ig) {
            psi_d[ib * nbs + ig]  *= inv_norm;
            hpsi_d[ib * nbs + ig] *= inv_norm;
        }
        for (int ig = nvalid; ig < nbs; ++ig) {
            psi_d[ib * nbs + ig] = static_cast<T>(0.0);
            hpsi_d[ib * nbs + ig] = static_cast<T>(0.0);
        }
    }
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
    if (ethr_band.size() < static_cast<size_t>(this->n_band_l)) {
        std::ostringstream oss;
        oss << "LOBPCG ethr_band size mismatch: size=" << ethr_band.size()
            << ", required=" << this->n_band_l;
        if (!this->diag_context.empty()) {
            oss << ", context={" << this->diag_context << "}";
        }
        throw std::invalid_argument(oss.str());
    }
}

template <typename T, typename Device>
bool DiagoLobpcg<T, Device>::test_error(
    const ct::Tensor& err_in, const std::vector<double>& ethr_band)
{
    Real* _err_st = err_in.data<Real>();
    bool  not_conv = false;
    std::vector<Real> tmp_cpu;
    if (err_in.device_type() == ct::DeviceType::GpuDevice) {
        tmp_cpu.resize(this->n_band_l);
        _err_st = tmp_cpu.data();
        syncmem_var_d2h_op()(_err_st, err_in.data<Real>(), this->n_band_l);
    }
    for (int ii = 0; ii < this->n_band_l; ii++)
        if (_err_st[ii] > ethr_band[ii]) not_conv = true;
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &not_conv, 1, MPI_C_BOOL, MPI_LOR, BP_WORLD);
#endif
    return not_conv;
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::report_not_converged(
    const char* problem_type,
    const int used_iter,
    const int max_iter,
    const std::vector<double>& ethr_band) const
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
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &notconv, 1, MPI_INT, MPI_SUM, BP_WORLD);
    const MPI_Datatype real_type = std::is_same<Real, double>::value ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Allreduce(MPI_IN_PLACE, &max_residual, 1, real_type, MPI_MAX, BP_WORLD);
#endif
    if (notconv > 0) {
        std::ostringstream msg;
        msg << "DiagoLobpcg::diag(" << problem_type << ") reached max_iter="
            << max_iter
            << " after " << used_iter
            << " iterations; notconv=" << notconv
            << ", max_residual=" << max_residual;
        if (!this->diag_context.empty()) {
            msg << ", context={" << this->diag_context << "}";
        }
        std::cout << "\n " << msg.str() << std::endl;
        if (this->notconv_max >= 0 && notconv > this->notconv_max) {
            throw std::runtime_error(msg.str());
        }
    }
}

template <typename T, typename Device>
bool DiagoLobpcg<T, Device>::profile_enabled() const
{
    return static_cast<long long>(this->n_basis) * this->n_band_l >= 200000;
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::profile_log(const char* problem_type,
                                         const char* stage,
                                         const int iter,
                                         const double seconds) const
{
    if (!this->profile_enabled()) {
        return;
    }
    std::ostringstream oss;
    oss << "LOBPCG_PROFILE " << problem_type
        << " iter=" << iter
        << " stage=" << stage
        << " seconds=" << std::fixed << std::setprecision(6) << seconds;
    if (!this->diag_context.empty()) {
        oss << " context={" << this->diag_context << "}";
    }
    std::cout << "\n " << oss.str() << std::endl;
    if (GlobalV::ofs_running.good()) {
        GlobalV::ofs_running << " " << oss.str() << std::endl;
        GlobalV::ofs_running.flush();
    }
}

// ============================================================================
// lobpcg_update — generalized R-R on subspace W = [X, Z, P]
//
// H_sub = <W|H|W>,  S_sub = <W|W>
// H_sub C = S_sub C Lambda  (hegvd)
// X_new = V_XX*X + V_ZX*Z + V_PX*P
// P_new = V_ZX*Z + V_PX*P         (soft restart)
// ============================================================================

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::lobpcg_update(
    ct::Tensor& psi, ct::Tensor& hpsi,
    ct::Tensor& grad, ct::Tensor& hgrad,
    ct::Tensor& pdir, ct::Tensor& hpdir,
    ct::Tensor& eigen)
{
    const int n    = this->n_band;
    const int nbs  = this->n_basis;
    const int nvalid = this->n_dim;
    const int local_sz = this->n_band_l * nbs;
    const Real eps = static_cast<Real>(100)
                   * std::numeric_limits<Real>::epsilon();

    if (this->n_band_l != this->n_band) {
        int block_count = this->has_pdir ? 3 : 2;
        int m = block_count * n;
        const T* basis_blocks[3] = {psi.data<T>(), grad.data<T>(), pdir.data<T>()};
        const T* hbasis_blocks[3] = {hpsi.data<T>(), hgrad.data<T>(), hpdir.data<T>()};

        T* hsub_d = this->hsub.data<T>();
        T* ssub_d = this->ssub.data<T>();

        auto store_block = [=](const T* src, T* dst, const int iblock, const int jblock) {
            for (int jc = 0; jc < n; ++jc) {
                std::copy(src + jc * n,
                          src + jc * n + n,
                          dst + (jblock * n + jc) * this->nsub + iblock * n);
            }
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
                    this->pmmcn.multiply(this->one_,
                                         basis_blocks[ib],
                                         hbasis_blocks[jb],
                                         this->zero_,
                                         this->tmp_hsub.data<T>());
                    store_hermitian_block(this->tmp_hsub.data<T>(), hsub_d, ib, jb);

                    this->pmmcn.multiply(this->one_,
                                         basis_blocks[ib],
                                         basis_blocks[jb],
                                         this->zero_,
                                         this->tmp_ssub.data<T>());
                    store_hermitian_block(this->tmp_ssub.data<T>(), ssub_d, ib, jb);
                }
            }
            hermitize(hsub_d, this->nsub, active_blocks * n);
            hermitize(ssub_d, this->nsub, active_blocks * n);
        };

        std::vector<Real> inv_subspace_norm;
        std::string scale_error;
        build_subspace(block_count);
        if (!scale_subspace_by_overlap_diag(hsub_d, ssub_d, this->nsub, m,
                                            inv_subspace_norm, scale_error)) {
            if (block_count == 3) {
                block_count = 2;
                m = block_count * n;
                build_subspace(block_count);
                scale_error.clear();
            }
            if (!scale_subspace_by_overlap_diag(hsub_d, ssub_d, this->nsub, m,
                                                inv_subspace_norm, scale_error)) {
                this->diag_log("lobpcg_update failed to normalize S_sub before hegvd",
                               hermitian_matrix_diagnostics("S_sub", ssub_d, this->nsub, m),
                               s_overlap_diagnostics(ssub_d, this->nsub, m),
                               scale_error);
                throw std::runtime_error("LOBPCG subspace overlap has invalid diagonal before hegvd");
            }
        }
        auto spd_check = check_subspace_spd<T, ct_Device>(ssub_d, this->nsub, m);
        if (!spd_check.ok && block_count == 3) {
            if (this->profile_enabled()) {
                this->diag_log("lobpcg_update dropping P due to ill-conditioned S_sub",
                               hermitian_matrix_diagnostics("S_sub", ssub_d, this->nsub, m),
                               s_overlap_diagnostics(ssub_d, this->nsub, m),
                               subspace_spd_diagnostics(spd_check));
            }
            block_count = 2;
            m = block_count * n;
            build_subspace(block_count);
            scale_error.clear();
            if (!scale_subspace_by_overlap_diag(hsub_d, ssub_d, this->nsub, m,
                                                inv_subspace_norm, scale_error)) {
                this->diag_log("lobpcg_update failed to normalize restarted S_sub before hegvd",
                               hermitian_matrix_diagnostics("S_sub", ssub_d, this->nsub, m),
                               s_overlap_diagnostics(ssub_d, this->nsub, m),
                               scale_error);
                throw std::runtime_error("LOBPCG restarted subspace overlap has invalid diagonal before hegvd");
            }
            spd_check = check_subspace_spd<T, ct_Device>(ssub_d, this->nsub, m);
        }

        if (!spd_check.ok) {
            this->diag_log("lobpcg_update rejected ill-conditioned S_sub before hegvd",
                           hermitian_matrix_diagnostics("S_sub", ssub_d, this->nsub, m),
                           s_overlap_diagnostics(ssub_d, this->nsub, m),
                           subspace_spd_diagnostics(spd_check));
            throw std::runtime_error("LOBPCG subspace overlap is ill-conditioned before hegvd");
        }

        try {
            ct::kernels::lapack_hegvd<T, ct_Device>()(
                m, this->nsub, hsub_d, ssub_d,
                this->sub_eigen.data<Real>(), hsub_d);
        } catch (const std::exception& e) {
            this->diag_log("lobpcg_update hegvd failed: " + std::string(e.what()),
                            hermitian_matrix_diagnostics("H_sub", hsub_d, this->nsub, m),
                            hermitian_matrix_diagnostics("S_sub", ssub_d, this->nsub, m),
                            s_overlap_diagnostics(ssub_d, this->nsub, m));
            throw;
        }

        if (!finite_real_block(this->sub_eigen.data<Real>(), m)
            || !finite_scalar_block(hsub_d, m * this->nsub)) {
            this->diag_log("lobpcg_update hegvd produced non-finite values",
                           hermitian_matrix_diagnostics("H_sub", hsub_d, this->nsub, m),
                           hermitian_matrix_diagnostics("S_sub", ssub_d, this->nsub, m),
                           subspace_spd_diagnostics(spd_check));
            throw std::runtime_error("LOBPCG subspace diagonalization produced non-finite values");
        }
        for (int jc = 0; jc < n; ++jc) {
            for (int ir = 0; ir < m; ++ir) {
                hsub_d[jc * this->nsub + ir] *= inv_subspace_norm[ir];
            }
        }

        const Real* sub = this->sub_eigen.data<Real>();
        Real* eig = eigen.data<Real>();
        for (int ib = 0; ib < n; ++ib) {
            eig[ib] = sub[ib];
        }

        T* x_new = this->work.data<T>();
        T* hx_new = this->hwork.data<T>();
        T* p_new = this->pwork.data<T>();
        T* hp_new = this->hpwork.data<T>();
        setmem_complex_op()(x_new, static_cast<T>(0.0), local_sz);
        setmem_complex_op()(hx_new, static_cast<T>(0.0), local_sz);
        setmem_complex_op()(p_new, static_cast<T>(0.0), local_sz);
        setmem_complex_op()(hp_new, static_cast<T>(0.0), local_sz);

        auto copy_coeff_block = [=](const int block, const int first_col, const int ncol, T* coeff) {
            setmem_complex_op()(coeff, static_cast<T>(0.0), n * n);
            for (int jc = 0; jc < ncol; ++jc) {
                std::copy(hsub_d + (first_col + jc) * this->nsub + block * n,
                          hsub_d + (first_col + jc) * this->nsub + block * n + n,
                          coeff + (first_col + jc) * n);
            }
        };

        for (int ib = 0; ib < block_count; ++ib) {
            copy_coeff_block(ib, 0, n, this->tmp_hsub.data<T>());
            this->plintrans.act(this->one_, basis_blocks[ib], this->tmp_hsub.data<T>(),
                                ib == 0 ? this->zero_ : this->one_, x_new);
            this->plintrans.act(this->one_, hbasis_blocks[ib], this->tmp_hsub.data<T>(),
                                ib == 0 ? this->zero_ : this->one_, hx_new);
        }

        const int tail_cols = m - n;
        if (tail_cols > 0) {
            for (int ib = 1; ib < block_count; ++ib) {
                copy_coeff_block(ib, 0, n, this->tmp_hsub.data<T>());
                this->plintrans.act(this->one_, basis_blocks[ib], this->tmp_hsub.data<T>(),
                                    ib == 1 ? this->zero_ : this->one_, p_new);
                this->plintrans.act(this->one_, hbasis_blocks[ib], this->tmp_hsub.data<T>(),
                                    ib == 1 ? this->zero_ : this->one_, hp_new);
            }
        }

        syncmem_complex_op()(psi.data<T>(),   x_new,  local_sz);
        syncmem_complex_op()(hpsi.data<T>(),  hx_new, local_sz);
        syncmem_complex_op()(pdir.data<T>(),  p_new,  local_sz);
        syncmem_complex_op()(hpdir.data<T>(), hp_new, local_sz);

        if (!finite_scalar_block(psi.data<T>(), local_sz)
            || !finite_scalar_block(hpsi.data<T>(), local_sz)
            || !finite_real_block(eigen.data<Real>(), n)) {
            throw std::runtime_error("LOBPCG band-parallel update produced non-finite values");
        }

        this->has_pdir = true;
        return;
    }

    std::vector<T> basis;
    std::vector<T> hbasis;
    basis.reserve((this->has_pdir ? 3 : 2) * local_sz);
    hbasis.reserve((this->has_pdir ? 3 : 2) * local_sz);

    append_orthonormal_block<T>(n, nbs, nvalid, eps,
                                psi.data<T>(), hpsi.data<T>(),
                                basis, hbasis);
    append_orthonormal_block<T>(n, nbs, nvalid, eps,
                                grad.data<T>(), hgrad.data<T>(),
                                basis, hbasis);
    if (this->has_pdir) {
        append_orthonormal_block<T>(n, nbs, nvalid, eps,
                                    pdir.data<T>(), hpdir.data<T>(),
                                    basis, hbasis);
    }

    const int m = static_cast<int>(basis.size() / nbs);
    if (m < n) {
        throw std::runtime_error("LOBPCG standard subspace lost rank");
    }

    T* hsub_d = this->hsub.data<T>();
    setmem_complex_op()(hsub_d, static_cast<T>(0.0), this->nsub * this->nsub);

    ModuleBase::gemm_op<T, Device>()('C', 'N',
                                     m, m, nvalid,
                                     this->one,
                                     basis.data(), nbs,
                                     hbasis.data(), nbs,
                                     this->zero,
                                     hsub_d, this->nsub);
#ifdef __MPI
    for (int jc = 0; jc < m; ++jc) {
        Parallel_Reduce::reduce_pool(hsub_d + jc * this->nsub, m);
    }
#endif
    hermitize(hsub_d, this->nsub, m);

    try {
        ct::kernels::lapack_heevd<T, ct_Device>()(
            m, hsub_d, this->nsub, this->sub_eigen.data<Real>());
    } catch (const std::exception& e) {
        this->diag_log("lobpcg_update heevd failed: " + std::string(e.what()),
                        hermitian_matrix_diagnostics("H_sub", hsub_d, this->nsub, m),
                        "S_sub unavailable after explicit orthonormalization");
        throw;
    }

    if (!finite_real_block(this->sub_eigen.data<Real>(), m)
        || !finite_scalar_block(hsub_d, m * this->nsub)) {
        throw std::runtime_error("LOBPCG subspace diagonalization produced non-finite values");
    }

    const Real* sub = this->sub_eigen.data<Real>();
    Real* eig = eigen.data<Real>();
    for (int ib = 0; ib < n; ++ib) {
        eig[ib] = sub[ib];
    }

    T* x_new = this->work.data<T>();
    T* hx_new = this->hwork.data<T>();
    T* p_new = this->pwork.data<T>();
    T* hp_new = this->hpwork.data<T>();
    setmem_complex_op()(x_new, static_cast<T>(0.0), local_sz);
    setmem_complex_op()(hx_new, static_cast<T>(0.0), local_sz);
    setmem_complex_op()(p_new, static_cast<T>(0.0), local_sz);
    setmem_complex_op()(hp_new, static_cast<T>(0.0), local_sz);

    ModuleBase::gemm_op<T, Device>()('N', 'N',
                                     nvalid, n, m,
                                     this->one,
                                     basis.data(), nbs,
                                     hsub_d, this->nsub,
                                     this->zero,
                                     x_new, nbs);
    ModuleBase::gemm_op<T, Device>()('N', 'N',
                                     nvalid, n, m,
                                     this->one,
                                     hbasis.data(), nbs,
                                     hsub_d, this->nsub,
                                     this->zero,
                                     hx_new, nbs);

    const int tail_cols = m - n;
    if (tail_cols > 0) {
        ModuleBase::gemm_op<T, Device>()('N', 'N',
                                         nvalid, n, tail_cols,
                                         this->one,
                                         basis.data() + n * nbs, nbs,
                                         hsub_d + n, this->nsub,
                                         this->zero,
                                         p_new, nbs);
        ModuleBase::gemm_op<T, Device>()('N', 'N',
                                         nvalid, n, tail_cols,
                                         this->one,
                                         hbasis.data() + n * nbs, nbs,
                                         hsub_d + n, this->nsub,
                                         this->zero,
                                         hp_new, nbs);
    }

    syncmem_complex_op()(psi.data<T>(),   x_new,  local_sz);
    syncmem_complex_op()(hpsi.data<T>(),  hx_new, local_sz);
    syncmem_complex_op()(pdir.data<T>(),  p_new,  local_sz);
    syncmem_complex_op()(hpdir.data<T>(), hp_new, local_sz);

    if (!finite_scalar_block(psi.data<T>(), local_sz)
        || !finite_scalar_block(hpsi.data<T>(), local_sz)
        || !finite_real_block(eigen.data<Real>(), n)) {
        throw std::runtime_error("LOBPCG standard update produced non-finite values");
    }

    this->has_pdir = true;
}

// ============================================================================
// lobpcg_update_s — generalized R-R on S-orthonormalized W = [X, Z, P]
// ============================================================================

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::lobpcg_update_s_parallel(
    ct::Tensor& psi, ct::Tensor& hpsi, ct::Tensor& spsi,
    ct::Tensor& grad, ct::Tensor& hgrad, ct::Tensor& sgrad,
    ct::Tensor& pdir, ct::Tensor& hpdir, ct::Tensor& spdir,
    ct::Tensor& eigen)
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
        for (int jc = 0; jc < n; ++jc) {
            std::copy(src + jc * n,
                      src + jc * n + n,
                      dst + (jblock * n + jc) * this->nsub + iblock * n);
        }
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
                this->pmmcn.multiply(this->one_,
                                     basis_blocks[ib],
                                     hbasis_blocks[jb],
                                     this->zero_,
                                     this->tmp_hsub.data<T>());
                store_hermitian_block(this->tmp_hsub.data<T>(), hsub_d, ib, jb);

                this->pmmcn.multiply(this->one_,
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

    std::vector<Real> inv_subspace_norm;
    std::string scale_error;
    build_subspace(block_count);
    if (!scale_subspace_by_overlap_diag(hsub_d, ssub_d, this->nsub, m,
                                        inv_subspace_norm, scale_error)) {
        if (block_count == 3) {
            block_count = 2;
            m = block_count * n;
            build_subspace(block_count);
            scale_error.clear();
        }
        if (!scale_subspace_by_overlap_diag(hsub_d, ssub_d, this->nsub, m,
                                            inv_subspace_norm, scale_error)) {
            this->diag_log("lobpcg_update_s_parallel failed to normalize S_sub before hegvd",
                           hermitian_matrix_diagnostics("S_sub", ssub_d, this->nsub, m),
                           s_overlap_diagnostics(ssub_d, this->nsub, m),
                           scale_error);
            throw std::runtime_error("LOBPCG generalized parallel subspace overlap has invalid diagonal");
        }
    }

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
    try {
        ct::kernels::lapack_heevd<T, ct_Device>()(m, s_evec.data(), m, s_eval.data());
    } catch (const std::exception& e) {
        this->diag_log("lobpcg_update_s_parallel S_sub heevd failed: " + std::string(e.what()),
                       hermitian_matrix_diagnostics("S_sub", ssub_d, this->nsub, m),
                       s_overlap_diagnostics(ssub_d, this->nsub, m),
                       scale_error);
        throw;
    }
    if (!finite_real_block(s_eval.data(), m) || !finite_scalar_block(s_evec.data(), m * m)) {
        this->diag_log("lobpcg_update_s_parallel S_sub heevd produced non-finite values",
                       hermitian_matrix_diagnostics("S_sub", ssub_d, this->nsub, m),
                       s_overlap_diagnostics(ssub_d, this->nsub, m),
                       "rank compression unavailable");
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
    int rank = m - first_kept;
    if (rank < n && block_count == 3) {
        block_count = 2;
        m = block_count * n;
        build_subspace(block_count);
        scale_error.clear();
        if (!scale_subspace_by_overlap_diag(hsub_d, ssub_d, this->nsub, m,
                                            inv_subspace_norm, scale_error)) {
            this->diag_log("lobpcg_update_s_parallel failed to normalize restarted S_sub before compression",
                           hermitian_matrix_diagnostics("S_sub", ssub_d, this->nsub, m),
                           s_overlap_diagnostics(ssub_d, this->nsub, m),
                           scale_error);
            throw std::runtime_error("LOBPCG generalized restarted parallel subspace overlap has invalid diagonal");
        }
        h_scaled.assign(m * m, static_cast<T>(0.0));
        s_scaled.assign(m * m, static_cast<T>(0.0));
        for (int jc = 0; jc < m; ++jc) {
            for (int ir = 0; ir < m; ++ir) {
                h_scaled[jc * m + ir] = hsub_d[jc * this->nsub + ir];
                s_scaled[jc * m + ir] = ssub_d[jc * this->nsub + ir];
            }
        }
        s_evec = s_scaled;
        s_eval.assign(m, static_cast<Real>(0.0));
        ct::kernels::lapack_heevd<T, ct_Device>()(m, s_evec.data(), m, s_eval.data());
        const Real restart_s_max = s_eval.empty() ? static_cast<Real>(0.0)
                                                  : *std::max_element(s_eval.begin(), s_eval.end());
        const Real restart_rank_floor = std::max(std::abs(restart_s_max) * static_cast<Real>(1.0e-10),
                                                 static_cast<Real>(100.0) * eps * std::max(static_cast<Real>(1.0),
                                                                                           std::abs(restart_s_max)));
        first_kept = 0;
        while (first_kept < m && s_eval[first_kept] <= restart_rank_floor) {
            ++first_kept;
        }
        rank = m - first_kept;
    }
    if (rank < n) {
        this->diag_log("lobpcg_update_s_parallel compressed subspace rank is too small",
                       hermitian_matrix_diagnostics("S_sub", ssub_d, this->nsub, m),
                       s_overlap_diagnostics(ssub_d, this->nsub, m),
                       "rank=" + std::to_string(rank) + " required=" + std::to_string(n));
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

    try {
        ct::kernels::lapack_heevd<T, ct_Device>()(
            rank, h_comp.data(), rank, this->sub_eigen.data<Real>());
    } catch (const std::exception& e) {
        this->diag_log("lobpcg_update_s_parallel compressed heevd failed: " + std::string(e.what()),
                       hermitian_matrix_diagnostics("H_sub", hsub_d, this->nsub, m),
                       hermitian_matrix_diagnostics("S_sub", ssub_d, this->nsub, m),
                       "rank=" + std::to_string(rank));
        throw;
    }
    if (!finite_real_block(this->sub_eigen.data<Real>(), rank)
        || !finite_scalar_block(h_comp.data(), rank * rank)) {
        this->diag_log("lobpcg_update_s_parallel compressed heevd produced non-finite values",
                       hermitian_matrix_diagnostics("H_sub", hsub_d, this->nsub, m),
                       hermitian_matrix_diagnostics("S_sub", ssub_d, this->nsub, m),
                       "rank=" + std::to_string(rank));
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
    for (int jc = 0; jc < n; ++jc) {
        for (int ir = 0; ir < m; ++ir) {
            hsub_d[jc * this->nsub + ir] = coeff_scaled[jc * m + ir] * inv_subspace_norm[ir];
        }
    }

    const Real* sub = this->sub_eigen.data<Real>();
    Real* eig = eigen.data<Real>();
    for (int ib = 0; ib < n; ++ib) {
        eig[ib] = sub[ib];
    }

    T* x_new = this->work.data<T>();
    T* hx_new = this->hwork.data<T>();
    T* sx_new = this->swork.data<T>();
    T* p_new = this->pwork.data<T>();
    T* hp_new = this->hpwork.data<T>();
    T* sp_new = this->spwork.data<T>();
    setmem_complex_op()(x_new, static_cast<T>(0.0), local_sz);
    setmem_complex_op()(hx_new, static_cast<T>(0.0), local_sz);
    setmem_complex_op()(sx_new, static_cast<T>(0.0), local_sz);
    setmem_complex_op()(p_new, static_cast<T>(0.0), local_sz);
    setmem_complex_op()(hp_new, static_cast<T>(0.0), local_sz);
    setmem_complex_op()(sp_new, static_cast<T>(0.0), local_sz);

    auto copy_coeff_block = [=](const int block, T* coeff) {
        setmem_complex_op()(coeff, static_cast<T>(0.0), n * n);
        for (int jc = 0; jc < n; ++jc) {
            std::copy(hsub_d + jc * this->nsub + block * n,
                      hsub_d + jc * this->nsub + block * n + n,
                      coeff + jc * n);
        }
    };

    for (int ib = 0; ib < block_count; ++ib) {
        copy_coeff_block(ib, this->tmp_hsub.data<T>());
        this->plintrans.act(this->one_, basis_blocks[ib], this->tmp_hsub.data<T>(),
                            ib == 0 ? this->zero_ : this->one_, x_new);
        this->plintrans.act(this->one_, hbasis_blocks[ib], this->tmp_hsub.data<T>(),
                            ib == 0 ? this->zero_ : this->one_, hx_new);
        this->plintrans.act(this->one_, sbasis_blocks[ib], this->tmp_hsub.data<T>(),
                            ib == 0 ? this->zero_ : this->one_, sx_new);
    }

    if (m > n) {
        for (int ib = 1; ib < block_count; ++ib) {
            copy_coeff_block(ib, this->tmp_hsub.data<T>());
            this->plintrans.act(this->one_, basis_blocks[ib], this->tmp_hsub.data<T>(),
                                ib == 1 ? this->zero_ : this->one_, p_new);
            this->plintrans.act(this->one_, hbasis_blocks[ib], this->tmp_hsub.data<T>(),
                                ib == 1 ? this->zero_ : this->one_, hp_new);
            this->plintrans.act(this->one_, sbasis_blocks[ib], this->tmp_hsub.data<T>(),
                                ib == 1 ? this->zero_ : this->one_, sp_new);
        }
    }

    syncmem_complex_op()(psi.data<T>(),   x_new,  local_sz);
    syncmem_complex_op()(hpsi.data<T>(),  hx_new, local_sz);
    syncmem_complex_op()(spsi.data<T>(),  sx_new, local_sz);
    syncmem_complex_op()(pdir.data<T>(),  p_new,  local_sz);
    syncmem_complex_op()(hpdir.data<T>(), hp_new, local_sz);
    syncmem_complex_op()(spdir.data<T>(), sp_new, local_sz);

    if (!finite_scalar_block(psi.data<T>(), local_sz)
        || !finite_scalar_block(hpsi.data<T>(), local_sz)
        || !finite_scalar_block(spsi.data<T>(), local_sz)
        || !finite_real_block(eigen.data<Real>(), n)) {
        throw std::runtime_error("LOBPCG generalized parallel update produced non-finite values");
    }

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
    basis.reserve((this->has_pdir ? 3 : 2) * local_sz);
    hbasis.reserve((this->has_pdir ? 3 : 2) * local_sz);
    sbasis.reserve((this->has_pdir ? 3 : 2) * local_sz);

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

    const Real* sub = this->sub_eigen.data<Real>();
    Real* eig = eigen.data<Real>();
    for (int ib = 0; ib < n; ++ib) {
        eig[ib] = sub[ib];
    }

    T* x_new = this->work.data<T>();
    T* hx_new = this->hwork.data<T>();
    T* sx_new = this->swork.data<T>();
    T* p_new = this->pwork.data<T>();
    T* hp_new = this->hpwork.data<T>();
    T* sp_new = this->spwork.data<T>();
    setmem_complex_op()(x_new, static_cast<T>(0.0), local_sz);
    setmem_complex_op()(hx_new, static_cast<T>(0.0), local_sz);
    setmem_complex_op()(sx_new, static_cast<T>(0.0), local_sz);
    setmem_complex_op()(p_new, static_cast<T>(0.0), local_sz);
    setmem_complex_op()(hp_new, static_cast<T>(0.0), local_sz);
    setmem_complex_op()(sp_new, static_cast<T>(0.0), local_sz);

    ModuleBase::gemm_op<T, Device>()('N', 'N',
                                     nvalid, n, m,
                                     this->one,
                                     basis.data(), nbs,
                                     hsub_d, this->nsub,
                                     this->zero,
                                     x_new, nbs);
    ModuleBase::gemm_op<T, Device>()('N', 'N',
                                     nvalid, n, m,
                                     this->one,
                                     hbasis.data(), nbs,
                                     hsub_d, this->nsub,
                                     this->zero,
                                     hx_new, nbs);
    ModuleBase::gemm_op<T, Device>()('N', 'N',
                                     nvalid, n, m,
                                     this->one,
                                     sbasis.data(), nbs,
                                     hsub_d, this->nsub,
                                     this->zero,
                                     sx_new, nbs);

    const int tail_cols = m - n;
    if (tail_cols > 0) {
        ModuleBase::gemm_op<T, Device>()('N', 'N',
                                         nvalid, n, tail_cols,
                                         this->one,
                                         basis.data() + n * nbs, nbs,
                                         hsub_d + n, this->nsub,
                                         this->zero,
                                         p_new, nbs);
        ModuleBase::gemm_op<T, Device>()('N', 'N',
                                         nvalid, n, tail_cols,
                                         this->one,
                                         hbasis.data() + n * nbs, nbs,
                                         hsub_d + n, this->nsub,
                                         this->zero,
                                         hp_new, nbs);
        ModuleBase::gemm_op<T, Device>()('N', 'N',
                                         nvalid, n, tail_cols,
                                         this->one,
                                         sbasis.data() + n * nbs, nbs,
                                         hsub_d + n, this->nsub,
                                         this->zero,
                                         sp_new, nbs);
    }

    syncmem_complex_op()(psi.data<T>(),   x_new,  local_sz);
    syncmem_complex_op()(hpsi.data<T>(),  hx_new, local_sz);
    syncmem_complex_op()(spsi.data<T>(),  sx_new, local_sz);
    syncmem_complex_op()(pdir.data<T>(),  p_new,  local_sz);
    syncmem_complex_op()(hpdir.data<T>(), hp_new, local_sz);
    syncmem_complex_op()(spdir.data<T>(), sp_new, local_sz);

    if (!finite_scalar_block(psi.data<T>(), local_sz)
        || !finite_scalar_block(hpsi.data<T>(), local_sz)
        || !finite_scalar_block(spsi.data<T>(), local_sz)
        || !finite_real_block(eigen.data<Real>(), n)) {
        throw std::runtime_error("LOBPCG generalized update produced non-finite values");
    }

    this->has_pdir = true;
}

// ============================================================================
// diag — main LOBPCG loop (NC, S=I)
// ============================================================================

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::diag(
    const HPsiFunc& hpsi_func, T* psi_in,
    Real* eigenvalue_in, const std::vector<double>& ethr_band)
{
    this->validate_ethr_band(ethr_band);

    this->has_pdir = false;
    const int scf_iter = DiagoIterAssist<T, Device>::SCF_ITER;

    this->psi = ct::TensorMap(psi_in, t_type, dev_type,
                               {this->n_band_l, this->n_basis});

    auto t0 = LobpcgClock::now();
    this->calc_prec();
    this->profile_log("S=I", "initial_calc_prec", 0,
                      std::chrono::duration<double>(LobpcgClock::now() - t0).count());

    t0 = LobpcgClock::now();
    this->calc_hpsi_with_block(hpsi_func, psi_in, this->hpsi);
    this->profile_log("S=I", "initial_hpsi", 0,
                      std::chrono::duration<double>(LobpcgClock::now() - t0).count());
    // Re-orthonormalize before initial R-R so H_sub is well-conditioned
    t0 = LobpcgClock::now();
    this->orth_cholesky(this->work, this->psi, this->hpsi, this->tmp_hsub);
    this->profile_log("S=I", "initial_orth", 0,
                      std::chrono::duration<double>(LobpcgClock::now() - t0).count());
    t0 = LobpcgClock::now();
    this->rayleigh_ritz(this->psi, this->hpsi, this->eigen);
    this->profile_log("S=I", "initial_rr", 0,
                      std::chrono::duration<double>(LobpcgClock::now() - t0).count());

    setmem_complex_op()(this->pdir.data<T>(),  static_cast<T>(0.0),
                         this->n_basis * this->n_band_l);
    setmem_complex_op()(this->hpdir.data<T>(), static_cast<T>(0.0),
                         this->n_basis * this->n_band_l);

    const int default_max_iter = (scf_iter > 1) ? this->nline : (this->nline * 20);
    const int max_iter = (this->max_iter > 0) ? this->max_iter : default_max_iter;
    int used_iter = 0;

    for (int ntry = 0; ntry < max_iter; ++ntry) {
        used_iter = ntry + 1;
        t0 = LobpcgClock::now();
        this->compute_residual(this->psi, this->hpsi, this->eigen,
                               this->prec, this->grad, this->err_st);
        this->profile_log("S=I", "residual", used_iter,
                          std::chrono::duration<double>(LobpcgClock::now() - t0).count());
        if (!this->test_error(this->err_st, ethr_band))
            break;

        const int psi_sz = this->n_basis * this->n_band_l;
        const int eig_sz = this->n_band;

        t0 = LobpcgClock::now();
        this->orth_projection(this->psi, this->tmp_hsub, this->grad);
        this->profile_log("S=I", "grad_projection", used_iter,
                          std::chrono::duration<double>(LobpcgClock::now() - t0).count());

        t0 = LobpcgClock::now();
        this->calc_hpsi_with_block(hpsi_func, this->grad.data<T>(), this->hgrad);
        this->profile_log("S=I", "grad_hpsi", used_iter,
                          std::chrono::duration<double>(LobpcgClock::now() - t0).count());

        // Backup stable state in case lobpcg_update corrupts psi/hpsi
        std::vector<T> psi_bak(psi_sz), hpsi_bak(psi_sz);
        std::vector<Real> eigen_bak(eig_sz);
        std::copy(this->psi.data<T>(),  this->psi.data<T>()  + psi_sz, psi_bak.data());
        std::copy(this->hpsi.data<T>(), this->hpsi.data<T>() + psi_sz, hpsi_bak.data());
        std::copy(this->eigen.data<Real>(), this->eigen.data<Real>() + eig_sz, eigen_bak.data());

        try {
            t0 = LobpcgClock::now();
            this->lobpcg_update(this->psi, this->hpsi,
                                 this->grad, this->hgrad,
                                 this->pdir, this->hpdir,
                                 this->eigen);
            this->profile_log("S=I", "lobpcg_update", used_iter,
                              std::chrono::duration<double>(LobpcgClock::now() - t0).count());
        } catch (const std::exception& e1) {
            std::copy(psi_bak.data(), psi_bak.data() + psi_sz, this->psi.data<T>());
            std::copy(hpsi_bak.data(), hpsi_bak.data() + psi_sz, this->hpsi.data<T>());
            std::copy(eigen_bak.data(), eigen_bak.data() + eig_sz, this->eigen.data<Real>());

            setmem_complex_op()(this->pdir.data<T>(),  static_cast<T>(0.0), psi_sz);
            setmem_complex_op()(this->hpdir.data<T>(), static_cast<T>(0.0), psi_sz);
            this->has_pdir = false;

            try {
                t0 = LobpcgClock::now();
                this->lobpcg_update(this->psi, this->hpsi,
                                     this->grad, this->hgrad,
                                     this->pdir, this->hpdir,
                                     this->eigen);
                this->profile_log("S=I", "lobpcg_update_retry", used_iter,
                                  std::chrono::duration<double>(LobpcgClock::now() - t0).count());
            } catch (const std::exception& e2) {
                std::copy(psi_bak.data(), psi_bak.data() + psi_sz, this->psi.data<T>());
                std::copy(hpsi_bak.data(), hpsi_bak.data() + psi_sz, this->hpsi.data<T>());
                std::copy(eigen_bak.data(), eigen_bak.data() + eig_sz, this->eigen.data<Real>());

                t0 = LobpcgClock::now();
                this->calc_hpsi_with_block(hpsi_func, this->psi.data<T>(), this->hpsi);
                this->orth_cholesky(this->work, this->psi, this->hpsi, this->tmp_hsub);
                this->rayleigh_ritz(this->psi, this->hpsi, this->eigen);
                this->profile_log("S=I", "fallback_rr", used_iter,
                                  std::chrono::duration<double>(LobpcgClock::now() - t0).count());
            }
        }

        const bool has_next_iteration = (ntry + 1) < max_iter;
        const bool restart_next = has_next_iteration && scf_iter == 1 && ((ntry + 1) % this->nline == 0);
        if (restart_next) {
            setmem_complex_op()(this->pdir.data<T>(),  static_cast<T>(0.0),
                                 this->n_basis * this->n_band_l);
            setmem_complex_op()(this->hpdir.data<T>(), static_cast<T>(0.0),
                                 this->n_basis * this->n_band_l);
            this->has_pdir = false;
        }
    }

    t0 = LobpcgClock::now();
    this->compute_residual(this->psi, this->hpsi, this->eigen,
                           this->prec, this->grad, this->err_st);
    this->profile_log("S=I", "final_residual", used_iter,
                      std::chrono::duration<double>(LobpcgClock::now() - t0).count());
    this->report_not_converged("S=I", used_iter, max_iter, ethr_band);
    DiagoIterAssist<T, Device>::avg_iter += static_cast<Real>(used_iter);

    syncmem_var_d2h_op()(eigenvalue_in,
                          this->eigen.data<Real>() + this->local_band_start(),
                          this->n_band_l);
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::diag(
    const HPsiFunc& hpsi_func, const SPsiFunc& spsi_func, T* psi_in,
    Real* eigenvalue_in, const std::vector<double>& ethr_band)
{
    this->validate_ethr_band(ethr_band);

    this->has_pdir = false;
    this->psi = ct::TensorMap(psi_in, t_type, dev_type,
                               {this->n_band_l, this->n_basis});

    this->calc_spsi_with_block(spsi_func, psi_in, this->spsi);
    {
        const T* spsi_d = this->spsi.data<T>();
        Real max_diff = static_cast<Real>(0.0);
        Real max_ref = static_cast<Real>(0.0);
        for (int ib = 0; ib < this->n_band_l; ++ib) {
            const int ioff = ib * this->n_basis;
            for (int ig = 0; ig < this->n_dim; ++ig) {
                const int idx = ioff + ig;
                max_diff = std::max(max_diff,
                                    static_cast<Real>(std::abs(spsi_d[idx] - psi_in[idx])));
                max_ref = std::max(max_ref,
                                   static_cast<Real>(std::abs(psi_in[idx])));
            }
        }
#ifdef __MPI
        if (this->n_band_l != this->n_band) {
            const MPI_Datatype real_type = std::is_same<Real, double>::value ? MPI_DOUBLE : MPI_FLOAT;
            MPI_Allreduce(MPI_IN_PLACE, &max_diff, 1, real_type, MPI_MAX, BP_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &max_ref, 1, real_type, MPI_MAX, BP_WORLD);
        }
#endif
        const Real tol = static_cast<Real>(1.0e-12)
                       * std::max(static_cast<Real>(1.0), max_ref);
        if (max_diff <= tol) {
            this->diag(hpsi_func, psi_in, eigenvalue_in, ethr_band);
            return;
        }
    }

    const int scf_iter = DiagoIterAssist<T, Device>::SCF_ITER;

    this->calc_prec();
    if (this->n_band_l != this->n_band) {
        this->calc_hpsi_with_block(hpsi_func, this->psi.data<T>(), this->hpsi);
        this->generalized_rayleigh_ritz_parallel(this->psi, this->hpsi, this->spsi, this->eigen);
    } else {
        this->repair_initial_subspace_s(hpsi_func, spsi_func);
        this->generalized_rayleigh_ritz(this->psi, this->hpsi, this->spsi, this->eigen);
    }

    setmem_complex_op()(this->pdir.data<T>(),  static_cast<T>(0.0),
                         this->n_basis * this->n_band_l);
    setmem_complex_op()(this->hpdir.data<T>(), static_cast<T>(0.0),
                         this->n_basis * this->n_band_l);
    setmem_complex_op()(this->spdir.data<T>(), static_cast<T>(0.0),
                         this->n_basis * this->n_band_l);

    const int default_max_iter = (scf_iter > 1) ? this->nline : (this->nline * 20);
    const int max_iter = (this->max_iter > 0) ? this->max_iter : default_max_iter;
    int used_iter = 0;
    std::vector<double> effective_ethr_band = ethr_band;
    if (this->notconv_max < 0) {
        // SCF can refine the density across outer iterations; avoid chasing a tiny diagonalization threshold.
        constexpr double scf_generalized_residual_floor = 1.0e-5;
        for (double& ethr : effective_ethr_band) {
            ethr = std::max(ethr, scf_generalized_residual_floor);
        }
    }

    for (int ntry = 0; ntry < max_iter; ++ntry) {
        used_iter = ntry + 1;
        this->compute_residual_s(this->psi, this->hpsi, this->spsi, this->eigen,
                                 this->prec, this->grad, this->err_st);
        if (!this->test_error(this->err_st, effective_ethr_band))
            break;

        const int psi_sz = this->n_basis * this->n_band_l;
        const int eig_sz = this->n_band;

        this->calc_spsi_with_block(spsi_func, this->grad.data<T>(), this->sgrad);
        this->orth_projection_s(this->psi, this->spsi, this->tmp_hsub,
                                this->sgrad, this->grad);
        this->calc_hpsi_with_block(hpsi_func, this->grad.data<T>(), this->hgrad);

        std::vector<T> psi_bak(psi_sz), hpsi_bak(psi_sz), spsi_bak(psi_sz);
        std::vector<Real> eigen_bak(eig_sz);
        std::copy(this->psi.data<T>(),  this->psi.data<T>()  + psi_sz, psi_bak.data());
        std::copy(this->hpsi.data<T>(), this->hpsi.data<T>() + psi_sz, hpsi_bak.data());
        std::copy(this->spsi.data<T>(), this->spsi.data<T>() + psi_sz, spsi_bak.data());
        std::copy(this->eigen.data<Real>(), this->eigen.data<Real>() + eig_sz, eigen_bak.data());

        try {
            if (this->n_band_l != this->n_band) {
                this->lobpcg_update_s_parallel(this->psi, this->hpsi, this->spsi,
                                               this->grad, this->hgrad, this->sgrad,
                                               this->pdir, this->hpdir, this->spdir,
                                               this->eigen);
            } else {
                this->lobpcg_update_s(this->psi, this->hpsi, this->spsi,
                                      this->grad, this->hgrad, this->sgrad,
                                      this->pdir, this->hpdir, this->spdir,
                                      this->eigen);
            }
        } catch (const std::exception&) {
            std::copy(psi_bak.data(), psi_bak.data() + psi_sz, this->psi.data<T>());
            std::copy(hpsi_bak.data(), hpsi_bak.data() + psi_sz, this->hpsi.data<T>());
            std::copy(spsi_bak.data(), spsi_bak.data() + psi_sz, this->spsi.data<T>());
            std::copy(eigen_bak.data(), eigen_bak.data() + eig_sz, this->eigen.data<Real>());

            setmem_complex_op()(this->pdir.data<T>(),  static_cast<T>(0.0), psi_sz);
            setmem_complex_op()(this->hpdir.data<T>(), static_cast<T>(0.0), psi_sz);
            setmem_complex_op()(this->spdir.data<T>(), static_cast<T>(0.0), psi_sz);
            this->has_pdir = false;

            try {
                if (this->n_band_l != this->n_band) {
                    this->lobpcg_update_s_parallel(this->psi, this->hpsi, this->spsi,
                                                   this->grad, this->hgrad, this->sgrad,
                                                   this->pdir, this->hpdir, this->spdir,
                                                   this->eigen);
                } else {
                    this->lobpcg_update_s(this->psi, this->hpsi, this->spsi,
                                          this->grad, this->hgrad, this->sgrad,
                                          this->pdir, this->hpdir, this->spdir,
                                          this->eigen);
                }
            } catch (const std::exception&) {
                std::copy(psi_bak.data(), psi_bak.data() + psi_sz, this->psi.data<T>());
                std::copy(hpsi_bak.data(), hpsi_bak.data() + psi_sz, this->hpsi.data<T>());
                std::copy(spsi_bak.data(), spsi_bak.data() + psi_sz, this->spsi.data<T>());
                std::copy(eigen_bak.data(), eigen_bak.data() + eig_sz, this->eigen.data<Real>());

                this->calc_hpsi_with_block(hpsi_func, this->psi.data<T>(), this->hpsi);
                this->calc_spsi_with_block(spsi_func, this->psi.data<T>(), this->spsi);
                if (this->n_band_l != this->n_band) {
                    this->generalized_rayleigh_ritz_parallel(this->psi, this->hpsi, this->spsi, this->eigen);
                } else {
                    this->generalized_rayleigh_ritz(this->psi, this->hpsi, this->spsi, this->eigen);
                }
            }
        }

        const bool has_next_iteration = (ntry + 1) < max_iter;
        const bool restart_next = has_next_iteration && scf_iter == 1 && ((ntry + 1) % this->nline == 0);
        if (has_next_iteration && !restart_next) {
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
                std::copy(psi_bak.data(), psi_bak.data() + psi_sz, this->psi.data<T>());
                std::copy(hpsi_bak.data(), hpsi_bak.data() + psi_sz, this->hpsi.data<T>());
                std::copy(spsi_bak.data(), spsi_bak.data() + psi_sz, this->spsi.data<T>());
                std::copy(eigen_bak.data(), eigen_bak.data() + eig_sz, this->eigen.data<Real>());
                setmem_complex_op()(this->pdir.data<T>(),  static_cast<T>(0.0), psi_sz);
                setmem_complex_op()(this->hpdir.data<T>(), static_cast<T>(0.0), psi_sz);
                setmem_complex_op()(this->spdir.data<T>(), static_cast<T>(0.0), psi_sz);
                this->has_pdir = false;
            }
        }
        if (restart_next) {
            setmem_complex_op()(this->pdir.data<T>(),  static_cast<T>(0.0), psi_sz);
            setmem_complex_op()(this->hpdir.data<T>(), static_cast<T>(0.0), psi_sz);
            setmem_complex_op()(this->spdir.data<T>(), static_cast<T>(0.0), psi_sz);
            this->has_pdir = false;
        }
    }

    this->compute_residual_s(this->psi, this->hpsi, this->spsi, this->eigen,
                             this->prec, this->grad, this->err_st);
    this->report_not_converged("S!=I", used_iter, max_iter, effective_ethr_band);
    DiagoIterAssist<T, Device>::avg_iter += static_cast<Real>(used_iter);

    syncmem_var_d2h_op()(eigenvalue_in,
                          this->eigen.data<Real>() + this->local_band_start(),
                          this->n_band_l);
}

template class DiagoLobpcg<std::complex<double>, base_device::DEVICE_CPU>;

} // namespace hsolver
