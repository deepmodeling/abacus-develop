#include "source_hsolver/diago_lobpcg.h"

#include "diago_iter_assist.h"
#include "source_base/global_function.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_comm.h"

#include <ATen/kernels/lapack.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace hsolver {

// ============================================================================
// Band-major explicit-loop helpers (CPU, n_band_l == n_band required)
//
// Psi is stored band-major:  psi_data[ib * n_basis + ig],
// shape [n_band_l, n_basis].  ig >= n_dim must be zero-padded.
//
// Subspace matrices (C, V) are column-major for direct LAPACK use:
//   C[col * ld + row]  =  C[j * nb + i]   (nb = leading dimension).
// ============================================================================

/// C(i,j) = sum_ig  conj( A(i,ig) ) * B(j,ig)    standard inner-product <A|B>
template <typename T>
static void inner_product_loop(int nb, int nbs,
                               T alpha, const T* A, const T* B,
                               T beta, T* C)
{
    for (int j = 0; j < nb; ++j) {
        for (int i = 0; i < nb; ++i) {
            T sum = static_cast<T>(0.0);
            for (int ig = 0; ig < nbs; ++ig) {
                sum += std::conj(A[i * nbs + ig]) * B[j * nbs + ig];
            }
            C[j * nb + i] = alpha * sum + beta * C[j * nb + i];
        }
    }
}

/// newRow_i = sum_k  V(k,i) * oldRow_k
/// V col-major:  V(k,i) = V[i * nb + k]
template <typename T>
static void rotate_loop(int nb, int nbs,
                        T alpha, const T* V, const T* A,
                        T beta, T* C)
{
    for (int i = 0; i < nb; ++i) {
        for (int ig = 0; ig < nbs; ++ig) {
            T sum = static_cast<T>(0.0);
            for (int k = 0; k < nb; ++k)
                sum += V[i * nb + k] * A[k * nbs + ig];
            C[i * nbs + ig] = alpha * sum + beta * C[i * nbs + ig];
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
static bool append_orthonormal_block(
    const int nvec, const int nbs, const typename GetTypeReal<T>::type thresh,
    const T* block, const T* hblock,
    std::vector<T>& basis, std::vector<T>& hbasis)
{
    using Real = typename GetTypeReal<T>::type;
    bool appended = false;
    const Real thresh2 = thresh * thresh;

    for (int ib = 0; ib < nvec; ++ib) {
        std::vector<T> q(nbs);
        std::vector<T> hq(nbs);
        for (int ig = 0; ig < nbs; ++ig) {
            q[ig] = block[ib * nbs + ig];
            hq[ig] = hblock[ib * nbs + ig];
        }

        const int nold = static_cast<int>(basis.size() / nbs);
        for (int jq = 0; jq < nold; ++jq) {
            T dot = static_cast<T>(0.0);
            for (int ig = 0; ig < nbs; ++ig) {
                dot += std::conj(basis[jq * nbs + ig]) * q[ig];
            }
#ifdef __MPI
            Parallel_Reduce::reduce_pool(&dot, 1);
#endif
            for (int ig = 0; ig < nbs; ++ig) {
                q[ig] -= dot * basis[jq * nbs + ig];
                hq[ig] -= dot * hbasis[jq * nbs + ig];
            }
        }

        Real norm2 = static_cast<Real>(0.0);
        for (int ig = 0; ig < nbs; ++ig) {
            norm2 += std::norm(q[ig]);
        }
#ifdef __MPI
        Parallel_Reduce::reduce_pool(&norm2, 1);
#endif
        if (norm2 <= thresh2) {
            continue;
        }

        const Real inv_norm = static_cast<Real>(1.0) / std::sqrt(norm2);
        for (int ig = 0; ig < nbs; ++ig) {
            basis.push_back(q[ig] * inv_norm);
            hbasis.push_back(hq[ig] * inv_norm);
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
    this->tmp_hsub = ct::Tensor(t_type, dev_type, {this->n_band, this->n_band});

    this->hpsi  = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});
    this->spsi  = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});
    this->grad  = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});
    this->hgrad = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});
    this->pdir  = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});
    this->hpdir = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});

    this->work   = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});
    this->hwork  = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});
    this->pwork  = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});
    this->hpwork = ct::Tensor(t_type, dev_type, {this->n_band_l, this->n_basis});

    this->prec = ct::Tensor(r_type, dev_type, {this->n_basis});
    this->h_prec = ct::TensorMap(
        (void*)this->h_prec_ptr, this->r_type,
        ct::DeviceType::CpuDevice, {this->n_basis});
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
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::calc_spsi_with_block(
    const SPsiFunc& spsi_func, T* psi_in, ct::Tensor& spsi_out)
{
    spsi_func(psi_in, spsi_out.data<T>(), this->n_basis, this->n_band_l);
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
    const int local_sz = this->n_band_l * nbs;

    // Zero-fill H_sub before inner_product to eliminate stale data
    for (int ii = 0; ii < nb * nb; ++ii)
        this->tmp_hsub.data<T>()[ii] = static_cast<T>(0.0);

    inner_product_loop<T>(nb, nbs, this->one_,
                          psi_inout.data<T>(), hpsi_inout.data<T>(),
                          this->zero_, this->tmp_hsub.data<T>());
#ifdef __MPI
    Parallel_Reduce::reduce_pool(this->tmp_hsub.data<T>(), nb * nb);
#endif
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

    ct::kernels::lapack_heevd<T, ct_Device>()(
        nb, this->tmp_hsub.data<T>(), nb, eigen_out.data<Real>());

    rotate_loop<T>(nb, nbs, this->one_,
                   this->tmp_hsub.data<T>(), psi_inout.data<T>(),
                   this->zero_, this->work.data<T>());
    syncmem_complex_op()(psi_inout.data<T>(), this->work.data<T>(), local_sz);

    rotate_loop<T>(nb, nbs, this->one_,
                   this->tmp_hsub.data<T>(), hpsi_inout.data<T>(),
                   this->zero_, this->work.data<T>());
    syncmem_complex_op()(hpsi_inout.data<T>(), this->work.data<T>(), local_sz);
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::generalized_rayleigh_ritz(
    ct::Tensor&, ct::Tensor&, ct::Tensor&, ct::Tensor&)
{
    ModuleBase::WARNING_QUIT("DiagoLobpcg",
        "Generalized R-R (USPP) is not implemented yet.");
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

    for (int ib = 0; ib < this->n_band_l; ib++) {
        const int  ioff   = ib * this->n_basis;
        const Real lambda = _eigen[ib];
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
    const ct::Tensor&, const ct::Tensor&, const ct::Tensor&,
    const ct::Tensor&, const ct::Tensor&, ct::Tensor&, ct::Tensor&)
{
    ModuleBase::WARNING_QUIT("DiagoLobpcg",
        "compute_residual_s (USPP) is not implemented yet.");
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
    const int nb  = this->n_band;
    const int nbs = this->n_basis;

    inner_product_loop<T>(nb, nbs, this->one_,
                          psi_in.data<T>(), grad_out.data<T>(),
                          this->zero_, hsub_work.data<T>());
#ifdef __MPI
    Parallel_Reduce::reduce_pool(hsub_work.data<T>(), nb * nb);
#endif

    const T* inner = hsub_work.data<T>();
    const T* psi   = psi_in.data<T>();
    T*       grad  = grad_out.data<T>();
    for (int jb = 0; jb < this->n_band_l; jb++) {
        for (int ig = 0; ig < nbs; ig++) {
            T sum = static_cast<T>(0.0);
            for (int ib = 0; ib < nb; ib++)
                sum += psi[ib * nbs + ig] * inner[jb * nb + ib];
            grad[jb * nbs + ig] -= sum;
        }
    }
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::orth_projection_s(
    const ct::Tensor&, const ct::Tensor&, ct::Tensor&, ct::Tensor&)
{
    ModuleBase::WARNING_QUIT("DiagoLobpcg",
        "orth_projection_s (USPP) is not implemented yet.");
}

// ============================================================================
// rotate_wf
// ============================================================================

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::rotate_wf(
    const ct::Tensor& hsub_in, ct::Tensor& psi_out, ct::Tensor& workspace_in)
{
    const int nb = this->n_band;
    const int nbs = this->n_basis;
    rotate_loop<T>(nb, nbs, this->one_,
                   hsub_in.data<T>(), psi_out.data<T>(),
                   this->zero_, workspace_in.data<T>());
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

    T* psi_d  = psi_out.data<T>();
    T* hpsi_d = hpsi_out.data<T>();
    const Real eps = static_cast<Real>(100)
                   * std::numeric_limits<Real>::epsilon();

    for (int ib = 0; ib < nb; ++ib) {
        for (int jb = 0; jb < ib; ++jb) {
            T dot = static_cast<T>(0.0);
            for (int ig = 0; ig < nbs; ++ig) {
                dot += std::conj(psi_d[jb * nbs + ig])
                     * psi_d[ib * nbs + ig];
            }
#ifdef __MPI
            Parallel_Reduce::reduce_pool(&dot, 1);
#endif
            for (int ig = 0; ig < nbs; ++ig) {
                psi_d[ib * nbs + ig]  -= dot * psi_d[jb * nbs + ig];
                hpsi_d[ib * nbs + ig] -= dot * hpsi_d[jb * nbs + ig];
            }
        }

        Real norm2 = static_cast<Real>(0.0);
        for (int ig = 0; ig < nbs; ++ig) {
            norm2 += std::norm(psi_d[ib * nbs + ig]);
        }
#ifdef __MPI
        Parallel_Reduce::reduce_pool(&norm2, 1);
#endif
        if (!(norm2 > eps)) {
            throw std::runtime_error("orth_cholesky failed: dependent vectors");
        }
        const Real inv_norm = static_cast<Real>(1.0) / std::sqrt(norm2);
        for (int ig = 0; ig < nbs; ++ig) {
            psi_d[ib * nbs + ig]  *= inv_norm;
            hpsi_d[ib * nbs + ig] *= inv_norm;
        }
    }
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::orth_cholesky_s(
    ct::Tensor&, ct::Tensor&, ct::Tensor&, ct::Tensor&, ct::Tensor&)
{
    ModuleBase::WARNING_QUIT("DiagoLobpcg",
        "orth_cholesky_s (USPP) is not implemented yet.");
}

// ============================================================================
// test_error
// ============================================================================

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
    const int local_sz = this->n_band_l * nbs;
    const Real eps = static_cast<Real>(100)
                   * std::numeric_limits<Real>::epsilon();

    std::vector<T> basis;
    std::vector<T> hbasis;
    basis.reserve((this->has_pdir ? 3 : 2) * local_sz);
    hbasis.reserve((this->has_pdir ? 3 : 2) * local_sz);

    append_orthonormal_block<T>(n, nbs, eps, psi.data<T>(), hpsi.data<T>(),
                                basis, hbasis);
    append_orthonormal_block<T>(n, nbs, eps, grad.data<T>(), hgrad.data<T>(),
                                basis, hbasis);
    if (this->has_pdir) {
        append_orthonormal_block<T>(n, nbs, eps, pdir.data<T>(), hpdir.data<T>(),
                                    basis, hbasis);
    }

    const int m = static_cast<int>(basis.size() / nbs);
    if (m < n) {
        throw std::runtime_error("LOBPCG subspace lost rank");
    }

    T* hsub_d = this->hsub.data<T>();
    setmem_complex_op()(hsub_d, static_cast<T>(0.0), this->nsub * this->nsub);

    for (int jc = 0; jc < m; ++jc) {
        for (int ir = 0; ir < m; ++ir) {
            T sum = static_cast<T>(0.0);
            for (int ig = 0; ig < nbs; ++ig) {
                sum += std::conj(basis[ir * nbs + ig])
                     * hbasis[jc * nbs + ig];
            }
            hsub_d[jc * this->nsub + ir] = sum;
        }
    }
#ifdef __MPI
    for (int jc = 0; jc < m; ++jc) {
        Parallel_Reduce::reduce_pool(hsub_d + jc * this->nsub, m);
    }
#endif
    clean_hermitian_diag(hsub_d, this->nsub, m);
    for (int jc = 0; jc < m; ++jc) {
        for (int ir = jc + 1; ir < m; ++ir) {
            const T avg = static_cast<T>(0.5)
                        * (hsub_d[jc * this->nsub + ir]
                        +  std::conj(hsub_d[ir * this->nsub + jc]));
            hsub_d[jc * this->nsub + ir] = avg;
            hsub_d[ir * this->nsub + jc] = std::conj(avg);
        }
    }

    ct::kernels::lapack_heevd<T, ct_Device>()(
        m, hsub_d, this->nsub, this->sub_eigen.data<Real>());

    const Real* sub = this->sub_eigen.data<Real>();
    Real* eig = eigen.data<Real>();
    for (int ib = 0; ib < n; ++ib) {
        eig[ib] = sub[ib];
    }

    T* x_new = this->work.data<T>();
    T* hx_new = this->hwork.data<T>();
    T* p_new = this->pwork.data<T>();
    T* hp_new = this->hpwork.data<T>();
    for (int i = 0; i < local_sz; ++i) {
        x_new[i] = static_cast<T>(0.0);
        hx_new[i] = static_cast<T>(0.0);
        p_new[i] = static_cast<T>(0.0);
        hp_new[i] = static_cast<T>(0.0);
    }

    for (int ib = 0; ib < n; ++ib) {
        for (int iq = 0; iq < m; ++iq) {
            const T coeff = hsub_d[ib * this->nsub + iq];
            for (int ig = 0; ig < nbs; ++ig) {
                const int dst = ib * nbs + ig;
                const int src = iq * nbs + ig;
                x_new[dst] += coeff * basis[src];
                hx_new[dst] += coeff * hbasis[src];
            }
            if (iq >= n) {
                for (int ig = 0; ig < nbs; ++ig) {
                    const int dst = ib * nbs + ig;
                    const int src = iq * nbs + ig;
                    p_new[dst] += coeff * basis[src];
                    hp_new[dst] += coeff * hbasis[src];
                }
            }
        }
    }

    syncmem_complex_op()(psi.data<T>(),   x_new,  local_sz);
    syncmem_complex_op()(hpsi.data<T>(),  hx_new, local_sz);
    syncmem_complex_op()(pdir.data<T>(),  p_new,  local_sz);
    syncmem_complex_op()(hpdir.data<T>(), hp_new, local_sz);

    this->has_pdir = true;
}

// ============================================================================
// diag — main LOBPCG loop (NC, S=I; n_band_l == n_band required)
// ============================================================================

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::diag(
    const HPsiFunc& hpsi_func, T* psi_in,
    Real* eigenvalue_in, const std::vector<double>& ethr_band)
{
    // ---- runtime guard -----------------------------------------------------
    if (this->n_band_l != this->n_band) {
        ModuleBase::WARNING_QUIT("DiagoLobpcg",
            "Explicit-loop LOBPCG requires n_band_l == n_band. "
            "Band-parallel MPI is not supported in this implementation.");
    }

    this->has_pdir = false;
    const int scf_iter = DiagoIterAssist<T, Device>::SCF_ITER;

    this->psi = ct::TensorMap(psi_in, t_type, dev_type,
                               {this->n_band_l, this->n_basis});

    this->calc_prec();

    this->calc_hpsi_with_block(hpsi_func, psi_in, this->hpsi);
    // Re-orthonormalize before initial R-R so H_sub is well-conditioned
    this->orth_cholesky(this->work, this->psi, this->hpsi, this->tmp_hsub);
    this->rayleigh_ritz(this->psi, this->hpsi, this->eigen);

    setmem_complex_op()(this->pdir.data<T>(),  static_cast<T>(0.0),
                         this->n_basis * this->n_band_l);
    setmem_complex_op()(this->hpdir.data<T>(), static_cast<T>(0.0),
                         this->n_basis * this->n_band_l);

    const int max_iter = (scf_iter > 1) ? this->nline : (this->nline * 20);

    for (int ntry = 0; ntry < max_iter; ++ntry) {
        this->compute_residual(this->psi, this->hpsi, this->eigen,
                               this->prec, this->grad, this->err_st);
        if (!this->test_error(this->err_st, ethr_band))
            break;

        const int psi_sz = this->n_basis * this->n_band_l;
        const int eig_sz = this->n_band_l;

        this->orth_projection(this->psi, this->tmp_hsub, this->grad);

        this->calc_hpsi_with_block(hpsi_func, this->grad.data<T>(), this->hgrad);

        // Backup stable state in case lobpcg_update corrupts psi/hpsi
        std::vector<T> psi_bak(psi_sz), hpsi_bak(psi_sz);
        std::vector<Real> eigen_bak(eig_sz);
        std::copy(this->psi.data<T>(),  this->psi.data<T>()  + psi_sz, psi_bak.data());
        std::copy(this->hpsi.data<T>(), this->hpsi.data<T>() + psi_sz, hpsi_bak.data());
        std::copy(this->eigen.data<Real>(), this->eigen.data<Real>() + eig_sz, eigen_bak.data());

        try {
            this->lobpcg_update(this->psi, this->hpsi,
                                 this->grad, this->hgrad,
                                 this->pdir, this->hpdir,
                                 this->eigen);
        } catch (const std::exception& e1) {
            std::copy(psi_bak.data(), psi_bak.data() + psi_sz, this->psi.data<T>());
            std::copy(hpsi_bak.data(), hpsi_bak.data() + psi_sz, this->hpsi.data<T>());
            std::copy(eigen_bak.data(), eigen_bak.data() + eig_sz, this->eigen.data<Real>());

            setmem_complex_op()(this->pdir.data<T>(),  static_cast<T>(0.0), psi_sz);
            setmem_complex_op()(this->hpdir.data<T>(), static_cast<T>(0.0), psi_sz);
            this->has_pdir = false;

            try {
                this->lobpcg_update(this->psi, this->hpsi,
                                     this->grad, this->hgrad,
                                     this->pdir, this->hpdir,
                                     this->eigen);
            } catch (const std::exception& e2) {
                std::copy(psi_bak.data(), psi_bak.data() + psi_sz, this->psi.data<T>());
                std::copy(hpsi_bak.data(), hpsi_bak.data() + psi_sz, this->hpsi.data<T>());
                std::copy(eigen_bak.data(), eigen_bak.data() + eig_sz, this->eigen.data<Real>());

                this->calc_hpsi_with_block(hpsi_func, this->psi.data<T>(), this->hpsi);
                this->orth_cholesky(this->work, this->psi, this->hpsi, this->tmp_hsub);
                this->rayleigh_ritz(this->psi, this->hpsi, this->eigen);
            }
        }

        this->calc_hpsi_with_block(hpsi_func, this->psi.data<T>(), this->hpsi);
        this->orth_cholesky(this->work, this->psi, this->hpsi, this->tmp_hsub);
        this->rayleigh_ritz(this->psi, this->hpsi, this->eigen);
        this->orth_projection(this->psi, this->tmp_hsub, this->pdir);
        this->calc_hpsi_with_block(hpsi_func, this->pdir.data<T>(), this->hpdir);

        if (scf_iter == 1 && ((ntry + 1) % this->nline == 0)) {
            this->calc_hpsi_with_block(hpsi_func, this->psi.data<T>(), this->hpsi);
            this->orth_cholesky(this->work, this->psi, this->hpsi, this->tmp_hsub);
            this->rayleigh_ritz(this->psi, this->hpsi, this->eigen);
            setmem_complex_op()(this->pdir.data<T>(),  static_cast<T>(0.0),
                                 this->n_basis * this->n_band_l);
            setmem_complex_op()(this->hpdir.data<T>(), static_cast<T>(0.0),
                                 this->n_basis * this->n_band_l);
            this->has_pdir = false;
        }
    }

    this->calc_hpsi_with_block(hpsi_func, this->psi.data<T>(), this->hpsi);
    this->orth_cholesky(this->work, this->psi, this->hpsi, this->tmp_hsub);
    this->rayleigh_ritz(this->psi, this->hpsi, this->eigen);
    // Ensure hpsi = H*psi is consistent for external residual checks
    this->calc_hpsi_with_block(hpsi_func, this->psi.data<T>(), this->hpsi);

    syncmem_var_d2h_op()(eigenvalue_in,
                          this->eigen.data<Real>(),
                          this->n_band_l);
}

template <typename T, typename Device>
void DiagoLobpcg<T, Device>::diag(
    const HPsiFunc&, const SPsiFunc&, T*, Real*, const std::vector<double>&)
{
    ModuleBase::WARNING_QUIT("DiagoLobpcg",
        "USPP/generalized LOBPCG (S != I) is not implemented yet.");
}

template class DiagoLobpcg<std::complex<double>, base_device::DEVICE_CPU>;

} // namespace hsolver
