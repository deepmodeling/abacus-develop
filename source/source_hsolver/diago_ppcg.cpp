#include "source_hsolver/diago_ppcg.h"

#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_base/tool_quit.h"
#include "source_hsolver/diago_iter_assist.h"
#include "source_hsolver/module_diag/diago_trace.h"

#include <ATen/kernels/lapack.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace hsolver
{

// ---- tiny helpers -----------------------------------------------------------
template <typename T>
static const T* p_one()
{
    static const T o = static_cast<T>(1.0);
    return &o;
}
template <typename T>
static const T* p_zero()
{
    static const T z = static_cast<T>(0.0);
    return &z;
}

// ---- constructor / destructor / init_iter -----------------------------------

template <typename T, typename Device>
DiagoPPCG<T, Device>::DiagoPPCG(const Real* precondition_in) : precondition(precondition_in)
{
    this->device = base_device::get_device_type(this->ctx);
}

template <typename T, typename Device>
DiagoPPCG<T, Device>::~DiagoPPCG()
{
    delmem_op()(hpsi);
    delmem_op()(w);
    delmem_op()(hw);
    delmem_op()(p);
    delmem_op()(hp);
    delmem_op()(p_new);
    delmem_op()(hp_new);
    delmem_op()(hpsi_new);
    delmem_op()(work);
    delmem_op()(d_bv_cache);
    delmem_op()(d_tmp_cache);
    delmem_op()(d_pack_basis);
    delmem_op()(d_pack_hprod);
    delmem_op()(d_block_h);
    delmem_op()(d_block_s);
    delmem_real_op()(d_eigen);
    delmem_real_op()(d_err);
    delmem_real_h()(h_eigen);
    delmem_real_h()(h_err);
#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
        delmem_real_op()(d_precondition);
#endif
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::init_iter(const int nband,
                                     const int nband_l,
                                     const int nbasis,
                                     const int ndim)
{
    this->n_band   = nband;
    this->n_band_l = nband_l;
    this->n_basis  = nbasis;
    this->n_dim    = ndim;
    this->n_work   = this->n_band_l + this->n_extra;

    const int bs = this->n_work * this->n_basis;

    // free any previous allocation
    delmem_op()(hpsi);     delmem_op()(w);      delmem_op()(hw);
    delmem_op()(p);        delmem_op()(hp);     delmem_op()(p_new);
    delmem_op()(hp_new);   delmem_op()(hpsi_new); delmem_op()(work);
    delmem_real_op()(d_eigen);  delmem_real_op()(d_err);
    delmem_real_h()(h_eigen);  delmem_real_h()(h_err);

    // allocate & zero device buffers
    resmem_op()(hpsi, bs);     setmem_op()(hpsi, 0, bs);
    resmem_op()(w, bs);        setmem_op()(w, 0, bs);
    resmem_op()(hw, bs);       setmem_op()(hw, 0, bs);
    resmem_op()(p, bs);        setmem_op()(p, 0, bs);
    resmem_op()(hp, bs);       setmem_op()(hp, 0, bs);
    resmem_op()(p_new, bs);    setmem_op()(p_new, 0, bs);
    resmem_op()(hp_new, bs);   setmem_op()(hp_new, 0, bs);
    resmem_op()(hpsi_new, bs); setmem_op()(hpsi_new, 0, bs);
    resmem_op()(work, bs);     setmem_op()(work, 0, bs);

    resmem_real_op()(d_eigen, this->n_work);
    setmem_real_op()(d_eigen, 0, this->n_work);
    resmem_real_op()(d_err, this->n_work);
    setmem_real_op()(d_err, 0, this->n_work);

    resmem_real_h()(h_eigen, this->n_work);
    resmem_real_h()(h_err, this->n_work);
    std::fill_n(h_eigen, this->n_work, Real(0));
    std::fill_n(h_err,   this->n_work, Real(0));

    // pre-allocate per-band subspace caches (B1: avoid alloc/free in inner loop)
    resmem_op()(d_bv_cache, 3 * this->n_basis);
    setmem_op()(d_bv_cache, 0, 3 * this->n_basis);
    resmem_op()(d_tmp_cache, 3);
    setmem_op()(d_tmp_cache, 0, 3);

    // pre-allocate blocked-mode pack buffers
    constexpr int k_max = 10;
    resmem_op()(d_pack_basis, 3 * k_max * this->n_basis);
    setmem_op()(d_pack_basis, 0, 3 * k_max * this->n_basis);
    resmem_op()(d_pack_hprod, 3 * k_max * this->n_basis);
    setmem_op()(d_pack_hprod, 0, 3 * k_max * this->n_basis);
    // pre-allocate Hsub/Ssub for blocked solves (max ns = 3*k_max = 30, ns2 = 900)
    resmem_op()(d_block_h, k_max * k_max * 9);
    resmem_op()(d_block_s, k_max * k_max * 9);

    this->is_locked.assign(this->n_work, 0);
    this->converge_count.assign(this->n_work, 0);
    this->ppcg_update_count = 0;

    // preconditioner: upload to device when running on GPU
#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        delmem_real_op()(d_precondition);
        resmem_real_op()(d_precondition, this->n_basis);
        syncmem_real_h2d()(d_precondition, this->precondition, this->n_basis);
    }
#endif
}

// ---- low-level vector operations --------------------------------------------

template <typename T, typename Device>
T DiagoPPCG<T, Device>::inner_product(const T* lhs, const T* rhs) const
{
    T* d_res = nullptr;
    resmem_op()(d_res, 1);
    setmem_op()(d_res, 0, 1);
    ModuleBase::gemv_op<T, Device>()('C', this->n_dim, 1,
                                     p_one<T>(), lhs, this->n_dim,
                                     rhs, 1,
                                     p_zero<T>(), d_res, 1);
    T result;
    syncmem_d2h()(&result, d_res, 1);
    delmem_op()(d_res);
    Parallel_Reduce::reduce_pool(&result, 1);
    return result;
}

template <typename T, typename Device>
typename DiagoPPCG<T, Device>::Real DiagoPPCG<T, Device>::vector_norm(const T* vec) const
{
    const Real n2 = std::max(Real(0),
                             ModuleBase::dot_real_op<T, Device>()(this->n_dim, vec, vec));
    return std::sqrt(n2);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::scale_vector(T* vec, const Real alpha) const
{
    ModuleBase::vector_mul_real_op<T, Device>()(this->n_dim, vec, vec, alpha);
    setmem_op()(vec + this->n_dim, 0, this->n_basis - this->n_dim);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::axpy_vector(T* y, const T* x, const T alpha) const
{
    T a = alpha;
    ModuleBase::axpy_op<T, Device>()(this->n_dim, &a, x, 1, y, 1);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::copy_vector(T* dst, const T* src) const
{
    syncmem_op()(dst, src, this->n_basis);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::zero_vector(T* vec) const
{
    setmem_op()(vec, 0, this->n_basis);
}

// ---- convergence test -------------------------------------------------------

template <typename T, typename Device>
bool DiagoPPCG<T, Device>::test_error(const std::vector<double>& ethr_band) const
{
    syncmem_real_d2h()(this->h_err, this->d_err, this->n_band_l);

    bool not_conv = false;
    for (int ib = 0; ib < this->n_band_l; ++ib)
        if (this->h_err[ib] > ethr_band[ib]) { not_conv = true; break; }
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &not_conv, 1, MPI_C_BOOL, MPI_LOR, BP_WORLD);
#endif
    return not_conv;
}

// ---- Hamiltonian application ------------------------------------------------

template <typename T, typename Device>
void DiagoPPCG<T, Device>::calc_hpsi(const HPsiFunc& hpsi_func,
                                     T* psi_in, T* hpsi_out) const
{
    hpsi_func(psi_in, hpsi_out, this->n_basis, this->n_work);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::calc_hpsi(const HPsiFunc& hpsi_func,
                                     T* psi_in, T* hpsi_out, int ncol) const
{
    hpsi_func(psi_in, hpsi_out, this->n_basis, ncol);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::apply_hpsi_to_active(const HPsiFunc& hpsi_func,
                                                T* vec_in, T* vec_out)
{
    // Apply H only to active (unlocked) columns.
    // Pack unlocked columns into work, apply H, scatter back, zero locked cols.
    std::vector<int> unlocked;
    unlocked.reserve(this->n_work);
    for (int ib = 0; ib < this->n_work; ++ib)
        if (!this->is_locked[ib]) unlocked.push_back(ib);

    const int nu = static_cast<int>(unlocked.size());
    if (nu == 0) return;

    // Pack -> work (reuse work buffer as temp; it will be overwritten later)
    for (int j = 0; j < nu; ++j)
    {
        const int ib = unlocked[j];
        syncmem_op()(this->work + j * this->n_basis,
                     vec_in + ib * this->n_basis, this->n_basis);
    }

    // H|work> -> hpsi_new (reused as output temp)
    setmem_op()(this->hpsi_new, 0, nu * this->n_basis);
    hpsi_func(this->work, this->hpsi_new, this->n_basis, nu);

    // Scatter back to vec_out at unlocked positions
    for (int j = 0; j < nu; ++j)
    {
        const int ib = unlocked[j];
        syncmem_op()(vec_out + ib * this->n_basis,
                     this->hpsi_new + j * this->n_basis, this->n_basis);
    }

    // Zero locked columns in output
    for (int ib = 0; ib < this->n_work; ++ib)
    {
        if (this->is_locked[ib])
            setmem_op()(vec_out + ib * this->n_basis, 0, this->n_basis);
    }
}

// ---- orthogonalization ------------------------------------------------------

template <typename T, typename Device>
void DiagoPPCG<T, Device>::modified_gram_schmidt(T* psi_in, T* hpsi_in) const
{
    for (int ib = 0; ib < this->n_work; ++ib)
    {
        T* xi  = psi_in  + ib * this->n_basis;
        T* hxi = hpsi_in + ib * this->n_basis;

        if (ib > 0)
        {
            // lagrange = psi[:,0:ib)^H * xi  -> device -> host
            T* d_lag = nullptr;
            resmem_op()(d_lag, ib);
            setmem_op()(d_lag, 0, ib);
            ModuleBase::gemv_op<T, Device>()('C', this->n_dim, ib,
                                             p_one<T>(), psi_in, this->n_basis,
                                             xi, 1, p_zero<T>(), d_lag, 1);
            std::vector<T> lag(ib);
            syncmem_d2h()(lag.data(), d_lag, ib);
            delmem_op()(d_lag);
            Parallel_Reduce::reduce_pool(lag.data(), ib);

            // upload to device for gemv input
            T* d_lag2 = nullptr;
            resmem_op()(d_lag2, ib);
            syncmem_h2d()(d_lag2, lag.data(), ib);

            T neg1 = static_cast<T>(-1.0);
            ModuleBase::gemv_op<T, Device>()('N', this->n_dim, ib,
                                             &neg1, psi_in,  this->n_basis,
                                             d_lag2, 1, p_one<T>(), xi, 1);
            ModuleBase::gemv_op<T, Device>()('N', this->n_dim, ib,
                                             &neg1, hpsi_in, this->n_basis,
                                             d_lag2, 1, p_one<T>(), hxi, 1);
            delmem_op()(d_lag2);
        }

        const Real nrm = this->vector_norm(xi);
        if (nrm <= Real(1.0e-14))
            ModuleBase::WARNING_QUIT("DiagoPPCG::modified_gram_schmidt",
                                     "linear dependent wavefunctions");
        this->scale_vector(xi,  Real(1) / nrm);
        this->scale_vector(hxi, Real(1) / nrm);
    }
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::orth_cholesky(T* psi_in, T* hpsi_in)
{
    // Only orthonormalize active (unlocked) bands.
    // Locked (converged) bands must be kept exactly as-is -- rotating
    // them together with active bands would slowly drift converged
    // eigenpairs and introduce ghost eigenvalues.
    std::vector<int> unlocked;
    unlocked.reserve(this->n_work);
    for (int ib = 0; ib < this->n_work; ++ib)
        if (!this->is_locked[ib]) unlocked.push_back(ib);

    const int nu = static_cast<int>(unlocked.size());
    if (nu <= 1) return;

    const int nw = this->n_work;
    const int nl = nw - nu;  // number of locked bands

    if (nl == 0)
    {
        // ---- fast path: no locked bands, operate on all columns ----
        T* d_s = nullptr;
        resmem_op()(d_s, nw * nw);
        setmem_op()(d_s, 0, nw * nw);
        ModuleBase::gemm_op<T, Device>()('C', 'N', nw, nw, this->n_dim,
                                         p_one<T>(), psi_in, this->n_basis,
                                         psi_in, this->n_basis,
                                         p_zero<T>(), d_s, nw);
        std::vector<T> s(nw * nw);
        syncmem_d2h()(s.data(), d_s, nw * nw);
        delmem_op()(d_s);
#ifdef __MPI
        Parallel_Reduce::reduce_pool(s.data(), nw * nw);
#endif
        // Regularise S to prevent potrf failure with badly conditioned psi.
        {
            Real s_max_diag = Real(0);
            for (int i = 0; i < nw; ++i) s_max_diag = std::max(s_max_diag, std::abs(std::real(s[i + i * nw])));
            Real s_reg = std::max(Real(1e-14), s_max_diag * Real(1e-12));
            for (int i = 0; i < nw; ++i) s[i + i * nw] += T(s_reg);
        }
        ct::kernels::lapack_potrf<T, ct::DEVICE_CPU>()('U', nw, s.data(), nw);
        for (int col = 0; col < nw; ++col)
            for (int row = col + 1; row < nw; ++row)
                s[row + col * nw] = T(0);
        ct::kernels::lapack_trtri<T, ct::DEVICE_CPU>()('U', 'N', nw, s.data(), nw);
        this->rotate_block(psi_in,  s.data(), this->work);
        this->rotate_block(hpsi_in, s.data(), this->work);
    }
    else
    {
        // ---- general path: locked bands present -- only orthonormalize unlocked ones,
        //      after projecting out locked-band components ----
        // 1. Pack unlocked psi -> this->work (columns 0..nu-1)
        for (int j = 0; j < nu; ++j) {
            const int ib = unlocked[j];
            syncmem_op()(this->work + j * this->n_basis,
                         psi_in + ib * this->n_basis, this->n_basis);
        }

        // 2. Orthogonalise unlocked psi against locked psi:
        //    C = psi_locked^H * psi_unlocked  (nl x nu)
        //    psi_unlocked -= psi_locked * C
        if (nl > 0) {
            T* d_c = nullptr;
            resmem_op()(d_c, nl * nu);
            setmem_op()(d_c, 0, nl * nu);
            // Compute C using a packed locked-psi view.  Locked columns are at
            // positions 0..nw-1 that are NOT in the unlocked list.
            // For simplicity we pack locked columns into hpsi_new as scratch.
            int lj = 0;
            for (int ib = 0; ib < nw; ++ib)
                if (this->is_locked[ib])
                    syncmem_op()(this->hpsi_new + (lj++) * this->n_basis,
                                 psi_in + ib * this->n_basis, this->n_basis);
            ModuleBase::gemm_op<T, Device>()('C', 'N', nl, nu, this->n_dim,
                                             p_one<T>(), this->hpsi_new, this->n_basis,
                                             this->work, this->n_basis,
                                             p_zero<T>(), d_c, nl);
            std::vector<T> c(nl * nu);
            syncmem_d2h()(c.data(), d_c, nl * nu);
            delmem_op()(d_c);
#ifdef __MPI
            Parallel_Reduce::reduce_pool(c.data(), nl * nu);
#endif
            // psi_unlocked -= psi_locked * C   AND also correct hpsi
            T* d_c2 = nullptr;
            resmem_op()(d_c2, nl * nu);
            syncmem_h2d()(d_c2, c.data(), nl * nu);
            T neg1 = static_cast<T>(-1.0);
            // 1) psi_u -= psi_l * C   (via GEMM into work)
            ModuleBase::gemm_op<T, Device>()('N', 'N', this->n_dim, nu, nl,
                                             &neg1, this->hpsi_new, this->n_basis,
                                             d_c2, nl,
                                             p_one<T>(), this->work, this->n_basis);
            // 2) hpsi_u -= hpsi_l * C -- critical: psi correction implies hpsi
            //    must also be corrected, otherwise hpsi != H*psi after projection.
            //    hpsi_new still holds psi_l, overwrite with hpsi_l, use p_new as scratch.
            lj = 0;
            for (int ib = 0; ib < nw; ++ib)
                if (this->is_locked[ib])
                    syncmem_op()(this->hpsi_new + (lj++) * this->n_basis,
                                 hpsi_in + ib * this->n_basis, this->n_basis);
            for (int j = 0; j < nu; ++j) {
                const int ib = unlocked[j];
                syncmem_op()(this->p_new + j * this->n_basis,
                             hpsi_in + ib * this->n_basis, this->n_basis);
            }
            ModuleBase::gemm_op<T, Device>()('N', 'N', this->n_dim, nu, nl,
                                             &neg1, this->hpsi_new, this->n_basis,
                                             d_c2, nl,
                                             p_one<T>(), this->p_new, this->n_basis);
            for (int j = 0; j < nu; ++j) {
                const int ib = unlocked[j];
                syncmem_op()(hpsi_in + ib * this->n_basis,
                             this->p_new + j * this->n_basis, this->n_basis);
            }
            delmem_op()(d_c2);
        }

        // 3. S = psi_u^H * psi_u  (nu x nu)
        T* d_s = nullptr;
        resmem_op()(d_s, nu * nu);
        setmem_op()(d_s, 0, nu * nu);
        ModuleBase::gemm_op<T, Device>()('C', 'N', nu, nu, this->n_dim,
                                         p_one<T>(), this->work, this->n_basis,
                                         this->work, this->n_basis,
                                         p_zero<T>(), d_s, nu);
        std::vector<T> s(nu * nu);
        syncmem_d2h()(s.data(), d_s, nu * nu);
        delmem_op()(d_s);
#ifdef __MPI
        Parallel_Reduce::reduce_pool(s.data(), nu * nu);
#endif
        // Regularise S to prevent potrf failure with badly conditioned psi.
        {
            Real s_max_diag = Real(0);
            for (int i = 0; i < nu; ++i) s_max_diag = std::max(s_max_diag, std::abs(std::real(s[i + i * nu])));
            Real s_reg = std::max(Real(1e-14), s_max_diag * Real(1e-12));
            for (int i = 0; i < nu; ++i) s[i + i * nu] += T(s_reg);
        }

        // 4. Cholesky: R = chol(S), then R^{-1}
        ct::kernels::lapack_potrf<T, ct::DEVICE_CPU>()('U', nu, s.data(), nu);
        for (int col = 0; col < nu; ++col)
            for (int row = col + 1; row < nu; ++row)
                s[row + col * nu] = T(0);
        ct::kernels::lapack_trtri<T, ct::DEVICE_CPU>()('U', 'N', nu, s.data(), nu);

        // 5. Rotate unlocked psi: psi_u = psi_u * R^{-1}
        //    Use hpsi_new as output workspace
        {
            T* d_c = nullptr;
            resmem_op()(d_c, nu * nu);
            syncmem_h2d()(d_c, s.data(), nu * nu);
            ModuleBase::gemm_op<T, Device>()('N', 'N',
                this->n_dim, nu, nu,
                p_one<T>(), this->work, this->n_basis,
                d_c, nu,
                p_zero<T>(), this->hpsi_new, this->n_basis);
            delmem_op()(d_c);
        }
        for (int j = 0; j < nu; ++j) {
            const int ib = unlocked[j];
            syncmem_op()(psi_in + ib * this->n_basis,
                         this->hpsi_new + j * this->n_basis, this->n_basis);
        }

        // 6. Pack unlocked hpsi, rotate, scatter
        for (int j = 0; j < nu; ++j) {
            const int ib = unlocked[j];
            syncmem_op()(this->work + j * this->n_basis,
                         hpsi_in + ib * this->n_basis, this->n_basis);
        }
        {
            // Re-use s (still holds R^{-1}) -> upload again
            T* d_c = nullptr;
            resmem_op()(d_c, nu * nu);
            syncmem_h2d()(d_c, s.data(), nu * nu);
            ModuleBase::gemm_op<T, Device>()('N', 'N',
                this->n_dim, nu, nu,
                p_one<T>(), this->work, this->n_basis,
                d_c, nu,
                p_zero<T>(), this->hpsi_new, this->n_basis);
            delmem_op()(d_c);
        }
        for (int j = 0; j < nu; ++j) {
            const int ib = unlocked[j];
            syncmem_op()(hpsi_in + ib * this->n_basis,
                         this->hpsi_new + j * this->n_basis, this->n_basis);
        }
    }
}

template <typename T, typename Device>
bool DiagoPPCG<T, Device>::check_orthonormality(T* psi_in, Real ortho_thr) const
{
    const int nw = this->n_work;

    T* d_s = nullptr;
    resmem_op()(d_s, nw * nw);
    setmem_op()(d_s, 0, nw * nw);
    ModuleBase::gemm_op<T, Device>()('C', 'N', nw, nw, this->n_dim,
                                     p_one<T>(), psi_in, this->n_basis,
                                     psi_in, this->n_basis,
                                     p_zero<T>(), d_s, nw);
    std::vector<T> s(nw * nw);
    syncmem_d2h()(s.data(), d_s, nw * nw);
    delmem_op()(d_s);
#ifdef __MPI
    Parallel_Reduce::reduce_pool(s.data(), nw * nw);
#endif

    Real frob2 = 0;
    for (int col = 0; col < nw; ++col)
        for (int row = 0; row < nw; ++row)
        {
            const T delta = s[row + col * nw]
                            - static_cast<T>(row == col ? 1.0 : 0.0);
            frob2 += std::norm(delta);
        }
    return std::sqrt(frob2) < ortho_thr;
}

// ---- rotation ---------------------------------------------------------------

template <typename T, typename Device>
void DiagoPPCG<T, Device>::rotate_block(T* block, const T* coeff,
                                        T* workspace) const
{
    // GEMM writes only n_dim rows; padding (n_dim..n_basis-1) is untouched.
    // workspace (this->work) is reused across calls -- zero it first so stale
    // padding from previous operations doesn't pollute psi/hpsi after syncmem.
    setmem_op()(workspace, 0, this->n_work * this->n_basis);

    // coeff is on host (small); upload -> gemm -> copy result back
    T* d_c = nullptr;
    resmem_op()(d_c, this->n_work * this->n_work);
    syncmem_h2d()(d_c, coeff, this->n_work * this->n_work);

    ModuleBase::gemm_op<T, Device>()('N', 'N',
                                     this->n_dim, this->n_work, this->n_work,
                                     p_one<T>(), block, this->n_basis,
                                     d_c, this->n_work,
                                     p_zero<T>(), workspace, this->n_basis);
    delmem_op()(d_c);
    syncmem_op()(block, workspace, this->n_work * this->n_basis);
}

// ---- Rayleigh-Ritz ----------------------------------------------------------

template <typename T, typename Device>
void DiagoPPCG<T, Device>::rayleigh_ritz(T* psi_in, T* hpsi_in)
{
    if (this->n_work == 0) return;
    const int nw = this->n_work;

    // Hsub = psi^H (H psi) -> device -> host
    T* d_h = nullptr;
    resmem_op()(d_h, nw * nw);
    setmem_op()(d_h, 0, nw * nw);
    ModuleBase::gemm_op<T, Device>()('C', 'N', nw, nw, this->n_dim,
                                     p_one<T>(), psi_in,  this->n_basis,
                                     hpsi_in, this->n_basis,
                                     p_zero<T>(), d_h, nw);
    std::vector<T> hsub(nw * nw);
    syncmem_d2h()(hsub.data(), d_h, nw * nw);
    delmem_op()(d_h);
#ifdef __MPI
    Parallel_Reduce::reduce_pool(hsub.data(), nw * nw);
#endif

    ct::kernels::lapack_heevd<T, ct::DEVICE_CPU>()(nw, hsub.data(), nw, this->h_eigen);
    syncmem_real_h2d()(this->d_eigen, this->h_eigen, nw);

    this->rotate_block(psi_in,  hsub.data(), this->work);
    this->rotate_block(hpsi_in, hsub.data(), this->work);
}

// ---- subspace residual -------------------------------------------------------

template <typename T, typename Device>
void DiagoPPCG<T, Device>::compute_subspace_residual(T* psi_in)
{
    // Post-Cholesky / post-RR: subspace residual only for ACTIVE
    // (unlocked) bands -- G_u = psi_u^H * hpsi_u,  W_u = hpsi_u - psi_u * G_u.
    // Computing the residual against ALL columns (including locked) strips away
    // smooth locked-band components, leaving rough high-frequency noise that the
    // preconditioner amplifies, eventually making S = psi^H*psi near-singular.
    const int nw = this->n_work;
    if (nw == 0) return;

    // --- collect unlocked columns ------------------------------------------
    std::vector<int> unlocked;
    unlocked.reserve(nw);
    for (int ib = 0; ib < nw; ++ib)
        if (!this->is_locked[ib]) unlocked.push_back(ib);
    const int nu = static_cast<int>(unlocked.size());

    // zero locked W columns
    for (int ib = 0; ib < nw; ++ib) {
        if (this->is_locked[ib])
            setmem_op()(this->w + ib * this->n_basis, 0, this->n_basis);
    }
    if (nu == 0) return;

    // --- pack unlocked psi -> work, unlocked hpsi -> hpsi_new (temp) ---------
    for (int j = 0; j < nu; ++j) {
        const int ib = unlocked[j];
        syncmem_op()(this->work     + j * this->n_basis,
                     psi_in         + ib * this->n_basis, this->n_basis);
        syncmem_op()(this->hpsi_new + j * this->n_basis,
                     this->hpsi     + ib * this->n_basis, this->n_basis);
    }

    // 1. G_u = psi_u^H * hpsi_u  (nu x nu) -> device -> host -> MPI reduce
    T* d_g = nullptr;
    resmem_op()(d_g, nu * nu);
    setmem_op()(d_g, 0, nu * nu);
    ModuleBase::gemm_op<T, Device>()('C', 'N', nu, nu, this->n_dim,
                                     p_one<T>(), this->work, this->n_basis,
                                     this->hpsi_new, this->n_basis,
                                     p_zero<T>(), d_g, nu);
    std::vector<T> g(nu * nu);
    syncmem_d2h()(g.data(), d_g, nu * nu);
    delmem_op()(d_g);
#ifdef __MPI
    Parallel_Reduce::reduce_pool(g.data(), nu * nu);
#endif

    // 2. h_eigen from G diagonal
    for (int j = 0; j < nu; ++j) {
        const int ib = unlocked[j];
        this->h_eigen[ib] = std::real(g[j + j * nu]);
    }

    // 3. W_u = 1.0 * hpsi_u  -  psi_u * G_u   (write into p_new, scatter back)
    setmem_op()(this->p_new, 0, nu * this->n_basis);
    syncmem_op()(this->p_new, this->hpsi_new, nu * this->n_basis);

    T* d_g2 = nullptr;
    resmem_op()(d_g2, nu * nu);
    syncmem_h2d()(d_g2, g.data(), nu * nu);
    T neg1 = static_cast<T>(-1.0);
    ModuleBase::gemm_op<T, Device>()('N', 'N', this->n_dim, nu, nu,
                                     &neg1, this->work, this->n_basis,
                                     d_g2, nu,
                                     p_one<T>(), this->p_new, this->n_basis);
    delmem_op()(d_g2);

    // 4. Scatter W_u -> w, zero padding
    for (int j = 0; j < nu; ++j) {
        const int ib = unlocked[j];
        syncmem_op()(this->w + ib * this->n_basis,
                     this->p_new + j * this->n_basis, this->n_basis);
        setmem_op()(this->w + ib * this->n_basis + this->n_dim, 0,
                    this->n_basis - this->n_dim);
    }

}

// ---- preconditioned residual ------------------------------------------------

template <typename T, typename Device>
void DiagoPPCG<T, Device>::calc_preconditioned_residual(T* psi_in, bool skip_residual)
{
    const Real* prec = (this->device == base_device::GpuDevice)
                           ? this->d_precondition
                           : this->precondition;

    // Compute subspace residual W = hpsi - psi*(psi^H*hpsi)
    // before applying the preconditioner.  This guarantees W perp span(psi).
    // When skip_residual is true (post-RR), W was already computed in the
    // RR step, so we only need error norms + preconditioner application.
    if (!skip_residual)
        this->compute_subspace_residual(psi_in);

    // Apply preconditioner and compute per-band error norms.
    // h_err is computed from the TRUE residual (before preconditioner flips the sign).
    for (int ib = 0; ib < this->n_work; ++ib)
    {
        T* wi = this->w + ib * this->n_basis;

        if (this->is_locked[ib]) { this->zero_vector(wi); continue; }

        // err = ||wi||  (true residual, before preconditioning)
        Real e2 = ModuleBase::dot_real_op<T, Device>()(this->n_dim, wi, wi);
        Parallel_Reduce::reduce_pool(e2);
        this->h_err[ib] = std::sqrt(std::max(Real(0), e2));

        // wi = -wi / prec
        ModuleBase::vector_mul_real_op<T, Device>()(this->n_dim, wi, wi, Real(-1));
        ModuleBase::vector_div_vector_op<T, Device>()(this->n_dim, wi, wi, prec);
        setmem_op()(wi + this->n_dim, 0, this->n_basis - this->n_dim);
    }

    syncmem_real_h2d()(this->d_eigen, this->h_eigen, this->n_work);
    syncmem_real_h2d()(this->d_err,   this->h_err,   this->n_work);
}

// ---- projection -------------------------------------------------------------

template <typename T, typename Device>
void DiagoPPCG<T, Device>::project_to_orthogonal_complement(T* psi_in,
                                                            T* block) const
{
    const int nw = this->n_work;

    // C = psi^H * block -> device -> host
    T* d_c = nullptr;
    resmem_op()(d_c, nw * nw);
    setmem_op()(d_c, 0, nw * nw);
    ModuleBase::gemm_op<T, Device>()('C', 'N', nw, nw, this->n_dim,
                                     p_one<T>(), psi_in, this->n_basis,
                                     block, this->n_basis,
                                     p_zero<T>(), d_c, nw);
    std::vector<T> coeff(nw * nw);
    syncmem_d2h()(coeff.data(), d_c, nw * nw);
    delmem_op()(d_c);
#ifdef __MPI
    Parallel_Reduce::reduce_pool(coeff.data(), nw * nw);
#endif

    // block = block - psi * coeff
    T* d_c2 = nullptr;
    resmem_op()(d_c2, nw * nw);
    syncmem_h2d()(d_c2, coeff.data(), nw * nw);
    T neg1 = static_cast<T>(-1.0);
    ModuleBase::gemm_op<T, Device>()('N', 'N', this->n_dim, nw, nw,
                                     &neg1, psi_in, this->n_basis,
                                     d_c2, nw,
                                     p_one<T>(), block, this->n_basis);
    delmem_op()(d_c2);
}

// ---- small generalized eigenproblem -----------------------------------------

template <typename T, typename Device>
bool DiagoPPCG<T, Device>::solve_small_problem(const int adim,
                                               T* hsmall, T* ssmall,
                                               T* coeff, Real* eval) const
{
    std::fill(coeff, coeff + 9, T(0));
    std::fill(eval,  eval + 3,  Real(0));
    if (adim <= 1) { coeff[0] = T(1); eval[0] = std::real(hsmall[0]); return true; }

    for (int i = 0; i < adim; ++i) ssmall[i + i * adim] += T(1.0e-12);

    try {
        ct::kernels::lapack_hegvd<T, ct::DEVICE_CPU>()(adim, adim, hsmall, ssmall, eval, coeff);
    } catch (const std::exception&) {
        coeff[0] = T(1); eval[0] = std::real(hsmall[0]); return false;
    }
    return true;
}

// ---- per-band PPCG subspace update ------------------------------------------

template <typename T, typename Device>
void DiagoPPCG<T, Device>::update_vectors_from_ppcg_subspace(T* psi_in)
{
    if (!this->block_sizes.empty()) { this->update_vectors_blocked(psi_in); this->ppcg_update_count++; return; }

    setmem_op()(this->p_new,    0, this->n_work * this->n_basis);
    setmem_op()(this->hp_new,   0, this->n_work * this->n_basis);
    setmem_op()(this->hpsi_new, 0, this->n_work * this->n_basis);
    setmem_op()(this->work,     0, this->n_work * this->n_basis);  // MPI: zero padding

#ifdef __MPI
    int my_rank = 0, n_ranks = 1;
    MPI_Comm_rank(BP_WORLD, &my_rank);
    MPI_Comm_size(BP_WORLD, &n_ranks);
#endif

    // Band-group distribution: locked bands on root, unlocked bands distributed.
    for (int ib = 0; ib < this->n_work; ++ib)
    {
        T* xnew   = this->work     + ib * this->n_basis;
        T* hxnew  = this->hpsi_new + ib * this->n_basis;
        T* pnext  = this->p_new    + ib * this->n_basis;
        T* hpnext = this->hp_new   + ib * this->n_basis;

        if (this->is_locked[ib])
        {
#ifdef __MPI
            if (my_rank != 0) continue;  // only root preserves locked bands
#endif
            T* xi  = psi_in      + ib * this->n_basis;
            T* hxi = this->hpsi  + ib * this->n_basis;
            this->copy_vector(xnew, xi);
            this->copy_vector(hxnew, hxi);
            this->zero_vector(pnext);
            this->zero_vector(hpnext);
            continue;
        }

#ifdef __MPI
        // Round-robin distribution of unlocked bands
        if (ib % n_ranks != my_rank) continue;
#endif

        T* xi  = psi_in      + ib * this->n_basis;
        T* hxi = this->hpsi  + ib * this->n_basis;
        T* wi  = this->w     + ib * this->n_basis;
        T* hwi = this->hw    + ib * this->n_basis;
        T* pi  = this->p     + ib * this->n_basis;
        T* hpi = this->hp    + ib * this->n_basis;

        const Real pnrm = this->vector_norm(pi);
        const int adim = (pnrm > Real(1.0e-12)) ? 3 : 2;

        const T* bv[3]  = {xi, wi, pi};
        const T* hbv[3] = {hxi, hwi, hpi};

        T hsmall[9] = {}, ssmall[9] = {}, coeff[9] = {};
        Real eval[3] = {};

        // Pack bv into pre-allocated cache so gemv sees contiguous columns.
        setmem_op()(this->d_bv_cache, 0, adim * this->n_basis);
        for (int j = 0; j < adim; ++j)
            syncmem_op()(this->d_bv_cache + j * this->n_basis, bv[j], this->n_basis);

        for (int col = 0; col < adim; ++col)
        {
            setmem_op()(this->d_tmp_cache, 0, adim);

            // hsmall[:,col] = bv^H * hbv[col]
            ModuleBase::gemv_op<T, Device>()('C', this->n_dim, adim,
                                             p_one<T>(), this->d_bv_cache, this->n_basis,
                                             hbv[col], 1,
                                             p_zero<T>(), this->d_tmp_cache, 1);
            T hc[3]; syncmem_d2h()(hc, this->d_tmp_cache, adim);
            for (int r = 0; r < adim; ++r) hsmall[r + col * adim] = hc[r];

            // ssmall[:,col] = bv^H * bv[col]
            setmem_op()(this->d_tmp_cache, 0, adim);
            ModuleBase::gemv_op<T, Device>()('C', this->n_dim, adim,
                                             p_one<T>(), this->d_bv_cache, this->n_basis,
                                             bv[col], 1,
                                             p_zero<T>(), this->d_tmp_cache, 1);
            syncmem_d2h()(hc, this->d_tmp_cache, adim);
            for (int r = 0; r < adim; ++r) ssmall[r + col * adim] = hc[r];
        }

        this->solve_small_problem(adim, hsmall, ssmall, coeff, eval);
        this->h_eigen[ib] = eval[0];

        this->zero_vector(xnew);   this->zero_vector(hxnew);
        this->zero_vector(pnext);  this->zero_vector(hpnext);

        for (int j = 0; j < adim; ++j)
        {
            this->axpy_vector(xnew,  bv[j],  coeff[j]);
            this->axpy_vector(hxnew, hbv[j], coeff[j]);
        }
        if (adim >= 2)
        {
            this->axpy_vector(pnext,  wi,  coeff[1]);
            this->axpy_vector(hpnext, hwi, coeff[1]);
        }
        if (adim == 3)
        {
            this->axpy_vector(pnext,  pi,  coeff[2]);
            this->axpy_vector(hpnext, hpi, coeff[2]);
        }
    }

#ifdef __MPI
    // Collect partial results from all MPI ranks.
    {
        const int count = this->n_work * this->n_basis;
        MPI_Allreduce(MPI_IN_PLACE, this->work,     count, MPI_DOUBLE_COMPLEX, MPI_SUM, BP_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, this->hpsi_new, count, MPI_DOUBLE_COMPLEX, MPI_SUM, BP_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, this->p_new,    count, MPI_DOUBLE_COMPLEX, MPI_SUM, BP_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, this->hp_new,   count, MPI_DOUBLE_COMPLEX, MPI_SUM, BP_WORLD);
    }
#endif

    syncmem_op()(psi_in,  this->work,     this->n_work * this->n_basis);
    syncmem_op()(this->hpsi, this->hpsi_new, this->n_work * this->n_basis);
    syncmem_op()(this->p,    this->p_new,    this->n_work * this->n_basis);
    syncmem_op()(this->hp,   this->hp_new,   this->n_work * this->n_basis);

    syncmem_real_h2d()(this->d_eigen, this->h_eigen, this->n_work);
}

// ---- block-diagonal PPCG subspace update ------------------------------------

template <typename T, typename Device>
void DiagoPPCG<T, Device>::update_vectors_blocked(T* psi_in)
{
    setmem_op()(this->p_new,    0, this->n_work * this->n_basis);
    setmem_op()(this->hp_new,   0, this->n_work * this->n_basis);
    setmem_op()(this->hpsi_new, 0, this->n_work * this->n_basis);
    setmem_op()(this->work,     0, this->n_work * this->n_basis);  // MPI: zero padding

    const int ldb = this->n_basis;
    const int target_bs = this->block_sizes.empty()
                          ? 10
                          : std::max(1, this->block_sizes[0]);

    // ---- Phase 1: collect all unlocked bands ----
    // Subspace dimension: 2*n_band for first iteration, 3*n_band thereafter.
    std::vector<int> all_unlocked;
    all_unlocked.reserve(this->n_work);
    for (int ib = 0; ib < this->n_work; ++ib)
        if (!this->is_locked[ib]) all_unlocked.push_back(ib);

    // 2D on first call (P=0), 3D thereafter.
    const int ndim_global = (this->ppcg_update_count == 0) ? 2 : 3;

    // ---- Phase 2: shared lambda -- pack, solve, scatter one block ------------
    auto process_block = [&](const std::vector<int>& indices, int ndim_eff)
    {
        const int k = static_cast<int>(indices.size());
        if (k == 0) return;
        const int ns = ndim_eff * k, ns2 = ns * ns;

        // Check if indices are contiguous -- skip pack when possible.
        bool contiguous = true;
        for (int i = 1; i < k; ++i) {
            if (indices[i] != indices[i-1] + 1) { contiguous = false; break; }
        }

        const T* X_ptr, *W_ptr, *P_ptr, *HX_ptr, *HW_ptr, *HP_ptr;
        if (contiguous) {
            const int off = indices[0];
            X_ptr  = psi_in    + off * ldb;
            W_ptr  = this->w   + off * ldb;
            P_ptr  = this->p   + off * ldb;
            HX_ptr = this->hpsi + off * ldb;
            HW_ptr = this->hw   + off * ldb;
            HP_ptr = this->hp   + off * ldb;
        } else {
            const T* src_basis[3] = { psi_in, this->w, this->p };
            const T* src_hprod[3] = { this->hpsi, this->hw, this->hp };
            for (int dim = 0; dim < ndim_eff; ++dim) {
                for (int i = 0; i < k; ++i) {
                    int ib = indices[i];
                    syncmem_op()(d_pack_basis + (dim * k + i) * ldb,
                                 src_basis[dim] + ib * ldb, ldb);
                    syncmem_op()(d_pack_hprod + (dim * k + i) * ldb,
                                 src_hprod[dim] + ib * ldb, ldb);
                }
            }
            X_ptr  = d_pack_basis + 0*k*ldb;
            W_ptr  = d_pack_basis + 1*k*ldb;
            P_ptr  = d_pack_basis + 2*k*ldb;
            HX_ptr = d_pack_hprod + 0*k*ldb;
            HW_ptr = d_pack_hprod + 1*k*ldb;
            HP_ptr = d_pack_hprod + 2*k*ldb;
        }

        T* d_h = this->d_block_h;  setmem_op()(d_h, 0, ns2);
        T* d_s = this->d_block_s;  setmem_op()(d_s, 0, ns2);

        // Hsub upper triangle
        // (0,0): X^H HX    (0,1): X^H HW
        ModuleBase::gemm_op<T, Device>()('C','N',k,k,this->n_dim,
            p_one<T>(), X_ptr, ldb, HX_ptr, ldb,
            p_zero<T>(), d_h+0*k+0*k*ns, ns);
        ModuleBase::gemm_op<T, Device>()('C','N',k,k,this->n_dim,
            p_one<T>(), X_ptr, ldb, HW_ptr, ldb,
            p_zero<T>(), d_h+1*k*ns+0*k, ns);
        // (1,1): W^H HW
        ModuleBase::gemm_op<T, Device>()('C','N',k,k,this->n_dim,
            p_one<T>(), W_ptr, ldb, HW_ptr, ldb,
            p_zero<T>(), d_h+1*k+1*k*ns, ns);

        // Ssub upper triangle
        // (0,0): X^H X     (0,1): X^H W
        ModuleBase::gemm_op<T, Device>()('C','N',k,k,this->n_dim,
            p_one<T>(), X_ptr, ldb, X_ptr, ldb,
            p_zero<T>(), d_s+0*k+0*k*ns, ns);
        ModuleBase::gemm_op<T, Device>()('C','N',k,k,this->n_dim,
            p_one<T>(), X_ptr, ldb, W_ptr, ldb,
            p_zero<T>(), d_s+1*k*ns+0*k, ns);
        // (1,1): W^H W
        ModuleBase::gemm_op<T, Device>()('C','N',k,k,this->n_dim,
            p_one<T>(), W_ptr, ldb, W_ptr, ldb,
            p_zero<T>(), d_s+1*k+1*k*ns, ns);

        if (ndim_eff >= 3) {
            // (0,2): X^H HP    (1,2): W^H HP    (2,2): P^H HP
            ModuleBase::gemm_op<T, Device>()('C','N',k,k,this->n_dim,
                p_one<T>(), X_ptr, ldb, HP_ptr, ldb,
                p_zero<T>(), d_h+2*k*ns+0*k, ns);
            ModuleBase::gemm_op<T, Device>()('C','N',k,k,this->n_dim,
                p_one<T>(), W_ptr, ldb, HP_ptr, ldb,
                p_zero<T>(), d_h+1*k+2*k*ns, ns);
            ModuleBase::gemm_op<T, Device>()('C','N',k,k,this->n_dim,
                p_one<T>(), P_ptr, ldb, HP_ptr, ldb,
                p_zero<T>(), d_h+2*k+2*k*ns, ns);
            // (0,2): X^H P     (1,2): W^H P     (2,2): P^H P
            ModuleBase::gemm_op<T, Device>()('C','N',k,k,this->n_dim,
                p_one<T>(), X_ptr, ldb, P_ptr, ldb,
                p_zero<T>(), d_s+2*k*ns+0*k, ns);
            ModuleBase::gemm_op<T, Device>()('C','N',k,k,this->n_dim,
                p_one<T>(), W_ptr, ldb, P_ptr, ldb,
                p_zero<T>(), d_s+1*k+2*k*ns, ns);
            ModuleBase::gemm_op<T, Device>()('C','N',k,k,this->n_dim,
                p_one<T>(), P_ptr, ldb, P_ptr, ldb,
                p_zero<T>(), d_s+2*k+2*k*ns, ns);
        }

        // D2H
        std::vector<T> hv(ns2), sv(ns2);
        syncmem_d2h()(hv.data(), d_h, ns2);
        syncmem_d2h()(sv.data(), d_s, ns2);
#ifdef __MPI
        Parallel_Reduce::reduce_pool(hv.data(), ns2);
        Parallel_Reduce::reduce_pool(sv.data(), ns2);
#endif

        // Fill lower triangle by Hermitian symmetry
        for (int c = 0; c < k; ++c)
            for (int r = 0; r < k; ++r) {
                hv[(1*k+r)+(0*k+c)*ns] = std::conj(hv[(0*k+c)+(1*k+r)*ns]);
                sv[(1*k+r)+(0*k+c)*ns] = std::conj(sv[(0*k+c)+(1*k+r)*ns]);
            }
        if (ndim_eff >= 3) {
            for (int c = 0; c < k; ++c)
                for (int r = 0; r < k; ++r) {
                    hv[(2*k+r)+(0*k+c)*ns] = std::conj(hv[(0*k+c)+(2*k+r)*ns]);
                    sv[(2*k+r)+(0*k+c)*ns] = std::conj(sv[(0*k+c)+(2*k+r)*ns]);
                    hv[(2*k+r)+(1*k+c)*ns] = std::conj(hv[(1*k+c)+(2*k+r)*ns]);
                    sv[(2*k+r)+(1*k+c)*ns] = std::conj(sv[(1*k+c)+(2*k+r)*ns]);
                }
        }

        // Scale regularization by max |S_ii| to handle near-singular S
        // from P~=0 blocks. s_max ~= 1 for orthonormal X; 1e-8 relative
        // regularization prevents Cholesky failure without affecting accuracy.
        Real s_max = Real(0);
        for (int i = 0; i < ns; ++i)
            s_max = std::max(s_max, std::abs(std::real(sv[i + i * ns])));
        Real s_reg = std::max(Real(1e-11), s_max * Real(1e-9));
        for (int i = 0; i < ns; ++i) sv[i + i * ns] += T(s_reg);

        std::vector<T>   ev(ns2, T(0));
        std::vector<Real> el(ns, Real(0));
        try {
            ct::kernels::lapack_hegvd<T, ct::DEVICE_CPU>()(ns, ns, hv.data(), sv.data(),
                                                            el.data(), ev.data());
        } catch (const std::exception&) {
            for (int i = 0; i < k; ++i) {
                int ib = indices[i];
                this->copy_vector(this->work     + ib * ldb, psi_in    + ib * ldb);
                this->copy_vector(this->hpsi_new + ib * ldb, this->hpsi + ib * ldb);
            }
            return;
        }

        // Scatter updated vectors back to their original positions
        for (int i = 0; i < k; ++i)
        {
            const int ig = indices[i];
            T* xn  = this->work     + ig * ldb;
            T* hn  = this->hpsi_new + ig * ldb;
            T* pn  = this->p_new    + ig * ldb;
            T* hpn = this->hp_new   + ig * ldb;
            this->zero_vector(xn);  this->zero_vector(hn);
            this->zero_vector(pn);  this->zero_vector(hpn);

            // When contiguous, bands are is = off + cs; avoid indices[] lookup.
            if (contiguous) {
                const int off = indices[0];
                for (int col = 0; col < ns; ++col) {
                    const int cs = col % k, cb = col / k, is = off + cs;
                    const T c = ev[col + i * ns];

                    const T *vs = nullptr, *hs = nullptr;
                    if (cb == 0)       { vs = psi_in + is * ldb; hs = this->hpsi + is * ldb; }
                    else if (cb == 1)  { vs = this->w + is * ldb; hs = this->hw   + is * ldb; }
                    else               { vs = this->p + is * ldb; hs = this->hp   + is * ldb; }

                    this->axpy_vector(xn, vs, c);
                    this->axpy_vector(hn, hs, c);
                    if (cb >= 1) { this->axpy_vector(pn, vs, c); this->axpy_vector(hpn, hs, c); }
                }
            } else {
                for (int col = 0; col < ns; ++col) {
                    const int cs = col % k, cb = col / k, is = indices[cs];
                    const T c = ev[col + i * ns];

                    const T *vs = nullptr, *hs = nullptr;
                    if (cb == 0)       { vs = psi_in + is * ldb; hs = this->hpsi + is * ldb; }
                    else if (cb == 1)  { vs = this->w + is * ldb; hs = this->hw   + is * ldb; }
                    else               { vs = this->p + is * ldb; hs = this->hp   + is * ldb; }

                    this->axpy_vector(xn, vs, c);
                    this->axpy_vector(hn, hs, c);
                    if (cb >= 1) { this->axpy_vector(pn, vs, c); this->axpy_vector(hpn, hs, c); }
                }
            }
        }
    };  // end process_block

    // ---- Phase 3: distribute blocks across MPI ranks ----
    // Build the full block list, then each rank processes a round-robin subset.
    {
        std::vector<std::vector<int>> all_blocks;
        for (size_t start = 0; start < all_unlocked.size(); start += target_bs) {
            size_t end = std::min(start + target_bs, all_unlocked.size());
            all_blocks.emplace_back(all_unlocked.begin() + start,
                                    all_unlocked.begin() + end);
        }

#ifdef __MPI
        int my_rank = 0, n_ranks = 1;
        MPI_Comm_rank(BP_WORLD, &my_rank);
        MPI_Comm_size(BP_WORLD, &n_ranks);

        for (size_t bi = my_rank; bi < all_blocks.size(); bi += n_ranks)
            process_block(all_blocks[bi], ndim_global);
#else
        for (auto& block : all_blocks)
            process_block(block, ndim_global);
#endif
    }

    // ---- Phase 4: locked bands -- only root rank keeps old values -----------
    // After MPI reduction, locked values come exclusively from root.
#ifdef __MPI
    int my_rank = 0;
    MPI_Comm_rank(BP_WORLD, &my_rank);
    if (my_rank == 0)
#endif
    {
        for (int ib = 0; ib < this->n_band_l; ++ib) {
            if (!this->is_locked[ib]) continue;
            this->copy_vector(this->work     + ib * ldb, psi_in    + ib * ldb);
            this->copy_vector(this->hpsi_new + ib * ldb, this->hpsi + ib * ldb);
        }
    }

#ifdef __MPI
    // Collect partial results from all MPI ranks..
    // Only processed columns are non-zero on each rank, so SUM is correct.
    const int count = this->n_work * ldb;
    MPI_Allreduce(MPI_IN_PLACE, this->work,     count, MPI_DOUBLE_COMPLEX, MPI_SUM, BP_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, this->hpsi_new, count, MPI_DOUBLE_COMPLEX, MPI_SUM, BP_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, this->p_new,    count, MPI_DOUBLE_COMPLEX, MPI_SUM, BP_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, this->hp_new,   count, MPI_DOUBLE_COMPLEX, MPI_SUM, BP_WORLD);
#endif

    syncmem_op()(psi_in,  this->work,     this->n_work * ldb);
    syncmem_op()(this->hpsi, this->hpsi_new, this->n_work * ldb);
    syncmem_op()(this->p,    this->p_new,    this->n_work * ldb);
    syncmem_op()(this->hp,   this->hp_new,   this->n_work * ldb);

    syncmem_real_h2d()(this->d_eigen, this->h_eigen, this->n_work);
}


// ---- main diagonalization entry point ---------------------------------------

template <typename T, typename Device>
int DiagoPPCG<T, Device>::diag(const HPsiFunc& hpsi_func,
                               T* psi_in,
                               Real* eigenvalue_in,
                               const std::vector<double>& ethr_band)
{
    ModuleBase::TITLE("DiagoPPCG", "diag");
    ModuleBase::timer::start("DiagoPPCG", "diag");
    DiagoTrace trace("PPCG");

    // ---- initial orthonormalization + Rayleigh-Ritz ----
    this->calc_hpsi(hpsi_func, psi_in, this->hpsi);
    this->modified_gram_schmidt(psi_in, this->hpsi);
    this->rayleigh_ritz(psi_in, this->hpsi);

    // ---- Compute post-RR residual W = H*Psi - Psi*diag(eigenvalues) ----
    // RR has globally rotated the subspace.  We must recompute the true
    // residual from the freshly rotated Psi before any convergence decision.
    for (int ib = 0; ib < this->n_work; ++ib) {
        T* wi  = this->w + ib * this->n_basis;
        T* xi  = psi_in + ib * this->n_basis;
        T* hxi = this->hpsi + ib * this->n_basis;
        syncmem_op()(wi, hxi, this->n_dim);
        T neg_e = static_cast<T>(-this->h_eigen[ib]);
        ModuleBase::axpy_op<T, Device>()(this->n_dim, &neg_e, xi, 1, wi, 1);
        setmem_op()(wi + this->n_dim, 0, this->n_basis - this->n_dim);
    }

    // Compute h_err from post-RR W and lock converged physical bands.
    for (int ib = 0; ib < this->n_work; ++ib) {
        if (this->is_locked[ib]) { this->zero_vector(this->w + ib * this->n_basis); continue; }
        Real e2 = ModuleBase::dot_real_op<T, Device>()(this->n_dim,
                        this->w + ib * this->n_basis, this->w + ib * this->n_basis);
        Parallel_Reduce::reduce_pool(e2);
        this->h_err[ib] = std::sqrt(std::max(Real(0), e2));
    }
    syncmem_real_h2d()(this->d_err, this->h_err, this->n_work);

    // Initial locking tolerance: sqrt(ethr).
    for (int ib = 0; ib < this->n_band_l; ++ib) {
        if (this->h_err[ib] <= std::sqrt(ethr_band[ib]))
            this->is_locked[ib] = 1;
    }

    // ---- Trace convergence init ----
    // trG = Sigma e_i for active (unlocked) physical bands after initial RR.
    Real trG = 0;
    int n_act = 0;
    for (int ib = 0; ib < this->n_band_l; ++ib) {
        if (!this->is_locked[ib]) { trG += this->h_eigen[ib]; n_act++; }
    }
    // Trace convergence tolerance: trtol = ethr * sqrt(nact).
    Real trtol = (n_act > 0) ? ethr_band[0] * std::sqrt(Real(n_act)) : Real(0);
    Real trdif = Real(-1);  // -1 = "undefined", always trigger at least one more iter

    int iter = 0;
    const int max_iter = std::max(1, DiagoIterAssist<T, Device>::PW_DIAG_NMAX);
    const int rr_period = 20;

    // did_rr: true when the previous iteration ended with an RR step.
    bool did_rr = false;

    for (; iter < max_iter; ++iter)
    {
        // ---- 1. preconditioned residuals ----
        this->calc_preconditioned_residual(psi_in, /*skip_residual=*/did_rr);
        did_rr = false;

        if (trace.enabled())
        {
            Real max_residual = Real(0);
            Real avg_residual = Real(0);
            int n_converged = 0;
            for (int ib = 0; ib < this->n_band_l; ++ib)
            {
                max_residual = std::max(max_residual, this->h_err[ib]);
                avg_residual += this->h_err[ib];
                if (this->is_locked[ib])
                {
                    ++n_converged;
                }
            }
            if (this->n_band_l > 0)
            {
                avg_residual /= this->n_band_l;
            }
            trace.record_iteration(iter,
                                   this->n_band_l,
                                   max_residual,
                                   avg_residual,
                                   n_converged,
                                   Real(-1),
                                   std::string("trdif=") + std::to_string(static_cast<double>(trdif))
                                       + " trtol=" + std::to_string(static_cast<double>(trtol))
                                       + (!this->block_sizes.empty() ? " blocked" : ""));
        }

        // ---- 2. convergence: per-band residual OR trace stabilised ----
        if (!this->test_error(ethr_band)) break;
        if (trdif >= Real(0) && trdif <= trtol) {
            break;
        }

        // ---- 3. project W, P to orthogonal complement ----
        this->project_to_orthogonal_complement(psi_in, this->w);
        this->project_to_orthogonal_complement(psi_in, this->p);

        // ---- 4. H|w>, H|p> (only active/unlocked columns) ----
        this->apply_hpsi_to_active(hpsi_func, this->w, this->hw);
        this->apply_hpsi_to_active(hpsi_func, this->p, this->hp);

        // ---- 5. subspace update ----
        this->update_vectors_from_ppcg_subspace(psi_in);

        // ---- 6. periodic Rayleigh-Ritz + locking (paper sec.3.4) ----
        if ((iter + 1) % rr_period == 0)
        {
            this->orth_cholesky(psi_in, this->hpsi);
            this->rayleigh_ritz(psi_in, this->hpsi);

            // ---- Recompute W = HPsi - Psi*diag(eigenvalues) after RR ----
            for (int ib = 0; ib < this->n_work; ++ib) {
                T* wi  = this->w + ib * this->n_basis;
                T* xi  = psi_in + ib * this->n_basis;
                T* hxi = this->hpsi + ib * this->n_basis;
                syncmem_op()(wi, hxi, this->n_dim);
                T neg_e = static_cast<T>(-this->h_eigen[ib]);
                ModuleBase::axpy_op<T, Device>()(this->n_dim, &neg_e, xi, 1, wi, 1);
                setmem_op()(wi + this->n_dim, 0, this->n_basis - this->n_dim);
            }

            // ---- Lock converged physical bands based on post-RR residual ----
            // Use sqrt(ethr) as lock tolerance.
            std::fill(this->is_locked.begin(), this->is_locked.end(), 0);
            for (int ib = 0; ib < this->n_band_l; ++ib) {
                Real e2 = ModuleBase::dot_real_op<T, Device>()(this->n_dim,
                                this->w + ib * this->n_basis, this->w + ib * this->n_basis);
                Parallel_Reduce::reduce_pool(e2);
                this->h_err[ib] = std::sqrt(std::max(Real(0), e2));
                if (this->h_err[ib] <= std::sqrt(ethr_band[ib]))
                    this->is_locked[ib] = 1;
            }
            syncmem_real_h2d()(this->d_err, this->h_err, this->n_work);

            // ---- After RR, trdif = -1, trG = sum e_i(active) ----
            trdif = Real(-1);
            trG = 0; n_act = 0;
            for (int ib = 0; ib < this->n_band_l; ++ib) {
                if (!this->is_locked[ib]) { trG += this->h_eigen[ib]; n_act++; }
            }
            trtol = (n_act > 0) ? ethr_band[0] * std::sqrt(Real(n_act)) : Real(0);

            // P directions are NOT cleared after RR -- old P directions are
            // orthogonalized against the new psi in the next iteration.
            // Clearing P would force a 2D restart and lose search info.

            did_rr = true;
        }
        else
        {
            // ---- non-RR iteration: orthonormalize + recompute subspace residual ----
            this->orth_cholesky(psi_in, this->hpsi);
            this->compute_subspace_residual(psi_in);

            // ---- Trace convergence: trG1 = sum h_eigen(active) ----
            Real trG1 = 0; n_act = 0;
            for (int ib = 0; ib < this->n_band_l; ++ib) {
                if (!this->is_locked[ib]) { trG1 += this->h_eigen[ib]; n_act++; }
            }
            trtol = (n_act > 0) ? ethr_band[0] * std::sqrt(Real(n_act)) : Real(0);
            if (n_act > 0) {
                trdif = std::abs(trG1 - trG);
                trG = trG1;
            } else {
                trdif = Real(0);  // all bands converged
            }
        }
    }

    // ---- final Rayleigh-Ritz + output ----
    this->rayleigh_ritz(psi_in, this->hpsi);
    for (int ib = 0; ib < this->n_band_l; ++ib)
        eigenvalue_in[ib] = this->h_eigen[ib];

    ModuleBase::timer::end("DiagoPPCG", "diag");

    return std::min(iter + 1, max_iter);
}

// ---- explicit template instantiations ---------------------------------------

template class DiagoPPCG<std::complex<float>,  base_device::DEVICE_CPU>;
template class DiagoPPCG<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class DiagoPPCG<std::complex<float>,  base_device::DEVICE_GPU>;
template class DiagoPPCG<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace hsolver
