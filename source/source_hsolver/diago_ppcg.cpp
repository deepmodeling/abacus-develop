#include "source_hsolver/diago_ppcg.h"

#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_base/tool_quit.h"
#include "source_hsolver/diago_iter_assist.h"
#include "source_hsolver/module_diag/diag_orthogonalizer.h"

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

// ---- orthogonalization ------------------------------------------------------

template <typename T, typename Device>
void DiagoPPCG<T, Device>::modified_gram_schmidt(T* psi_in, T* hpsi_in) const
{
    DiagOrthogonalizer<T, Device>(this->n_dim, this->n_basis)
        .modified_gram_schmidt(psi_in, hpsi_in, this->n_work);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::orth_cholesky(T* psi_in, T* hpsi_in)
{
    DiagOrthogonalizer<T, Device>(this->n_dim, this->n_basis)
        .cholesky_orth(psi_in, hpsi_in, this->work, this->n_work);
}

template <typename T, typename Device>
bool DiagoPPCG<T, Device>::check_orthonormality(T* psi_in) const
{
    return DiagOrthogonalizer<T, Device>(this->n_dim, this->n_basis)
        .check_orthonormality(psi_in, this->n_work, Real(1e-1));
}

// ---- rotation ---------------------------------------------------------------

template <typename T, typename Device>
void DiagoPPCG<T, Device>::rotate_block(T* block, const T* coeff,
                                        T* workspace) const
{
    // coeff is on host (small); upload → gemm → copy result back
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

    // Hsub = psi^H (H psi) → device → host
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

// ---- preconditioned residual ------------------------------------------------

template <typename T, typename Device>
void DiagoPPCG<T, Device>::calc_preconditioned_residual(T* psi_in)
{
    const Real* prec = (this->device == base_device::GpuDevice)
                           ? this->d_precondition
                           : this->precondition;

    for (int ib = 0; ib < this->n_work; ++ib)
    {
        T* wi  = this->w + ib * this->n_basis;
        T* xi  = psi_in   + ib * this->n_basis;
        T* hxi = this->hpsi + ib * this->n_basis;

        if (this->is_locked[ib]) { this->zero_vector(wi); continue; }

        // lambda = Re <xi | H | xi>
        const Real lam = ModuleBase::dot_real_op<T, Device>()(this->n_dim, xi, hxi);
        this->h_eigen[ib] = lam;

        // wi = hxi - lam * xi
        syncmem_op()(wi, hxi, this->n_dim);
        T nlam = static_cast<T>(-lam);
        ModuleBase::axpy_op<T, Device>()(this->n_dim, &nlam, xi, 1, wi, 1);

        // err = ||wi||
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
    DiagOrthogonalizer<T, Device>(this->n_dim, this->n_basis)
        .project_out(psi_in, block, this->n_work, this->n_work);
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
    if (!this->block_sizes.empty()) { this->update_vectors_blocked(psi_in); return; }

    setmem_op()(this->p_new,    0, this->n_work * this->n_basis);
    setmem_op()(this->hp_new,   0, this->n_work * this->n_basis);
    setmem_op()(this->hpsi_new, 0, this->n_work * this->n_basis);

    for (int ib = 0; ib < this->n_work; ++ib)
    {
        T* xi  = psi_in      + ib * this->n_basis;
        T* hxi = this->hpsi  + ib * this->n_basis;
        T* wi  = this->w     + ib * this->n_basis;
        T* hwi = this->hw    + ib * this->n_basis;
        T* pi  = this->p     + ib * this->n_basis;
        T* hpi = this->hp    + ib * this->n_basis;

        T* xnew   = this->work     + ib * this->n_basis;
        T* hxnew  = this->hpsi_new + ib * this->n_basis;
        T* pnext  = this->p_new    + ib * this->n_basis;
        T* hpnext = this->hp_new   + ib * this->n_basis;

        if (this->is_locked[ib])
        {
            this->copy_vector(xnew, xi);
            this->copy_vector(hxnew, hxi);
            this->zero_vector(pnext);
            this->zero_vector(hpnext);
            continue;
        }

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

    const int ldb = this->n_basis;
    const int target_bs = this->block_sizes.empty()
                          ? 10
                          : std::max(1, this->block_sizes[0]);

    // ---- Phase 1: classify unlocked bands by P-norm (2D vs 3D subspace) ----
    std::vector<int> idx_2d, idx_3d;
    idx_2d.reserve(this->n_band_l);
    idx_3d.reserve(this->n_band_l);

    for (int ib = 0; ib < this->n_band_l; ++ib)
    {
        if (this->is_locked[ib]) continue;

        // Per-band P-norm check — same threshold as per-band solver (adim=2 vs 3).
        Real p_norm2 = 0;
        {
            const T* pi = this->p + ib * ldb;
            for (int ig = 0; ig < this->n_dim; ++ig) {
                const T& v = pi[ig];
                p_norm2 += std::real(v) * std::real(v) + std::imag(v) * std::imag(v);
            }
        }
#ifdef __MPI
        Parallel_Reduce::reduce_pool(p_norm2);
#endif
        if (p_norm2 < Real(1e-30))
            idx_2d.push_back(ib);
        else
            idx_3d.push_back(ib);
    }

    // ---- Phase 2: shared lambda — pack, solve, scatter one block ------------
    auto process_block = [&](const std::vector<int>& indices, int ndim_eff)
    {
        const int k = static_cast<int>(indices.size());
        if (k == 0) return;
        const int ns = ndim_eff * k, ns2 = ns * ns;

        // Check if indices are contiguous — skip pack when possible.
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

        for (int i = 0; i < ns; ++i) sv[i + i * ns] += T(1.0e-12);

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

    // ---- Phase 3: process 2D and 3D groups in blocks -----------------------
    for (size_t start = 0; start < idx_2d.size(); start += target_bs)
    {
        size_t end = std::min(start + target_bs, idx_2d.size());
        std::vector<int> block(idx_2d.begin() + start, idx_2d.begin() + end);
        process_block(block, 2);
    }
    for (size_t start = 0; start < idx_3d.size(); start += target_bs)
    {
        size_t end = std::min(start + target_bs, idx_3d.size());
        std::vector<int> block(idx_3d.begin() + start, idx_3d.begin() + end);
        process_block(block, 3);
    }

    // ---- Phase 4: locked bands — keep old values ---------------------------
    for (int ib = 0; ib < this->n_band_l; ++ib)
    {
        if (!this->is_locked[ib]) continue;
        this->copy_vector(this->work     + ib * ldb, psi_in    + ib * ldb);
        this->copy_vector(this->hpsi_new + ib * ldb, this->hpsi + ib * ldb);
    }

    // ---- Phase 5: extra (buffer) bands — per-band PPCG ---------------------
    for (int ib = this->n_band_l; ib < this->n_work; ++ib)
    {
        T* xi  = psi_in      + ib * ldb;
        T* hxi = this->hpsi  + ib * ldb;
        T* wi  = this->w     + ib * ldb;
        T* hwi = this->hw    + ib * ldb;
        T* pi  = this->p     + ib * ldb;
        T* hpi = this->hp    + ib * ldb;

        T* xnew   = this->work     + ib * ldb;
        T* hxnew  = this->hpsi_new + ib * ldb;
        T* pnext  = this->p_new    + ib * ldb;
        T* hpnext = this->hp_new   + ib * ldb;

        if (this->is_locked[ib]) {
            this->copy_vector(xnew, xi);
            this->copy_vector(hxnew, hxi);
            continue;
        }

        T* bv[3]  = { xi,  wi,  pi };
        T* hbv[3] = { hxi, hwi, hpi };

        Real p_norm = this->vector_norm(pi);
        int  adim = (p_norm > Real(1e-15)) ? 3 : 2;

        setmem_op()(this->d_bv_cache, 0, adim * ldb);
        for (int j = 0; j < adim; ++j)
            syncmem_op()(this->d_bv_cache + j * ldb, bv[j], ldb);

        T hsmall[9], ssmall[9], coeff[9];
        setmem_op()(this->d_tmp_cache, 0, 3);
        for (int col = 0; col < adim; ++col) {
            ModuleBase::gemv_op<T, Device>()('C', this->n_dim, adim,
                p_one<T>(), this->d_bv_cache, ldb, hbv[col], 1,
                p_zero<T>(), this->d_tmp_cache, 1);
            T hc[3]; syncmem_d2h()(hc, this->d_tmp_cache, adim);
            for (int r = 0; r < adim; ++r) hsmall[r + col * adim] = hc[r];

            setmem_op()(this->d_tmp_cache, 0, 3);
            ModuleBase::gemv_op<T, Device>()('C', this->n_dim, adim,
                p_one<T>(), this->d_bv_cache, ldb, bv[col], 1,
                p_zero<T>(), this->d_tmp_cache, 1);
            syncmem_d2h()(hc, this->d_tmp_cache, adim);
            for (int r = 0; r < adim; ++r) ssmall[r + col * adim] = hc[r];
        }

        Real eval[3];
        this->solve_small_problem(adim, hsmall, ssmall, coeff, eval);
        this->h_eigen[ib] = eval[0];

        this->zero_vector(xnew);   this->zero_vector(hxnew);
        this->zero_vector(pnext);  this->zero_vector(hpnext);

        for (int j = 0; j < adim; ++j) {
            this->axpy_vector(xnew,  bv[j],  coeff[j]);
            this->axpy_vector(hxnew, hbv[j], coeff[j]);
        }
        if (adim >= 2) {
            this->axpy_vector(pnext,  wi,  coeff[1]);
            this->axpy_vector(hpnext, hwi, coeff[1]);
        }
        if (adim == 3) {
            this->axpy_vector(pnext,  pi,  coeff[2]);
            this->axpy_vector(hpnext, hpi, coeff[2]);
        }
    }

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

    // ---- initial orthonormalization + Rayleigh-Ritz ----
    this->calc_hpsi(hpsi_func, psi_in, this->hpsi);
    this->modified_gram_schmidt(psi_in, this->hpsi);
    this->rayleigh_ritz(psi_in, this->hpsi);

    int iter = 0;
    const int max_iter = std::max(1, DiagoIterAssist<T, Device>::PW_DIAG_NMAX);
    for (; iter < max_iter; ++iter)
    {
        // 1. preconditioned residuals
        this->calc_preconditioned_residual(psi_in);

        // diagnostics
        if (iter % 10 == 0 || iter == max_iter - 1)
        {
            int nl = 0;
            for (int ib = 0; ib < this->n_band_l; ++ib)
                if (this->is_locked[ib]) nl++;
            std::cerr << "[PPCG] iter=" << iter
                      << " err[0]=" << this->h_err[0]
                      << " err[end]=" << this->h_err[this->n_band_l - 1]
                      << " ethr=" << ethr_band[0]
                      << " locked=" << nl << "/" << this->n_band_l
                      << " blocked=" << (!this->block_sizes.empty() ? "yes" : "no")
                      << " dev=" << (this->device == base_device::GpuDevice ? "GPU" : "CPU")
                      << std::endl;
        }

        // 2. lock converged bands
        for (int ib = 0; ib < this->n_band_l; ++ib)
        {
            if (this->is_locked[ib]) continue;
            if (this->h_err[ib] <= ethr_band[ib])
            {
                if (++this->converge_count[ib] >= 2)
                {
                    this->is_locked[ib] = 1;
                    this->h_err[ib] = Real(0);
                }
            }
            else this->converge_count[ib] = 0;
        }

        // 3. global convergence
        if (!this->test_error(ethr_band)) break;

        // 4. project W, P to orthogonal complement
        this->project_to_orthogonal_complement(psi_in, this->w);
        this->project_to_orthogonal_complement(psi_in, this->p);

        // 5. H|w>, H|p>
        this->calc_hpsi(hpsi_func, this->w, this->hw);
        this->calc_hpsi(hpsi_func, this->p, this->hp);

        // 6. subspace update
        this->update_vectors_from_ppcg_subspace(psi_in);

        // 7. periodic re-orthonormalization
        if ((iter + 1) % 15 == 0)
        {
            this->orth_cholesky(psi_in, this->hpsi);
            this->rayleigh_ritz(psi_in, this->hpsi);
        }
        else if (!this->check_orthonormality(psi_in))
        {
            this->orth_cholesky(psi_in, this->hpsi);
        }
    }

    // final Rayleigh-Ritz + output
    this->rayleigh_ritz(psi_in, this->hpsi);
    for (int ib = 0; ib < this->n_band_l; ++ib)
        eigenvalue_in[ib] = this->h_eigen[ib];

    ModuleBase::timer::end("DiagoPPCG", "diag");

    std::cerr << "[PPCG] done: niter=" << std::min(iter + 1, max_iter)
              << " final_err[0]=" << this->h_err[0]
              << " final_err[end]=" << this->h_err[this->n_band_l - 1]
              << " eigen[0]=" << eigenvalue_in[0] << std::endl;

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
