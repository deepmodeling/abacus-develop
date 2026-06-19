#include "source_hsolver/diago_ppcg.h"

#include "diago_iter_assist.h"
#include "source_base/global_variable.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_comm.h" // POOL_WORLD/BP_WORLD
#include "source_base/parallel_reduce.h"
#include "source_hsolver/kernels/bpcg_kernel_op.h" // reuse normalize_op / apply_eigenvalues_op / precondition_op

#include "source_base/module_container/base/third_party/lapack.h"

#ifdef __MPI
#include <mpi.h>
#endif

namespace hsolver {

namespace lapackConnector = container::lapackConnector;

template <typename T, typename Device>
DiagoPPCG<T, Device>::DiagoPPCG(const Real* precondition)
{
    this->device_type = ct::DeviceTypeToEnum<Device>::value;
    this->ctx = {}; // default device context
    this->r_type = ct::DataTypeToEnum<Real>::value;
    this->t_type = ct::DataTypeToEnum<T>::value;

    this->h_prec = std::move(ct::TensorMap((void*)precondition, this->r_type, this->device_type, {this->n_basis}));
}

template <typename T, typename Device>
DiagoPPCG<T, Device>::~DiagoPPCG() = default;

template <typename T, typename Device>
void DiagoPPCG<T, Device>::init_iter(const int nband, const int nband_l, const int nbasis, const int ndim)
{
    this->n_band = nband;
    this->n_band_l = nband_l;
    this->n_basis = nbasis;
    this->n_dim = ndim;

    this->prec = ct::Tensor(this->r_type, this->device_type, {this->n_basis});

    this->HX = ct::Tensor(this->t_type, this->device_type, {this->n_band_l, this->n_basis});
    this->R = ct::Tensor(this->t_type, this->device_type, {this->n_band_l, this->n_basis});
    this->W = ct::Tensor(this->t_type, this->device_type, {this->n_band_l, this->n_basis});
    this->HW = ct::Tensor(this->t_type, this->device_type, {this->n_band_l, this->n_basis});
    this->P = ct::Tensor(this->t_type, this->device_type, {this->n_band_l, this->n_basis});
    this->HP = ct::Tensor(this->t_type, this->device_type, {this->n_band_l, this->n_basis});

    const int max_cols = 3 * this->n_band_l;
    this->V = ct::Tensor(this->t_type, this->device_type, {max_cols, this->n_basis});
    this->HV = ct::Tensor(this->t_type, this->device_type, {max_cols, this->n_basis});

    const int max_small = 3 * this->n_band;
    this->hcc = ct::Tensor(this->t_type, this->device_type, {max_small, max_small});
    this->scc = ct::Tensor(this->t_type, this->device_type, {max_small, max_small});
    this->vcc = ct::Tensor(this->t_type, this->device_type, {max_small, max_small});
    this->eval = ct::Tensor(this->r_type, this->device_type, {max_small});
    // Zero-initialize so that uninitialised entries are immediately visible
    // as exact 0.0 rather than denormal garbage (e.g. 4.68e-310).
    Parallel_Reduce::ZEROS(this->eval.data<Real>(), max_small);

    this->work = ct::Tensor(this->t_type, this->device_type, {max_cols, this->n_basis});

    this->calc_prec();
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::calc_prec()
{
    syncmem_var_h2d_op()(this->prec.data<Real>(), this->h_prec.data<Real>(), this->n_basis);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::apply_h(const HPsiFunc& hpsi_func, const ct::Tensor& in_vecs, ct::Tensor& out_vecs, const int nvec)
{
    // hpsi_func(psi_in, hpsi_out, ld_psi, nvec)
    hpsi_func(in_vecs.data<T>(), out_vecs.data<T>(), this->n_basis, nvec);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::project_out(const ct::Tensor& basis,
                                      const int ncols_basis,
                                      ct::Tensor& vecs,
                                      const int ncols_vecs)
{
    if (ncols_basis <= 0 || ncols_vecs <= 0)
    {
        return;
    }

    // coeff = basis^H * vecs  (ncols_basis x ncols_vecs)
    const int ldh = ncols_basis;
    ct::Tensor coeff(this->t_type, this->device_type, {ldh, ncols_vecs});

#ifdef __MPI
    this->pmmcn.set_dimension(BP_WORLD,
                              POOL_WORLD,
                              ncols_basis,
                              this->n_basis,
                              ncols_vecs,
                              this->n_basis,
                              this->n_dim,
                              ldh);
#else
    this->pmmcn.set_dimension(ncols_basis,
                              this->n_basis,
                              ncols_vecs,
                              this->n_basis,
                              this->n_dim,
                              ldh);
#endif

    this->pmmcn.multiply(1.0, basis.data<T>(), vecs.data<T>(), 0.0, coeff.data<T>());

    // vecs -= basis * coeff
    ModuleBase::gemm_op<T, Device>()('N',
                                    'N',
                                    this->n_dim,
                                    ncols_vecs,
                                    ncols_basis,
                                    this->neg_one,
                                    basis.data<T>(),
                                    this->n_basis,
                                    coeff.data<T>(),
                                    ldh,
                                    this->one,
                                    vecs.data<T>(),
                                    this->n_basis);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::orthonormalize_block(ct::Tensor& A, ct::Tensor* HA, const int ncols)
{
    if (ncols <= 0)
    {
        return;
    }

    // gram = A^H A  (ncols x ncols)
    ct::Tensor gram(this->t_type, this->device_type, {ncols, ncols});
#ifdef __MPI
    this->pmmcn.set_dimension(BP_WORLD,
                              POOL_WORLD,
                              ncols,
                              this->n_basis,
                              ncols,
                              this->n_basis,
                              this->n_dim,
                              ncols);
#else
    this->pmmcn.set_dimension(ncols,
                              this->n_basis,
                              ncols,
                              this->n_basis,
                              this->n_dim,
                              ncols);
#endif
    this->pmmcn.multiply(static_cast<T>(1.0), A.data<T>(), A.data<T>(), static_cast<T>(0.0), gram.data<T>());

    // Cholesky: gram = U^H U (upper), then invert U in-place -> gram holds inv(U) in upper triangle
    int info = 0;
    lapackConnector::potrf('U', ncols, gram.data<T>(), ncols, info);
    assert(info == 0);
    lapackConnector::trtri('U', 'N', ncols, gram.data<T>(), ncols, info);
    assert(info == 0);

    // Zero out lower triangle so a dense GEMM applies only the upper-triangular factor.
    T* g = gram.data<T>();
    for (int j = 0; j < ncols; ++j)
    {
        for (int i = j + 1; i < ncols; ++i)
        {
            g[i + j * ncols] = static_cast<T>(0.0);
        }
    }

    // A <- A * inv(U)
    ModuleBase::gemm_op<T, Device>()('N',
                                    'N',
                                    this->n_dim,
                                    ncols,
                                    ncols,
                                    this->one,
                                    A.data<T>(),
                                    this->n_basis,
                                    gram.data<T>(),
                                    ncols,
                                    this->zero,
                                    this->work.data<T>(),
                                    this->n_basis);
    syncmem_complex_op()(A.data<T>(), this->work.data<T>(), static_cast<size_t>(ncols) * this->n_basis);

    if (HA)
    {
        ModuleBase::gemm_op<T, Device>()('N',
                                        'N',
                                        this->n_dim,
                                        ncols,
                                        ncols,
                                        this->one,
                                        HA->data<T>(),
                                        this->n_basis,
                                        gram.data<T>(),
                                        ncols,
                                        this->zero,
                                        this->work.data<T>(),
                                        this->n_basis);
        syncmem_complex_op()(HA->data<T>(), this->work.data<T>(), static_cast<size_t>(ncols) * this->n_basis);
    }
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::pack_basis(const int ncols, const bool has_p)
{
    // V columns: [X, W, P?]
    // Copy X
    syncmem_complex_op()(this->V.data<T>(), this->X.data<T>(), this->n_band_l * this->n_basis);
    // Copy W
    syncmem_complex_op()(this->V.data<T>() + this->n_band_l * this->n_basis,
                         this->W.data<T>(),
                         this->n_band_l * this->n_basis);

    if (has_p)
    {
        syncmem_complex_op()(this->V.data<T>() + 2 * this->n_band_l * this->n_basis,
                             this->P.data<T>(),
                             this->n_band_l * this->n_basis);
    }

    // HV: [HX, HW, HP?]
    syncmem_complex_op()(this->HV.data<T>(), this->HX.data<T>(), this->n_band_l * this->n_basis);

    syncmem_complex_op()(this->HV.data<T>() + this->n_band_l * this->n_basis,
                         this->HW.data<T>(),
                         this->n_band_l * this->n_basis);

    if (has_p)
    {
        syncmem_complex_op()(this->HV.data<T>() + 2 * this->n_band_l * this->n_basis,
                             this->HP.data<T>(),
                             this->n_band_l * this->n_basis);
    }

    (void)ncols;
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::compute_projected_mats(const int ncols)
{
    // hcc = V^H HV, scc = V^H V (col-major, ldh = ncols)
    const int ld_small = 3 * this->n_band;
#ifdef __MPI
    this->pmmcn.set_dimension(BP_WORLD,
                              POOL_WORLD,
                              ncols,
                              this->n_basis,
                              ncols,
                              this->n_basis,
                              this->n_dim,
                              ld_small);
#else
    this->pmmcn.set_dimension(ncols,
                              this->n_basis,
                              ncols,
                              this->n_basis,
                              this->n_dim,
                              ld_small);
#endif

    this->pmmcn.multiply(1.0, this->V.data<T>(), this->HV.data<T>(), 0.0, this->hcc.data<T>());
    this->pmmcn.multiply(1.0, this->V.data<T>(), this->V.data<T>(), 0.0, this->scc.data<T>());
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::solve_projected(const int ncols)
{
    // Solve (hcc) c = lambda (scc) c, eigenvectors in vcc
    const int ld_small = 3 * this->n_band;
    hsolver::hegvd_op<T, Device>()(this->ctx,
                                   ncols,
                                   ld_small,
                                   this->hcc.data<T>(),
                                   this->scc.data<T>(),
                                   this->eval.data<Real>(),
                                   this->vcc.data<T>());
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::update_from_projected(const int ncols, const bool has_p, const bool update_p)
{
    // Update X, HX from V, HV using the first n_band eigenvectors.
    // X_new = V * vcc(:, 1:nband)
    // HX_new = HV * vcc(:, 1:nband)
    const T* coeff = this->vcc.data<T>();
    const int ld_small = 3 * this->n_band;

    // work (n_basis x n_band_l)
    ModuleBase::gemm_op<T, Device>()('N',
                                    'N',
                                    this->n_dim,
                                    this->n_band_l,
                                    ncols,
                                    this->one,
                                    this->V.data<T>(),
                                    this->n_basis,
                                    coeff,
                                    ld_small,
                                    this->zero,
                                    this->work.data<T>(),
                                    this->n_basis);
    syncmem_complex_op()(this->X.data<T>(), this->work.data<T>(), this->n_band_l * this->n_basis);

    ModuleBase::gemm_op<T, Device>()('N',
                                    'N',
                                    this->n_dim,
                                    this->n_band_l,
                                    ncols,
                                    this->one,
                                    this->HV.data<T>(),
                                    this->n_basis,
                                    coeff,
                                    ld_small,
                                    this->zero,
                                    this->work.data<T>(),
                                    this->n_basis);
    syncmem_complex_op()(this->HX.data<T>(), this->work.data<T>(), this->n_band_l * this->n_basis);

    if (!update_p)
    {
        return;
    }

    // Update P (search directions) from blocks W and P (exclude X block to keep meaning)
    // P_new = W * Cw + P * Cp, where Cw = coeff(rows b..2b-1, cols 0..b-1)
    // and Cp = coeff(rows 2b..3b-1, cols 0..b-1)
    const int b = this->n_band_l;
    const int ld = ld_small;

    // When the subspace is smaller than the full [X,W,P] block (ncols < 3b),
    // only a prefix of W and/or P participates.  Keep the inner dimensions
    // consistent so we never read garbage rows from vcc.
    const int ncols_W = std::max(0, std::min(b, ncols - b));
    if (ncols_W > 0)
    {
        const T* Cw = coeff + b; // row offset b in vcc

        ModuleBase::gemm_op<T, Device>()('N',
                                        'N',
                                        this->n_dim,
                                        b,
                                        ncols_W,
                                        this->one,
                                        this->W.data<T>(),
                                        this->n_basis,
                                        Cw,
                                        ld,
                                        this->zero,
                                        this->P.data<T>(),
                                        this->n_basis);

        ModuleBase::gemm_op<T, Device>()('N',
                                        'N',
                                        this->n_dim,
                                        b,
                                        ncols_W,
                                        this->one,
                                        this->HW.data<T>(),
                                        this->n_basis,
                                        Cw,
                                        ld,
                                        this->zero,
                                        this->HP.data<T>(),
                                        this->n_basis);
    }

    if (has_p)
    {
        const int ncols_P = std::max(0, std::min(b, ncols - 2 * b));
        if (ncols_P > 0)
        {
            const T* Cp = coeff + 2 * b;
            // The P block inside V / HV always stores b columns (pack_basis
            // writes the full block).  Only the first ncols_P of them
            // correspond to valid rows in Cp.
            ModuleBase::gemm_op<T, Device>()('N',
                                            'N',
                                            this->n_dim,
                                            b,
                                            ncols_P,
                                            this->one,
                                            this->V.data<T>() + 2 * b * this->n_basis,
                                            this->n_basis,
                                            Cp,
                                            ld,
                                            this->one,
                                            this->P.data<T>(),
                                            this->n_basis);
            ModuleBase::gemm_op<T, Device>()('N',
                                            'N',
                                            this->n_dim,
                                            b,
                                            ncols_P,
                                            this->one,
                                            this->HV.data<T>() + 2 * b * this->n_basis,
                                            this->n_basis,
                                            Cp,
                                            ld,
                                            this->one,
                                            this->HP.data<T>(),
                                            this->n_basis);
        }
    }

        // Keep P orthogonal to X and keep HP consistent with P.
        // If we do: P <- P - X * (X^H P), then we must also do: HP <- HP - HX * (X^H P)
        // to preserve the relation HP = H * P inside the subspace.
        {
        const int bproj = this->n_band_l;
        const int ld_coef = bproj;
        ct::Tensor coef_xp(this->t_type, this->device_type, {ld_coef, bproj});

    #ifdef __MPI
        this->pmmcn.set_dimension(BP_WORLD,
                      POOL_WORLD,
                      bproj,
                      this->n_basis,
                      bproj,
                      this->n_basis,
                      this->n_dim,
                      ld_coef);
    #else
        this->pmmcn.set_dimension(bproj,
                      this->n_basis,
                      bproj,
                      this->n_basis,
                      this->n_dim,
                      ld_coef);
    #endif
        // coef_xp = X^H * P
        this->pmmcn.multiply(1.0, this->X.data<T>(), this->P.data<T>(), 0.0, coef_xp.data<T>());

        // P  -= X  * coef_xp
        ModuleBase::gemm_op<T, Device>()('N',
                        'N',
                        this->n_dim,
                        bproj,
                        bproj,
                        this->neg_one,
                        this->X.data<T>(),
                        this->n_basis,
                        coef_xp.data<T>(),
                        ld_coef,
                        this->one,
                        this->P.data<T>(),
                        this->n_basis);

        // HP -= HX * coef_xp
        ModuleBase::gemm_op<T, Device>()('N',
                        'N',
                        this->n_dim,
                        bproj,
                        bproj,
                        this->neg_one,
                        this->HX.data<T>(),
                        this->n_basis,
                        coef_xp.data<T>(),
                        ld_coef,
                        this->one,
                        this->HP.data<T>(),
                        this->n_basis);
        }

        // Block-orthonormalize P and apply the same transformation to HP.
        // (Avoid calling normalize_op(P) alone, which would desynchronize HP.)
        this->orthonormalize_block(this->P, &this->HP, this->n_band_l);
}

template <typename T, typename Device>
bool DiagoPPCG<T, Device>::check_convergence(const ct::Tensor& residual, const std::vector<double>& ethr_band)
{
    // Check ||r_i|| <= ethr_band[i] for all local bands.
    bool not_conv = false;
    for (int ib = 0; ib < this->n_band_l; ++ib)
    {
        const T* ri = residual.data<T>() + ib * this->n_basis;
        const Real nrm2 = std::sqrt(ModuleBase::dot_real_op<T, Device>()(this->n_dim, ri, ri, true));
        if (ib < static_cast<int>(ethr_band.size()) && nrm2 > static_cast<Real>(ethr_band[ib]))
        {
            not_conv = true;
            break;
        }
    }

#ifdef __MPI
    // Any rank not converged means global not converged.
    int local = not_conv ? 1 : 0;
    int global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, BP_WORLD);
    not_conv = (global != 0);
#endif

    return !not_conv;
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::compute_residual_and_precond(const std::vector<double>& ethr_band, bool& not_conv)
{
    // Residual R = HX - X * diag(eig)
    // First, R <- HX
    syncmem_complex_op()(this->R.data<T>(), this->HX.data<T>(), this->n_band_l * this->n_basis);

    // tmp = X * diag(e)
    apply_eigenvalues_op<T, Device>()(this->n_dim,
                                      this->n_basis,
                                      this->n_band_l,
                                      this->work.data<T>(),
                                      this->X.data<T>(),
                                      this->eval.data<Real>());
    // R -= tmp
    for (int ib = 0; ib < this->n_band_l; ++ib)
    {
        const T alpha = static_cast<T>(-1.0);
        ModuleBase::axpy_op<T, Device>()(this->n_dim,
                                         &alpha,
                                         this->work.data<T>() + ib * this->n_basis,
                                         1,
                                         this->R.data<T>() + ib * this->n_basis,
                                         1);
    }

    // not_conv if any band residual above threshold
    not_conv = !this->check_convergence(this->R, ethr_band);
    if (!not_conv)
    {
        return;
    }

    // W = - M^{-1} R
    syncmem_complex_op()(this->W.data<T>(), this->R.data<T>(), this->n_band_l * this->n_basis);
    precondition_op<T, Device>()(this->n_dim,
                                 this->W.data<T>(),
                                 0,
                                 this->n_band_l,
                                 this->prec.data<Real>(),
                                 this->eval.data<Real>());
    for (int ib = 0; ib < this->n_band_l; ++ib)
    {
        ModuleBase::vector_mul_real_op<T, Device>()(this->n_dim,
                                                    this->W.data<T>() + ib * this->n_basis,
                                                    this->W.data<T>() + ib * this->n_basis,
                                                    static_cast<Real>(-1.0));
    }

    // Project W out of X and P (if P contains previous directions)
    this->project_out(this->X, this->n_band_l, this->W, this->n_band_l);

    // Also project out of P to reduce near-dependencies
    // (P may be all-zero at the first iter, projection is harmless)
    this->project_out(this->P, this->n_band_l, this->W, this->n_band_l);

    // Normalize each vector in W
    normalize_op<T, Device>()(this->n_dim, this->W.data<T>(), 0, this->n_band_l, nullptr);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::diag(const HPsiFunc& hpsi_func,
                               T* psi_in,
                               Real* eigenvalue_in,
                               const std::vector<double>& ethr_band)
{
    // Map external psi to X
    this->X = ct::TensorMap((void*)psi_in, this->t_type, this->device_type, {this->n_band_l, this->n_basis});

    // Normalize initial X
    normalize_op<T, Device>()(this->n_dim, this->X.data<T>(), 0, this->n_band_l, nullptr);

    // HX = H X
    this->apply_h(hpsi_func, this->X, this->HX, this->n_band_l);

    // Make X block-orthonormal (and keep HX consistent) before any projection/RR.
    this->orthonormalize_block(this->X, &this->HX, this->n_band_l);

    // Initial Rayleigh-Ritz on X alone: solve (X^H H X) c = (X^H X) c Λ
    {
        const int ncols = this->n_band;
        ct::Tensor hxx(this->t_type, this->device_type, {ncols, ncols});
        ct::Tensor sxx(this->t_type, this->device_type, {ncols, ncols});
        ct::Tensor vxx(this->t_type, this->device_type, {ncols, ncols});
        ct::Tensor exx(this->r_type, this->device_type, {ncols});

    #ifdef __MPI
        this->pmmcn.set_dimension(BP_WORLD,
                      POOL_WORLD,
                      this->n_band_l,
                      this->n_basis,
                      this->n_band_l,
                      this->n_basis,
                      this->n_dim,
                      ncols);
    #else
        this->pmmcn.set_dimension(this->n_band_l,
                      this->n_basis,
                      this->n_band_l,
                      this->n_basis,
                      this->n_dim,
                      ncols);
    #endif
        this->pmmcn.multiply(1.0, this->X.data<T>(), this->HX.data<T>(), 0.0, hxx.data<T>());
        this->pmmcn.multiply(1.0, this->X.data<T>(), this->X.data<T>(), 0.0, sxx.data<T>());

        hsolver::hegvd_op<T, Device>()(this->ctx,
                                       ncols,
                                       ncols,
                                       hxx.data<T>(),
                                       sxx.data<T>(),
                                       exx.data<Real>(),
                                       vxx.data<T>());

        // Rotate X, HX: X <- X * vxx, HX <- HX * vxx
        ModuleBase::gemm_op<T, Device>()('N',
                                        'N',
                                        this->n_dim,
                                        this->n_band_l,
                                        ncols,
                                        this->one,
                                        this->X.data<T>(),
                                        this->n_basis,
                                        vxx.data<T>(),
                                        ncols,
                                        this->zero,
                                        this->work.data<T>(),
                                        this->n_basis);
        syncmem_complex_op()(this->X.data<T>(), this->work.data<T>(), this->n_band_l * this->n_basis);

        ModuleBase::gemm_op<T, Device>()('N',
                                        'N',
                                        this->n_dim,
                                        this->n_band_l,
                                        ncols,
                                        this->one,
                                        this->HX.data<T>(),
                                        this->n_basis,
                                        vxx.data<T>(),
                                        ncols,
                                        this->zero,
                                        this->work.data<T>(),
                                        this->n_basis);

        syncmem_complex_op()(this->HX.data<T>(), this->work.data<T>(), this->n_band_l * this->n_basis);

        syncmem_var_op()(this->eval.data<Real>(), exx.data<Real>(), this->n_band);
    }

    // Clear P/HP to zero for the first outer iteration
    Parallel_Reduce::ZEROS(this->P.data<T>(), this->n_band_l * this->n_basis);
    Parallel_Reduce::ZEROS(this->HP.data<T>(), this->n_band_l * this->n_basis);

    // Compute initial residual and preconditioned direction W.
    bool not_conv = true;
    this->compute_residual_and_precond(ethr_band, not_conv);

    // Determine how many inner iterations to allow.
    // When 3*n_band fits in the ambient space the P block is safe and
    // 2-3 iterations accelerate convergence.  Otherwise stick to 1 to
    // avoid near-singular overlap matrices.
    const bool p_safe = (3 * this->n_band <= this->n_dim - this->p_safe_margin_);
    const int max_iter = p_safe ? this->max_inner_iter_ : 1;
    const int max_w_cols = std::max(0, std::min(this->n_band_l, this->n_dim - this->n_band_l));
    int active_w_cols = not_conv ? max_w_cols : 0;
    if (not_conv && active_w_cols > 0)
    {
        // HW = H W
        this->apply_h(hpsi_func, this->W, this->HW, active_w_cols);

        // Keep W and HW consistent while improving conditioning.
        this->orthonormalize_block(this->W, &this->HW, active_w_cols);
    }
    for (int iter = 0; iter < max_iter && not_conv; ++iter)
    {
        const bool has_p = (iter > 0) && p_safe;
        const int raw_ncols = this->n_band + active_w_cols + (has_p ? this->n_band : 0);
        const int ncols = std::min(raw_ncols, this->n_dim);

        // Pack basis V/HV
        this->pack_basis(ncols, has_p);

        // Solve projected generalized eigenproblem
        this->compute_projected_mats(ncols);
        this->solve_projected(ncols);

        // Update X/HX and P/HP
        const bool update_p = (iter + 1 < max_iter);
        this->update_from_projected(ncols, has_p, update_p);

        if (iter + 1 >= max_iter)
        {
            break;
        }

        // Residual for next convergence check
        this->compute_residual_and_precond(ethr_band, not_conv);

        if (!not_conv)
        {
            break;
        }

        active_w_cols = max_w_cols;
        // Update HW for the next iteration
        this->apply_h(hpsi_func, this->W, this->HW, active_w_cols);

        // Keep W and HW consistent
        this->orthonormalize_block(this->W, &this->HW, active_w_cols);
    }

    // Final Rayleigh-Ritz on the current X subspace to ensure (X, eval) consistency.
    // This mirrors BPCG's exit behavior (subspace diagonalization before returning).
    {
        const int ncols = this->n_band;
        ct::Tensor hxx(this->t_type, this->device_type, {ncols, ncols});
        ct::Tensor sxx(this->t_type, this->device_type, {ncols, ncols});
        ct::Tensor vxx(this->t_type, this->device_type, {ncols, ncols});
        ct::Tensor exx(this->r_type, this->device_type, {ncols});

#ifdef __MPI
        this->pmmcn.set_dimension(BP_WORLD,
                                  POOL_WORLD,
                                  this->n_band_l,
                                  this->n_basis,
                                  this->n_band_l,
                                  this->n_basis,
                                  this->n_dim,
                                  ncols);
#else
        this->pmmcn.set_dimension(this->n_band_l,
                                  this->n_basis,
                                  this->n_band_l,
                                  this->n_basis,
                                  this->n_dim,
                                  ncols);
#endif
        this->pmmcn.multiply(1.0, this->X.data<T>(), this->HX.data<T>(), 0.0, hxx.data<T>());
        this->pmmcn.multiply(1.0, this->X.data<T>(), this->X.data<T>(), 0.0, sxx.data<T>());

        hsolver::hegvd_op<T, Device>()(this->ctx,
                                       ncols,
                                       ncols,
                                       hxx.data<T>(),
                                       sxx.data<T>(),
                                       exx.data<Real>(),
                                       vxx.data<T>());

        // Rotate X, HX: X <- X * vxx, HX <- HX * vxx
        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         this->n_dim,
                                         this->n_band_l,
                                         ncols,
                                         this->one,
                                         this->X.data<T>(),
                                         this->n_basis,
                                         vxx.data<T>(),
                                         ncols,
                                         this->zero,
                                         this->work.data<T>(),
                                         this->n_basis);
        syncmem_complex_op()(this->X.data<T>(), this->work.data<T>(), this->n_band_l * this->n_basis);

        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         this->n_dim,
                                         this->n_band_l,
                                         ncols,
                                         this->one,
                                         this->HX.data<T>(),
                                         this->n_basis,
                                         vxx.data<T>(),
                                         ncols,
                                         this->zero,
                                         this->work.data<T>(),
                                         this->n_basis);
        syncmem_complex_op()(this->HX.data<T>(), this->work.data<T>(), this->n_band_l * this->n_basis);

        syncmem_var_op()(this->eval.data<Real>(), exx.data<Real>(), this->n_band);
    }

    // Copy eigenvalues out
    syncmem_var_d2h_op()(eigenvalue_in, this->eval.data<Real>(), this->n_band);
}

// explicit instantiation
template class DiagoPPCG<std::complex<double>, base_device::DEVICE_CPU>;
template class DiagoPPCG<std::complex<float>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class DiagoPPCG<std::complex<double>, base_device::DEVICE_GPU>;
template class DiagoPPCG<std::complex<float>, base_device::DEVICE_GPU>;
#endif

} // namespace hsolver
