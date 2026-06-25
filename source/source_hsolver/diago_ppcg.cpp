#include "source_hsolver/diago_ppcg.h"

#include "diago_iter_assist.h"
#include "source_base/global_function.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_reduce.h"
#include "source_hsolver/kernels/bpcg_kernel_op.h"
#include "para_linear_transform.h"

#include <ATen/kernels/blas.h>
#include <ATen/kernels/lapack.h>
#include <ATen/ops/einsum_op.h>
#include <algorithm>
#include <cmath>
#include <limits>

namespace hsolver {

template <typename T, typename Device>
DiagoPPCG<T, Device>::DiagoPPCG(const Real* precondition_in)
{
    this->r_type   = ct::DataTypeToEnum<Real>::value;
    this->t_type   = ct::DataTypeToEnum<T>::value;
    this->device_type    = ct::DeviceTypeToEnum<Device>::value;

    this->h_prec  = std::move(ct::TensorMap((void *) precondition_in, r_type, device_type, {this->n_basis}));
}

template <typename T, typename Device>
DiagoPPCG<T, Device>::~DiagoPPCG() {
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::init_iter(const int nband, const int nband_l, const int nbasis, const int ndim) {
    this->n_band        = nband;
    this->n_band_l      = nband_l;
    this->n_basis       = nbasis;
    this->n_dim         = ndim;

    this->beta          = std::move(ct::Tensor(r_type, device_type, {this->n_band_l}));
    this->eigen         = std::move(ct::Tensor(r_type, device_type, {this->n_band}));
    this->err_st        = std::move(ct::Tensor(r_type, device_type, {this->n_band_l}));

    this->hsub          = std::move(ct::Tensor(t_type, device_type, {this->n_band, this->n_band}));

    this->hpsi          = std::move(ct::Tensor(t_type, device_type, {this->n_band_l, this->n_basis}));
    this->work          = std::move(ct::Tensor(t_type, device_type, {this->n_band_l, this->n_basis}));
    this->hgrad         = std::move(ct::Tensor(t_type, device_type, {this->n_band_l, this->n_basis}));
    this->grad_old      = std::move(ct::Tensor(t_type, device_type, {this->n_band_l, this->n_basis}));
    this->z_old         = std::move(ct::Tensor(t_type, device_type, {this->n_band_l, this->n_basis}));

    this->prec          = std::move(ct::Tensor(r_type, device_type, {this->n_basis}));

    this->grad          = std::move(ct::Tensor(t_type, device_type, {this->n_band_l, this->n_basis}));

    this->conv_mask.assign(this->n_band_l, 0);

#ifdef __MPI
    this->pmmcn.set_dimension(BP_WORLD, POOL_WORLD, n_band_l, n_basis, n_band_l, n_basis, n_dim, n_band);
    this->plintrans.set_dimension(n_dim, nband_l, n_band_l, n_basis, BP_WORLD, false);
#else
    this->pmmcn.set_dimension(n_band_l, n_basis, n_band_l, n_basis, n_dim, n_band);
    this->plintrans.set_dimension(n_dim, nband_l, n_band_l, n_basis, false);
#endif
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::update_convergence(const std::vector<double>& ethr_band)
{
    Real* err = this->err_st.data<Real>();
    for (int ib = 0; ib < this->n_band_l; ib++) {
        if (err[ib] <= static_cast<Real>(ethr_band[ib])) {
            this->conv_mask[ib] = 1;
        }
    }
}

template <typename T, typename Device>
bool DiagoPPCG<T, Device>::test_error(const ct::Tensor& err_in, const std::vector<double>& ethr_band)
{
    Real* _err_st = err_in.data<Real>();
    bool not_conv = false;
    for (int ii = 0; ii < this->n_band_l; ii++) {
        if (_err_st[ii] > ethr_band[ii]) {
            not_conv = true;
            break;
        }
    }
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &not_conv, 1, MPI_C_BOOL, MPI_LOR, BP_WORLD);
#endif
    return not_conv;
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::line_minimize(
    ct::Tensor& grad_in,
    ct::Tensor& hgrad_in,
    ct::Tensor& psi_out,
    ct::Tensor& hpsi_out)
{
    line_minimize_with_block_op<T, Device>()(grad_in.data<T>(),
                                             hgrad_in.data<T>(),
                                             psi_out.data<T>(),
                                             hpsi_out.data<T>(),
                                             this->n_dim,
                                             this->n_basis,
                                             this->n_band_l);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::orth_cholesky(
        ct::Tensor& workspace_in,
        ct::Tensor& psi_out,
        ct::Tensor& hpsi_out,
        ct::Tensor& hsub_out)
{
    this->pmmcn.multiply(1.0, psi_out.data<T>(), psi_out.data<T>(), 0.0, hsub_out.data<T>());

    ct::kernels::set_matrix<T, ct_Device>()(
        'L', hsub_out.data<T>(), this->n_band);

    ct::kernels::lapack_potrf<T, ct_Device>()(
        'U', this->n_band, hsub_out.data<T>(), this->n_band);
    ct::kernels::lapack_trtri<T, ct_Device>()(
        'U', 'N', this->n_band, hsub_out.data<T>(), this->n_band);

    this->rotate_wf(hsub_out, psi_out, workspace_in);
    this->rotate_wf(hsub_out, hpsi_out, workspace_in);

    // Rotate conjugate history so that z_old and grad_old remain
    // consistent with the rotated basis.  Without this, the PR beta
    // formula would mix vectors from different bases.
    this->rotate_wf(hsub_out, this->z_old, workspace_in);
    this->rotate_wf(hsub_out, this->grad_old, workspace_in);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::calc_grad_with_block(
        const ct::Tensor& prec_in,
        ct::Tensor& err_out,
        ct::Tensor& beta_out,
        ct::Tensor& psi_in,
        ct::Tensor& hpsi_in,
        ct::Tensor& grad_out,
        ct::Tensor& grad_old_out)
{
    // PPCG gradient computation with Polak-Ribiere beta and soft-locking.
    //
    // Uses batched MPI reductions: scalars are collected into arrays and
    // reduced in a single call per quantity, instead of per-band calls.
    //
    // For each active band i:
    //   1. Normalize psi_i
    //   2. Compute eps_i = <psi_i|H|psi_i>
    //   3. Raw residual: r_i = H|psi_i> - eps_i * psi_i
    //   4. Preconditioned residual: z_i = -P^{-1} r_i
    //   5. PR beta: beta_i = max(0, <z_i, r_i - r_old_i> / <z_old_i, r_old_i>)
    //   6. Search direction: d_i = z_i + beta_i * d_old_i
    //
    // Converged bands (conv_mask[ib] == 1) are soft-locked:
    // their search direction is set to zero and error is set to 0.

    Real* prec = prec_in.data<Real>();
    Real* err_arr = err_out.data<Real>();
    Real* beta_arr = beta_out.data<Real>();
    T* psi = psi_in.data<T>();
    T* hpsi = hpsi_in.data<T>();
    T* grad_new = grad_out.data<T>();
    T* grad_old = grad_old_out.data<T>();
    T* z_old_p = this->z_old.template data<T>();

    const int nb = this->n_band_l;
    const int ns = this->n_basis;

    // Temporary arrays for batching MPI reductions
    std::vector<Real> norms(nb, 0.0);
    std::vector<Real> epsilos(nb, 0.0);
    std::vector<Real> beta_num_zr(nb, 0.0);
    std::vector<Real> beta_num_zo(nb, 0.0);
    std::vector<Real> err_loc(nb, 0.0);
    std::vector<Real> zPz_cur(nb, 0.0);  // <z_old, P * z_old> from current (possibly rotated) z_old

    // ---- Phase 1: compute psi norms for active bands ----
    for (int ib = 0; ib < nb; ib++) {
        if (this->conv_mask[ib]) continue;
        T* p_psi = psi + ib * ns;
        auto A = reinterpret_cast<const Real*>(p_psi);
        norms[ib] = BlasConnector::dot(2 * ns, A, 1, A, 1);
    }
    Parallel_Reduce::reduce_pool(norms.data(), nb);

    // ---- Phase 2: normalize psi, hpsi and compute epsilos ----
    for (int ib = 0; ib < nb; ib++) {
        if (this->conv_mask[ib]) continue;
        T* p_psi  = psi + ib * ns;
        T* p_hpsi = hpsi + ib * ns;
        Real nrm = 1.0 / std::sqrt(norms[ib]);
        Real eps = 0.0;
        for (int i = 0; i < ns; i++) {
            p_psi[i]  *= nrm;
            p_hpsi[i] *= nrm;
            eps += std::real(p_hpsi[i] * std::conj(p_psi[i]));
        }
        epsilos[ib] = eps;
    }
    Parallel_Reduce::reduce_pool(epsilos.data(), nb);

    // ---- Phase 3: residuals, preconditioned z, PR beta numerators, error ----
    for (int ib = 0; ib < nb; ib++) {
        if (this->conv_mask[ib]) continue;
        T* p_psi   = psi   + ib * ns;
        T* p_hpsi  = hpsi  + ib * ns;
        T* p_d_new = grad_new + ib * ns;
        T* p_z_old = z_old_p  + ib * ns;
        Real eps   = epsilos[ib];

        Real bzr = 0.0, bzo = 0.0, eloc = 0.0, zpz = 0.0;
        for (int i = 0; i < ns; i++) {
            T r_new = p_hpsi[i] - T(eps) * p_psi[i];
            T z_new = T(-1.0) * r_new / T(prec[i]);
            T r_old = T(prec[i]) * T(-1.0) * p_z_old[i];

            bzr += std::real(z_new * std::conj(r_new));
            bzo += std::real(z_new * std::conj(r_old));
            eloc += std::norm(r_new);
            zpz  += prec[i] * std::norm(p_z_old[i]);

            p_d_new[i] = z_new;
        }
        beta_num_zr[ib] = bzr;
        beta_num_zo[ib] = bzo;
        err_loc[ib] = eloc;
        zPz_cur[ib] = zpz;
    }
    Parallel_Reduce::reduce_pool(beta_num_zr.data(), nb);
    Parallel_Reduce::reduce_pool(beta_num_zo.data(), nb);
    Parallel_Reduce::reduce_pool(err_loc.data(), nb);
    Parallel_Reduce::reduce_pool(zPz_cur.data(), nb);

    // ---- Phase 4: compute beta, mix directions, save state ----
    for (int ib = 0; ib < nb; ib++) {
        if (this->conv_mask[ib]) {
            // Soft-locked: zero search direction, no error contribution.
            T* p_d_new = grad_new + ib * ns;
            std::fill(p_d_new, p_d_new + ns, T(0));
            err_arr[ib]  = 0;
            beta_arr[ib] = 0;
            continue;
        }

        T* p_d_new = grad_new + ib * ns;
        T* p_d_old = grad_old + ib * ns;
        T* p_z_old = z_old_p  + ib * ns;

        // PR denominator: <z_old, r_old> = -<z_old, P*z_old> = -zPz_cur.
        // Recomputing from current z_old ensures consistency after
        // Cholesky rotation of z_old in orth_cholesky.
        Real denom = -zPz_cur[ib];
        Real beta_pr = 0.0;
        if (std::abs(denom) > 1e-30) {
            beta_pr = (beta_num_zr[ib] - beta_num_zo[ib]) / denom;
            if (beta_pr < 0.0) { beta_pr = 0.0; }
        }

        beta_arr[ib] = beta_pr;
        err_arr[ib] = std::sqrt(err_loc[ib]);

        // Mix: d_new = z_new + beta_pr * d_old
        for (int i = 0; i < ns; i++) {
            p_d_new[i] += T(beta_pr) * p_d_old[i];
        }
        // Save z_new (before mixing) for next iteration's PR beta.
        // Recover: z_new = d_new - beta_pr * d_old
        for (int i = 0; i < ns; i++) {
            p_z_old[i] = p_d_new[i] - T(beta_pr) * p_d_old[i];
        }
    }
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::calc_prec()
{
    syncmem_var_h2d_op()(this->prec.template data<Real>(), this->h_prec.template data<Real>(), this->n_basis);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::orth_projection(
        const ct::Tensor& psi_in,
        ct::Tensor& hsub_in,
        ct::Tensor& grad_out)
{
    this->pmmcn.multiply(1.0, psi_in.data<T>(), grad_out.data<T>(), 0.0, hsub_in.data<T>());
    this->plintrans.act(-1.0, psi_in.data<T>(), hsub_in.data<T>(), 1.0, grad_out.data<T>());
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::rotate_wf(
        const ct::Tensor& hsub_in,
        ct::Tensor& psi_out,
        ct::Tensor& workspace_in)
{
    this->plintrans.act(1.0, psi_out.data<T>(), hsub_in.data<T>(), 0.0, workspace_in.data<T>());
    syncmem_complex_op()(psi_out.template data<T>(), workspace_in.template data<T>(), this->n_band_l * this->n_basis);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::calc_hpsi_with_block(
        const HPsiFunc& hpsi_func,
        T *psi_in,
        ct::Tensor& hpsi_out)
{
    hpsi_func(psi_in, hpsi_out.data<T>(), this->n_basis, this->n_band_l);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::diag_hsub(
        const ct::Tensor& psi_in,
        const ct::Tensor& hpsi_in,
        ct::Tensor& hsub_out,
        ct::Tensor& eigenvalue_out)
{
    this->pmmcn.multiply(1.0, hpsi_in.data<T>(), psi_in.data<T>(), 0.0, hsub_out.data<T>());

    const int n = this->n_band;
    int lwork = 2 * n + n * n;
    if (lwork < 1 + 6 * n + 2 * n * n) {
        lwork = 1 + 6 * n + 2 * n * n;
    }
    ct::Tensor work(t_type, device_type, {lwork});
    work.zero();

    int lrwork = 1 + 5 * n + 2 * n * n;
    ct::Tensor rwork(r_type, device_type, {lrwork});
    rwork.zero();

    int liwork = 3 + 5 * n;
    ct::Tensor iwork(ct::DataType::DT_INT, device_type, {liwork});
    iwork.zero();

    int info = 0;
    container::lapackConnector::dnevd('V', 'U', n,
                                       hsub_out.data<T>(), n,
                                       eigenvalue_out.data<Real>(),
                                       work.data<T>(), lwork,
                                       rwork.data<Real>(), lrwork,
                                       iwork.data<int>(), liwork,
                                       info);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::calc_hsub_with_block(
        const HPsiFunc& hpsi_func,
        T *psi_in,
        ct::Tensor& psi_out,
        ct::Tensor& hpsi_out,
        ct::Tensor& hsub_out,
        ct::Tensor& workspace_in,
        ct::Tensor& eigenvalue_out)
{
    this->calc_hpsi_with_block(hpsi_func, psi_in, hpsi_out);
    this->diag_hsub(psi_out, hpsi_out, hsub_out, eigenvalue_out);
    this->rotate_wf(hsub_out, psi_out, workspace_in);
    this->rotate_wf(hsub_out, hpsi_out, workspace_in);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::calc_hsub_with_block_exit(
        ct::Tensor& psi_out,
        ct::Tensor& hpsi_out,
        ct::Tensor& hsub_out,
        ct::Tensor& workspace_in,
        ct::Tensor& eigenvalue_out)
{
    this->diag_hsub(psi_out, hpsi_out, hsub_out, eigenvalue_out);
    this->rotate_wf(hsub_out, psi_out, workspace_in);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::diag(const HPsiFunc& hpsi_func,
                                T* psi_in,
                                Real* eigenvalue_in,
                                const std::vector<double>& ethr_band)
{
    const int current_scf_iter = hsolver::DiagoIterAssist<T, Device>::SCF_ITER;

    this->psi = std::move(ct::TensorMap(psi_in, t_type, device_type, {this->n_band_l, this->n_basis}));

    this->calc_prec();

    // Reset convergence mask for this diag call.
    this->conv_mask.assign(this->n_band_l, 0);

    // Initial subspace diagonalization.
    this->calc_hsub_with_block(hpsi_func, psi_in, this->psi, this->hpsi,
                               this->hsub, this->work, this->eigen);

    setmem_complex_op()(this->grad_old.template data<T>(), 0, this->n_basis * this->n_band_l);
    setmem_complex_op()(this->z_old.template data<T>(),    0, this->n_basis * this->n_band_l);

    int ntry = 0;
    int max_iter = current_scf_iter > 1 ?
                   this->nline :
                   this->nline * 6;
    do
    {
        ++ntry;

        // Soft-lock bands that converged in the previous iteration.
        if (ntry > 1) {
            this->update_convergence(ethr_band);
        }

        this->calc_grad_with_block(this->prec, this->err_st, this->beta,
                                   this->psi, this->hpsi, this->grad, this->grad_old);

        this->orth_projection(this->psi, this->hsub, this->grad);

        syncmem_complex_op()(this->grad_old.template data<T>(),
                             this->grad.template data<T>(),
                             n_basis * n_band_l);

        this->calc_hpsi_with_block(hpsi_func, this->grad.template data<T>(), this->hgrad);

        this->line_minimize(this->grad, this->hgrad, this->psi, this->hpsi);

        this->orth_cholesky(this->work, this->psi, this->hpsi, this->hsub);

        if (current_scf_iter == 1 && ntry % this->nline == 0) {
            this->calc_hsub_with_block(hpsi_func, psi_in, this->psi, this->hpsi,
                                       this->hsub, this->work, this->eigen);
            // Reset convergence mask after subspace re-diagonalization
            // because psi vectors have been rotated.
            this->conv_mask.assign(this->n_band_l, 0);
        }
    } while (ntry < max_iter && this->test_error(this->err_st, ethr_band));

    this->calc_hsub_with_block_exit(this->psi, this->hpsi, this->hsub, this->work, this->eigen);

    int start_nband = 0;
#ifdef __MPI
    if (this->plintrans.nproc_col > 1)
    {
        start_nband = this->plintrans.start_colB[GlobalV::MY_BNDGROUP];
    }
#endif
    syncmem_var_d2h_op()(eigenvalue_in, this->eigen.template data<Real>() + start_nband, this->n_band_l);
}

template class DiagoPPCG<std::complex<float>, base_device::DEVICE_CPU>;
template class DiagoPPCG<std::complex<double>, base_device::DEVICE_CPU>;

} // namespace hsolver
