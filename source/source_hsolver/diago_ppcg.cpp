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

    this->one = &one_;
    this->zero = &zero_;
    this->neg_one = &neg_one_;
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
    this->beta_denom    = std::move(ct::Tensor(r_type, device_type, {this->n_band_l}));

    this->hsub          = std::move(ct::Tensor(t_type, device_type, {this->n_band, this->n_band}));

    this->hpsi          = std::move(ct::Tensor(t_type, device_type, {this->n_band_l, this->n_basis}));
    this->work          = std::move(ct::Tensor(t_type, device_type, {this->n_band_l, this->n_basis}));
    this->hgrad         = std::move(ct::Tensor(t_type, device_type, {this->n_band_l, this->n_basis}));
    this->grad_old      = std::move(ct::Tensor(t_type, device_type, {this->n_band_l, this->n_basis}));
    this->z_old         = std::move(ct::Tensor(t_type, device_type, {this->n_band_l, this->n_basis}));

    this->prec          = std::move(ct::Tensor(r_type, device_type, {this->n_basis}));

    this->grad          = std::move(ct::Tensor(t_type, device_type, {this->n_band_l, this->n_basis}));
#ifdef __MPI
    this->pmmcn.set_dimension(BP_WORLD, POOL_WORLD, n_band_l, n_basis, n_band_l, n_basis, n_dim, n_band);
    this->plintrans.set_dimension(n_dim, nband_l, n_band_l, n_basis, BP_WORLD, false);
#else
    this->pmmcn.set_dimension(n_band_l, n_basis, n_band_l, n_basis, n_dim, n_band);
    this->plintrans.set_dimension(n_dim, nband_l, n_band_l, n_basis, false);
#endif
}

template <typename T, typename Device>
bool DiagoPPCG<T, Device>::test_error(const ct::Tensor& err_in, const std::vector<double>& ethr_band)
{
    Real* _err_st = err_in.data<Real>();
    bool not_conv = false;
    std::vector<Real> tmp_cpu;
    if (err_in.device_type() == ct::DeviceType::GpuDevice) {
        tmp_cpu.resize(this->n_band_l);
        _err_st = tmp_cpu.data();
        syncmem_var_d2h_op()(_err_st, err_in.data<Real>(), this->n_band_l);
    }
    for (int ii = 0; ii < this->n_band_l; ii++) {
        if (_err_st[ii] > ethr_band[ii]) {
            not_conv = true;
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
    // PPCG-specific gradient computation using the Polak-Ribiere beta formula.
    //
    // For each band i:
    //   1. Normalize psi_i
    //   2. Compute eps_i = <psi_i|H|psi_i>
    //   3. Raw residual: r_i = H|psi_i> - eps_i * psi_i
    //   4. Preconditioned residual: z_i = -P^{-1} r_i
    //   5. PR beta: beta_i = max(0, <z_i, r_i - r_old_i> / <z_old_i, r_old_i>)
    //   6. Search direction: d_i = z_i + beta_i * d_old_i
    //
    // Key difference from BPCG: step 5 uses the Polak-Ribiere formula
    // (projection-corrected) instead of the Fletcher-Reeves formula.

    Real* prec = prec_in.data<Real>();
    Real* err_arr = err_out.data<Real>();
    Real* beta_arr = beta_out.data<Real>();
    T* psi = psi_in.data<T>();
    T* hpsi = hpsi_in.data<T>();
    T* grad_new = grad_out.data<T>();
    T* grad_old = grad_old_out.data<T>();
    T* z_old_arr = this->z_old.template data<T>();
    Real* beta_denom_arr = this->beta_denom.template data<Real>();

    for (int ib = 0; ib < this->n_band_l; ib++)
    {
        const int offset = ib * this->n_basis;
        T* p_psi   = psi + offset;
        T* p_hpsi  = hpsi + offset;
        T* p_d_new = grad_new + offset;
        T* p_d_old = grad_old + offset;
        T* p_z_old = z_old_arr + offset;

        // Step 1: normalize psi
        auto A = reinterpret_cast<const Real*>(p_psi);
        Real norm = BlasConnector::dot(2 * this->n_basis, A, 1, A, 1);
        Parallel_Reduce::reduce_pool(norm);
        norm = 1.0 / sqrt(norm);

        // Step 2: compute epsilo = <psi|H|psi>
        Real epsilo = 0.0;
        for (int i = 0; i < this->n_basis; i++)
        {
            p_psi[i] *= norm;
            p_hpsi[i] *= norm;
            epsilo += std::real(p_hpsi[i] * std::conj(p_psi[i]));
        }
        Parallel_Reduce::reduce_pool(epsilo);

        // Step 3 & 4: compute raw residual r and preconditioned z = -P^{-1} r
        // Simultaneously accumulate PR beta numerator:
        //   num = <z_new, r_new - r_old> = <z_new, r_new> - <z_new, r_old>
        Real beta_num_zr = 0.0;   // <z_new, r_new>
        Real beta_num_zo = 0.0;   // <z_new, r_old>
        Real err = 0.0;

        for (int i = 0; i < this->n_basis; i++)
        {
            T r_new = p_hpsi[i] - T(epsilo) * p_psi[i];      // unpreconditioned residual
            T z_new = T(-1.0) * r_new / T(prec[i]);           // preconditioned residual
            T r_old = T(prec[i]) * T(-1.0) * p_z_old[i];     // recover old raw residual from old preconditioned z

            beta_num_zr += std::real(z_new * std::conj(r_new));
            beta_num_zo += std::real(z_new * std::conj(r_old));
            err += std::norm(r_new);

            p_d_new[i] = z_new;
        }
        Parallel_Reduce::reduce_pool(beta_num_zr);
        Parallel_Reduce::reduce_pool(beta_num_zo);
        Parallel_Reduce::reduce_pool(err);

        // Step 5: Polak-Ribiere beta (clamped to [0, inf) for stability)
        Real beta_pr = 0.0;
        Real denom = beta_denom_arr[ib];
        if (std::abs(denom) > 1e-30)
        {
            beta_pr = (beta_num_zr - beta_num_zo) / denom;
            if (beta_pr < 0.0) { beta_pr = 0.0; }
        }

        // Step 6: mix with old search direction   d_new = z_new + beta_pr * d_old
        for (int i = 0; i < this->n_basis; i++)
        {
            p_d_new[i] += T(beta_pr) * p_d_old[i];
        }

        // Store PR beta and denominator for next iteration.
        // beta_denom_arr[ib] stores <z_cur, r_cur> which will be the
        // denominator <z_old, r_old> in the next iteration.
        beta_arr[ib] = beta_pr;
        beta_denom_arr[ib] = beta_num_zr + 1e-30;
        err_arr[ib] = sqrt(err);

        // Save preconditioned z (before mixing) for next iteration's PR beta.
        // We stored z_new in p_d_new before mixing; now copy back to z_old.
        // Recover z_new from: z_new = d_new - beta_pr * d_old
        // Since d_new = p_d_new and d_old = p_d_old, and beta_pr may be zero:
        for (int i = 0; i < this->n_basis; i++)
        {
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
    return;
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::rotate_wf(
        const ct::Tensor& hsub_in,
        ct::Tensor& psi_out,
        ct::Tensor& workspace_in)
{
    this->plintrans.act(1.0, psi_out.data<T>(), hsub_in.data<T>(), 0.0, workspace_in.data<T>());
    syncmem_complex_op()(psi_out.template data<T>(), workspace_in.template data<T>(), this->n_band_l * this->n_basis);
    return;
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
    ct::kernels::lapack_heevd<T, ct_Device>()(this->n_band, hsub_out.data<T>(), this->n_band, eigenvalue_out.data<Real>());
    return;
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
    return;
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
    return;
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

    // Initial subspace diagonalization to get a good starting guess.
    this->calc_hsub_with_block(hpsi_func, psi_in, this->psi, this->hpsi,
                               this->hsub, this->work, this->eigen);

    // Initialize PPCG state: first search direction is the steepest descent.
    // grad_old and z_old start as zero; beta_denom starts as infinity so
    // the first beta is computed as 0 (pure steepest descent).
    setmem_complex_op()(this->grad_old.template data<T>(), 0, this->n_basis * this->n_band_l);
    setmem_complex_op()(this->z_old.template data<T>(),    0, this->n_basis * this->n_band_l);
    setmem_var_op()(this->beta_denom.template data<Real>(), std::numeric_limits<Real>::infinity(), this->n_band_l);

    int ntry = 0;
    int max_iter = current_scf_iter > 1 ?
                   this->nline :
                   this->nline * 6;
    do
    {
        ++ntry;

        // PPCG step 1: compute the preconditioned conjugate direction.
        //   Uses the Polak-Ribiere beta formula (projection-corrected):
        //     beta_i = max(0, <z_i, r_i - r_old_i> / <z_old_i, r_old_i>)
        //     d_i = z_i + beta_i * d_old_i
        //   where z_i = -P^{-1} r_i is the preconditioned residual.
        this->calc_grad_with_block(this->prec, this->err_st, this->beta,
                                   this->psi, this->hpsi, this->grad, this->grad_old);

        // PPCG step 2: explicit projection.
        //   Project the search direction d to be orthogonal to all current
        //   approximate eigenvectors: d_i = d_i - sum_j <psi_j|d_i> psi_j.
        //   This keeps the search within the relevant subspace.
        this->orth_projection(this->psi, this->hsub, this->grad);

        // PPCG step 3: save the current search direction for the next iteration.
        //   z_old is already saved inside calc_grad_with_block.
        syncmem_complex_op()(this->grad_old.template data<T>(),
                             this->grad.template data<T>(),
                             n_basis * n_band_l);

        // PPCG step 4: apply H to the search direction.
        this->calc_hpsi_with_block(hpsi_func, this->grad.template data<T>(), this->hgrad);

        // PPCG step 5: line minimization along the search direction.
        //   psi_new = psi * cos(theta) + d * sin(theta)
        //   where theta minimizes the Rayleigh quotient.
        this->line_minimize(this->grad, this->hgrad, this->psi, this->hpsi);

        // PPCG step 6: Cholesky orthonormalization.
        //   Ensures all bands remain orthonormal after the line minimization.
        this->orth_cholesky(this->work, this->psi, this->hpsi, this->hsub);

        // Periodic subspace re-alignment in the first SCF iteration.
        if (current_scf_iter == 1 && ntry % this->nline == 0) {
            this->calc_hsub_with_block(hpsi_func, psi_in, this->psi, this->hpsi,
                                       this->hsub, this->work, this->eigen);
        }
    } while (ntry < max_iter && this->test_error(this->err_st, ethr_band));

    // Final subspace diagonalization and rotation to get accurate eigenvalues.
    this->calc_hsub_with_block_exit(this->psi, this->hpsi, this->hsub, this->work, this->eigen);

    int start_nband = 0;
#ifdef __MPI
    if (this->plintrans.nproc_col > 1)
    {
        start_nband = this->plintrans.start_colB[GlobalV::MY_BNDGROUP];
    }
#endif
    syncmem_var_d2h_op()(eigenvalue_in, this->eigen.template data<Real>() + start_nband, this->n_band_l);

    return;
}

template class DiagoPPCG<std::complex<float>, base_device::DEVICE_CPU>;
template class DiagoPPCG<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class DiagoPPCG<std::complex<float>, base_device::DEVICE_GPU>;
template class DiagoPPCG<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace hsolver
