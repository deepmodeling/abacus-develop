#include "source_hsolver/diago_bpcg.h"

#include "diago_iter_assist.h"
#include "source_base/global_function.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_comm.h"
#include "source_hsolver/kernels/bpcg_kernel_op.h"
#include "para_linear_transform.h"

#include <ATen/kernels/blas.h>
#include <ATen/kernels/lapack.h>
#include <ATen/ops/einsum_op.h>
#include <limits>
#include <algorithm>   // for std::swap

namespace hsolver {

template<typename T, typename Device>
DiagoBPCG<T, Device>::DiagoBPCG(const Real* precondition_in)
{
    this->r_type   = ct::DataTypeToEnum<Real>::value;
    this->t_type   = ct::DataTypeToEnum<T>::value;
    this->device_type    = ct::DeviceTypeToEnum<Device>::value;

    this->h_prec  = std::move(ct::TensorMap((void *) precondition_in, r_type, device_type, {this->n_basis}));

    this->one = &one_;
    this->zero = &zero_;
    this->neg_one = &neg_one_;
}

template<typename T, typename Device>
DiagoBPCG<T, Device>::~DiagoBPCG() {
}

template<typename T, typename Device>
void DiagoBPCG<T, Device>::init_iter(const int nband, const int nband_l, const int nbasis, const int ndim) {
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

template<typename T, typename Device>
bool DiagoBPCG<T, Device>::test_error(const ct::Tensor& err_in, const std::vector<double>& ethr_band)
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

template<typename T, typename Device>
void DiagoBPCG<T, Device>::line_minimize(
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

template<typename T, typename Device>
void DiagoBPCG<T, Device>::orth_cholesky(
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

template<typename T, typename Device>
void DiagoBPCG<T, Device>::calc_grad_with_block(
        const ct::Tensor& prec_in,
        ct::Tensor& err_out,
        ct::Tensor& beta_out,
        ct::Tensor& psi_in,
        ct::Tensor& hpsi_in,
        ct::Tensor& grad_out,
        ct::Tensor& grad_old_out)
{
    calc_grad_with_block_op<T, Device>()(prec_in.data<Real>(),
                                         err_out.data<Real>(),
                                         beta_out.data<Real>(),
                                         psi_in.data<T>(),
                                         hpsi_in.data<T>(),
                                         grad_out.data<T>(),
                                         grad_old_out.data<T>(),
                                         this->n_dim,
                                         this->n_basis,
                                         this->n_band_l);
}

template<typename T, typename Device>
void DiagoBPCG<T, Device>::calc_prec()
{
    syncmem_var_h2d_op()(this->prec.template data<Real>(), this->h_prec.template data<Real>(), this->n_basis);
}

template<typename T, typename Device>
void DiagoBPCG<T, Device>::orth_projection(
        const ct::Tensor& psi_in,
        ct::Tensor& hsub_in,
        ct::Tensor& grad_out)
{
    this->pmmcn.multiply(1.0, psi_in.data<T>(), grad_out.data<T>(), 0.0, hsub_in.data<T>());
    this->plintrans.act(-1.0, psi_in.data<T>(), hsub_in.data<T>(), 1.0, grad_out.data<T>());
    return;
}

template<typename T, typename Device>
void DiagoBPCG<T, Device>::rotate_wf(
        const ct::Tensor& hsub_in,
        ct::Tensor& psi_out,
        ct::Tensor& workspace_in)
{
    this->plintrans.act(1.0, psi_out.data<T>(), hsub_in.data<T>(), 0.0, workspace_in.data<T>());
    syncmem_complex_op()(psi_out.template data<T>(), workspace_in.template data<T>(), this->n_band_l * this->n_basis);
    return;
}

template<typename T, typename Device>
void DiagoBPCG<T, Device>::calc_hpsi_with_block(
        const HPsiFunc& hpsi_func,
        T *psi_in,
        ct::Tensor& hpsi_out)
{
    hpsi_func(psi_in, hpsi_out.data<T>(), this->n_basis, this->n_band_l);
}

template<typename T, typename Device>
void DiagoBPCG<T, Device>::diag_hsub(
        const ct::Tensor& psi_in,
        const ct::Tensor& hpsi_in,
        ct::Tensor& hsub_out,
        ct::Tensor& eigenvalue_out)
{
    this->pmmcn.multiply(1.0, hpsi_in.data<T>(), psi_in.data<T>(), 0.0, hsub_out.data<T>());
    ct::kernels::lapack_heevd<T, ct_Device>()(this->n_band, hsub_out.data<T>(), this->n_band, eigenvalue_out.data<Real>());
    return;
}

template<typename T, typename Device>
void DiagoBPCG<T, Device>::calc_hsub_with_block(
        const HPsiFunc& hpsi_func,
        T *psi_in,
        ct::Tensor& psi_out,
        ct::Tensor& hpsi_out,
        ct::Tensor& hsub_out,
        ct::Tensor& workspace_in,
        ct::Tensor& eigenvalue_out)
{
    this->calc_hpsi_with_block(hpsi_func, psi_in, hpsi_out);
    this->diag_hsub(psi_out,hpsi_out, hsub_out, eigenvalue_out);
    this->rotate_wf(hsub_out, psi_out, workspace_in);
    this->rotate_wf(hsub_out, hpsi_out, workspace_in);
    return;
}

template<typename T, typename Device>
void DiagoBPCG<T, Device>::calc_hsub_with_block_exit(
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

// ===== 核心修改：grad_old 通过 swap 而非拷贝来更新 =====
template <typename T, typename Device>
void DiagoBPCG<T, Device>::diag(const HPsiFunc& hpsi_func,
                                T* psi_in,
                                Real* eigenvalue_in,
                                const std::vector<double>& ethr_band)
{
    const int current_scf_iter = hsolver::DiagoIterAssist<T, Device>::SCF_ITER;
    this->psi = std::move(ct::TensorMap(psi_in, t_type, device_type, {this->n_band_l, this->n_basis}));

    this->calc_prec();

    this->calc_hsub_with_block(hpsi_func, psi_in, this->psi, this->hpsi, this->hsub, this->work, this->eigen);

    setmem_complex_op()(this->grad_old.template data<T>(), 0, this->n_basis * this->n_band_l);
    setmem_var_op()(this->beta.template data<Real>(), std::numeric_limits<Real>::infinity(), this->n_band_l);

    int ntry = 0;
    int max_iter = current_scf_iter > 1 ? this->nline : this->nline * 6;

    // 标志位：第一次迭代需要显式拷贝，之后通过 swap 实现
    bool first_iter = true;

    do
    {
        // 除第一次外，将上一轮的 grad 通过交换传给 grad_old，避免深拷贝
        if (!first_iter)
        {
            std::swap(this->grad, this->grad_old);
        }

        ++ntry;
        this->calc_grad_with_block(this->prec, this->err_st, this->beta,
                                 this->psi, this->hpsi, this->grad, this->grad_old);

        this->orth_projection(this->psi, this->hsub, this->grad);

        // 第一次迭代后需要把当前 grad 赋值给 grad_old，供下一轮使用
        if (first_iter)
        {
            syncmem_complex_op()(this->grad_old.template data<T>(), this->grad.template data<T>(), n_basis * n_band_l);
            first_iter = false;
        }

        this->calc_hpsi_with_block(hpsi_func, this->grad.template data<T>(), this->hgrad);

        this->line_minimize(this->grad, this->hgrad, this->psi, this->hpsi);

        this->orth_cholesky(this->work, this->psi, this->hpsi, this->hsub);

        if (current_scf_iter == 1 && ntry % this->nline == 0) {
            this->calc_hsub_with_block(hpsi_func, psi_in, this->psi, this->hpsi, this->hsub, this->work, this->eigen);
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

    return;
}

// 显式模板实例化
template class DiagoBPCG<std::complex<float>, base_device::DEVICE_CPU>;
template class DiagoBPCG<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class DiagoBPCG<std::complex<float>, base_device::DEVICE_GPU>;
template class DiagoBPCG<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace hsolver