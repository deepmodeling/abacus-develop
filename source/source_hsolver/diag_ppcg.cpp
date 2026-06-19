#include "diago_iter_assist.h"
#include "para_linear_transform.h"
#include "source_base/global_function.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_comm.h"
#include "source_hsolver/diago_ppcg.h"

#include <ATen/kernels/blas.h>
#include <ATen/kernels/lapack.h>
#include <cstring>
#include <limits>

namespace hsolver
{

template <typename T, typename Device>
DiagoPPCG<T, Device>::DiagoPPCG(const Real* precondition_in)
{
    this->r_type = ct::DataTypeToEnum<Real>::value;
    this->t_type = ct::DataTypeToEnum<T>::value;
    this->device_type = ct::DeviceTypeToEnum<Device>::value;

    this->h_prec = std::move(ct::TensorMap((void*)precondition_in, r_type, ct::DeviceType::CpuDevice, {this->n_basis}));

    this->one = &one_;
    this->zero = &zero_;
    this->neg_one = &neg_one_;
}

template <typename T, typename Device>
DiagoPPCG<T, Device>::~DiagoPPCG()
{
    // h_prec is a ref to outside data, do not free.
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::init_iter(const int nband, const int nband_l, const int nbasis, const int ndim)
{
    this->n_band = nband;
    this->n_band_l = nband_l;
    this->n_basis = nbasis;
    this->n_dim = ndim;

    this->eigen = std::move(ct::Tensor(r_type, device_type, {this->n_band}));
    this->err_st = std::move(ct::Tensor(r_type, device_type, {this->n_band_l}));

    this->psi = std::move(ct::Tensor(t_type, device_type, {this->n_band_l, this->n_basis}));
    this->hpsi = std::move(ct::Tensor(t_type, device_type, {this->n_band_l, this->n_basis}));
    this->w = std::move(ct::Tensor(t_type, device_type, {this->n_band_l, this->n_basis}));
    this->hw = std::move(ct::Tensor(t_type, device_type, {this->n_band_l, this->n_basis}));
    this->p = std::move(ct::Tensor(t_type, device_type, {this->n_band_l, this->n_basis}));
    this->hp = std::move(ct::Tensor(t_type, device_type, {this->n_band_l, this->n_basis}));
    this->work = std::move(ct::Tensor(t_type, device_type, {this->n_band_l, this->n_basis}));

    this->prec = std::move(ct::Tensor(r_type, device_type, {this->n_basis}));

    this->nlocked = 0;
    this->eigen_locked.resize(this->n_band, static_cast<Real>(0.0));

#ifdef __MPI
    this->pmmcn.set_dimension(BP_WORLD, POOL_WORLD, n_band_l, n_basis, n_band_l, n_basis, n_dim, n_band);
    this->plintrans.set_dimension(n_dim, nband_l, n_band_l, n_basis, BP_WORLD, false);

    this->all_n_band_l.resize(this->plintrans.nproc_col);
    MPI_Allgather(&this->n_band_l, 1, MPI_INT, this->all_n_band_l.data(), 1, MPI_INT, BP_WORLD);
    this->band_displs.resize(this->plintrans.nproc_col);
    this->band_displs[0] = 0;
    for (int i = 1; i < this->plintrans.nproc_col; ++i)
    {
        this->band_displs[i] = this->band_displs[i - 1] + this->all_n_band_l[i - 1];
    }
#else
    this->pmmcn.set_dimension(n_band_l, n_basis, n_band_l, n_basis, n_dim, n_band);
    this->plintrans.set_dimension(n_dim, nband_l, n_band_l, n_basis, false);
    this->all_n_band_l = {this->n_band_l};
    this->band_displs = {0};
#endif
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::calc_prec()
{
    syncmem_var_h2d_op()(this->prec.template data<Real>(), this->h_prec.template data<Real>(), this->n_basis);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::calc_hpsi(const HPsiFunc& hpsi_func, T* psi_in, ct::Tensor& hpsi_out)
{
    hpsi_func(psi_in, hpsi_out.data<T>(), this->n_basis, this->n_band_l);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::calc_grad(const ct::Tensor& prec_in,
                                     ct::Tensor& err_out,
                                     ct::Tensor& psi_in,
                                     ct::Tensor& hpsi_in,
                                     ct::Tensor& grad_out,
                                     const int nlocked_in)
{
    int start_nband = 0;
#ifdef __MPI
    if (this->plintrans.nproc_col > 1)
    {
        start_nband = this->plintrans.start_colB[GlobalV::MY_BNDGROUP];
    }
#endif
    int local_nlocked = std::max(0, nlocked_in - start_nband);
    local_nlocked = std::min(local_nlocked, this->n_band_l);

    // Zero out locked bands
    for (int ib = 0; ib < local_nlocked; ++ib)
    {
        setmem_complex_op()(grad_out.data<T>() + ib * this->n_basis, 0, this->n_basis);
        err_out.data<Real>()[ib] = static_cast<Real>(0.0);
    }

    for (int ib = local_nlocked; ib < this->n_band_l; ++ib)
    {
        T* psi_ptr = psi_in.data<T>() + ib * this->n_basis;
        T* hpsi_ptr = hpsi_in.data<T>() + ib * this->n_basis;
        T* grad_ptr = grad_out.data<T>() + ib * this->n_basis;

        // 1. Normalize psi (and hpsi consistently)
        Real norm = ModuleBase::dot_real_op<T, Device>()(this->n_dim, psi_ptr, psi_ptr, true);
        norm = 1.0 / sqrt(norm);
        ModuleBase::vector_div_constant_op<T, Device>()(this->n_dim, psi_ptr, psi_ptr, norm);
        ModuleBase::vector_div_constant_op<T, Device>()(this->n_dim, hpsi_ptr, hpsi_ptr, norm);

        // 2. Rayleigh quotient: epsilo = <psi|hpsi>
        Real epsilo = ModuleBase::dot_real_op<T, Device>()(this->n_dim, psi_ptr, hpsi_ptr, true);

        // 3. Residual: grad = hpsi - epsilo * psi
        ModuleBase::vector_add_vector_op<T, Device>()(this->n_dim, grad_ptr, hpsi_ptr, 1.0, psi_ptr, -epsilo);

        // 4. Error = ||raw residual||
        Real err = ModuleBase::dot_real_op<T, Device>()(this->n_dim, grad_ptr, grad_ptr, true);
        err_out.data<Real>()[ib] = sqrt(err);

        // 5. Apply preconditioner: grad = grad / prec
        ModuleBase::vector_div_vector_op<T, Device>()(this->n_dim, grad_ptr, grad_ptr, prec_in.data<Real>());
    }
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::update_locking(const ct::Tensor& err_in, const std::vector<double>& ethr_band)
{
    // Gather local errors to global array
    std::vector<Real> local_err(this->n_band_l);
    if (err_in.device_type() == ct::DeviceType::GpuDevice)
    {
        syncmem_var_d2h_op()(local_err.data(), err_in.data<Real>(), this->n_band_l);
    }
    else
    {
        std::memcpy(local_err.data(), err_in.data<Real>(), this->n_band_l * sizeof(Real));
    }

    std::vector<Real> global_err(this->n_band, static_cast<Real>(0.0));
    std::vector<double> global_ethr(this->n_band, 0.0);

#ifdef __MPI
    MPI_Datatype mpi_real_type = (sizeof(Real) == sizeof(float)) ? MPI_FLOAT : MPI_DOUBLE;
    MPI_Allgatherv(local_err.data(),
                   this->n_band_l,
                   mpi_real_type,
                   global_err.data(),
                   this->all_n_band_l.data(),
                   this->band_displs.data(),
                   mpi_real_type,
                   BP_WORLD);

    std::vector<double> local_ethr_double(ethr_band.begin(), ethr_band.end());
    MPI_Allgatherv(local_ethr_double.data(),
                   this->n_band_l,
                   MPI_DOUBLE,
                   global_ethr.data(),
                   this->all_n_band_l.data(),
                   this->band_displs.data(),
                   MPI_DOUBLE,
                   BP_WORLD);
#else
    for (int i = 0; i < this->n_band_l; ++i)
    {
        global_err[i] = local_err[i];
        global_ethr[i] = ethr_band[i];
    }
#endif

    // Gather current eigenvalues from device
    std::vector<Real> current_eigen(this->n_band, static_cast<Real>(0.0));
    syncmem_var_d2h_op()(current_eigen.data(), this->eigen.data<Real>(), this->n_band);

    // Scan from current nlocked forward and lock converged bands
    while (this->nlocked < this->n_band)
    {
        if (global_err[this->nlocked] <= global_ethr[this->nlocked])
        {
            this->eigen_locked[this->nlocked] = current_eigen[this->nlocked];
            this->nlocked++;
        }
        else
        {
            break;
        }
    }
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::orth_projection(const ct::Tensor& psi_in, ct::Tensor& hsub_tmp, ct::Tensor& grad_out)
{
    // hsub_tmp = psi^H * grad (n_band x n_band global)
    this->pmmcn.multiply(1.0, psi_in.data<T>(), grad_out.data<T>(), 0.0, hsub_tmp.data<T>());

    // grad = grad - psi * hsub_tmp
    this->plintrans.act(-1.0, psi_in.data<T>(), hsub_tmp.data<T>(), 1.0, grad_out.data<T>());
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::orth_cholesky(ct::Tensor& workspace_in,
                                         ct::Tensor& psi_out,
                                         ct::Tensor& hpsi_out,
                                         ct::Tensor& hsub_out)
{
    // hsub_out = psi_out^H * psi_out
    this->pmmcn.multiply(1.0, psi_out.data<T>(), psi_out.data<T>(), 0.0, hsub_out.data<T>());

    ct::kernels::set_matrix<T, ct_Device>()('L', hsub_out.data<T>(), this->n_band);

    ct::kernels::lapack_potrf<T, ct_Device>()('U', this->n_band, hsub_out.data<T>(), this->n_band);
    ct::kernels::lapack_trtri<T, ct_Device>()('U', 'N', this->n_band, hsub_out.data<T>(), this->n_band);

    // Rotate psi and hpsi
    this->plintrans.act(1.0, psi_out.data<T>(), hsub_out.data<T>(), 0.0, workspace_in.data<T>());
    syncmem_complex_op()(psi_out.data<T>(), workspace_in.data<T>(), this->n_band_l * this->n_basis);

    this->plintrans.act(1.0, hpsi_out.data<T>(), hsub_out.data<T>(), 0.0, workspace_in.data<T>());
    syncmem_complex_op()(hpsi_out.data<T>(), workspace_in.data<T>(), this->n_band_l * this->n_basis);
}

template <typename T, typename Device>
bool DiagoPPCG<T, Device>::test_error(const ct::Tensor& err_in, const std::vector<double>& ethr_band)
{
    Real* _err_st = err_in.data<Real>();
    bool not_conv = false;
    std::vector<Real> tmp_cpu;
    if (err_in.device_type() == ct::DeviceType::GpuDevice)
    {
        tmp_cpu.resize(this->n_band_l);
        _err_st = tmp_cpu.data();
        syncmem_var_d2h_op()(_err_st, err_in.data<Real>(), this->n_band_l);
    }
    for (int ii = 0; ii < this->n_band_l; ++ii)
    {
        if (_err_st[ii] > ethr_band[ii])
        {
            not_conv = true;
        }
    }
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &not_conv, 1, MPI_C_BOOL, MPI_LOR, BP_WORLD);
#endif
    return not_conv;
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::diag_subspace(const HPsiFunc& hpsi_func,
                                         const bool has_p,
                                         ct::Tensor& psi_out,
                                         ct::Tensor& hpsi_out,
                                         ct::Tensor& p_out,
                                         ct::Tensor& hp_out,
                                         const int nlocked_in)
{
    const int n_sub = has_p ? 3 * this->n_band : 2 * this->n_band;

    // 1. Compute H|w>
    this->calc_hpsi(hpsi_func, this->w.data<T>(), this->hw);

    // 2. Compute overlap (S) and Hamiltonian (H) projection blocks.
    //    Only upper-triangular blocks are computed explicitly;
    //    lower-triangular parts are filled by Hermitian conjugate.
    ct::Tensor b_00(t_type, device_type, {this->n_band, this->n_band});
    ct::Tensor b_01(t_type, device_type, {this->n_band, this->n_band});
    ct::Tensor b_11(t_type, device_type, {this->n_band, this->n_band});

    this->pmmcn.multiply(one_, psi_out.data<T>(), psi_out.data<T>(), zero_, b_00.data<T>());
    this->pmmcn.multiply(one_, psi_out.data<T>(), this->w.data<T>(), zero_, b_01.data<T>());
    this->pmmcn.multiply(one_, this->w.data<T>(), this->w.data<T>(), zero_, b_11.data<T>());

    ct::Tensor bh_00(t_type, device_type, {this->n_band, this->n_band});
    ct::Tensor bh_01(t_type, device_type, {this->n_band, this->n_band});
    ct::Tensor bh_11(t_type, device_type, {this->n_band, this->n_band});

    this->pmmcn.multiply(one_, psi_out.data<T>(), hpsi_out.data<T>(), zero_, bh_00.data<T>());
    this->pmmcn.multiply(one_, psi_out.data<T>(), this->hw.data<T>(), zero_, bh_01.data<T>());
    this->pmmcn.multiply(one_, this->w.data<T>(), this->hw.data<T>(), zero_, bh_11.data<T>());

    ct::Tensor b_02, b_12, b_22, bh_02, bh_12, bh_22;
    if (has_p)
    {
        b_02 = ct::Tensor(t_type, device_type, {this->n_band, this->n_band});
        b_12 = ct::Tensor(t_type, device_type, {this->n_band, this->n_band});
        b_22 = ct::Tensor(t_type, device_type, {this->n_band, this->n_band});
        bh_02 = ct::Tensor(t_type, device_type, {this->n_band, this->n_band});
        bh_12 = ct::Tensor(t_type, device_type, {this->n_band, this->n_band});
        bh_22 = ct::Tensor(t_type, device_type, {this->n_band, this->n_band});

        this->pmmcn.multiply(one_, psi_out.data<T>(), p_out.data<T>(), zero_, b_02.data<T>());
        this->pmmcn.multiply(one_, this->w.data<T>(), p_out.data<T>(), zero_, b_12.data<T>());
        this->pmmcn.multiply(one_, p_out.data<T>(), p_out.data<T>(), zero_, b_22.data<T>());

        this->pmmcn.multiply(one_, psi_out.data<T>(), hp_out.data<T>(), zero_, bh_02.data<T>());
        this->pmmcn.multiply(one_, this->w.data<T>(), hp_out.data<T>(), zero_, bh_12.data<T>());
        this->pmmcn.multiply(one_, p_out.data<T>(), hp_out.data<T>(), zero_, bh_22.data<T>());
    }

    // 3. Assemble projected matrices on CPU
    ct::Tensor hsub_cpu(t_type, ct::DeviceType::CpuDevice, {n_sub, n_sub});
    ct::Tensor ssub_cpu(t_type, ct::DeviceType::CpuDevice, {n_sub, n_sub});
    ct::Tensor vcc_cpu(t_type, ct::DeviceType::CpuDevice, {n_sub, n_sub});
    ct::Tensor eigen_cpu(r_type, ct::DeviceType::CpuDevice, {n_sub});

    // Helper to copy block and optionally Hermitian-conjugate transpose
    auto copy_block = [&](const ct::Tensor& dev_block, int row_off, int col_off, bool to_h, bool hc) {
        std::vector<T> tmp(this->n_band * this->n_band);
        syncmem_complex_d2h_op()(tmp.data(), dev_block.data<T>(), this->n_band * this->n_band);
        T* dest = to_h ? hsub_cpu.data<T>() : ssub_cpu.data<T>();
        for (int j = 0; j < this->n_band; ++j)
        {
            for (int i = 0; i < this->n_band; ++i)
            {
                T val = hc ? std::conj(tmp[j + i * this->n_band]) : tmp[i + j * this->n_band];
                dest[(row_off + i) + (col_off + j) * n_sub] = val;
            }
        }
    };

    // S_sub assembly
    copy_block(b_00, 0, 0, false, false);
    copy_block(b_01, 0, this->n_band, false, false);
    copy_block(b_01, this->n_band, 0, false, true); // b_10 = b_01^H
    copy_block(b_11, this->n_band, this->n_band, false, false);

    // H_sub assembly
    copy_block(bh_00, 0, 0, true, false);
    copy_block(bh_01, 0, this->n_band, true, false);
    copy_block(bh_01, this->n_band, 0, true, true); // bh_10 = bh_01^H
    copy_block(bh_11, this->n_band, this->n_band, true, false);

    if (has_p)
    {
        copy_block(b_02, 0, 2 * this->n_band, false, false);
        copy_block(b_02, 2 * this->n_band, 0, false, true);
        copy_block(b_12, this->n_band, 2 * this->n_band, false, false);
        copy_block(b_12, 2 * this->n_band, this->n_band, false, true);
        copy_block(b_22, 2 * this->n_band, 2 * this->n_band, false, false);

        copy_block(bh_02, 0, 2 * this->n_band, true, false);
        copy_block(bh_02, 2 * this->n_band, 0, true, true);
        copy_block(bh_12, this->n_band, 2 * this->n_band, true, false);
        copy_block(bh_12, 2 * this->n_band, this->n_band, true, true);
        copy_block(bh_22, 2 * this->n_band, 2 * this->n_band, true, false);
    }

    // 4. Freeze locked bands: force their rows/columns to diagonal standard basis
    if (nlocked_in > 0)
    {
        for (int i = 0; i < nlocked_in; ++i)
        {
            for (int j = 0; j < n_sub; ++j)
            {
                T s_val = (j == i) ? one_ : zero_;
                T h_val = (j == i) ? static_cast<T>(this->eigen_locked[i]) : zero_;
                hsub_cpu.data<T>()[i + j * n_sub] = h_val;
                hsub_cpu.data<T>()[j + i * n_sub] = h_val;
                ssub_cpu.data<T>()[i + j * n_sub] = s_val;
                ssub_cpu.data<T>()[j + i * n_sub] = s_val;
            }
        }
    }

    // 5. Solve generalized eigenvalue problem H_sub * v = lambda * S_sub * v
    hsolver::hegvd_op<T, base_device::DEVICE_CPU>()(nullptr,
                                                    n_sub,
                                                    n_sub,
                                                    hsub_cpu.data<T>(),
                                                    ssub_cpu.data<T>(),
                                                    eigen_cpu.data<Real>(),
                                                    vcc_cpu.data<T>());

    // Ensure locked eigenvalues remain unchanged (overwrite in case of numerical drift)
    for (int i = 0; i < nlocked_in && i < this->n_band; ++i)
    {
        eigen_cpu.data<Real>()[i] = this->eigen_locked[i];
    }

    // 6. Move eigenvectors back to device
    ct::Tensor vcc_dev = vcc_cpu.to_device<ct_Device>();

    // 7. Update psi = X*C_X + W*C_W + (P*C_P)
    setmem_complex_op()(this->work.data<T>(), 0, this->n_band_l * this->n_basis);
    this->plintrans.act(1.0, psi_out.data<T>(), vcc_dev.data<T>() + 0, 0.0, this->work.data<T>());
    this->plintrans.act(1.0, this->w.data<T>(), vcc_dev.data<T>() + this->n_band, 1.0, this->work.data<T>());
    if (has_p)
    {
        this->plintrans.act(1.0, p_out.data<T>(), vcc_dev.data<T>() + 2 * this->n_band, 1.0, this->work.data<T>());
    }
    syncmem_complex_op()(psi_out.data<T>(), this->work.data<T>(), this->n_band_l * this->n_basis);

    // 8. Update hpsi = HX*C_X + HW*C_W + (HP*C_P)
    setmem_complex_op()(this->work.data<T>(), 0, this->n_band_l * this->n_basis);
    this->plintrans.act(1.0, hpsi_out.data<T>(), vcc_dev.data<T>() + 0, 0.0, this->work.data<T>());
    this->plintrans.act(1.0, this->hw.data<T>(), vcc_dev.data<T>() + this->n_band, 1.0, this->work.data<T>());
    if (has_p)
    {
        this->plintrans.act(1.0, hp_out.data<T>(), vcc_dev.data<T>() + 2 * this->n_band, 1.0, this->work.data<T>());
    }
    syncmem_complex_op()(hpsi_out.data<T>(), this->work.data<T>(), this->n_band_l * this->n_basis);

    // 9. Update p = W*C_W + (P*C_P)  -- LOBPCG style, no X component
    setmem_complex_op()(this->work.data<T>(), 0, this->n_band_l * this->n_basis);
    this->plintrans.act(1.0, this->w.data<T>(), vcc_dev.data<T>() + this->n_band, 0.0, this->work.data<T>());
    if (has_p)
    {
        this->plintrans.act(1.0, p_out.data<T>(), vcc_dev.data<T>() + 2 * this->n_band, 1.0, this->work.data<T>());
    }
    syncmem_complex_op()(p_out.data<T>(), this->work.data<T>(), this->n_band_l * this->n_basis);

    // 10. Update hp = HW*C_W + (HP*C_P)
    setmem_complex_op()(this->work.data<T>(), 0, this->n_band_l * this->n_basis);
    this->plintrans.act(1.0, this->hw.data<T>(), vcc_dev.data<T>() + this->n_band, 0.0, this->work.data<T>());
    if (has_p)
    {
        this->plintrans.act(1.0, hp_out.data<T>(), vcc_dev.data<T>() + 2 * this->n_band, 1.0, this->work.data<T>());
    }
    syncmem_complex_op()(hp_out.data<T>(), this->work.data<T>(), this->n_band_l * this->n_basis);

    // 11. Update eigenvalues with the lowest n_band eigenvalues from subspace diagonalization
    syncmem_var_h2d_op()(this->eigen.data<Real>(), eigen_cpu.data<Real>(), this->n_band);
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::diag(const HPsiFunc& hpsi_func,
                                T* psi_in,
                                Real* eigenvalue_in,
                                const std::vector<double>& ethr_band)
{
    const int current_scf_iter = hsolver::DiagoIterAssist<T, Device>::SCF_ITER;

    // Map the input psi pointer
    this->psi = std::move(ct::TensorMap(psi_in, t_type, device_type, {this->n_band_l, this->n_basis}));

    // Update precondition array
    this->calc_prec();

    // Initial subspace diagonalization to improve the initial guess
    this->calc_hpsi(hpsi_func, psi_in, this->hpsi);

    // Build and diagonalize the subspace Hamiltonian in the psi basis
    ct::Tensor hsub_init(t_type, device_type, {this->n_band, this->n_band});
    this->pmmcn.multiply(one_, this->hpsi.data<T>(), this->psi.data<T>(), zero_, hsub_init.data<T>());
    ct::kernels::lapack_heevd<T, ct_Device>()(this->n_band,
                                              hsub_init.data<T>(),
                                              this->n_band,
                                              this->eigen.data<Real>());

    // Rotate psi and hpsi with the eigenvectors of the subspace problem
    this->plintrans.act(1.0, this->psi.data<T>(), hsub_init.data<T>(), 0.0, this->work.data<T>());
    syncmem_complex_op()(this->psi.data<T>(), this->work.data<T>(), this->n_band_l * this->n_basis);
    this->plintrans.act(1.0, this->hpsi.data<T>(), hsub_init.data<T>(), 0.0, this->work.data<T>());
    syncmem_complex_op()(this->hpsi.data<T>(), this->work.data<T>(), this->n_band_l * this->n_basis);

    // Initialize search direction to zero
    setmem_complex_op()(this->p.data<T>(), 0, this->n_band_l * this->n_basis);
    setmem_complex_op()(this->hp.data<T>(), 0, this->n_band_l * this->n_basis);

    // Allocate a reusable tensor for projection overlap
    ct::Tensor hsub_orth(t_type, device_type, {this->n_band, this->n_band});

    int ntry = 0;
    int max_iter = current_scf_iter > 1 ? this->nline : this->nline * 6;
    this->nlocked = 0;

    do
    {
        ++ntry;

        // 1. Calculate preconditioned residual w and error for active bands only
        this->calc_grad(this->prec, this->err_st, this->psi, this->hpsi, this->w, this->nlocked);

        // 2. Update locking status: scan from current nlocked forward
        this->update_locking(this->err_st, ethr_band);

        // 3. Exit if all bands have converged
        if (this->nlocked >= this->n_band)
        {
            break;
        }

        // 4. Project active residual to orthogonal complement of psi
        this->orth_projection(this->psi, hsub_orth, this->w);

        // 5. Expanded subspace diagonalization with locking
        //    Locked bands are frozen in the subspace problem
        this->diag_subspace(hpsi_func, ntry > 1, this->psi, this->hpsi, this->p, this->hp, this->nlocked);

        // Note: orth_cholesky is intentionally skipped here.
        // The Rayleigh-Ritz step already provides orthonormal vectors (within numerical precision).
        // Global Cholesky would destroy the locking by remixing all bands.

    } while (ntry < max_iter && this->nlocked < this->n_band);

    // Final subspace diagonalization to obtain accurate eigenvalues
    this->pmmcn.multiply(one_, this->hpsi.data<T>(), this->psi.data<T>(), zero_, hsub_orth.data<T>());
    ct::kernels::lapack_heevd<T, ct_Device>()(this->n_band,
                                              hsub_orth.data<T>(),
                                              this->n_band,
                                              this->eigen.data<Real>());
    this->plintrans.act(1.0, this->psi.data<T>(), hsub_orth.data<T>(), 0.0, this->work.data<T>());
    syncmem_complex_op()(this->psi.data<T>(), this->work.data<T>(), this->n_band_l * this->n_basis);

    // Copy eigenvalues to output
    int start_nband = 0;
#ifdef __MPI
    if (this->plintrans.nproc_col > 1)
    {
        start_nband = this->plintrans.start_colB[GlobalV::MY_BNDGROUP];
    }
#endif
    syncmem_var_d2h_op()(eigenvalue_in, this->eigen.data<Real>() + start_nband, this->n_band_l);
}

// Explicit template instantiations
template class DiagoPPCG<std::complex<float>, base_device::DEVICE_CPU>;
template class DiagoPPCG<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class DiagoPPCG<std::complex<float>, base_device::DEVICE_GPU>;
template class DiagoPPCG<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace hsolver
