#include "gint_rho_gpu.h"
#include "gint_common.h"
#include "gint_helper.h"
#include "batch_biggrid.h"
#include "kernel/phi_operator_gpu.h"
#include "source_base/module_device/device_check.h"

#include <algorithm>

namespace ModuleGint
{

void Gint_rho_gpu::cal_gint()
{
    ModuleBase::TITLE("Gint", "cal_gint_rho");
    ModuleBase::timer::start("Gint", "cal_gint_rho");
    switch (gint_info_->get_exec_precision())
    {
    case GintPrecision::fp32:
        cal_gint_impl_<float>();
        break;
    case GintPrecision::fp64:
    default:
        cal_gint_impl_<double>();
        break;
    }
    ModuleBase::timer::end("Gint", "cal_gint_rho");
}

template<typename Real>
void Gint_rho_gpu::cal_gint_impl_()
{
    // 1. Initialize dm_gint as HContainer<Real>
    std::vector<HContainer<Real>> dm_gint_vec(nspin_);
    for (int is = 0; is < nspin_; is++)
    {
        dm_gint_vec[is] = gint_info_->get_hr<Real>();
    }

    // 2. Transfer dm from 2D parallel distribution to gint serial distribution
    transfer_dm_2d_to_gint(*gint_info_, dm_vec_, dm_gint_vec);

    // 3. Transfer dm to GPU
    std::vector<CudaMemWrapper<Real>> dm_gint_d_vec(nspin_);
    std::vector<CudaMemWrapper<Real>> rho_d_vec(nspin_);
    for (int is = 0; is < nspin_; is++)
    {
        dm_gint_d_vec[is] = CudaMemWrapper<Real>(dm_gint_vec[is].get_nnr(), 0, false);
        rho_d_vec[is] = CudaMemWrapper<Real>(gint_info_->get_local_mgrid_num(), 0, false);
        CHECK_CUDA(cudaMemcpy(dm_gint_d_vec[is].get_device_ptr(), dm_gint_vec[is].get_wrapper(),
            dm_gint_vec[is].get_nnr() * sizeof(Real), cudaMemcpyHostToDevice));
    }

    // 4. Calculate rho on GPU
#pragma omp parallel num_threads(gint_info_->get_streams_num())
    {
        // 20240620 Note that it must be set again here because
        // cuda's device is not safe in a multi-threaded environment.
        CHECK_CUDA(cudaSetDevice(gint_info_->get_dev_id()));
        cudaStream_t stream;
        CHECK_CUDA(cudaStreamCreate(&stream));
        PhiOperatorGpu<Real> phi_op(gint_info_->get_gpu_vars(), stream);
        CudaMemWrapper<Real> phi(BatchBigGrid::get_max_phi_len(), stream, false);
        CudaMemWrapper<Real> phi_dm(BatchBigGrid::get_max_phi_len(), stream, false);
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < gint_info_->get_bgrid_batches_num(); ++i)
        {
            const auto& bgrid_batch = gint_info_->get_bgrid_batches()[i];
            if(bgrid_batch->empty())
            {
                continue;
            }
            phi_op.set_bgrid_batch(bgrid_batch);
            phi_op.set_phi(phi.get_device_ptr());
            for(int is = 0; is < nspin_; is++)
            {
                phi_op.phi_mul_dm(phi.get_device_ptr(), dm_gint_d_vec[is].get_device_ptr(), dm_gint_vec[is],
                                  is_dm_symm_, phi_dm.get_device_ptr());
                phi_op.phi_dot_phi(phi.get_device_ptr(), phi_dm.get_device_ptr(), rho_d_vec[is].get_device_ptr());
            }
       }
       CHECK_CUDA(cudaStreamSynchronize(stream));
       CHECK_CUDA(cudaStreamDestroy(stream));
    }

    // 5. Transfer rho back to CPU and convert to double if needed
    const int local_mgrid_num = gint_info_->get_local_mgrid_num();
    for (int is = 0; is < nspin_; is++)
    {
        std::vector<Real> rho_tmp(local_mgrid_num);
        CHECK_CUDA(cudaMemcpy(rho_tmp.data(), rho_d_vec[is].get_device_ptr(),
            local_mgrid_num * sizeof(Real), cudaMemcpyDeviceToHost));
        for (int ir = 0; ir < local_mgrid_num; ++ir)
        {
            rho_[is][ir] = static_cast<double>(rho_tmp[ir]);
        }
    }
}

}  // namespace ModuleGint