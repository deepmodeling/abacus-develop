#include "gint_vl_gpu.h"
#include "gint_common.h"
#include "gint_helper.h"
#include "batch_biggrid.h"
#include "kernel/phi_operator_gpu.h"
#include "source_base/module_device/device_check.h"

#include <algorithm>
#include <type_traits>

namespace ModuleGint
{

void Gint_vl_gpu::cal_gint()
{
    ModuleBase::TITLE("Gint", "cal_gint_vl");
    ModuleBase::timer::start("Gint", "cal_gint_vl");
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
    ModuleBase::timer::end("Gint", "cal_gint_vl");
}

// Helper: finalize hr_gint (double path — no cast needed)
inline void finalize_hr_gint_gpu_(HContainer<double>& hr_gint, HContainer<double>* hR)
{
    compose_hr_gint(hr_gint);
    transfer_hr_gint_to_hR(hr_gint, *hR);
}

// Helper: finalize hr_gint (non-double path — cast to double first)
template<typename Real>
void finalize_hr_gint_gpu_(HContainer<Real>& hr_gint, HContainer<double>* hR)
{
    HContainer<double> hr_gint_dp = make_cast_hcontainer<double>(hr_gint);
    compose_hr_gint(hr_gint_dp);
    transfer_hr_gint_to_hR(hr_gint_dp, *hR);
}

template<typename Real>
void Gint_vl_gpu::cal_gint_impl_()
{
    // 1. Initialize hr_gint as HContainer<Real>
    HContainer<Real> hr_gint = gint_info_->get_hr<Real>();

    // 2. Convert vr_eff to Real and transfer to GPU
    const int local_mgrid_num = gint_info_->get_local_mgrid_num();
    CudaMemWrapper<Real> vr_eff_d(local_mgrid_num, 0, false);
    CudaMemWrapper<Real> hr_gint_d(hr_gint.get_nnr(), 0, false);

    if (std::is_same<Real, double>::value)
    {
        // No conversion needed
        CHECK_CUDA(cudaMemcpy(vr_eff_d.get_device_ptr(), reinterpret_cast<const Real*>(vr_eff_),
            local_mgrid_num * sizeof(Real), cudaMemcpyHostToDevice));
    }
    else
    {
        // Convert double vr_eff to Real (float)
        std::vector<Real> vr_eff_buffer(local_mgrid_num);
        std::transform(vr_eff_, vr_eff_ + local_mgrid_num, vr_eff_buffer.begin(),
            [](const double v) { return static_cast<Real>(v); });
        CHECK_CUDA(cudaMemcpy(vr_eff_d.get_device_ptr(), vr_eff_buffer.data(),
            local_mgrid_num * sizeof(Real), cudaMemcpyHostToDevice));
    }

    // 3. Calculate hr_gint on GPU
#pragma omp parallel num_threads(gint_info_->get_streams_num())
    {
        // 20240620 Note that it must be set again here because
        // cuda's device is not safe in a multi-threaded environment.
        CHECK_CUDA(cudaSetDevice(gint_info_->get_dev_id()));
        cudaStream_t stream;
        CHECK_CUDA(cudaStreamCreate(&stream));
        PhiOperatorGpu<Real> phi_op(gint_info_->get_gpu_vars(), stream);
        CudaMemWrapper<Real> phi(BatchBigGrid::get_max_phi_len(), stream, false);
        CudaMemWrapper<Real> phi_vldr3(BatchBigGrid::get_max_phi_len(), stream, false);
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
            phi_op.phi_mul_vldr3(vr_eff_d.get_device_ptr(), static_cast<Real>(dr3_),
                 phi.get_device_ptr(), phi_vldr3.get_device_ptr());
            phi_op.phi_mul_phi(phi.get_device_ptr(), phi_vldr3.get_device_ptr(),
                 hr_gint, hr_gint_d.get_device_ptr());
        }
        CHECK_CUDA(cudaStreamSynchronize(stream));
        CHECK_CUDA(cudaStreamDestroy(stream));
    }

    // 4. Transfer hr_gint back to CPU
    CHECK_CUDA(cudaMemcpy(hr_gint.get_wrapper(), hr_gint_d.get_device_ptr(),
        hr_gint.get_nnr() * sizeof(Real), cudaMemcpyDeviceToHost));

    // 5. Compose and transfer to hR (with cast if needed)
    finalize_hr_gint_gpu_(hr_gint, hR_);
}

}