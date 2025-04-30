#include "gint_vl_gpu.h"
#include "gint_common.h"
#include "gint_helper.h"
#include "batch_biggrid.h"
#include "kernel/phi_operator_gpu.h"

namespace ModuleGint
{

void Gint_vl_gpu::cal_gint()
{
    init_hr_gint_();
    transfer_cpu_to_gpu_();
    cal_hr_gint_();
    transfer_gpu_to_cpu_();
    compose_hr_gint(hr_gint_);
    transfer_hr_gint_to_hR(toConstSharedPtr(hr_gint_), hR_);
}

void Gint_vl_gpu::init_hr_gint_()
{
    hr_gint_ = gint_info_->get_hr<double>();
}

void Gint_vl_gpu::transfer_cpu_to_gpu_()
{
    hr_gint_d_ = CudaMemWrapper<double>(hr_gint_->get_nnr(), 0, false);
    vr_eff_d_ = CudaMemWrapper<double>(gint_info_->get_local_mgrid_num(), 0, false);
    checkCuda(cudaMemcpy(vr_eff_d_.get_device_ptr(), vr_eff_,
        gint_info_->get_local_mgrid_num() * sizeof(double), cudaMemcpyHostToDevice));
}

void Gint_vl_gpu::transfer_gpu_to_cpu_()
{
    checkCuda(cudaMemcpy(hr_gint_->get_wrapper(), hr_gint_d_.get_device_ptr(), 
        hr_gint_->get_nnr() * sizeof(double), cudaMemcpyDeviceToHost));
}

void Gint_vl_gpu::cal_hr_gint_()
{
#pragma omp parallel 
    {
        // 20240620 Note that it must be set again here because 
        // cuda's device is not safe in a multi-threaded environment.
        checkCuda(cudaSetDevice(gint_info_->get_dev_id()));
        cudaStream_t stream;
        checkCuda(cudaStreamCreate(&stream));
        PhiOperatorGpu phi_op(gint_info_->get_gpu_vars(), stream);
        CudaMemWrapper<double> phi(BatchBigGrid::get_max_phi_len(), stream, false);
        CudaMemWrapper<double> phi_vldr3(BatchBigGrid::get_max_phi_len(), stream, false);
        #pragma omp for schedule(dynamic)
        for(const auto& bgrid_batch: gint_info_->get_bgrid_batches())
        {
            if(bgrid_batch->empty())
            {
                continue;
            }
            phi_op.set_bgrid_batch(bgrid_batch);
            phi_op.set_phi(phi.get_device_ptr());
            phi_op.phi_mul_vldr3(vr_eff_d_.get_device_ptr(), dr3_,
                 phi.get_device_ptr(), phi_vldr3.get_device_ptr());
            phi_op.phi_mul_phi_vldr3(phi.get_device_ptr(), phi_vldr3.get_device_ptr(),
                 hr_gint_, hr_gint_d_.get_device_ptr());
        }
        checkCuda(cudaStreamSynchronize(stream));
        checkCuda(cudaStreamDestroy(stream));
    }
}

}