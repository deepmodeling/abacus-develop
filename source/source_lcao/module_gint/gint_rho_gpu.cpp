#include "gint_rho_gpu.h"
#include "gint_common.h"
#include "gint_helper.h"
#include "batch_biggrid.h"
#include "kernel/phi_operator_gpu.h"

namespace ModuleGint
{

void Gint_rho_gpu::cal_gint()
{
    ModuleBase::TITLE("Gint", "cal_gint_rho");
    ModuleBase::timer::tick("Gint", "cal_gint_rho");
    init_dm_gint_();
    transfer_dm_2d_to_gint(*gint_info_, dm_vec_, dm_gint_vec_);
    cal_rho_();
    ModuleBase::timer::tick("Gint", "cal_gint_rho");
}

void Gint_rho_gpu::init_dm_gint_()
{
    dm_gint_vec_.resize(nspin_);
    for (int is = 0; is < nspin_; is++)
    {
        dm_gint_vec_[is] = gint_info_->get_hr<double>();
    }
}

void Gint_rho_gpu::transfer_cpu_to_gpu_()
{
    dm_gint_d_vec_.resize(nspin_);
    rho_d_vec_.resize(nspin_);
    for (int is = 0; is < nspin_; is++)
    {
        dm_gint_d_vec_[is] = CudaMemWrapper<double>(dm_gint_vec_[is].get_nnr(), 0, false);
        rho_d_vec_[is] = CudaMemWrapper<double>(gint_info_->get_local_mgrid_num(), 0, false);
        checkCuda(cudaMemcpy(dm_gint_d_vec_[is].get_device_ptr(), dm_gint_vec_[is].get_wrapper(), 
            dm_gint_vec_[is].get_nnr() * sizeof(double), cudaMemcpyHostToDevice));
    }
}

void Gint_rho_gpu::transfer_gpu_to_cpu_()
{
    for (int is = 0; is < nspin_; is++)
    {
        checkCuda(cudaMemcpy(rho_[is], rho_d_vec_[is].get_device_ptr(), 
            gint_info_->get_local_mgrid_num() * sizeof(double), cudaMemcpyDeviceToHost));
    }
}

void Gint_rho_gpu::cal_rho_()
{
    transfer_cpu_to_gpu_();
#pragma omp parallel num_threads(gint_info_->get_streams_num())
    {
        // 20240620 Note that it must be set again here because 
        // cuda's device is not safe in a multi-threaded environment.
        checkCuda(cudaSetDevice(gint_info_->get_dev_id()));
        cudaStream_t stream;
        checkCuda(cudaStreamCreate(&stream));
        PhiOperatorGpu phi_op(gint_info_->get_gpu_vars(), stream);
        CudaMemWrapper<double> phi(BatchBigGrid::get_max_phi_len(), stream, false);
        CudaMemWrapper<double> phi_dm(BatchBigGrid::get_max_phi_len(), stream, false);
        #pragma omp for schedule(dynamic)
        for(const auto& bgrid_batch: gint_info_->get_bgrid_batches())
        {
            if(bgrid_batch->empty())
            {
                continue;
            }
            phi_op.set_bgrid_batch(bgrid_batch);
            phi_op.set_phi(phi.get_device_ptr());
            for(int is = 0; is < nspin_; is++)
            {
                phi_op.phi_mul_dm(phi.get_device_ptr(), dm_gint_d_vec_[is].get_device_ptr(), dm_gint_vec_[is],
                                  is_dm_symm_, phi_dm.get_device_ptr());
                phi_op.phi_dot_phi(phi.get_device_ptr(), phi_dm.get_device_ptr(), rho_d_vec_[is].get_device_ptr());
            }
       }
       checkCuda(cudaStreamSynchronize(stream));
       checkCuda(cudaStreamDestroy(stream));
    }
    transfer_gpu_to_cpu_();
}

}  // namespace ModuleGint