#include "gint_tau_gpu.h"
#include "gint_common.h"
#include "gint_helper.h"
#include "batch_biggrid.h"
#include "kernel/phi_operator_gpu.h"

namespace ModuleGint
{

void Gint_tau_gpu::cal_gint()
{
    ModuleBase::TITLE("Gint", "cal_gint_tau");
    ModuleBase::timer::tick("Gint", "cal_gint_tau");
    init_dm_gint_();
    transfer_dm_2d_to_gint(*gint_info_, dm_vec_, dm_gint_vec_);
    cal_tau_();
    ModuleBase::timer::tick("Gint", "cal_gint_tau");
}

void Gint_tau_gpu::init_dm_gint_()
{
    dm_gint_vec_.resize(nspin_);
    for (int is = 0; is < nspin_; is++)
    {
        dm_gint_vec_[is] = gint_info_->get_hr<double>();
    }
}

void Gint_tau_gpu::transfer_cpu_to_gpu_()
{
    dm_gint_d_vec_.resize(nspin_);
    kin_d_vec_.resize(nspin_);
    for (int is = 0; is < nspin_; is++)
    {
        dm_gint_d_vec_[is] = CudaMemWrapper<double>(dm_gint_vec_[is].get_nnr(), 0, false);
        kin_d_vec_[is] = CudaMemWrapper<double>(gint_info_->get_local_mgrid_num(), 0, false);
        checkCuda(cudaMemcpy(dm_gint_d_vec_[is].get_device_ptr(), dm_gint_vec_[is].get_wrapper(), 
            dm_gint_vec_[is].get_nnr() * sizeof(double), cudaMemcpyHostToDevice));
    }
}

void Gint_tau_gpu::transfer_gpu_to_cpu_()
{
    for (int is = 0; is < nspin_; is++)
    {
        checkCuda(cudaMemcpy(kin_[is], kin_d_vec_[is].get_device_ptr(), 
            gint_info_->get_local_mgrid_num() * sizeof(double), cudaMemcpyDeviceToHost));
    }
}

void Gint_tau_gpu::cal_tau_()
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
        CudaMemWrapper<double> dphi_x(BatchBigGrid::get_max_phi_len(), stream, false);
        CudaMemWrapper<double> dphi_y(BatchBigGrid::get_max_phi_len(), stream, false);
        CudaMemWrapper<double> dphi_z(BatchBigGrid::get_max_phi_len(), stream, false);
        CudaMemWrapper<double> dphi_x_dm(BatchBigGrid::get_max_phi_len(), stream, false);
        CudaMemWrapper<double> dphi_y_dm(BatchBigGrid::get_max_phi_len(), stream, false);
        CudaMemWrapper<double> dphi_z_dm(BatchBigGrid::get_max_phi_len(), stream, false);
        #pragma omp for schedule(dynamic)
        for(const auto& bgrid_batch: gint_info_->get_bgrid_batches())
        {
            if(bgrid_batch->empty())
            {
                continue;
            }
            phi_op.set_bgrid_batch(bgrid_batch);
            phi_op.set_phi_dphi(nullptr,
                                dphi_x.get_device_ptr(), dphi_y.get_device_ptr(), dphi_z.get_device_ptr());
            for(int is = 0; is < nspin_; is++)
            {
                const bool is_symm = true;
                phi_op.phi_mul_dm(dphi_x.get_device_ptr(), dm_gint_d_vec_[is].get_device_ptr(), dm_gint_vec_[is],
                                  is_symm, dphi_x_dm.get_device_ptr());
                phi_op.phi_mul_dm(dphi_y.get_device_ptr(), dm_gint_d_vec_[is].get_device_ptr(), dm_gint_vec_[is],
                                  is_symm, dphi_y_dm.get_device_ptr());
                phi_op.phi_mul_dm(dphi_z.get_device_ptr(), dm_gint_d_vec_[is].get_device_ptr(), dm_gint_vec_[is],
                                  is_symm, dphi_z_dm.get_device_ptr());
                phi_op.phi_dot_phi(dphi_x.get_device_ptr(), dphi_x_dm.get_device_ptr(), kin_d_vec_[is].get_device_ptr());
                phi_op.phi_dot_phi(dphi_y.get_device_ptr(), dphi_y_dm.get_device_ptr(), kin_d_vec_[is].get_device_ptr());
                phi_op.phi_dot_phi(dphi_z.get_device_ptr(), dphi_z_dm.get_device_ptr(), kin_d_vec_[is].get_device_ptr());
            }
       }
       checkCuda(cudaStreamSynchronize(stream));
       checkCuda(cudaStreamDestroy(stream));
    }
    transfer_gpu_to_cpu_();
}

}
