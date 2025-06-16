#include "gint_fvl_gpu.h"
#include "gint_common.h"
#include "gint_helper.h"
#include "batch_biggrid.h"
#include "kernel/phi_operator_gpu.h"

namespace ModuleGint
{

void Gint_fvl_gpu::cal_gint()
{
    ModuleBase::TITLE("Gint", "cal_gint_fvl");
    ModuleBase::timer::tick("Gint", "cal_gint_fvl");
    init_dm_gint_();
    transfer_dm_2d_to_gint(*gint_info_, dm_vec_, dm_gint_vec_);
    cal_fvl_svl_();
    ModuleBase::timer::tick("Gint", "cal_gint_fvl");
}

void Gint_fvl_gpu::init_dm_gint_()
{
    dm_gint_vec_.resize(nspin_);
    for (int is = 0; is < nspin_; is++)
    {
        dm_gint_vec_[is] = gint_info_->get_hr<double>();
    }
}

void Gint_fvl_gpu::transfer_cpu_to_gpu_()
{
    dm_gint_d_vec_.resize(nspin_);
    vr_eff_d_vec_.resize(nspin_);
    for (int is = 0; is < nspin_; is++)
    {
        dm_gint_d_vec_[is] = CudaMemWrapper<double>(dm_gint_vec_[is].get_nnr(), 0, false);
        checkCuda(cudaMemcpy(dm_gint_d_vec_[is].get_device_ptr(), dm_gint_vec_[is].get_wrapper(), 
                             dm_gint_vec_[is].get_nnr() * sizeof(double), cudaMemcpyHostToDevice));
        vr_eff_d_vec_[is] = CudaMemWrapper<double>(gint_info_->get_local_mgrid_num(), 0, false);
        checkCuda(cudaMemcpy(vr_eff_d_vec_[is].get_device_ptr(), vr_eff_[is],
                             gint_info_->get_local_mgrid_num() * sizeof(double), cudaMemcpyHostToDevice));
    }
    if (isforce_)
    {
        fvl_d_ = CudaMemWrapper<double>(gint_info_->get_nat() * 3, 0, true);
    }
    if (isstress_)
    {
        svl_d_ = CudaMemWrapper<double>(6, 0, true);
    }
}

void Gint_fvl_gpu::transfer_gpu_to_cpu_()
{
    if (isforce_)
    {
        fvl_d_.copy_device_to_host_sync();
        for (int iat = 0; iat < gint_info_->get_nat(); iat++)
        {
            for (int j = 0; j < 3; j++)
            {
                fvl_[0](iat, j) += fvl_d_.get_host_ptr()[iat * 3 + j];
            }
        }
    }
    if (isstress_)
    {
        svl_d_.copy_device_to_host_sync();
        svl_[0](0, 0) += svl_d_.get_host_ptr()[0];
        svl_[0](0, 1) += svl_d_.get_host_ptr()[1];
        svl_[0](0, 2) += svl_d_.get_host_ptr()[2];
        svl_[0](1, 1) += svl_d_.get_host_ptr()[3];
        svl_[0](1, 2) += svl_d_.get_host_ptr()[4];
        svl_[0](2, 2) += svl_d_.get_host_ptr()[5];
    }
}

void Gint_fvl_gpu::cal_fvl_svl_()
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
        CudaMemWrapper<double> phi_vldr3(BatchBigGrid::get_max_phi_len(), stream, false);
        CudaMemWrapper<double> phi_vldr3_dm(BatchBigGrid::get_max_phi_len(), stream, false);
        CudaMemWrapper<double> dphi_x(BatchBigGrid::get_max_phi_len(), stream, false);
        CudaMemWrapper<double> dphi_y(BatchBigGrid::get_max_phi_len(), stream, false);
        CudaMemWrapper<double> dphi_z(BatchBigGrid::get_max_phi_len(), stream, false);

        #pragma omp for schedule(dynamic)
        for(const auto& bgrid_batch: gint_info_->get_bgrid_batches())
        {
            if(bgrid_batch->empty())
            {
                continue;
            }
            phi_op.set_bgrid_batch(bgrid_batch);
            phi_op.set_phi_dphi(phi.get_device_ptr(),
                                dphi_x.get_device_ptr(),
                                dphi_y.get_device_ptr(),
                                dphi_z.get_device_ptr());
            for(int is = 0; is < nspin_; is++)
            {
                const bool is_symm = false;
                phi_op.phi_mul_vldr3(vr_eff_d_vec_[is].get_device_ptr(), dr3_,
                                     phi.get_device_ptr(), phi_vldr3.get_device_ptr());
                phi_op.phi_mul_dm(phi_vldr3.get_device_ptr(), dm_gint_d_vec_[is].get_device_ptr(),
                                  dm_gint_vec_[is], is_symm, phi_vldr3_dm.get_device_ptr());
                if (isforce_)
                {
                    phi_op.phi_dot_dphi(phi_vldr3_dm.get_device_ptr(),
                                        dphi_x.get_device_ptr(), dphi_y.get_device_ptr(),
                                        dphi_z.get_device_ptr(), fvl_d_.get_device_ptr());
                }
                if (isstress_)
                {
                    phi_op.phi_dot_dphi_r(phi_vldr3_dm.get_device_ptr(),
                                          dphi_x.get_device_ptr(), dphi_y.get_device_ptr(),
                                          dphi_z.get_device_ptr(), svl_d_.get_device_ptr());
                }
            }
       }
       checkCuda(cudaStreamSynchronize(stream));
       checkCuda(cudaStreamDestroy(stream));
    }
    transfer_gpu_to_cpu_();
}

}