#include "gint_vl_metagga_nspin4_gpu.h"
#include "gint_common.h"
#include "gint_helper.h"
#include "batch_biggrid.h"
#include "kernel/phi_operator_gpu.h"

namespace ModuleGint
{

void Gint_vl_metagga_nspin4_gpu::cal_gint()
{
    ModuleBase::TITLE("Gint", "cal_gint_vl");
    ModuleBase::timer::tick("Gint", "cal_gint_vl");
    init_hr_gint_();
    cal_hr_gint_();
    compose_hr_gint(hr_gint_part_, hr_gint_full_);
    transfer_hr_gint_to_hR(hr_gint_full_, *hR_);
    ModuleBase::timer::tick("Gint", "cal_gint_vl");
}

void Gint_vl_metagga_nspin4_gpu::init_hr_gint_()
{
    hr_gint_part_.resize(nspin_);
    for(int i = 0; i < nspin_; i++)
    {
        hr_gint_part_[i] = gint_info_->get_hr<double>();
    }
    const int npol = 2;
    hr_gint_full_ = gint_info_->get_hr<std::complex<double>>(npol);
}

void Gint_vl_metagga_nspin4_gpu::transfer_cpu_to_gpu_()
{
    vr_eff_d_.resize(nspin_);
    vofk_d_.resize(nspin_);
    hr_gint_part_d_.resize(nspin_);
    for(int i = 0; i < nspin_; i++)
    {
        hr_gint_part_d_[i] = CudaMemWrapper<double>(hr_gint_part_[i].get_nnr(), 0, false);
        vr_eff_d_[i] = CudaMemWrapper<double>(gint_info_->get_local_mgrid_num(), 0, false);
        vofk_d_[i] = CudaMemWrapper<double>(gint_info_->get_local_mgrid_num(), 0, false);
        checkCuda(cudaMemcpy(vr_eff_d_[i].get_device_ptr(), vr_eff_[i],
                  gint_info_->get_local_mgrid_num() * sizeof(double), cudaMemcpyHostToDevice));
        checkCuda(cudaMemcpy(vofk_d_[i].get_device_ptr(), vofk_[i],
                  gint_info_->get_local_mgrid_num() * sizeof(double), cudaMemcpyHostToDevice));
    }
}

void Gint_vl_metagga_nspin4_gpu::transfer_gpu_to_cpu_()
{
    for(int i = 0; i < nspin_; i++)
    {
        checkCuda(cudaMemcpy(hr_gint_part_[i].get_wrapper(), hr_gint_part_d_[i].get_device_ptr(), 
                             hr_gint_part_[i].get_nnr() * sizeof(double), cudaMemcpyDeviceToHost));
    }
}

void Gint_vl_metagga_nspin4_gpu::cal_hr_gint_()
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
        CudaMemWrapper<double> dphi_x(BatchBigGrid::get_max_phi_len(), stream, false);
        CudaMemWrapper<double> dphi_y(BatchBigGrid::get_max_phi_len(), stream, false);
        CudaMemWrapper<double> dphi_z(BatchBigGrid::get_max_phi_len(), stream, false);
        CudaMemWrapper<double> dphi_x_vldr3(BatchBigGrid::get_max_phi_len(), stream, false);
        CudaMemWrapper<double> dphi_y_vldr3(BatchBigGrid::get_max_phi_len(), stream, false);
        CudaMemWrapper<double> dphi_z_vldr3(BatchBigGrid::get_max_phi_len(), stream, false);
        #pragma omp for schedule(dynamic)
        for(const auto& bgrid_batch: gint_info_->get_bgrid_batches())
        {
            if(bgrid_batch->empty())
            {
                continue;
            }
            phi_op.set_bgrid_batch(bgrid_batch);
            phi_op.set_phi_dphi(phi.get_device_ptr(),
                                dphi_x.get_device_ptr(), dphi_y.get_device_ptr(), dphi_z.get_device_ptr());
            for(int is = 0; is < nspin_; is++)
            {
                phi_op.phi_mul_vldr3(vr_eff_d_[is].get_device_ptr(), dr3_,
                                     phi.get_device_ptr(), phi_vldr3.get_device_ptr());
                phi_op.phi_mul_vldr3(vofk_d_[is].get_device_ptr(), dr3_,
                                     dphi_x.get_device_ptr(), dphi_x_vldr3.get_device_ptr());
                phi_op.phi_mul_vldr3(vofk_d_[is].get_device_ptr(), dr3_,
                                     dphi_y.get_device_ptr(), dphi_y_vldr3.get_device_ptr());
                phi_op.phi_mul_vldr3(vofk_d_[is].get_device_ptr(), dr3_,
                                     dphi_z.get_device_ptr(), dphi_z_vldr3.get_device_ptr());
                phi_op.phi_mul_phi(phi.get_device_ptr(), phi_vldr3.get_device_ptr(),
                                   hr_gint_part_[is], hr_gint_part_d_[is].get_device_ptr());
                phi_op.phi_mul_phi(dphi_x.get_device_ptr(), dphi_x_vldr3.get_device_ptr(),
                                   hr_gint_part_[is], hr_gint_part_d_[is].get_device_ptr());
                phi_op.phi_mul_phi(dphi_y.get_device_ptr(), dphi_y_vldr3.get_device_ptr(),
                                   hr_gint_part_[is], hr_gint_part_d_[is].get_device_ptr());
                phi_op.phi_mul_phi(dphi_z.get_device_ptr(), dphi_z_vldr3.get_device_ptr(),
                                   hr_gint_part_[is], hr_gint_part_d_[is].get_device_ptr());
            }
        }
        checkCuda(cudaStreamSynchronize(stream));
        checkCuda(cudaStreamDestroy(stream));
    }
    transfer_gpu_to_cpu_();
}

} // namespace ModuleGint