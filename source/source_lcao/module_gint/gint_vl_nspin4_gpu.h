#pragma once

#include <memory>
#include <vector>
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "gint.h"
#include "gint_info.h"
#include "source_lcao/module_gint/temp_gint/kernel/cuda_mem_wrapper.h"

namespace ModuleGint
{

class Gint_vl_nspin4_gpu : public Gint
{
    public:
    Gint_vl_nspin4_gpu(
        std::vector<const double*> vr_eff,
        HContainer<std::complex<double>>* hR)
        : vr_eff_(vr_eff), hR_(hR), dr3_(gint_info_->get_mgrid_volume()) {}
    
    void cal_gint();

    private:

    void init_hr_gint_();
    
    void transfer_cpu_to_gpu_();

    void transfer_gpu_to_cpu_();
    
    // note that only the upper triangle matrix of hR is calculated
    // that's why we need compose_hr_gint() to fill the lower triangle matrix.
    void cal_hr_gint_();

    // input
    std::vector<const double*> vr_eff_;

    // output
    HContainer<std::complex<double>>* hR_;

    // Intermediate variables
    const double dr3_;

    const int nspin_ = 4;

    std::vector<HContainer<double>> hr_gint_part_;
    
    std::vector<CudaMemWrapper<double>> vr_eff_d_;
    std::vector<CudaMemWrapper<double>> hr_gint_part_d_;
};

} // namespace ModuleGint