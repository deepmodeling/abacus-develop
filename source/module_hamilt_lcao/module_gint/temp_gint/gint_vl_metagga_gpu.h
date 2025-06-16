#pragma once

#include <memory>
#include <vector>
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "gint.h"
#include "gint_info.h"
#include "module_hamilt_lcao/module_gint/temp_gint/kernel/cuda_mem_wrapper.h"

namespace ModuleGint
{

class Gint_vl_metagga_gpu : public Gint
{
    public:
    Gint_vl_metagga_gpu(
        const double* vr_eff,
        const double* vofk,
        HContainer<double>* hR)
        : vr_eff_(vr_eff), vofk_(vofk), hR_(hR), dr3_(gint_info_->get_mgrid_volume()) {}
    
    void cal_gint();

    private:

    void init_hr_gint_();

    void transfer_cpu_to_gpu_();

    void transfer_gpu_to_cpu_();
    
    // note that only the upper triangle matrix of hR is calculated
    // that's why we need compose_hr_gint() to fill the lower triangle matrix.
    void cal_hr_gint_();

    // input
    const double* vr_eff_;
    const double* vofk_;

    // output
    HContainer<double>* hR_;

    // Intermediate variables
    double dr3_;

    HContainer<double> hr_gint_;
    
    CudaMemWrapper<double> hr_gint_d_;
    CudaMemWrapper<double> vr_eff_d_;
    CudaMemWrapper<double> vofk_d_;
};

}