#pragma once

#include <memory>
#include <vector>
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "gint.h"
#include "gint_info.h"
#include "source_lcao/module_gint/temp_gint/kernel/cuda_mem_wrapper.h"

namespace ModuleGint
{

class Gint_vl_gpu : public Gint
{
    public:
    Gint_vl_gpu(
        const double* vr_eff,
        HContainer<double>* hR)
        : vr_eff_(vr_eff), hR_(hR), dr3_(gint_info_->get_mgrid_volume()) {}
    
    void cal_gint();

    private:

    void init_hr_gint_();

    void transfer_cpu_to_gpu_();

    void transfer_gpu_to_cpu_();

    void cal_hr_gint_();

    // input
    const double* vr_eff_;

        
    // output
    HContainer<double>* hR_;

    // Intermediate variables
    double dr3_;

    HContainer<double> hr_gint_;
    
    CudaMemWrapper<double> hr_gint_d_;
    CudaMemWrapper<double> vr_eff_d_;
};

}