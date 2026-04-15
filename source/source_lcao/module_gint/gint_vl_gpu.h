#pragma once

#include <memory>
#include <vector>
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "gint.h"
#include "gint_info.h"
#include "source_lcao/module_gint/kernel/cuda_mem_wrapper.h"

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
    template<typename Real>
    void cal_gint_impl_();

    // input
    const double* vr_eff_ = nullptr;

        
    // output
    HContainer<double>* hR_;

    // Intermediate variables
    double dr3_;
};

}