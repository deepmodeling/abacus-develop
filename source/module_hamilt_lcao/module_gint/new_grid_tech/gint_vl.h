#pragma once

#include <memory>
#include <vector>
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "gint_info.h"

namespace Gint
{

class Gint_vl
{
    public:
    Gint_vl(
        std::shared_ptr<const GintInfo>gint_info,
        double* vr_eff,
        hamilt::HContainer<double>* hR)
        : gint_info_(gint_info), vr_eff_(vr_eff), hR_(hR){};
    
    void cal_gint();

    private:

    void init_hRGint_();
    
    // note that only the upper triangle matrix of hR is calculated
    // that's why we need compose_hRGint_() to fill the lower triangle matrix.
    void cal_hRGint_();

    void compose_hRGint_();

    void transfer_hRGint_to_hR_();

    std::shared_ptr<const GintInfo> gint_info_;

    // input
    double* vr_eff_;

    // output
    hamilt::HContainer<double>* hR_;

    //========================
    // Intermediate variables
    //========================
    std::shared_ptr<hamilt::HContainer<double>> hRGint_;
};

}