#pragma once

#include <memory>
#include <vector>
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "gint.h"
#include "gint_info.h"

namespace ModuleGint
{

class Gint_vl_metagga : public Gint
{
    public:
    Gint_vl_metagga(
        const double* vr_eff,
        const double* vofk,
        HContainer<double>* hR)
        : vr_eff_(vr_eff), vofk_(vofk), hR_(hR), dr3_(gint_info_->get_mgrid_volume()){};
    
    void cal_gint() override;

    private:

    void init_hr_gint_();
    
    // note that only the upper triangle matrix of hR is calculated
    // that's why we need compose_hr_gint() to fill the lower triangle matrix.
    void cal_hr_gint_();

    // input
    const double* vr_eff_;
    const double* vofk_;

    // output
    HContainer<double>* hR_;

    //========================
    // Intermediate variables
    //========================
    double dr3_;

    std::shared_ptr<HContainer<double>> hr_gint_;

};

}