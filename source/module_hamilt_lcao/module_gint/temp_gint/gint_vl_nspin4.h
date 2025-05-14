#pragma once

#include <memory>
#include <vector>
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "gint.h"
#include "gint_info.h"

namespace ModuleGint
{

class Gint_vl_nspin4 : public Gint
{
    public:
    Gint_vl_nspin4(
        std::vector<const double*> vr_eff,
        HContainer<std::complex<double>>* hR)
        : vr_eff_(vr_eff), hR_(hR), dr3_(gint_info_->get_mgrid_volume()){};
    
    void cal_gint() override;

    private:

    void init_hr_gint_();
    
    // note that only the upper triangle matrix of hR is calculated
    // that's why we need compose_hr_gint() to fill the lower triangle matrix.
    void cal_hr_gint_();

    // input
    std::vector<const double*> vr_eff_;

    // output
    HContainer<std::complex<double>>* hR_;

    //========================
    // Intermediate variables
    //========================
    const double dr3_;

    const int nspin_ = 4;

    std::vector<std::shared_ptr<HContainer<double>>> hr_gint_part_;
    std::shared_ptr<HContainer<std::complex<double>>> hr_gint_full_;

};

} // namespace ModuleGint