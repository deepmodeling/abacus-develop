#pragma once

#include <memory>
#include <vector>
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "gint.h"
#include "gint_info.h"

namespace ModuleGint
{

class Gint_vl_metagga_nspin4 : public Gint
{
    public:
    Gint_vl_metagga_nspin4(
        std::vector<const double*> vr_eff,
        std::vector<const double*> vofk,
        HContainer<std::complex<double>>* hR)
        : vr_eff_(vr_eff), vofk_(vofk), hR_(hR), dr3_(gint_info_->get_mgrid_volume()){};
    
    void cal_gint() override;

    private:
    void init_hr_gint_();

    void cal_hr_gint_();

    // input
    std::vector<const double*> vr_eff_;
    std::vector<const double*> vofk_;
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

}