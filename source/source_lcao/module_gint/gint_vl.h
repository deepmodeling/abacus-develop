#pragma once

#include <memory>
#include <utility>
#include <vector>
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "gint.h"
#include "gint_info.h"

namespace ModuleGint
{

class Gint_vl : public Gint
{
    public:
    Gint_vl(
        const double* vr_eff,
        HContainer<double>* hR)
        : vr_eff_(vr_eff), hR_(hR), dr3_(gint_info_ == nullptr ? 0.0 : gint_info_->get_mgrid_volume()) {}
    Gint_vl(Gint_vl&& other) noexcept
        : Gint(std::move(other)),
          vr_eff_(other.vr_eff_),
          hR_(other.hR_),
          dr3_(other.dr3_)
    {
        other.vr_eff_ = nullptr;
        other.hR_ = nullptr;
        other.dr3_ = 0.0;
    }
    
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
