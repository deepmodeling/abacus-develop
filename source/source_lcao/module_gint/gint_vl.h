#pragma once

#include <memory>
#include <vector>
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "gint.h"
#include "gint_info.h"
#include "gint_precision.h"

namespace ModuleGint
{

class Gint_vl : public Gint
{
    public:
    Gint_vl(
        const double* vr_eff,
        HContainer<double>* hR,
        const GintExecConfig& cfg = {})
        : vr_eff_(vr_eff), hR_(hR), cfg_(cfg), dr3_(gint_info_->get_mgrid_volume()) {}
    
    void cal_gint();

    private:

    template<typename Real>
    void cal_gint_impl_();

    template<typename Real>
    HContainer<Real> init_hr_gint_() const;

    template<typename Real>
    const Real* get_vr_eff_data_(std::vector<Real>& vr_eff_buffer) const;

    // input
    const double* vr_eff_ = nullptr;

    // output
    HContainer<double>* hR_;

    GintExecConfig cfg_;

    // Intermediate variables
    double dr3_;
};

}
