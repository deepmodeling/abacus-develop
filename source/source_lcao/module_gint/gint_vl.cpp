#include <algorithm>

#include "gint_common.h"
#include "gint_vl.h"
#include "phi_operator.h"
#include "gint_helper.h"

namespace ModuleGint
{

void Gint_vl::cal_gint()
{
    ModuleBase::TITLE("Gint", "cal_gint_vl");
    ModuleBase::timer::start("Gint", "cal_gint_vl");
    switch (gint_info_->get_exec_precision())
    {
    case GintPrecision::fp32:
        cal_gint_impl_<float>();
        break;
    case GintPrecision::fp64:
    default:
        cal_gint_impl_<double>();
        break;
    }
    ModuleBase::timer::end("Gint", "cal_gint_vl");
}

//========================
// Private functions
//========================

// Overloaded helpers (C++11-compatible alternative to if constexpr).
// The double overload is preferred by overload resolution when Real=double;
// the template overload handles all other types (e.g. float).

inline const double* get_vr_eff_data_(const double* vr_eff, int /*local_mgrid_num*/,
                                      std::vector<double>& /*vr_eff_buffer*/)
{
    return vr_eff;
}

template<typename Real>
const Real* get_vr_eff_data_(const double* vr_eff, int local_mgrid_num,
                             std::vector<Real>& vr_eff_buffer)
{
    vr_eff_buffer.resize(local_mgrid_num);
    std::transform(vr_eff, vr_eff + local_mgrid_num, vr_eff_buffer.begin(), [](const double value) {
        return static_cast<Real>(value);
    });
    return vr_eff_buffer.data();
}

template<typename Real>
void Gint_vl::cal_gint_impl_()
{
    // hr_gint is always allocated as HContainer<double>: when Real=float, the
    // fp32 multiplies feed into a fp64 accumulator inside phi_mul_phi to avoid
    // catastrophic precision loss in the global reduction.
    HContainer<double> hr_gint = gint_info_->get_hr<double>();
    std::vector<Real> vr_eff_buffer;
    const Real* vr_eff = get_vr_eff_data_(
        vr_eff_, gint_info_->get_local_mgrid_num(), vr_eff_buffer);

#pragma omp parallel
    {
        PhiOperator phi_op;
        std::vector<Real> phi;
        std::vector<Real> phi_vldr3;
#pragma omp for schedule(dynamic)
        for (int i = 0; i < gint_info_->get_bgrids_num(); i++)
        {
            const auto& biggrid = gint_info_->get_biggrids()[i];
            if (biggrid->get_atoms().empty())
            {
                continue;
            }
            phi_op.set_bgrid(biggrid);
            const int phi_len = phi_op.get_rows() * phi_op.get_cols();
            phi.resize(phi_len);
            phi_vldr3.resize(phi_len);
            phi_op.set_phi(phi.data());
            phi_op.phi_mul_vldr3(vr_eff, static_cast<Real>(dr3_), phi.data(), phi_vldr3.data());
            phi_op.phi_mul_phi(phi.data(), phi_vldr3.data(), hr_gint, PhiOperator::TriPart::Upper);
        }
    }

    compose_hr_gint(hr_gint);
    hr_gint_to_hR(hr_gint, *hR_);
}

} // namespace ModuleGint
