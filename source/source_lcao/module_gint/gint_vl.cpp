#include <algorithm>
#include <type_traits>

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
    case GintRealPrecision::fp32:
        cal_gint_impl_<float>();
        break;
    case GintRealPrecision::fp64:
    default:
        cal_gint_impl_<double>();
        break;
    }
    ModuleBase::timer::end("Gint", "cal_gint_vl");
}

//========================
// Private functions
//========================

template<typename Real>
HContainer<Real> Gint_vl::init_hr_gint_() const
{
    return gint_info_->get_hr<Real>();
}

template<typename Real>
const Real* Gint_vl::get_vr_eff_data_(std::vector<Real>& vr_eff_buffer) const
{
    if constexpr (std::is_same_v<Real, double>)
    {
        return vr_eff_;
    }

    const int local_mgrid_num = gint_info_->get_local_mgrid_num();
    vr_eff_buffer.resize(local_mgrid_num);
    std::transform(vr_eff_, vr_eff_ + local_mgrid_num, vr_eff_buffer.begin(), [](const double value) {
        return static_cast<Real>(value);
    });
    return vr_eff_buffer.data();
}

template<typename Real>
void Gint_vl::cal_gint_impl_()
{
    HContainer<Real> hr_gint = init_hr_gint_<Real>();
    std::vector<Real> vr_eff_buffer;
    const Real* vr_eff = get_vr_eff_data_<Real>(vr_eff_buffer);

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
            phi_op.phi_mul_phi(phi.data(), phi_vldr3.data(), hr_gint, PhiOperator::Triangular_Matrix::Upper);
        }
    }

    if constexpr (std::is_same_v<Real, double>)
    {
        compose_hr_gint(hr_gint);
        transfer_hr_gint_to_hR(hr_gint, *hR_);
    }
    else
    {
        HContainer<double> hr_gint_dp = make_cast_hcontainer<double>(hr_gint);
        compose_hr_gint(hr_gint_dp);
        transfer_hr_gint_to_hR(hr_gint_dp, *hR_);
    }
}

} // namespace ModuleGint
