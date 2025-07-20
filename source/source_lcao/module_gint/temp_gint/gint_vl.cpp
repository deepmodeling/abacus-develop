#include "gint_common.h"
#include "gint_vl.h"
#include "phi_operator.h"
#include "gint_helper.h"

namespace ModuleGint
{

void Gint_vl::cal_gint()
{
    ModuleBase::TITLE("Gint", "cal_gint_vl");
    ModuleBase::timer::tick("Gint", "cal_gint_vl");
    init_hr_gint_();
    cal_hr_gint_();
    compose_hr_gint(hr_gint_);
    transfer_hr_gint_to_hR(hr_gint_, *hR_);
    ModuleBase::timer::tick("Gint", "cal_gint_vl");
}

//========================
// Private functions
//========================

void Gint_vl::init_hr_gint_()
{
    hr_gint_ = gint_info_->get_hr<double>();
}

void Gint_vl::cal_hr_gint_()
{
#pragma omp parallel
    {
        PhiOperator phi_op;
        std::vector<double> phi;
        std::vector<double> phi_vldr3;
#pragma omp for schedule(dynamic)
        for(const auto& biggrid: gint_info_->get_biggrids())
        {
            if(biggrid->get_atoms().empty())
            {
                continue;
            }
            phi_op.set_bgrid(biggrid);
            const int phi_len = phi_op.get_rows() * phi_op.get_cols();
            phi.resize(phi_len);
            phi_vldr3.resize(phi_len);
            phi_op.set_phi(phi.data());
            phi_op.phi_mul_vldr3(vr_eff_, dr3_, phi.data(), phi_vldr3.data());
            phi_op.phi_mul_phi(phi.data(), phi_vldr3.data(), hr_gint_, PhiOperator::Triangular_Matrix::Upper);
        }
    }
}

}