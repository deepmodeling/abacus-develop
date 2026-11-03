#include "gint_dvlocal.h"
#include "phi_operator.h"

namespace ModuleGint
{

void Gint_dvlocal::cal_dvlocal()
{
    ModuleBase::TITLE("Gint", "cal_gint_dvlocal");
    ModuleBase::timer::start("Gint", "cal_gint_dvlocal");
    init_hr_gint_();
    cal_hr_gint_();
    ModuleBase::timer::end("Gint", "cal_gint_dvlocal");
}

void Gint_dvlocal::init_hr_gint_()
{
    pvdpRx = gint_info_->get_hr<double>();
    pvdpRy = gint_info_->get_hr<double>();
    pvdpRz = gint_info_->get_hr<double>();
}

void Gint_dvlocal::cal_hr_gint_()
{
#pragma omp parallel
    {
        PhiOperator phi_op;
        std::vector<double> phi;
        std::vector<double> phi_vldr3;
        std::vector<double> dphi_x;
        std::vector<double> dphi_y;
        std::vector<double> dphi_z;
#pragma omp for schedule(dynamic)
        for (int i = 0; i < gint_info_->get_bgrids_num(); i++)
        {
            const auto& biggrid = gint_info_->get_biggrids()[i];
            if(biggrid->get_atoms().size() == 0)
            {
                continue;
            }
            phi_op.set_bgrid(biggrid);
            const int phi_len = phi_op.get_rows() * phi_op.get_cols();
            phi.resize(phi_len);
            phi_vldr3.resize(phi_len);
            dphi_x.resize(phi_len);
            dphi_y.resize(phi_len);
            dphi_z.resize(phi_len);
            phi_op.set_phi_dphi(phi.data(), dphi_x.data(), dphi_y.data(), dphi_z.data());
            phi_op.phi_mul_vldr3(vr_eff_, dr3_, phi.data(), phi_vldr3.data());
            const PhiOperator::TriPart tri
                = full_triangle_ ? PhiOperator::TriPart::Full : PhiOperator::TriPart::Upper;
            phi_op.phi_mul_phi(phi_vldr3.data(), dphi_x.data(), pvdpRx, tri);
            phi_op.phi_mul_phi(phi_vldr3.data(), dphi_y.data(), pvdpRy, tri);
            phi_op.phi_mul_phi(phi_vldr3.data(), dphi_z.data(), pvdpRz, tri);
        }
    }
}

}
