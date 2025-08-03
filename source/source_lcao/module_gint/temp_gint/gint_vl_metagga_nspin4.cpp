#include "source_base/global_function.h"
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "gint_common.h"
#include "gint_vl_metagga_nspin4.h"
#include "phi_operator.h"
#include "gint_helper.h"

namespace ModuleGint
{

void Gint_vl_metagga_nspin4::cal_gint()
{
    ModuleBase::TITLE("Gint", "cal_gint_vl");
    ModuleBase::timer::tick("Gint", "cal_gint_vl");
    init_hr_gint_();
    cal_hr_gint_();
    compose_hr_gint(hr_gint_part_, hr_gint_full_);
    transfer_hr_gint_to_hR(hr_gint_full_, *hR_);
    ModuleBase::timer::tick("Gint", "cal_gint_vl");
}

void Gint_vl_metagga_nspin4::init_hr_gint_()
{
    hr_gint_part_.resize(nspin_);
    for(int i = 0; i < nspin_; i++)
    {
        hr_gint_part_[i] = gint_info_->get_hr<double>();
    }
    const int npol = 2;
    hr_gint_full_ = gint_info_->get_hr<std::complex<double>>(npol);
}

void Gint_vl_metagga_nspin4::cal_hr_gint_()
{
#pragma omp parallel
    {
        PhiOperator phi_op;
        std::vector<double> phi;
        std::vector<double> phi_vldr3;
        std::vector<double> dphi_x;
        std::vector<double> dphi_y;
        std::vector<double> dphi_z;
        std::vector<double> dphi_x_vldr3;
        std::vector<double> dphi_y_vldr3;
        std::vector<double> dphi_z_vldr3;
#pragma omp for schedule(dynamic)
        for(const auto& biggrid: gint_info_->get_biggrids())
        {
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
            dphi_x_vldr3.resize(phi_len);
            dphi_y_vldr3.resize(phi_len);
            dphi_z_vldr3.resize(phi_len);
            phi_op.set_phi_dphi(phi.data(), dphi_x.data(), dphi_y.data(), dphi_z.data());
            for(int is = 0; is < nspin_; is++)
            {
                phi_op.phi_mul_vldr3(vr_eff_[is], dr3_, phi.data(), phi_vldr3.data());
                phi_op.phi_mul_vldr3(vofk_[is], dr3_, dphi_x.data(), dphi_x_vldr3.data());
                phi_op.phi_mul_vldr3(vofk_[is], dr3_, dphi_y.data(), dphi_y_vldr3.data());
                phi_op.phi_mul_vldr3(vofk_[is], dr3_, dphi_z.data(), dphi_z_vldr3.data());
                phi_op.phi_mul_phi(phi.data(), phi_vldr3.data(), hr_gint_part_[is], PhiOperator::Triangular_Matrix::Upper);
                phi_op.phi_mul_phi(dphi_x.data(), dphi_x_vldr3.data(), hr_gint_part_[is], PhiOperator::Triangular_Matrix::Upper);
                phi_op.phi_mul_phi(dphi_y.data(), dphi_y_vldr3.data(), hr_gint_part_[is], PhiOperator::Triangular_Matrix::Upper);
                phi_op.phi_mul_phi(dphi_z.data(), dphi_z_vldr3.data(), hr_gint_part_[is], PhiOperator::Triangular_Matrix::Upper);
            }
        }
    }
}

} // namespace ModuleGint