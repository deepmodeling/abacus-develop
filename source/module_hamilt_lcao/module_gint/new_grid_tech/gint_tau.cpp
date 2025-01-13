#include "module_base/global_function.h"
#include "gint_tau.h"
#include "gint_common.h"
#include "phi_operator.h"

namespace ModuleGint
{

void Gint_tau::cal_gint()
{
    init_DMRGint_();
    transfer_DM_to_DMGint(gint_info_, DMR_vec_, DMRGint_vec_);
    cal_tau_();
}

void Gint_tau::init_DMRGint_()
{
    DMRGint_vec_.resize(nspin_);
    for (int is = 0; is < nspin_; is++)
    {
        DMRGint_vec_[is] = gint_info_->get_hr<double>();
    }
}

void Gint_tau::cal_tau_()
{
#pragma omp parallel
    {
        PhiOperator phi_op;
        std::vector<double> dphi_x;
        std::vector<double> dphi_y;
        std::vector<double> dphi_z;
        std::vector<double> dphi_x_DM;
        std::vector<double> dphi_y_DM;
        std::vector<double> dphi_z_DM;
#pragma omp for schedule(dynamic)
        for(const auto& biggrid: gint_info_->get_biggrids())
        {
            if(biggrid->get_atoms().size() == 0)
            {
                continue;
            }
            phi_op.set_bgrid(biggrid);
            const int phi_len = phi_op.get_rows() * phi_op.get_cols();
            dphi_x.resize(phi_len);
            dphi_y.resize(phi_len);
            dphi_z.resize(phi_len);
            dphi_x_DM.resize(phi_len);
            dphi_y_DM.resize(phi_len);
            dphi_z_DM.resize(phi_len);
            phi_op.set_phi_dphi(nullptr, dphi_x.data(), dphi_y.data(), dphi_z.data());
            for (int is = 0; is < nspin_; is++)
            {
                phi_op.phi_mul_dm(dphi_x.data(), *DMRGint_vec_[is], true, dphi_x_DM.data());
                phi_op.phi_mul_dm(dphi_y.data(), *DMRGint_vec_[is], true, dphi_y_DM.data());
                phi_op.phi_mul_dm(dphi_z.data(), *DMRGint_vec_[is], true, dphi_z_DM.data());
                phi_op.phi_dot_phi_dm(dphi_x.data(), dphi_x_DM.data(), kin_[is]);
                phi_op.phi_dot_phi_dm(dphi_y.data(), dphi_y_DM.data(), kin_[is]);
                phi_op.phi_dot_phi_dm(dphi_z.data(), dphi_z_DM.data(), kin_[is]);
            }
        }
    }
}

}