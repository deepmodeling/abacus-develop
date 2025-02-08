#include "module_base/global_function.h"
#include "gint_tau.h"
#include "gint_common.h"
#include "phi_operator.h"

namespace ModuleGint
{

void Gint_tau::cal_gint()
{
    init_dm_gint_();
    transfer_dm_2d_to_gint(gint_info_, dm_vec_, dm_gint_vec_);
    cal_tau_();
}

void Gint_tau::init_dm_gint_()
{
    dm_gint_vec_.resize(nspin_);
    for (int is = 0; is < nspin_; is++)
    {
        dm_gint_vec_[is] = gint_info_->get_hr<double>();
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
        std::vector<double> dphi_x_dm;
        std::vector<double> dphi_y_dm;
        std::vector<double> dphi_z_dm;
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
            dphi_x_dm.resize(phi_len);
            dphi_y_dm.resize(phi_len);
            dphi_z_dm.resize(phi_len);
            phi_op.set_phi_dphi(nullptr, dphi_x.data(), dphi_y.data(), dphi_z.data());
            for (int is = 0; is < nspin_; is++)
            {
                phi_op.phi_mul_dm(dphi_x.data(), *dm_gint_vec_[is], true, dphi_x_dm.data());
                phi_op.phi_mul_dm(dphi_y.data(), *dm_gint_vec_[is], true, dphi_y_dm.data());
                phi_op.phi_mul_dm(dphi_z.data(), *dm_gint_vec_[is], true, dphi_z_dm.data());
                phi_op.phi_dot_phi_dm(dphi_x.data(), dphi_x_dm.data(), kin_[is]);
                phi_op.phi_dot_phi_dm(dphi_y.data(), dphi_y_dm.data(), kin_[is]);
                phi_op.phi_dot_phi_dm(dphi_z.data(), dphi_z_dm.data(), kin_[is]);
            }
        }
    }
}

}