#include "module_base/global_function.h"
#include "gint_fvl_meta.h"
#include "gint_common.h"
#include "phi_operator.h"

namespace ModuleGint
{

void Gint_fvl_meta::cal_gint()
{
    init_dm_gint_();
    transfer_dm_2d_to_gint(gint_info_, dm_vec_, dm_gint_vec_);
    cal_fvl_svl_();
}

void Gint_fvl_meta::init_dm_gint_()
{
    dm_gint_vec_.resize(nspin_);
    for (int is = 0; is < nspin_; is++)
    {
        dm_gint_vec_[is] = gint_info_->get_hr<double>();
    }
}

void Gint_fvl_meta::cal_fvl_svl_()
{
#pragma omp parallel
    {
        PhiOperator phi_op;
        std::vector<double> phi;
        std::vector<double> phi_vldr3;
        std::vector<double> phi_vldr3_dm;
        std::vector<double> dphi_x;
        std::vector<double> dphi_y;
        std::vector<double> dphi_z;
        std::vector<double> dphi_x_vldr3;
        std::vector<double> dphi_y_vldr3;
        std::vector<double> dphi_z_vldr3;
        std::vector<double> dphi_x_vldr3_dm;
        std::vector<double> dphi_y_vldr3_dm;
        std::vector<double> dphi_z_vldr3_dm;
        std::vector<double> ddphi_xx;
        std::vector<double> ddphi_xy;
        std::vector<double> ddphi_xz;
        std::vector<double> ddphi_yy;
        std::vector<double> ddphi_yz;
        std::vector<double> ddphi_zz;
        ModuleBase::matrix* fvl_thread = nullptr;
        ModuleBase::matrix* svl_thread = nullptr;
        if(isforce_)
        {
            fvl_thread = new ModuleBase::matrix(*fvl_);
            fvl_thread->zero_out();
        }
        if(isstress_)
        {
            svl_thread = new ModuleBase::matrix(*svl_);
            svl_thread->zero_out();
        }
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
            phi_vldr3_dm.resize(phi_len);
            dphi_x.resize(phi_len);
            dphi_y.resize(phi_len);
            dphi_z.resize(phi_len);
            dphi_x_vldr3.resize(phi_len);
            dphi_y_vldr3.resize(phi_len);
            dphi_z_vldr3.resize(phi_len);
            dphi_x_vldr3_dm.resize(phi_len);
            dphi_y_vldr3_dm.resize(phi_len);
            dphi_z_vldr3_dm.resize(phi_len);
            ddphi_xx.resize(phi_len);
            ddphi_xy.resize(phi_len);
            ddphi_xz.resize(phi_len);
            ddphi_yy.resize(phi_len);
            ddphi_yz.resize(phi_len);
            ddphi_zz.resize(phi_len);
            phi_op.set_phi_dphi(phi.data(), dphi_x.data(), dphi_y.data(), dphi_z.data());
            phi_op.set_ddphi(ddphi_xx.data(), ddphi_xy.data(), ddphi_xz.data(),
                             ddphi_yy.data(), ddphi_yz.data(), ddphi_zz.data());
            for (int is = 0; is < nspin_; is++)
            {
                phi_op.phi_mul_vldr3(vr_eff_[is], dr3_, phi.data(), phi_vldr3.data());
                phi_op.phi_mul_vldr3(vofk_[is], dr3_, dphi_x.data(), dphi_x_vldr3.data());
                phi_op.phi_mul_vldr3(vofk_[is], dr3_, dphi_y.data(), dphi_y_vldr3.data());
                phi_op.phi_mul_vldr3(vofk_[is], dr3_, dphi_z.data(), dphi_z_vldr3.data());
                phi_op.phi_mul_dm(phi_vldr3.data(), *dm_gint_vec_[is], false, phi_vldr3_dm.data());
                phi_op.phi_mul_dm(dphi_x_vldr3.data(), *dm_gint_vec_[is], false, dphi_x_vldr3_dm.data());
                phi_op.phi_mul_dm(dphi_y_vldr3.data(), *dm_gint_vec_[is], false, dphi_y_vldr3_dm.data());
                phi_op.phi_mul_dm(dphi_z_vldr3.data(), *dm_gint_vec_[is], false, dphi_z_vldr3_dm.data());
                if(isforce_)
                {
                    phi_op.phi_dot_dphi(phi_vldr3_dm.data(), dphi_x.data(), dphi_y.data(), dphi_z.data(), fvl_thread);
                    phi_op.phi_dot_dphi(dphi_x_vldr3_dm.data(), ddphi_xx.data(), ddphi_xy.data(), ddphi_xz.data(), fvl_thread);
                    phi_op.phi_dot_dphi(dphi_y_vldr3_dm.data(), ddphi_xy.data(), ddphi_yy.data(), ddphi_yz.data(), fvl_thread);
                    phi_op.phi_dot_dphi(dphi_z_vldr3_dm.data(), ddphi_xz.data(), ddphi_yz.data(), ddphi_zz.data(), fvl_thread);
                }
                if(isstress_)
                {
                    phi_op.phi_dot_dphi_r(phi_vldr3_dm.data(), dphi_x.data(), dphi_y.data(), dphi_z.data(), svl_thread);
                    phi_op.phi_dot_dphi_r(dphi_x_vldr3_dm.data(), ddphi_xx.data(), ddphi_xy.data(), ddphi_xz.data(), svl_thread);
                    phi_op.phi_dot_dphi_r(dphi_y_vldr3_dm.data(), ddphi_xy.data(), ddphi_yy.data(), ddphi_yz.data(), svl_thread);
                    phi_op.phi_dot_dphi_r(dphi_z_vldr3_dm.data(), ddphi_xz.data(), ddphi_yz.data(), ddphi_zz.data(), svl_thread);
                }
            }
        }
#pragma omp critical
        {
            if(isforce_)
            {
                fvl_[0] += fvl_thread[0];
                delete fvl_thread;
            }
            if(isstress_)
            {
                svl_[0] += svl_thread[0];
                delete svl_thread;
            }
        }
    }
}

} // namespace ModuleGint