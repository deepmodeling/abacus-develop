#include "module_base/array_pool.h"
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
#pragma omp for schedule(dynamic)
        for(const auto& biggrid: gint_info_->get_biggrids())
        {
            if(biggrid->get_atoms().size() == 0)
            {
                continue;
            }
            phi_op.set_bgrid(biggrid);
            ModuleBase::Array_Pool<double> dphi_x(phi_op.get_rows(), phi_op.get_cols());
            ModuleBase::Array_Pool<double> dphi_y(phi_op.get_rows(), phi_op.get_cols());
            ModuleBase::Array_Pool<double> dphi_z(phi_op.get_rows(), phi_op.get_cols());
            ModuleBase::Array_Pool<double> dphi_x_DM(phi_op.get_rows(), phi_op.get_cols());
            ModuleBase::Array_Pool<double> dphi_y_DM(phi_op.get_rows(), phi_op.get_cols());
            ModuleBase::Array_Pool<double> dphi_z_DM(phi_op.get_rows(), phi_op.get_cols());
            phi_op.set_phi_dphi(nullptr, dphi_x.get_ptr_1D(), dphi_y.get_ptr_1D(), dphi_z.get_ptr_1D());
            for (int is = 0; is < nspin_; is++)
            {
                ModuleBase::zeros(dphi_x_DM.get_ptr_1D(), phi_op.get_rows()*phi_op.get_cols());
                ModuleBase::zeros(dphi_y_DM.get_ptr_1D(), phi_op.get_rows()*phi_op.get_cols());
                ModuleBase::zeros(dphi_z_DM.get_ptr_1D(), phi_op.get_rows()*phi_op.get_cols());
                phi_op.phi_mul_dm(dphi_x.get_ptr_1D(), *DMRGint_vec_[is], true, dphi_x_DM.get_ptr_1D());
                phi_op.phi_mul_dm(dphi_y.get_ptr_1D(), *DMRGint_vec_[is], true, dphi_y_DM.get_ptr_1D());
                phi_op.phi_mul_dm(dphi_z.get_ptr_1D(), *DMRGint_vec_[is], true, dphi_z_DM.get_ptr_1D());
                phi_op.phi_dot_phi_dm(dphi_x.get_ptr_1D(), dphi_x_DM.get_ptr_1D(), kin_[is]);
                phi_op.phi_dot_phi_dm(dphi_y.get_ptr_1D(), dphi_y_DM.get_ptr_1D(), kin_[is]);
                phi_op.phi_dot_phi_dm(dphi_z.get_ptr_1D(), dphi_z_DM.get_ptr_1D(), kin_[is]);
            }
        }
    }
}

}