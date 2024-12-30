#include "module_base/array_pool.h"
#include "module_base/global_function.h"
#include "gint_fvl.h"
#include "gint_common.h"
#include "phi_operator.h"

namespace ModuleGint
{

void Gint_fvl::cal_gint()
{
    init_DMRGint_();
    transfer_DM_to_DMGint(gint_info_, DMR_vec_, DMRGint_vec_);
    cal_fvl_svl_();
}

void Gint_fvl::init_DMRGint_()
{
    DMRGint_vec_.resize(nspin_);
    for (int is = 0; is < nspin_; is++)
    {
        DMRGint_vec_[is] = gint_info_->get_hr<double>();
    }
}

void Gint_fvl::cal_fvl_svl_()
{
#pragma omp parallel
    {
        PhiOperator phi_op;
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
            ModuleBase::Array_Pool<double> phi(phi_op.get_rows(), phi_op.get_cols());
            ModuleBase::Array_Pool<double> phi_vldr3(phi_op.get_rows(), phi_op.get_cols());
            ModuleBase::Array_Pool<double> phi_vldr3_DM(phi_op.get_rows(), phi_op.get_cols());
            ModuleBase::Array_Pool<double> dphi_x(phi_op.get_rows(), phi_op.get_cols());
            ModuleBase::Array_Pool<double> dphi_y(phi_op.get_rows(), phi_op.get_cols());
            ModuleBase::Array_Pool<double> dphi_z(phi_op.get_rows(), phi_op.get_cols());
            phi_op.set_phi_dphi(phi.get_ptr_1D(), dphi_x.get_ptr_1D(), dphi_y.get_ptr_1D(), dphi_z.get_ptr_1D());
            for (int is = 0; is < nspin_; is++)
            {
                ModuleBase::zeros(phi_vldr3.get_ptr_1D(), phi_op.get_rows()*phi_op.get_cols());
                ModuleBase::zeros(phi_vldr3_DM.get_ptr_1D(), phi_op.get_rows()*phi_op.get_cols());
                phi_op.phi_mul_vldr3(vr_eff_[is], dr3_, phi.get_ptr_2D(), phi_vldr3.get_ptr_2D());
                phi_op.phi_mul_dm(phi_vldr3.get_ptr_2D(), *DMRGint_vec_[is], false, phi_vldr3_DM.get_ptr_2D());
                if(isforce_)
                {
                    phi_op.phi_dot_dphi(phi_vldr3_DM.get_ptr_2D(), dphi_x.get_ptr_2D(), dphi_y.get_ptr_2D(), dphi_z.get_ptr_2D(), fvl_thread);
                }
                if(isstress_)
                {
                    phi_op.phi_dot_dphi_r(phi_vldr3_DM.get_ptr_2D(), dphi_x.get_ptr_2D(), dphi_y.get_ptr_2D(), dphi_z.get_ptr_2D(), svl_thread);
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


}