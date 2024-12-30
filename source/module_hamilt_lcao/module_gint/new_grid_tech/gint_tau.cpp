#include "module_base/array_pool.h"
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
        const int npol = (nspin_ == 4 ? 2 : 1);
        DMRGint_vec_[is] = gint_info_->get_hr<double>(npol);
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
            ModuleBase::Array_Pool<double> phi(phi_op.get_rows(), phi_op.get_cols());
            ModuleBase::Array_Pool<double> phi_DMR(phi_op.get_rows(), phi_op.get_cols());
            phi_op.set_phi(phi.get_ptr_1D());
            for (int is = 0; is < nspin_; is++)
            {
                phi_op.phi_mul_dm(phi.get_ptr_2D(), *DMRGint_vec_[is], true, phi_DMR.get_ptr_2D());
                phi_op.phi_dot_phi_dm(phi.get_ptr_2D(), phi_DMR.get_ptr_2D(), tau[is]);
            }
        }
    }
}

}