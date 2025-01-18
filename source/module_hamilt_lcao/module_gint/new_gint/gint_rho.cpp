#include "module_base/global_function.h"
#include "gint_rho.h"
#include "gint_common.h"
#include "phi_operator.h"

namespace ModuleGint
{

void Gint_rho::cal_gint()
{
    init_DMRGint_();
    transfer_DM_to_DMGint(gint_info_, DMR_vec_, DMRGint_vec_);
    cal_rho_();
}

void Gint_rho::init_DMRGint_()
{
    DMRGint_vec_.resize(nspin_);
    for (int is = 0; is < nspin_; is++)
    {
        DMRGint_vec_[is] = gint_info_->get_hr<double>();
    }
}

void Gint_rho::cal_rho_()
{
#pragma omp parallel
    {
        PhiOperator phi_op;
        std::vector<double> phi;
        std::vector<double> phi_DMR;
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
            phi_DMR.resize(phi_len);
            phi_op.set_phi(phi.data());
            for (int is = 0; is < nspin_; is++)
            {
                phi_op.phi_mul_dm(phi.data(), *DMRGint_vec_[is], true, phi_DMR.data());
                phi_op.phi_dot_phi_dm(phi.data(), phi_DMR.data(), rho_[is]);
            }
        }
    }
}


}