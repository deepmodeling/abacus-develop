#include "module_base/global_function.h"
#include "gint_rho.h"
#include "gint_common.h"
#include "phi_operator.h"

namespace ModuleGint
{

void Gint_rho::cal_gint()
{
    init_dm_gint_();
    transfer_dm_2d_to_gint(gint_info_, dm_vec_, dm_gint_vec_);
    cal_rho_();
}

void Gint_rho::init_dm_gint_()
{
    dm_gint_vec_.resize(nspin_);
    for (int is = 0; is < nspin_; is++)
    {
        dm_gint_vec_[is] = gint_info_->get_hr<double>();
    }
}

void Gint_rho::cal_rho_()
{
#pragma omp parallel
    {
        PhiOperator phi_op;
        std::vector<double> phi;
        std::vector<double> phi_dm;
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
            phi_dm.resize(phi_len);
            phi_op.set_phi(phi.data());
            for (int is = 0; is < nspin_; is++)
            {
                phi_op.phi_mul_dm(phi.data(), *dm_gint_vec_[is], true, phi_dm.data());
                phi_op.phi_dot_phi_dm(phi.data(), phi_dm.data(), rho_[is]);
            }
        }
    }
}


}