#include "source_base/global_function.h"
#include "gint_drho.h"
#include "gint_common.h"
#include "phi_operator.h"

namespace ModuleGint
{

void Gint_drho::cal_gint()
{
    ModuleBase::TITLE("Gint", "cal_gint_drho");
    ModuleBase::timer::start("Gint", "cal_gint_drho");
    std::vector<HContainer<double>> dm_gint_vec = init_dm_gint_();
    dm_2d_to_gint(*gint_info_, dm_vec_, dm_gint_vec);
    cal_drho_(dm_gint_vec);
    ModuleBase::timer::end("Gint", "cal_gint_drho");
}

std::vector<HContainer<double>> Gint_drho::init_dm_gint_() const
{
    std::vector<HContainer<double>> dm_gint_vec(nspin_);
    for (int is = 0; is < nspin_; is++)
    {
        dm_gint_vec[is] = gint_info_->get_hr<double>();
    }
    return dm_gint_vec;
}

void Gint_drho::cal_drho_(const std::vector<HContainer<double>>& dm_gint_vec) const
{
#pragma omp parallel
    {
        PhiOperator phi_op;
        std::vector<double> phi;
        std::vector<double> dphi_x;
        std::vector<double> dphi_y;
        std::vector<double> dphi_z;
        std::vector<double> phi_dm;
#pragma omp for schedule(dynamic)
        for (int i = 0; i < gint_info_->get_bgrids_num(); i++)
        {
            const auto& biggrid = gint_info_->get_biggrids()[i];
            if (biggrid->get_atoms().empty())
            {
                continue;
            }
            phi_op.set_bgrid(biggrid);
            const int phi_len = phi_op.get_rows() * phi_op.get_cols();
            phi.resize(phi_len);
            dphi_x.resize(phi_len);
            dphi_y.resize(phi_len);
            dphi_z.resize(phi_len);
            phi_dm.resize(phi_len);
            // phi and its gradient, exactly as in Gint_dvlocal
            phi_op.set_phi_dphi(phi.data(), dphi_x.data(), dphi_y.data(), dphi_z.data());
            for (int is = 0; is < nspin_; is++)
            {
                // contract the gradient orbital (the FIRST/row index of D) with D, then
                // dot with the value orbital phi (the second/column index):
                //   phi_dm[ir,L] = sum_K dphi^d[ir,K] D[K,L]
                //   drho^d[ir]  += sum_L phi[ir,L] phi_dm[ir,L]
                //               = sum_{K,L} D[K,L] dphi^d_K(ir) phi_L(ir)
                // is_symm must stay false here: the symmetric phi_mul_dm fast path folds
                // the contraction assuming phi_dot_phi reuses the SAME operand, which is
                // not the case once the value orbital (phi) differs from the gradient one.
                phi_op.phi_mul_dm(dphi_x.data(), dm_gint_vec[is], false, phi_dm.data());
                phi_op.phi_dot_phi(phi.data(), phi_dm.data(), drho_x_[is]);
                phi_op.phi_mul_dm(dphi_y.data(), dm_gint_vec[is], false, phi_dm.data());
                phi_op.phi_dot_phi(phi.data(), phi_dm.data(), drho_y_[is]);
                phi_op.phi_mul_dm(dphi_z.data(), dm_gint_vec[is], false, phi_dm.data());
                phi_op.phi_dot_phi(phi.data(), phi_dm.data(), drho_z_[is]);
            }
        }
    }
}

}
