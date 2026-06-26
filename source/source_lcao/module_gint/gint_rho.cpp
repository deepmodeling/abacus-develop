#include "source_base/global_function.h"
#include "gint_rho.h"
#include "gint_common.h"
#include "phi_operator.h"

namespace ModuleGint
{

void Gint_rho::cal_gint()
{
    ModuleBase::TITLE("Gint", "cal_gint_rho");
    ModuleBase::timer::start("Gint", "cal_gint_rho");
    switch (gint_info_->get_exec_precision())
    {
    case GintPrecision::fp32:
        cal_gint_impl_<float>();
        break;
    case GintPrecision::fp64:
    default:
        cal_gint_impl_<double>();
        break;
    }
    ModuleBase::timer::end("Gint", "cal_gint_rho");
}

template<typename Real>
void Gint_rho::cal_gint_impl_()
{
    std::vector<HContainer<Real>> dm_gint_vec = init_dm_gint_<Real>();
    // rho_[is] is always double; phi_dot_phi accumulates into it directly.
    std::vector<double*> rho_data(nspin_);
    for (int is = 0; is < nspin_; ++is)
    {
        rho_data[is] = rho_[is];
    }
    dm_2d_to_gint(*gint_info_, dm_vec_, dm_gint_vec);
    cal_rho_(dm_gint_vec, rho_data);
}

template<typename Real>
std::vector<HContainer<Real>> Gint_rho::init_dm_gint_() const
{
    std::vector<HContainer<Real>> dm_gint_vec(nspin_);
    for (int is = 0; is < nspin_; is++)
    {
        dm_gint_vec[is] = gint_info_->get_hr<Real>();
    }
    return dm_gint_vec;
}

template<typename Real>
void Gint_rho::cal_rho_(
    const std::vector<HContainer<Real>>& dm_gint_vec,
    const std::vector<double*>& rho_data) const
{
#pragma omp parallel
    {
        PhiOperator phi_op;
        std::vector<Real> phi;
        // phi_dm is always double: phi_mul_dm writes the cast-to-double result
        // into it, and phi_dot_phi reads it as fp64 (so the rho reduction's
        // right-hand side is uniformly fp64 even on the fp32 path).
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
            phi_dm.resize(phi_len);
            phi_op.set_phi(phi.data());
            for (int is = 0; is < nspin_; is++)
            {
                phi_op.phi_mul_dm(phi.data(), dm_gint_vec[is], is_dm_symm_, phi_dm.data());
                phi_op.phi_dot_phi(phi.data(), phi_dm.data(), rho_data[is]);
            }
        }
    }
}

}
