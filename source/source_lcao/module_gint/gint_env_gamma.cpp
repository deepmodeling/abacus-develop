#include "gint_env_gamma.h"
#include "gint_common.h"
#include "phi_operator.h"

namespace ModuleGint
{

Gint_env_gamma::Gint_env_gamma(
    const double* psid,
    const Parallel_Orbitals* pv,
    const int nbands,
    const int nlocal,
    double* rho)
    :rho_(rho)
{
    wfc_gint_.resize(nbands * gint_info_->get_lgd());
    wfc_2d_to_gint(psid, nbands, nlocal, *pv, wfc_gint_.data(), *gint_info_);
}

void Gint_env_gamma::cal_env_band(const int iband)
{
    ModuleBase::TITLE("Gint", "cal_gint_env");
    ModuleBase::timer::tick("Gint", "cal_gint_env");
    ModuleBase::GlobalFunc::ZEROS(rho_, gint_info_->get_local_mgrid_num());
    const double* wfc_gint_band = &wfc_gint_[iband * gint_info_->get_lgd()];
#pragma omp parallel
    {
        PhiOperator phi_op;
        std::vector<double> phi;
#pragma omp for schedule(dynamic)
        for(const auto& biggrid: gint_info_->get_biggrids())
        {
            if(biggrid->get_atoms().empty())
            {
                continue;
            }
            phi_op.set_bgrid(biggrid);
            const int phi_len = phi_op.get_rows() * phi_op.get_cols();
            phi.resize(phi_len);
            phi_op.set_phi(phi.data());
            phi_op.cal_env_gamma(phi.data(), wfc_gint_band, gint_info_->get_trace_lo(), rho_);
        }
    }
    ModuleBase::timer::tick("Gint", "cal_gint_env");
}


}