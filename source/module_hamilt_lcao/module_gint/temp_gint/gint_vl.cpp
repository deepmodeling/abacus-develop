#include "module_base/blas_connector.h"
#include "gint_common.h"
#include "gint_vl.h"
#include "phi_operator.h"
#include "gint_helper.h"

namespace ModuleGint
{

void Gint_vl::cal_gint()
{
    init_hr_gint_();
    cal_hr_gint_();
    compose_hr_gint(hr_gint_);
    transfer_hr_gint_to_hR(toConstSharedPtr(hr_gint_), hR_);
}

//========================
// Private functions
//========================

void Gint_vl::init_hr_gint_()
{
    hr_gint_ = gint_info_->get_hr<double>();
}

void Gint_vl::cal_hr_gint_()
{
// be careful!!
// each thread will have a copy of hr_gint_, this may cause a lot of memory usage
#pragma omp parallel
    {
        PhiOperator phi_op;
        std::vector<double> phi;
        std::vector<double> phi_vldr3;
        HContainer<double> hr_gint_local(*hr_gint_);
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
            phi_op.set_phi(phi.data());
            phi_op.phi_mul_vldr3(vr_eff_, dr3_, phi.data(), phi_vldr3.data());
            phi_op.phi_mul_phi_vldr3(phi.data(), phi_vldr3.data(), &hr_gint_local);
        }
#pragma omp critical
        {
            BlasConnector::axpy(hr_gint_local.get_nnr(), 1.0, hr_gint_local.get_wrapper(),
                                1, hr_gint_->get_wrapper(), 1);
        }
    }
}

}