#include "module_base/array_pool.h"
#include "module_base/timer.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
#include "module_base/blas_connector.h"
#include "gint_common.h"
#include "gint_vl.h"
#include "phi_operator.h"

namespace ModuleGint
{

void Gint_vl::cal_gint()
{
    ModuleBase::timer::tick("Gint_vl", "cal_gint");

    init_hRGint_();
    cal_hRGint_();
    compose_hRGint(hRGint_.get());
    transfer_hRGint_to_hR(hRGint_.get(), hR_);
    
    ModuleBase::timer::tick("Gint_vl", "cal_gint");
}

//========================
// Private functions
//========================

void Gint_vl::init_hRGint_()
{
    hRGint_ = gint_info_->get_hr<double>();
}

void Gint_vl::cal_hRGint_()
{
// be careful!!
// each thread will have a copy of hRGint_, this may cause a lot of memory usage
#pragma omp parallel
    {
        PhiOperator phi_op;
        hamilt::HContainer<double> hRGint_local(*hRGint_);
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
            phi_op.set_phi(phi.get_ptr_1D());
            phi_op.phi_mul_vldr3(vr_eff_, dr3_, phi.get_ptr_2D(), phi_vldr3.get_ptr_2D());
            phi_op.phi_mul_phi_vldr3(phi.get_ptr_2D(), phi_vldr3.get_ptr_2D(), &hRGint_local);
        }
#pragma omp critical
        {
            BlasConnector::axpy(hRGint_local.get_nnr(), 1.0, hRGint_local.get_wrapper(),
                                1, hRGint_->get_wrapper(), 1);
        }
    }
}

}