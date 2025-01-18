#include "module_base/global_function.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_base/blas_connector.h"
#include "gint_common.h"
#include "gint_vl_nspin4.h"
#include "phi_operator.h"
#include "gint_helper.h"

namespace ModuleGint
{
void Gint_vl_nspin4::cal_gint()
{
    init_hRGint_();
    cal_hRGint_();
    compose_hRGint(hRGint_part_, hRGint_full_);
    transfer_hRGint_to_hR(toConstSharedPtr(hRGint_full_), hR_);
}

void Gint_vl_nspin4::init_hRGint_()
{
    hRGint_part_.resize(nspin_);
    for(int i = 0; i < nspin_; i++)
    {
        hRGint_part_[i] = gint_info_->get_hr<double>();
    }
    const int npol = 2;
    hRGint_full_ = gint_info_->get_hr<std::complex<double>>(npol);
}

void Gint_vl_nspin4::cal_hRGint_()
{
#pragma omp parallel
    {
        PhiOperator phi_op;
        std::vector<double> phi;
        std::vector<double> phi_vldr3;
        std::vector<HContainer<double>> hRGint_part_thread(nspin_, *hRGint_part_[0]);
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
            for(int is = 0; is < nspin_; is++)
            {
                phi_op.phi_mul_vldr3(vr_eff_[is], dr3_, phi.data(), phi_vldr3.data());
                phi_op.phi_mul_phi_vldr3(phi.data(), phi_vldr3.data(), &hRGint_part_thread[is]);
            }
        }
#pragma omp critical
        {
            for(int is = 0; is < nspin_; is++)
            {
                {
                    BlasConnector::axpy(hRGint_part_thread[is].get_nnr(), 1.0, hRGint_part_thread[is].get_wrapper(),
                                        1, hRGint_part_[is]->get_wrapper(), 1);
                }
            }
        }
    }
}

} // namespace ModuleGint