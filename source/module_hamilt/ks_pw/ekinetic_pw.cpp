#include "ekinetic_pw.h"

#include "module_base/timer.h"
#include "module_base/tool_quit.h"

namespace hamilt
{

Ekinetic<OperatorPW>::Ekinetic(
    double tpiba2_in,
    const double* gk2_in,
    const int gk2_dim_in
)
{
    this->cal_type = 10;
    this->tpiba2 = tpiba2_in;
    this->gk2 = gk2_in;
    this->gk2_dim = gk2_dim_in;
    if( this->tpiba2 < 1e-10 || this->gk2 == nullptr
    )
    {
        ModuleBase::WARNING_QUIT("EkineticPW", "Constuctor of Operator::EkineticPW is failed, please check your code!");
    }
}

void Ekinetic<OperatorPW>::act
(
    const psi::Psi<std::complex<double>> *psi_in, 
    const int n_npwx, 
    const std::complex<double>* tmpsi_in, 
    std::complex<double>* tmhpsi
)const  
{
    ModuleBase::timer::tick("Operator", "EkineticPW");
    const int npw = psi_in->get_ngk(this->ik);
    this->max_npw = psi_in->get_nbasis() / psi_in->npol;


    const double *gk2_ik = &(this->gk2[this->ik * this->gk2_dim]);

    for (int ib = 0; ib < n_npwx; ++ib)
    {
        for (int ig = 0; ig < npw; ++ig)
        {
            tmhpsi[ig] += gk2_ik[ig] * tpiba2 * tmpsi_in[ig];
        }
        tmhpsi += this->max_npw;
        tmpsi_in += this->max_npw;
    }

    ModuleBase::timer::tick("Operator", "EkineticPW");
    return;
}

} // namespace hamilt