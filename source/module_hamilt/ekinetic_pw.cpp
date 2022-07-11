#include "ekinetic_pw.h"

#include "module_base/timer.h"
#include "module_base/tool_quit.h"

namespace hamilt
{
EkineticPW::EkineticPW(
    int max_npw_in,
    int npol_in,
    double tpiba2_in,
    const int* ngk_in, 
    const double* gk2_in
)
{
    this->max_npw = max_npw_in;
    this->npol = npol_in;
    this->tpiba2 = tpiba2_in;
    this->ngk = ngk_in;
    this->gk2 = gk2_in;
    if(this->max_npw == 0 || this->npol == 0 || this->tpiba2 < 1e-10  
    || this->ngk == nullptr || this->gk2 == nullptr
    )
    {
        ModuleBase::WARNING_QUIT("EkineticPW", "Constuctor of Operator::EkineticPW is failed, please check your code!");
    }
}

void EkineticPW::act(const std::complex<double> *psi_in, std::complex<double> *hpsi, const size_t size) const
{
    ModuleBase::timer::tick("Operator", "EkineticPW");
    int m = int(size / this->max_npw / this->npol);
    if (int(size - m * this->max_npw * this->npol) != 0)
        m++;
    const int npw = this->ngk[this->ik];
    const double *gk2_ik = &(this->gk2[this->ik * this->max_npw]);

    std::complex<double> *tmhpsi;
    const std::complex<double> *tmpsi_in;

    tmhpsi = hpsi;
    tmpsi_in = psi_in;
    for (int ib = 0; ib < m; ++ib)
    {
        for (int ig = 0; ig < npw; ++ig)
        {
            tmhpsi[ig] = gk2_ik[ig] * tpiba2 * tmpsi_in[ig];
        }
        if (this->npol == 2)
        {
            for (int ig = npw; ig < this->max_npw; ++ig)
            {
                tmhpsi[ig] = 0;
            }
            tmhpsi += this->max_npw;
            tmpsi_in += this->max_npw;
            for (int ig = 0; ig < npw; ++ig)
            {
                tmhpsi[ig] = gk2_ik[ig] * tpiba2 * tmpsi_in[ig];
            }
            for (int ig = npw; ig < this->max_npw; ++ig)
            {
                tmhpsi[ig] = 0;
            }
        }
        tmhpsi += this->max_npw;
        tmpsi_in += this->max_npw;
    }

    ModuleBase::timer::tick("Operator", "EkineticPW");
    return;
}

} // namespace hamilt