#include "meta_pw.h"

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "src_pw/global.h"
#include "module_xc/xc_functional.h"
#include "module_base/tool_quit.h"

namespace hamilt
{

MetaPW::MetaPW(
    int max_npw_in,
    int npol_in,
    double tpiba_in,
    const int* ngk_in, 
    const int* isk_in,
    const ModuleBase::matrix* vk_in,
    ModulePW::PW_Basis_K* wfcpw_in
)
{
    this->ngk = ngk_in;
    this->isk = isk_in;
    this->max_npw = max_npw_in;
    this->tpiba = tpiba_in;
    this->vk = vk_in;
    this->wfcpw = wfcpw_in;
    this->npol = npol_in;
    if(this->ngk == nullptr || this->isk == nullptr || this->max_npw == 0 
    || this->tpiba < 1e-10 || this->vk == nullptr || this->wfcpw == nullptr
    || this->npol == 0)
    {
        ModuleBase::WARNING_QUIT("MetaPW", "Constuctor of Operator::MetaPW is failed, please check your code!");
    }
}

void MetaPW::act(const std::complex<double> *psi_in, std::complex<double> *hpsi, const size_t size) const
{
    if (XC_Functional::get_func_type() != 3)
    {
        return;
    }

    ModuleBase::timer::tick("Operator", "MetaPW");
    int m = int(size / this->max_npw / this->npol);
    if (int(size - m * this->max_npw * this->npol) != 0)
    {
        m++;
    }
    const int npw = this->ngk[this->ik];
    const int current_spin = this->isk[this->ik];

    std::complex<double> *tmhpsi;
    const std::complex<double> *tmpsi_in;

    tmhpsi = hpsi;
    tmpsi_in = psi_in;
    std::complex<double> *porter = new std::complex<double>[wfcpw->nmaxgr];
    for (int ib = 0; ib < m; ++ib)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int ig = 0; ig < npw; ig++)
            {
                double fact = wfcpw->getgpluskcar(ik, ig)[j] * this->tpiba;
                porter[ig] = tmpsi_in[ig] * complex<double>(0.0, fact);
            }

            wfcpw->recip2real(porter, porter, ik);

            const double* pvk = &(this->vk[0](current_spin, 0));
            for (int ir = 0; ir < this->vk->nc; ir++)
            {
                porter[ir] *= pvk[ir];
            }
            wfcpw->real2recip(porter, porter, ik);

            for (int ig = 0; ig < npw; ig++)
            {
                double fact = wfcpw->getgpluskcar(ik, ig)[j] * this->tpiba;
                tmhpsi[ig] -= complex<double>(0.0, fact) * porter[ig];
            }
        } // x,y,z directions
    }
    delete[] porter;
    ModuleBase::timer::tick("Operator", "MetaPW");
}

} // namespace hamilt