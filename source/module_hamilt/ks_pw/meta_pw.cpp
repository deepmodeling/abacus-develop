#include "meta_pw.h"

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "src_pw/global.h"
#include "module_xc/xc_functional.h"
#include "module_base/tool_quit.h"

namespace hamilt
{

Meta<OperatorPW>::Meta(
    double tpiba_in, 
    const int* isk_in,
    const ModuleBase::matrix* vk_in,
    ModulePW::PW_Basis_K* wfcpw_in
)
{
    this->cal_type = 13;
    this->isk = isk_in;
    this->tpiba = tpiba_in;
    this->vk = vk_in;
    this->wfcpw = wfcpw_in;
    if(this->isk == nullptr || this->tpiba < 1e-10 || this->vk == nullptr || this->wfcpw == nullptr)
    {
        ModuleBase::WARNING_QUIT("MetaPW", "Constuctor of Operator::MetaPW is failed, please check your code!");
    }
}

void Meta<OperatorPW>::act
(
    const psi::Psi<std::complex<double>> *psi_in, 
    const int n_npwx, 
    const std::complex<double>* tmpsi_in, 
    std::complex<double>* tmhpsi
)const
{
    if (XC_Functional::get_func_type() != 3)
    {
        return;
    }

    ModuleBase::timer::tick("Operator", "MetaPW");

    const int npw = psi_in->get_ngk(this->ik);
    const int current_spin = this->isk[this->ik];
    this->max_npw = psi_in->get_nbasis() / psi_in->npol;
    //npol == 2 case has not been considered
    this->npol = psi_in->npol;

    std::complex<double> *porter = new std::complex<double>[wfcpw->nmaxgr];
    for (int ib = 0; ib < n_npwx; ++ib)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int ig = 0; ig < npw; ig++)
            {
                double fact = wfcpw->getgpluskcar(ik, ig)[j] * this->tpiba;
                porter[ig] = tmpsi_in[ig] * complex<double>(0.0, fact);
            }

            wfcpw->recip2real(porter, porter, ik);
            if(this->vk->nc != 0)
            {
                const double* pvk = &(this->vk[0](current_spin, 0));
                for (int ir = 0; ir < this->vk->nc; ir++)
                {
                    porter[ir] *= pvk[ir];
                }
            }
            wfcpw->real2recip(porter, porter, ik);

            for (int ig = 0; ig < npw; ig++)
            {
                double fact = wfcpw->getgpluskcar(ik, ig)[j] * this->tpiba;
                tmhpsi[ig] -= complex<double>(0.0, fact) * porter[ig];
            }
        } // x,y,z directions
        tmhpsi += this->max_npw;
        tmpsi_in += this->max_npw;
    }
    delete[] porter;
    ModuleBase::timer::tick("Operator", "MetaPW");
}

} // namespace hamilt