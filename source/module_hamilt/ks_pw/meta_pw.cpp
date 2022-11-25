#include "meta_pw.h"

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "src_pw/global.h"
#include "module_xc/xc_functional.h"
#include "module_base/tool_quit.h"

using hamilt::Meta;
using hamilt::OperatorPW;

template<typename FPTYPE, typename Device>
Meta<OperatorPW<FPTYPE, Device>>::Meta(
    FPTYPE tpiba_in, 
    const int* isk_in,
    const ModuleBase::matrix* vk_in,
    ModulePW::PW_Basis_K* wfcpw_in
)
{
    this->classname = "Meta";
    this->cal_type = pw_meta;
    this->isk = isk_in;
    this->tpiba = tpiba_in;
    this->vk = vk_in;
    this->wfcpw = wfcpw_in;
    if(this->isk == nullptr || this->tpiba < 1e-10 || this->vk == nullptr || this->wfcpw == nullptr)
    {
        ModuleBase::WARNING_QUIT("MetaPW", "Constuctor of Operator::MetaPW is failed, please check your code!");
    }
}

template<typename FPTYPE, typename Device>
void Meta<OperatorPW<FPTYPE, Device>>::act(
    const psi::Psi<std::complex<FPTYPE>, Device> *psi_in, 
    const int n_npwx, 
    const std::complex<FPTYPE>* tmpsi_in, 
    std::complex<FPTYPE>* tmhpsi
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

    std::complex<FPTYPE> *porter = new std::complex<FPTYPE>[wfcpw->nmaxgr];
    for (int ib = 0; ib < n_npwx; ++ib)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int ig = 0; ig < npw; ig++)
            {
                FPTYPE fact = wfcpw->getgpluskcar(this->ik, ig)[j] * this->tpiba;
                porter[ig] = tmpsi_in[ig] * complex<FPTYPE>(0.0, fact);
            }

            wfcpw->recip2real(porter, porter, this->ik);
            if(this->vk->nc != 0)
            {
                const FPTYPE* pvk = &(this->vk[0](current_spin, 0));
                for (int ir = 0; ir < this->vk->nc; ir++)
                {
                    porter[ir] *= pvk[ir];
                }
            }
            wfcpw->real2recip(porter, porter, this->ik);

            for (int ig = 0; ig < npw; ig++)
            {
                FPTYPE fact = wfcpw->getgpluskcar(this->ik, ig)[j] * this->tpiba;
                tmhpsi[ig] -= complex<FPTYPE>(0.0, fact) * porter[ig];
            }
        } // x,y,z directions
        tmhpsi += this->max_npw;
        tmpsi_in += this->max_npw;
    }
    delete[] porter;
    ModuleBase::timer::tick("Operator", "MetaPW");
}

template<typename FPTYPE, typename Device>
template<typename T_in, typename Device_in>
hamilt::Meta<OperatorPW<FPTYPE, Device>>::Meta(const Meta<OperatorPW<T_in, Device_in>> *meta) {
    this->classname = "Meta";
    this->cal_type = pw_meta;
    this->ik = meta->get_ik();
    this->isk = meta->get_isk();
    this->tpiba = meta->get_tpiba();
    this->vk = meta->get_vk();
    this->wfcpw = meta->get_wfcpw();
    if(this->isk == nullptr || this->tpiba < 1e-10 || this->vk == nullptr || this->wfcpw == nullptr)
    {
        ModuleBase::WARNING_QUIT("MetaPW", "Constuctor of Operator::MetaPW is failed, please check your code!");
    }
}

namespace hamilt{
template class Meta<OperatorPW<double, psi::DEVICE_CPU>>;
template Meta<OperatorPW<double, psi::DEVICE_CPU>>::Meta(const Meta<OperatorPW<double, psi::DEVICE_CPU>> *meta);
#if ((defined __CUDA) || (defined __ROCM))
template class Meta<OperatorPW<double, psi::DEVICE_GPU>>;
template Meta<OperatorPW<double, psi::DEVICE_CPU>>::Meta(const Meta<OperatorPW<double, psi::DEVICE_GPU>> *meta);
template Meta<OperatorPW<double, psi::DEVICE_GPU>>::Meta(const Meta<OperatorPW<double, psi::DEVICE_CPU>> *meta);
template Meta<OperatorPW<double, psi::DEVICE_GPU>>::Meta(const Meta<OperatorPW<double, psi::DEVICE_GPU>> *meta);
#endif
} // namespace hamilt