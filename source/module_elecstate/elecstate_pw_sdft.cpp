#include "./elecstate_pw_sdft.h"

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_parameter/parameter.h"
namespace elecstate
{

template <typename T, typename Device>
void ElecStatePW_SDFT<T, Device>::psiToRho(const psi::Psi<T>& psi)
{
    ModuleBase::TITLE(this->classname, "psiToRho");
    ModuleBase::timer::tick(this->classname, "psiToRho");
    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(this->charge->rho[is], this->charge->nrxx);
        if (XC_Functional::get_func_type() == 3)
        {
            ModuleBase::GlobalFunc::ZEROS(this->charge->kin_r[is], this->charge->nrxx);
        }
    }

    if (GlobalV::MY_STOGROUP == 0)
    {
        this->calEBand();

        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(this->charge->rho[is], this->charge->nrxx);
        }

        for (int ik = 0; ik < psi.get_nk(); ++ik)
        {
            psi.fix_k(ik);
            this->updateRhoK(psi);
        }
        this->parallelK();
    }
    ModuleBase::timer::tick(this->classname, "psiToRho");
    return;
}

// template class ElecStatePW_SDFT<std::complex<float>, base_device::DEVICE_CPU>;
template class ElecStatePW_SDFT<std::complex<double>, base_device::DEVICE_CPU>;
} // namespace elecstate