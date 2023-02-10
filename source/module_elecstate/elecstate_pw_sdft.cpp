#include "./elecstate_pw_sdft.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "module_base/global_function.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
namespace elecstate
{
    void ElecStatePW_SDFT::psiToRho(const psi::Psi<std::complex<double>>& psi)
    {
        ModuleBase::TITLE(this->classname, "psiToRho");
        ModuleBase::timer::tick(this->classname, "psiToRho");
        for(int is=0; is < GlobalV::NSPIN; is++)
		{
			ModuleBase::GlobalFunc::ZEROS(this->charge->rho[is], this->charge->nrxx);
            if (XC_Functional::get_func_type() == 3)
		    {
                ModuleBase::GlobalFunc::ZEROS(this->charge->kin_r[is], this->charge->nrxx);
            }
		}
        
        if(GlobalV::MY_STOGROUP == 0)
		{
            this->calEBand();

            for(int is=0; is<GlobalV::NSPIN; is++)
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
}