#include "pot_grad_xc.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_lr/potentials/xc_kernel.h"
#include "source_base/timer.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include <set>
namespace LR
{
    // constructor for exchange-correlation kernel
    PotGradXCLR::PotGradXCLR(const KernelXC& xc_kernel, const ModulePW::PW_Basis& rho_basis, const UnitCell& ucell, const int& nrxx)
        :xc_kernel_components_(xc_kernel),
        PotLRBase(rho_basis, (PARAM.inp.nspin == 1 || (PARAM.inp.nspin == 4 && !PARAM.globalv.domag && !PARAM.globalv.domag_z) ? 1 : 2), nrxx, ucell.tpiba)
    {}

    void PotGradXCLR::cal_v_eff(double** rho, const UnitCell& ucell, ModuleBase::matrix& v_eff, const std::vector<int>& ispin_op) const
    {
        ModuleBase::TITLE("PotGradXCLR", "cal_v_eff");
        ModuleBase::timer::start("PotGradXCLR", "cal_v_eff");
        const int nspin = v_eff.nr;
        if (XC_Functional::get_func_type() == 1 || XC_Functional::get_func_type() == 2 || XC_Functional::get_func_type() == 4)//LDA or GGA or HYBGGA
            if (nspin == 1)// for LDA-spin0, just fxc*rho where fxc=v2rho2; for GGA, v2rho2 has been replaced by the true fxc
                for (int ir = 0;ir < nrxx_;++ir)
                    v_eff(0, ir) += this->xc_kernel_components_.v3rho3.at(ir) * rho[0][ir] * rho[0][ir];
            else  //remain for spin 4
                throw std::domain_error("nspin =" + std::to_string(nspin)
                    + " unfinished in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
        else
            throw std::domain_error("GlobalV::XC_Functional::get_func_type() =" + std::to_string(XC_Functional::get_func_type())
                + " unfinished in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));

        ModuleBase::timer::end("PotGradXCLR", "cal_v_eff");
    }

}