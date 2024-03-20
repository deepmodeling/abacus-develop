#include "H_Hartree_pw.h"
#include "efield.h"
#include "gatefield.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
#include "module_base/tool_title.h"
#include "pot_local.h"
#include "pot_surchem.hpp"
#include "pot_xc.h"
#include "potential_new.h"
#include "pot_local_paw.h"
#ifdef __LCAO
#include "H_TDDFT_pw.h"
#endif

namespace elecstate
{

PotBase* Potential::get_pot_type(const std::string& pot_type)
{
    ModuleBase::TITLE("Potential", "get_pot_type");
    if (pot_type == "local")
    {
        if(!GlobalV::use_paw)
        {
            return new PotLocal(this->vloc_, &(this->structure_factors_->strucFac), this->rho_basis_);
        }
        else
        {
            return new PotLocal_PAW();
        }
    }
    else if (pot_type == "hartree")
    {
        return new PotHartree(this->rho_basis_);
    }
    else if (pot_type == "xc")
    {
        return new PotXC(this->rho_basis_, this->etxc_, this->vtxc_, &(this->vofk_effective));
    }
    else if (pot_type == "surchem")
    {
        return new PotSurChem(this->rho_basis_,
                              this->structure_factors_,
                              this->v_effective_fixed.data(),
                              &GlobalC::solvent_model);
    }
    else if (pot_type == "efield")
    {
        return new PotEfield(this->rho_basis_, this->ucell_, GlobalV::DIP_COR_FLAG);
    }
    else if (pot_type == "gatefield")
    {
        return new PotGate(this->rho_basis_, this->ucell_);
    }
#ifdef __LCAO
    else if (pot_type == "tddft")
    {
        return new H_TDDFT_pw(this->rho_basis_, this->ucell_);
    }
#endif
    else
    {
        ModuleBase::WARNING_QUIT("Potential::get_pot_type", "Please input correct component of potential!");
        __builtin_unreachable();
    }
}

} // namespace elecstate