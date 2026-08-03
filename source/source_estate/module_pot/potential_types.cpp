#include "H_Hartree_pw.h"
#include "efield.h"
#include "gatefield.h"
#include "pot_local.h"
#include "pot_sep.h"
#include "pot_surchem.hpp"
#include "pot_xc.h"
#include "potential_new.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"
#include "source_io/module_parameter/parameter.h"
#include "H_TDDFT_pw.h"
#ifdef __MLALGO
#include "pot_ml_exx.h"
#endif

namespace elecstate
{

PotBase* Potential::get_pot_type(const std::string& pot_type)
{
    ModuleBase::TITLE("Potential", "get_pot_type");
    if (pot_type == "local")
    {
        return new PotLocal(this->vloc_, &(this->structure_factors_->strucFac), this->rho_basis_, this->vl_of_0);
    }
    else if (pot_type == "hartree")
    {
        return new PotHartree(this->rho_basis_);
    }
    else if (pot_type == "xc")
    {
        return new PotXC(this->rho_basis_, this->etxc_, this->vtxc_, &(this->vofk_eff));
    }
    else if (pot_type == "surchem")
    {
        return new PotSurChem(this->rho_basis_, this->structure_factors_, this->v_eff_fixed.data(), this->solvent_);
    }
    else if (pot_type == "efield")
    {
        return new PotEfield(this->rho_basis_, this->ucell_, this->solvent_, PARAM.inp.dip_cor_flag);
    }
    else if (pot_type == "gatefield")
    {
        return new PotGate(this->rho_basis_, this->ucell_);
    }
    else if (pot_type == "tddft")
    {
        // The RT-TDDFT ESolver injects this manager before potential
        // registration so the potential does not own a separate time counter.
        return new H_TDDFT_pw(this->rho_basis_, this->ucell_, this->td_field_manager_);
    }
#ifdef __MLALGO
    else if (pot_type == "ml_exx")
    {
        return new PotML_EXX(this->rho_basis_, this->ucell_);
    }
#endif
    else if (pot_type == "dfthalf")
    {
        return new PotSep(&(this->structure_factors_->strucFac), this->rho_basis_, this->vsep_cell);
    }
    else
    {
        ModuleBase::WARNING_QUIT("Potential::get_pot_type", "Please input correct component of potential!");
        __builtin_unreachable();
    }
}

} // namespace elecstate
