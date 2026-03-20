#include "potential_factory.h"
#include "H_Hartree_pw.h"
#include "efield.h"
#include "gatefield.h"
#include "pot_local.h"
#include "pot_surchem.hpp"
#include "pot_xc.h"
#include "pot_sep.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/tool_quit.h"
#include "source_io/module_parameter/parameter.h"
#ifdef __LCAO
#include "H_TDDFT_pw.h"
#endif
#ifdef __MLALGO
#include "pot_ml_exx.h"
#endif

namespace elecstate
{

PotentialFactory::PotentialFactory(const ModuleBase::matrix* vloc_in,
                                   Structure_Factor* structure_factors_in,
                                   const ModulePW::PW_Basis* rho_basis_in,
                                   double* etxc_in,
                                   double* vtxc_in,
                                   ModuleBase::matrix& vofk_eff_out,
                                   surchem* solvent_in,
                                   VSep* vsep_cell_in,
                                   const UnitCell* ucell_in)
    : vloc_(vloc_in),
      structure_factors_(structure_factors_in),
      rho_basis_(rho_basis_in),
      etxc_(etxc_in),
      vtxc_(vtxc_in),
      vofk_eff_(vofk_eff_out),
      solvent_(solvent_in),
      vsep_cell_(vsep_cell_in),
      ucell_(ucell_in)
{
}

std::vector<PotBase*> PotentialFactory::create_components(const std::vector<std::string>& components_list)
{
    ModuleBase::TITLE("PotentialFactory", "create_components");
    std::vector<PotBase*> components;

    for (const auto& pot_type : components_list)
    {
        PotBase* tmp = nullptr;
        if (pot_type == "local")
        {
            tmp = create_pot_local();
        }
        else if (pot_type == "hartree")
        {
            tmp = create_pot_hartree();
        }
        else if (pot_type == "xc")
        {
            tmp = create_pot_xc();
        }
        else if (pot_type == "surchem")
        {
            tmp = create_pot_surchem();
        }
        else if (pot_type == "efield")
        {
            tmp = create_pot_efield();
        }
        else if (pot_type == "gatefield")
        {
            tmp = create_pot_gatefield();
        }
#ifdef __LCAO
        else if (pot_type == "tddft")
        {
            tmp = create_pot_tddft();
        }
#endif
#ifdef __MLALGO
        else if (pot_type == "ml_exx")
        {
            tmp = create_pot_ml_exx();
        }
#endif
        else if (pot_type == "dfthalf")
        {
            tmp = create_pot_dfthalf();
        }
        else
        {
            ModuleBase::WARNING_QUIT("PotentialFactory::create_components",
                                     "Unknown potential type: " + pot_type);
        }

        if (tmp != nullptr)
        {
            components.push_back(tmp);
        }
    }

    return components;
}

PotBase* PotentialFactory::create_pot_local()
{
    return new PotLocal(this->vloc_,
                        &(this->structure_factors_->strucFac),
                        this->rho_basis_,
                        nullptr);
}

PotBase* PotentialFactory::create_pot_hartree()
{
    return new PotHartree(this->rho_basis_);
}

PotBase* PotentialFactory::create_pot_xc()
{
    return new PotXC(this->rho_basis_, this->etxc_, this->vtxc_, &(this->vofk_eff_));
}

PotBase* PotentialFactory::create_pot_surchem()
{
    return new PotSurChem(this->rho_basis_,
                          this->structure_factors_,
                          nullptr,
                          this->solvent_);
}

PotBase* PotentialFactory::create_pot_efield()
{
    return new PotEfield(this->rho_basis_,
                         this->ucell_,
                         this->solvent_,
                         PARAM.inp.dip_cor_flag);
}

PotBase* PotentialFactory::create_pot_gatefield()
{
    return new PotGate(this->rho_basis_, this->ucell_);
}

#ifdef __LCAO
PotBase* PotentialFactory::create_pot_tddft()
{
    return new H_TDDFT_pw(this->rho_basis_, this->ucell_);
}
#endif

#ifdef __MLALGO
PotBase* PotentialFactory::create_pot_ml_exx()
{
    return new PotML_EXX(this->rho_basis_, this->ucell_);
}
#endif

PotBase* PotentialFactory::create_pot_dfthalf()
{
    return new PotSep(&(this->structure_factors_->strucFac), this->rho_basis_, this->vsep_cell_);
}

} // namespace elecstate
