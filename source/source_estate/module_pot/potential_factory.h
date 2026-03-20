#ifndef POTENTIAL_FACTORY_H
#define POTENTIAL_FACTORY_H

#include <string>
#include <vector>
#include "pot_base.h"
#include "source_pw/module_pwdft/structure_factor.h"

class UnitCell;
class surchem;
class VSep;

namespace elecstate
{

class PotentialFactory
{
  public:
    PotentialFactory(const ModuleBase::matrix* vloc_in,
                     Structure_Factor* structure_factors_in,
                     const ModulePW::PW_Basis* rho_basis_in,
                     double* etxc_in,
                     double* vtxc_in,
                     ModuleBase::matrix& vofk_eff_out,
                     surchem* solvent_in,
                     VSep* vsep_cell_in,
                     const UnitCell* ucell_in);

    std::vector<PotBase*> create_components(const std::vector<std::string>& components_list);

  private:
    PotBase* create_pot_local();
    PotBase* create_pot_hartree();
    PotBase* create_pot_xc();
    PotBase* create_pot_surchem();
    PotBase* create_pot_efield();
    PotBase* create_pot_gatefield();
#ifdef __LCAO
    PotBase* create_pot_tddft();
#endif
#ifdef __MLALGO
    PotBase* create_pot_ml_exx();
#endif
    PotBase* create_pot_dfthalf();

    const ModuleBase::matrix* vloc_ = nullptr;
    Structure_Factor* structure_factors_ = nullptr;
    const ModulePW::PW_Basis* rho_basis_ = nullptr;
    double* etxc_ = nullptr;
    double* vtxc_ = nullptr;
    ModuleBase::matrix& vofk_eff_;
    surchem* solvent_ = nullptr;
    VSep* vsep_cell_ = nullptr;
    const UnitCell* ucell_ = nullptr;
};

} // namespace elecstate

#endif
