#ifndef SETUP_ESTATE_PW_H
#define SETUP_ESTATE_PW_H

#include "source_base/module_device/device.h"
#include "source_cell/unitcell.h"
#include "source_cell/klist.h"
#include "source_pw/module_pwdft/structure_factor.h"
#include "source_estate/elecstate.h"
#include "source_pw/module_pwdft/vl_pw.h"
#include "source_pw/module_pwdft/vsep_pw.h"

namespace elecstate
{

void setup_estate_pw(
    UnitCell& ucell,
    K_Vectors& kv,
    Structure_Factor& sf,
    elecstate::ElecState*& pelec,
    Charge& chr,
    pseudopot_cell_vl& locpp,
    pseudopot_cell_vnl& ppcell,
    VSep*& vsep_cell,
    ModulePW::PW_Basis_K* pw_wfc,
    ModulePW::PW_Basis* pw_rho,
    ModulePW::PW_Basis* pw_rhod,
    ModulePW::PW_Basis_Big* pw_big,
    surchem& solvent,
    const Input_para& inp);

void teardown_estate_pw(elecstate::ElecState*& pelec, VSep*& vsep_cell);

template <typename T, typename Device>
void setup_estate_pw_impl(
    UnitCell& ucell,
    K_Vectors& kv,
    Structure_Factor& sf,
    elecstate::ElecState*& pelec,
    Charge& chr,
    pseudopot_cell_vl& locpp,
    pseudopot_cell_vnl& ppcell,
    VSep*& vsep_cell,
    ModulePW::PW_Basis_K* pw_wfc,
    ModulePW::PW_Basis* pw_rho,
    ModulePW::PW_Basis* pw_rhod,
    ModulePW::PW_Basis_Big* pw_big,
    surchem& solvent,
    const Input_para& inp);

template <typename T, typename Device>
void teardown_estate_pw_impl(elecstate::ElecState*& pelec, VSep*& vsep_cell);

}

#endif
