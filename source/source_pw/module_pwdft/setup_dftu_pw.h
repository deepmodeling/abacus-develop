#ifndef SETUP_DFTU_PW_H
#define SETUP_DFTU_PW_H

#include "source_cell/unitcell.h"
#include "source_base/matrix.h"
#include "source_estate/module_charge/charge_mixing.h"

struct Input_para;
class Plus_U_Base; // mohan add 2025-11-06

namespace DFTU_BASE
{

void iter_init_dftu_pw(const int iter,
                       const int istep,
                       Plus_U_Base& dftu, // mohan add 2025-11-06
                       const void* psi,
                       const ModuleBase::matrix& wg,
                       const UnitCell& ucell,
                       Charge_Mixing* p_chgmix,
                       const int* isk);

}

#endif
