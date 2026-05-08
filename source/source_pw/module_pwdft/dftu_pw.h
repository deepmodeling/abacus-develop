#ifndef DFTU_PW_H
#define DFTU_PW_H

#include "source_cell/unitcell.h"
#include "source_base/matrix.h"
#include "source_estate/module_charge/charge_mixing.h"

struct Input_para;
class Plus_U;

namespace pw
{

void iter_init_dftu_pw(const int iter,
                       const int istep,
                       Plus_U& dftu,
                       const void* psi,
                       const ModuleBase::matrix& wg,
                       const UnitCell& ucell,
                       Charge_Mixing* p_chgmix,
                       const int* isk);

}

#endif
