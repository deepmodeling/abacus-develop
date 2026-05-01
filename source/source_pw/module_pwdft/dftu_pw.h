#ifndef DFTU_PW_H
#define DFTU_PW_H

#include "source_io/module_parameter/parameter.h"
#include "source_cell/unitcell.h"
#include "source_base/matrix.h"
#include "source_estate/module_charge/charge_mixing.h"

class Plus_U;

namespace pw
{

void iter_init_dftu_pw(const int iter,
                       const int istep,
                       Plus_U& dftu,
                       const void* psi,
                       const ModuleBase::matrix& wg,
                       const UnitCell& ucell,
                       Charge_Mixing* p_chgmix);

}

#endif
