#ifndef DFTU_PW_H
#define DFTU_PW_H

#include "source_cell/unitcell.h"
#include "source_base/matrix.h"

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
                       const Input_para& inp);

}

#endif
