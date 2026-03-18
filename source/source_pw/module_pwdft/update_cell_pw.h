#ifndef UPDATE_CELL_PW_H
#define UPDATE_CELL_PW_H

#include "source_io/module_parameter/parameter.h"
#include "source_cell/unitcell.h"
#include "source_cell/klist.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_pw/module_pwdft/vnl_pw.h"

namespace pw
{

void update_cell_pw(const UnitCell& ucell,
                    pseudopot_cell_vnl& ppcell,
                    const K_Vectors& kv,
                    ModulePW::PW_Basis_K* pw_wfc,
                    const Input_para& inp);

}
#endif
