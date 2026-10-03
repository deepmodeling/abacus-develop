#ifndef CTRL_OUTPUT_FP_H
#define CTRL_OUTPUT_FP_H

#include "source_estate/elecstate_lcao.h"
#include "source_basis/module_pw/pw_basis_big.h"

struct Input_para;

namespace ModuleIO
{

void ctrl_output_fp(UnitCell& ucell,
                    const Input_para& inp,
                    elecstate::ElecState* pelec,
                    ModulePW::PW_Basis_Big* pw_big,
                    ModulePW::PW_Basis* pw_rhod,
                    Charge& chr,
                    surchem& solvent,
                    Parallel_Grid& para_grid,
                    const int istep);

}
#endif
