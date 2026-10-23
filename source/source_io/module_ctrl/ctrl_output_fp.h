#ifndef CTRL_OUTPUT_FP_H
#define CTRL_OUTPUT_FP_H

#include <fstream>

struct Input_para;
class UnitCell;
class Charge;
class surchem;
class Parallel_Grid;

namespace elecstate
{
class ElecState;
}

namespace ModulePW
{
class PW_Basis;
class PW_Basis_Big;
}

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
                    const int istep,
                    std::ofstream& ofs_running);

}
#endif
