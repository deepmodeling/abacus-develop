#ifndef GLOBAL_H
#define GLOBAL_H

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_estate/module_charge/charge_mixing.h"
#include "source_pw/module_pwdft/VNL_in_pw.h"
#include "source_io/restart.h"
#include "source_relax/relax_driver.h"
#ifdef __EXX
#include "source_hamilt/module_xc/exx_info.h"
#include "source_lcao/module_ri/exx_lip.h"
#endif
#include "source_estate/magnetism.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_base/module_device/device_check.h"

//==========================================================
// EXPLAIN : define "GLOBAL CLASS"
//==========================================================
namespace GlobalC
{
//#ifdef __EXX
    extern Exx_Info exx_info;
//#endif
} // namespace GlobalC

#include "source_cell/parallel_kpoints.h"
#include "source_cell/unitcell.h"
namespace GlobalC
{
extern Restart restart; // Peize Lin add 2020.04.04
} // namespace GlobalC

// extern Magnetism mag;

#endif
