#ifndef GLOBAL_H
#define GLOBAL_H

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_io/restart.h"

#include "source_relax/ions_move_methods.h"
#include "source_relax/lattice_change_methods.h"
#include "source_cell/unitcell.h"
#include "source_relax/bfgs.h"
#include "source_io/module_parameter/input_parameter.h"

#ifdef __EXX
#include "source_hamilt/module_xc/exx_info.h"
#include "source_lcao/module_ri/exx_lip.h"
#endif

#include "source_hamilt/module_xc/xc_functional.h"
#include "source_base/module_device/device_check.h"

#endif
