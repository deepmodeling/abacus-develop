#include "module_elecstate/elecstate_getters.h"

#include "module_cell/unitcell.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_parameter/parameter.h"
#include "module_hamilt_general/module_xc/xc_functional.h"

namespace elecstate
{

int get_xc_func_type()
{
    return XC_Functional::get_func_type();
}


} // namespace elecstate
