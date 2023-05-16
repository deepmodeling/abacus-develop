#include "module_elecstate/elecstate_getters.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace elecstate
{

double get_ucell_omega()
{
    return GlobalC::ucell.omega;
}

} // namespace elecstate
