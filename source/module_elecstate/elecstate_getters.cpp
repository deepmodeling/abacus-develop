#include "module_elecstate/elecstate_getters.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace elecstate
{

const double get_ucell_omega()
{
    return GlobalC::ucell.omega;
}

}
