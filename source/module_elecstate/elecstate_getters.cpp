#include "module_elecstate/elecstate_getters.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace elecstate
{

int get_rhopw_nrxx()
{
    return GlobalC::rhopw->nrxx;
}
int get_rhopw_nxyz()
{
    return GlobalC::rhopw->nxyz;
}
double get_ucell_omega()
{
    return GlobalC::ucell.omega;
}

}
