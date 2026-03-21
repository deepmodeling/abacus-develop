#include "source_estate/elecstate.h"
#include "source_hamilt/module_xc/exx_info.h" // use GlobalC::exx_info

namespace elecstate
{

/// @brief calculation if converged
/// @date Peize Lin add 2016-12-03
void ElecState::set_exx(const double& Eexx)
{
    ModuleBase::TITLE("energy", "set_exx");

    if (GlobalC::exx_info.info_global.cal_exx)
    {
        this->f_en.exx = GlobalC::exx_info.info_global.hybrid_alpha * Eexx;
    }
    return;
}

}
