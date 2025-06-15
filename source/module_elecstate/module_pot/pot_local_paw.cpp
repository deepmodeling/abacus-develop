#include "pot_local_paw.h"

#include "module_base/timer.h"
#include "module_base/tool_title.h"

#include <complex>

namespace elecstate
{

//==========================================================
// This routine computes the local potential in real space
//==========================================================
void PotLocal_PAW::cal_fixed_v(double *vl_pseudo // store the local pseudopotential
)
{
    ModuleBase::TITLE("PotLocal_PAW", "cal_fixed_v");
    ModuleBase::timer::tick("PotLocal_PAW", "cal_fixed_v");

    // GlobalV::ofs_running <<" set local pseudopotential done." << std::endl;
    ModuleBase::timer::tick("PotLocal_PAW", "cal_fixed_v");
    return;
}

} // namespace elecstate