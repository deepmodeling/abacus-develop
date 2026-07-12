#include "pot_sep.h"

#include "source_base/timer.h"
#include "source_base/tool_title.h"

namespace elecstate
{
void PotSep::cal_fixed_v(double* vl_pseudo)
{
    ModuleBase::TITLE("PotSep", "cal_fixed_v");
    ModuleBase::timer::start("PotSep", "cal_fixed_v");

    if (vsep_cell != nullptr)
    {
        for (int ir = 0; ir < this->rho_basis_->nrxx; ++ir)
        {
            vl_pseudo[ir] += vsep_cell->vsep_r[ir];
        }
    }

    ModuleBase::timer::end("PotSep", "cal_fixed_v");
    return;
}
} // namespace elecstate
