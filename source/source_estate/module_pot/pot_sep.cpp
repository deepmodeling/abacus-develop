#include "pot_sep.h"

#include "source_base/timer.h"
#include "source_base/tool_title.h"

namespace elecstate
{
void PotSep::cal_fixed_v(double* vl_pseudo)
{
    ModuleBase::TITLE("PotSep", "cal_fixed_v");
    ModuleBase::timer::tick("PotSep", "cal_fixed_v");

    // GlobalC::vsep_cell.generate_vsep_r(this->rho_basis_[0], this->sf_[0]);

    // const_cast<VSep*>(this->vsep_)->generate_vsep_r(this->rho_basis_[0], this->sf_[0]);

    if (vsep_cell != nullptr)
    {
        for (int ir = 0; ir < this->rho_basis_->nrxx; ++ir)
        {
            vl_pseudo[ir] += vsep_cell->vsep_r[ir];
        }
    }

    ModuleBase::timer::tick("PotSep", "cal_fixed_v");
    return;
}
} // namespace elecstate
