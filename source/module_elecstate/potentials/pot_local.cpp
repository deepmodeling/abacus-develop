#include "pot_local.h"

#include "module_base/timer.h"
#include "module_base/tool_title.h"

#include <complex>

namespace elecstate
{

//==========================================================
// This routine computes the local potential in real space
//==========================================================
void PotLocal::cal_fixed_v(double *vl_pseudo // store the local pseudopotential
)
{
    ModuleBase::TITLE("PotLocal", "cal_fixed_v");
    ModuleBase::timer::tick("PotLocal", "cal_fixed_v");

    std::complex<double> *vg = new std::complex<double>[this->rho_basis_->npw];

    ModuleBase::GlobalFunc::ZEROS(vg, this->rho_basis_->npw);

    for (int it = 0; it < this->ntype_; it++)
    {
        for (int ig = 0; ig < this->rho_basis_->npw; ig++)
        {
            vg[ig] += this->vloc_[0](it, this->rho_basis_->ig2igg[ig]) * this->sf_[0](it, ig);
        }
    }

    // recip2real should be a const function, but now it isn't
    // a dangerous usage appears here, which should be fix in the future.
    const_cast<ModulePW::PW_Basis *>(this->rho_basis_)->recip2real(vg, vl_pseudo);

    delete[] vg;

    // GlobalV::ofs_running <<" set local pseudopotential done." << std::endl;
    ModuleBase::timer::tick("PotLocal", "cal_fixed_v");
    return;
}

} // namespace elecstate