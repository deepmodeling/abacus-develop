#include "pot_rvv10.h"

#include "source_base/tool_quit.h"

#include <exception>
#include <vector>

namespace elecstate
{
PotRvv10::PotRvv10(const ModulePW::PW_Basis* density_basis, double* etxc, double* vtxc)
    : etxc_(etxc), vtxc_(vtxc), evaluator_(Rvv10::rvv10_b_default, Rvv10::rvv10_c_default)
{
    rho_basis_ = density_basis;
    dynamic_mode = true;
    fixed_mode = false;
}

void PotRvv10::cal_v_eff(const Charge* chg, const UnitCell* cell, ModuleBase::matrix& v_eff)
{
    if (!chg || !cell || !rho_basis_ || chg->rhopw != rho_basis_ || (chg->nspin != 1 && chg->nspin != 2)
        || chg->nrxx != rho_basis_->nrxx || !chg->rho || !chg->rho_core || v_eff.nr != chg->nspin
        || v_eff.nc != chg->nrxx)
    {
        ModuleBase::WARNING_QUIT("PotRvv10", "rVV10 requires a consistent nspin=1 or nspin=2 density basis");
    }
    std::vector<double> valence(chg->nrxx, 0.0);
    for (int is = 0; is < chg->nspin; ++is)
        for (int r = 0; r < chg->nrxx; ++r)
            valence[r] += chg->rho[is][r];
    std::vector<double> total = valence;
    for (int r = 0; r < chg->nrxx; ++r)
        total[r] += chg->rho_core[r];
    try
    {
        const auto nl = evaluator_.evaluate(*rho_basis_, total, valence);

        // PotXC is registered before this component and owns the semilocal
        // contribution. This component only adds the nonlocal correction.
        *etxc_ += nl.energy;
        *vtxc_ += nl.vtxc;
        for (int is = 0; is < chg->nspin; ++is)
            for (int r = 0; r < chg->nrxx; ++r)
                v_eff(is, r) += nl.potential[r];
    }
    catch (const std::exception& error)
    {
        ModuleBase::WARNING_QUIT("PotRvv10", error.what());
    }
}
} // namespace elecstate
