#include "pot_nonlocal_xc.h"

#include "source_hamilt/module_vdw/nonlocal_xc.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"

#include <vector>

namespace elecstate
{
PotNonlocalXC::PotNonlocalXC(const ModulePW::PW_Basis* density_basis,
                             double* etxc,
                             double* vtxc,
                             const std::string& method,
                             double b,
                             double c)
    : etxc_(etxc), vtxc_(vtxc)
{
    backend_ = vdw::make_nonlocal_xc(density_basis, method, b, c);
    rho_basis_ = density_basis;
    dynamic_mode = true;
    fixed_mode = false;
}

PotNonlocalXC::~PotNonlocalXC() = default;

void PotNonlocalXC::cal_v_eff(const Charge* chg, const UnitCell* cell, ModuleBase::matrix& v_eff)
{
    ModuleBase::TITLE("PotNonlocalXC", "cal_v_eff");
    ModuleBase::timer::start("PotNonlocalXC", "cal_v_eff");
    if (!chg || !cell || !rho_basis_ || chg->rhopw != rho_basis_ || chg->nspin != 1
        || chg->nrxx != rho_basis_->nrxx || !chg->rho || !chg->rho_core || v_eff.nr != chg->nspin
        || v_eff.nc != chg->nrxx)
    {
        ModuleBase::WARNING_QUIT("PotNonlocalXC", "nonlocal XC requires a consistent nspin=1 density basis");
    }
    // nspin=1 is enforced above, so the valence density is the single channel.
    const std::vector<double> valence(chg->rho[0], chg->rho[0] + chg->nrxx);
    std::vector<double> total(chg->nrxx, 0.0);
    for (int r = 0; r < chg->nrxx; ++r)
    {
        total[r] = valence[r] + chg->rho_core[r];
    }
    const vdw::NonlocalXCResult nl = backend_->evaluate(*rho_basis_, total, valence);
    if (nl.potential.size() != static_cast<std::size_t>(chg->nrxx))
    {
        ModuleBase::WARNING_QUIT("PotNonlocalXC", "nonlocal XC backend returned an invalid potential size");
    }

    // PotXC is registered before this component and owns the semilocal
    // contribution. This component only adds the selected nonlocal correction.
    *etxc_ += nl.energy;
    *vtxc_ += nl.vtxc;
    for (int r = 0; r < chg->nrxx; ++r)
    {
        v_eff(0, r) += nl.potential[r];
    }
    ModuleBase::timer::end("PotNonlocalXC", "cal_v_eff");
}
} // namespace elecstate
