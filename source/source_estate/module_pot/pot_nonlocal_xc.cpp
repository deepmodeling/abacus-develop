#include "pot_nonlocal_xc.h"

#include "source_hamilt/module_vdw/nonlocal_xc.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"

#include <exception>
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
    ModuleBase::TITLE("PotNonlocalXC", "cal_veff");
    ModuleBase::timer::start("PotNonlocalXC", "cal_v_eff");
    if (!chg || !cell || !rho_basis_ || chg->rhopw != rho_basis_ || chg->nspin != 1
        || chg->nrxx != rho_basis_->nrxx || !chg->rho || !chg->rho_core || v_eff.nr != chg->nspin
        || v_eff.nc != chg->nrxx)
    {
        ModuleBase::WARNING_QUIT("PotNonlocalXC", "the first rVV10 interface requires a consistent nspin=1 density basis");
    }
    std::vector<double> valence(chg->nrxx, 0.0);
    for (int is = 0; is < chg->nspin; ++is)
        for (int r = 0; r < chg->nrxx; ++r)
            valence[r] += chg->rho[is][r];
    std::vector<double> total = valence;
    for (int r = 0; r < chg->nrxx; ++r)
        total[r] += chg->rho_core[r];
    vdw::NonlocalXCResult nl = {0.0, 0.0, std::vector<double>()};
    try
    {
        nl = backend_->evaluate(*rho_basis_, total, valence);
    }
    catch (const std::exception& error)
    {
        ModuleBase::WARNING_QUIT("PotNonlocalXC", error.what());
    }
    if (nl.potential.size() != static_cast<std::size_t>(chg->nrxx))
    {
        ModuleBase::WARNING_QUIT("PotNonlocalXC", "nonlocal XC backend returned an invalid potential size");
    }

    // PotXC is registered before this component and owns the semilocal
    // contribution. This component only adds the selected nonlocal correction.
    *etxc_ += nl.energy;
    *vtxc_ += nl.vtxc;
    for (int is = 0; is < chg->nspin; ++is)
        for (int r = 0; r < chg->nrxx; ++r)
            v_eff(is, r) += nl.potential[r];
    ModuleBase::timer::end("PotNonlocalXC", "cal_v_eff");
}
} // namespace elecstate
