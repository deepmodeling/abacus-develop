#include "rho_tau_lcao.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_lcao/module_gint/gint_interface.h"

void LCAO_domain::dm2rho(std::vector<hamilt::HContainer<double>*> &dmr,
    const int nspin,
    Charge* chg)
{
    ModuleBase::TITLE("LCAO_domain", "dm2rho");
    ModuleBase::timer::tick("LCAO_domain", "dm2rho");

    for (int is = 0; is < nspin; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(chg->rho[is], chg->nrxx);
    }

    ModuleGint::cal_gint_rho(dmr, nspin, chg->rho);

    chg->renormalize_rho();

    // symmetrize of charge density should be here, mohan 20251023

    ModuleBase::timer::tick("LCAO_domain", "dm2rho");
    return;
}


void LCAO_domain::dm2tau(std::vector<hamilt::HContainer<double>*> &dmr,
    const int nspin,
    Charge* chg)
{
    ModuleBase::TITLE("LCAO_domain", "dm2tau");
    ModuleBase::timer::tick("LCAO_domain", "dm2tau");

    if (XC_Functional::get_ked_flag())
    {
        for (int is = 0; is < nspin; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(chg->kin_r[is], chg->nrxx);
        }
        ModuleGint::cal_gint_tau(dmr, nspin, chg->kin_r);
    }

    ModuleBase::timer::tick("LCAO_domain", "dm2tau");
}


