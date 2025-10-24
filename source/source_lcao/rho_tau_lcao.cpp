#include "rho_tau_lcao.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_lcao/module_gint/gint_interface.h"

void LCAO_domain::dm2rho(std::vector<hamilt::HContainer<double>*> &dmr,
    const int nspin,
    Charge* chr)
{
    ModuleBase::TITLE("LCAO_domain", "dm2rho");
    ModuleBase::timer::tick("LCAO_domain", "dm2rho");

    for (int is = 0; is < nspin; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(chr->rho[is], chr->nrxx);
    }

    ModuleGint::cal_gint_rho(dmr, nspin, chr->rho);

    chr->renormalize_rho();

    // should be moved somewhere else, mohan 20251024
	if (XC_Functional::get_ked_flag())
	{
		dm2tau(dmr, nspin, chr);
	}

    // symmetrize of charge density should be here, mohan 20251023

    ModuleBase::timer::tick("LCAO_domain", "dm2rho");
    return;
}


void LCAO_domain::dm2tau(std::vector<hamilt::HContainer<double>*> &dmr,
    const int nspin,
    Charge* chr)
{
    ModuleBase::TITLE("LCAO_domain", "dm2tau");
    ModuleBase::timer::tick("LCAO_domain", "dm2tau");

	for (int is = 0; is < nspin; is++)
	{
		ModuleBase::GlobalFunc::ZEROS(chr->kin_r[is], chr->nrxx);
	}
	ModuleGint::cal_gint_tau(dmr, nspin, chr->kin_r);

    ModuleBase::timer::tick("LCAO_domain", "dm2tau");
}
