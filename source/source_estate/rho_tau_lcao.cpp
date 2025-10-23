#include "rho_tau_lcao.h"
#include "source_estate/module_dm/density_matrix.h" // use density matrix

namespace elecstate
{

// gamma only case
void psi2tau_lcao(std::vector<hamilt::HContainer<double>*> &dmr,
    const int nspin,
    Charge* chg)
{
    if (XC_Functional::get_ked_flag())
	{
		for (int is = 0; is < nspin; is++)
		{
			ModuleBase::GlobalFunc::ZEROS(chg->kin_r[is], chg->nrxx);
		}
		ModuleGint::cal_gint_tau(dmr, nspin, chg->kin_r);
    }
}

// multi-k case
void psi2rho_lcao(std::vector<hamilt::HContainer<double>*> &dmr,
    const int nspin,
    Charge* chg)
{
    ModuleBase::TITLE("elecstate", "psi2rho_lcao");
    ModuleBase::timer::tick("elecstate", "psi2rho_lcao");

    for (int is = 0; is < nspin; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(chg->rho[is], chg->nrxx);
    }

    ModuleGint::cal_gint_rho(dmr, inp.nspin, chg->rho);

    charge->renormalize_rho();

    // symmetrize of charge density should be here, mohan 20251023

    ModuleBase::timer::tick("elecstate", "psi2rho_lcao");
    return;
}

}// end namespace
