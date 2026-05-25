#include "source_estate/module_dm/init_dm.h"
#include "source_estate/module_dm/cal_dm_psi.h"
#include "source_estate/elecstate_tools.h"
#include "source_estate/cal_ux.h"
#include "source_lcao/rho_tau_lcao.h" // mohan add 2025-11-12
#include "source_lcao/module_rt/td_info.h"

template <typename TK>
void elecstate::init_dm(UnitCell& ucell,
		elecstate::ElecState* pelec,
        LCAO_domain::Setup_DM<TK> &dmat,
        psi::Psi<TK>* psi,
		Charge &chr,
        const int iter,
        const int exx_two_level_step)
{
    ModuleBase::TITLE("elecstate", "init_dm");

	if (iter == 1 && exx_two_level_step == 0)
	{
		std::cout << " LCAO WAVEFUN -> CHARGE " << std::endl;

		elecstate::calEBand(pelec->ekb, pelec->wg, pelec->f_en);

		elecstate::cal_dm_psi(dmat.dm->get_paraV_pointer(), pelec->wg, *psi, *dmat.dm);
		if (PARAM.inp.esolver_type!="tddft" && PARAM.inp.td_stype == 2)
		{
			dmat.dm->cal_DMR_td(ucell, TD_info::cart_At);
		}
		else
		{
			dmat.dm->cal_DMR();
		}

        // mohan add 2025-11-12, use density matrix to calculate the charge density
        LCAO_domain::dm2rho(dmat.dm->get_DMR_vector(), PARAM.inp.nspin, &chr);

		elecstate::cal_ux(ucell);

		//! update the potentials by using new electron charge density
		pelec->pot->update_from_charge(&chr, &ucell);

		//! compute the correction energy for metals
		pelec->f_en.descf = pelec->cal_delta_escf();
	}

    return;
}


template void elecstate::init_dm<double>(UnitCell& ucell,
		elecstate::ElecState* pelec,
        LCAO_domain::Setup_DM<double> &dmat,
        psi::Psi<double>* psi,
		Charge &chr,
        const int iter,
        const int exx_two_level_step);

template void elecstate::init_dm<std::complex<double>>(UnitCell& ucell,
		elecstate::ElecState* pelec,
        LCAO_domain::Setup_DM<std::complex<double>> &dmat,
        psi::Psi<std::complex<double>>* psi,
		Charge &chr,
        const int iter,
        const int exx_two_level_step);

