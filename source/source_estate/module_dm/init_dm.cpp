#include "source_estate/module_dm/init_dm.h"
#include "source_estate/module_dm/cal_dm_psi.h"
#include "source_estate/elecstate_tools.h"
#include "source_estate/cal_ux.h"

template <typename TK>
void elecstate::init_dm(UnitCell& ucell,
		elecstate::ElecStateLCAO<TK>* pelec,
        psi::Psi<TK>* psi,
		Charge &chr,
        const int iter,
        const int exx_two_level_step)
{
    ModuleBase::TITLE("elecstate", "init_dm");

	if (iter == 1 && exx_two_level_step == 0)
	{
		std::cout << " WAVEFUN -> CHARGE " << std::endl;

		// calculate the density matrix using read in wave functions
		// and then calculate the charge density on grid.

		pelec->skip_weights = true;
		elecstate::calculate_weights(pelec->ekb,
				pelec->wg,
				pelec->klist,
				pelec->eferm,
				pelec->f_en,
				pelec->nelec_spin,
				pelec->skip_weights);

		elecstate::calEBand(pelec->ekb, pelec->wg, pelec->f_en);
		elecstate::cal_dm_psi(pelec->DM->get_paraV_pointer(), pelec->wg, *psi, *(pelec->DM));
		pelec->DM->cal_DMR();

		pelec->psiToRho(*psi);
		pelec->skip_weights = false;

		elecstate::cal_ux(ucell);

		//! update the potentials by using new electron charge density
		pelec->pot->update_from_charge(&chr, &ucell);

		//! compute the correction energy for metals
		pelec->f_en.descf = pelec->cal_delta_escf();
	}

    return;
}


template void elecstate::init_dm<double>(UnitCell& ucell,
		elecstate::ElecStateLCAO<double>* pelec,
        psi::Psi<double>* psi,
		Charge &chr,
        const int iter,
        const int exx_two_level_step);

template void elecstate::init_dm<std::complex<double>>(UnitCell& ucell,
		elecstate::ElecStateLCAO<std::complex<double>>* pelec,
        psi::Psi<std::complex<double>>* psi,
		Charge &chr,
        const int iter,
        const int exx_two_level_step);

