#include "source_estate/module_dm/dm_routine.h"
#include "source_estate/module_dm/dm_from_psi.h"
#include "source_estate/elecstate_tools.h"
#include "source_cell/cal_ux.h"

template <typename TK>
void module_dm::init_dm(UnitCell& ucell,
        elecstate::ElecState* pelec,
        module_dm::Setup_DM<TK> &dmat,
        psi::Psi<TK>* psi,
        Charge &chr,
        const Parallel_Orbitals& pv,
        const int iter,
        const int exx_two_level_step,
        const Init_DM_Config& cfg)
{
    ModuleBase::TITLE("elecstate", "init_dm");

    if (iter == 1 && exx_two_level_step == 0)
    {
        std::cout << " LCAO WAVEFUN -> CHARGE " << std::endl;

        elecstate::calEBand(pelec->ekb, pelec->wg, pelec->f_en);

        module_dm::dm_from_psi(&pv, pelec->wg, *psi, *dmat.dm);
        if (cfg.esolver_type != "tddft" && cfg.td_stype == 2)
        {
            dmat.dm->cal_dmr_td(*cfg.td_phase_hybrid, cfg.td_cart_At, -1);
        }
        else
        {
            dmat.dm->cal_dmr(-1);
        }

        // use density matrix to calculate the charge density
        cfg.dm2rho_func(dmat.dm->get_dmr_vec(), cfg.nspin, &chr, cfg.nelec, ucell.omega, false);

        unitcell::cal_ux(ucell, cfg.nspin);

        //! update the potentials by using new electron charge density
        pelec->pot->update_from_charge(&chr, &ucell);

        //! compute the correction energy for metals
        pelec->f_en.descf = pelec->cal_delta_escf();
    }

    return;
}


template void module_dm::init_dm<double>(UnitCell& ucell,
        elecstate::ElecState* pelec,
        module_dm::Setup_DM<double> &dmat,
        psi::Psi<double>* psi,
        Charge &chr,
        const Parallel_Orbitals& pv,
        const int iter,
        const int exx_two_level_step,
        const Init_DM_Config& cfg);

template void module_dm::init_dm<std::complex<double>>(UnitCell& ucell,
        elecstate::ElecState* pelec,
        module_dm::Setup_DM<std::complex<double>> &dmat,
        psi::Psi<std::complex<double>>* psi,
        Charge &chr,
        const Parallel_Orbitals& pv,
        const int iter,
        const int exx_two_level_step,
        const Init_DM_Config& cfg);

