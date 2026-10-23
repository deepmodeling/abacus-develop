#ifndef CHG_ROUTINE_H
#define CHG_ROUTINE_H

#include "source_estate/elecstate.h" // use pelec
#include "source_estate/module_charge/charge.h" // use chr
#include "source_estate/module_charge/chg_mix.h" // use p_chgmix
#include "source_basis/module_pw/pw_basis.h" // use PW_Basis
#include "source_io/module_parameter/input_parameter.h" // use Input_para
#include "source_cell/unitcell.h"

// Plus_U_Base forward declaration, full definition in source_pw/module_pwdft/dftu_base.h
class Plus_U_Base;

namespace module_charge
{

/// Aggregated SCF convergence thresholds and status flags for chgmixing_ks
struct ScfMixingCtx
{
    double hsolver_error;  ///< solver error from diagonalization
    double scf_thr;        ///< charge density convergence threshold
    double scf_ene_thr;    ///< energy convergence threshold
    bool converged_u;      ///< whether DFT+U has converged
    bool ks_run;           ///< whether the current run is a KS calculation (PARAM.globalv.ks_run)
    double drho;            ///< charge density deviation (in/out)
    bool oscillate_esolver; ///< whether esolver oscillates (out)
    bool conv_esolver;      ///< whether esolver converged (out)
};

void chgmixing_ks(const int iter,
        UnitCell& ucell,
        elecstate::ElecState* pelec,
        Charge &chr,
        const ModulePW::PW_Basis& rhopw,
        Charge_Mixing* p_chgmix,
        ScfMixingCtx& ctx,
        const Input_para& inp);

void chgmixing_ks_pw(const int iter,
        Charge_Mixing* p_chgmix,
        Plus_U_Base& dftu,
        const bool mag_converged, ///< whether DeltaSpin magnetization converged; pass true when sc_mag_switch is off
        const Input_para& inp); // input parameters

void chgmixing_ks_lcao(const int iter, // scf iteration number
        Charge_Mixing* p_chgmix, // charge mixing class
        Plus_U_Base& dftu,
        const int nnr, // dimension of density matrix
        const Input_para& inp); // input parameters

}


#endif
