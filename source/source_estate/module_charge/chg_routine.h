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
class K_Vectors;

namespace module_charge
{

/// Aggregated SCF convergence thresholds and status flags for chgmixing_ks
struct ScfMixingCtx
{
    int iter_limit;        ///< maximum iterations for the current SCF attempt
    double hsolver_error;  ///< solver error from diagonalization
    double scf_thr;        ///< charge density convergence threshold
    double scf_ene_thr;    ///< energy convergence threshold
    bool converged_u;      ///< whether DFT+U has converged
    bool ks_run;           ///< whether the current run is a KS calculation (PARAM.globalv.ks_run)
    bool model_changed = false; ///< a potential component changed its model in this iteration (SCCS delayed start); the SCF is not converged
    double drho;            ///< charge density deviation (in/out)
    bool oscillate_esolver; ///< whether esolver oscillates (out)
    bool conv_esolver;      ///< whether esolver converged (out)
};

/// Lets a potential component change its model before mixing (SCCS delayed
/// start); true only in the iteration in which it does, which then restarts the
/// mixing. The first rank of the first pool announces the change on the screen.
bool update_pot_model(const int iter,
        const double drho,
        const K_Vectors& kv,
        elecstate::Potential& pot,
        Charge_Mixing* p_chgmix);

/// sccs_debug report of the SCCS correction in this SCF iteration, in the
/// running log of the first rank of the first pool.
void print_corr_iter(const K_Vectors& kv,
        const int sccs_debug,
        const double drho,
        const elecstate::Potential& pot);

void chgmixing_ks(const int iter,
        UnitCell& ucell,
        elecstate::ElecState* pelec,
        Charge &chr,
        const ModulePW::PW_Basis& rhopw,
        Charge_Mixing* p_chgmix,
        ScfMixingCtx& ctx,
        const Input_para& inp);

void chgmixing_ks_pw(const int iter,
        const int iter_limit,
        Charge_Mixing* p_chgmix,
        Plus_U_Base& dftu,
        const bool mag_converged, ///< whether DeltaSpin magnetization converged; pass true when sc_mag_switch is off
        const Input_para& inp); // input parameters

void chgmixing_ks_lcao(const int iter, // scf iteration number
        const int iter_limit,
        Charge_Mixing* p_chgmix, // charge mixing class
        Plus_U_Base& dftu,
        const int nnr, // dimension of density matrix
        const Input_para& inp); // input parameters

}


#endif
