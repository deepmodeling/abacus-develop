/**
 * @file deltaspin_init.h
 * @brief Solver-independent initialization of the spin-constrained state.
 *
 * @par Purpose
 * Bridges UnitCell/STRU parsing results (target moments, initial lambda,
 * constraint flags) and INPUT parameters into ScState. Kept separate from
 * SpinConstrain so that the data initialization logic has no dependency
 * on solver-side objects (Hamiltonian, wavefunctions, k-points).
 */
#ifndef DELTASPIN_INIT_H
#define DELTASPIN_INIT_H

#include "deltaspin_state.h"

class UnitCell;

namespace spinconstrain
{

/// Input parameters for DeltaSpin initialization (same units as INPUT/STRU;
/// unit conversion to Ry happens inside init_sc_state / ScState).
struct ScInitParams {
    double sc_thr;        ///< RMS(Mi - M_target) convergence threshold (uB)
    int    nsc;           ///< Maximum inner lambda optimization steps
    int    nsc_min;       ///< Minimum steps before early exit checks
    double alpha_trial;   ///< Initial trial step size (eV/uB^2)
    double sccut;         ///< Maximum lambda change per step (eV/uB)
    double sc_drop_thr;   ///< Fraction of initial RMS for adaptive threshold
    bool   direction_only;///< Only optimize spin direction
    int    nspin;         ///< 2=collinear, 4=non-collinear
};

/**
 * @brief Populate ScState from UnitCell and input parameters.
 *
 * @details Performs the solver-independent portion of init_sc():
 *   1. Set input parameters (thresholds, step sizes; unit conversion to Ry)
 *   2. Get atom/orbital/lnchi counts from UnitCell for indexing
 *   3. Set nspin and npol (nspin=4 -> npol=2, nspin=2 -> npol=1)
 *   4. Load target_mag, lambda, constrain from UnitCell (parsed from STRU)
 *   5. For nspin=2: force x,y constraint flags to 0 (collinear: only z constrained)
 *   6. Set atom labels, direction_only, tpiba; zero-initialize decay_grad
 *
 * @param params Input parameters (see ScInitParams)
 * @param ucell Unit cell with STRU constraint data
 * @param state ScState to fill (in/out)
 */
void init_sc_state(const ScInitParams& params, const UnitCell& ucell, ScState& state);

} // namespace spinconstrain

#endif // DELTASPIN_INIT_H
