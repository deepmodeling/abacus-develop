/**
 * @file deltaspin_init.cpp
 * @brief Solver-independent initialization of the spin-constrained state,
 *        plus the thin SpinConstrain::init_sc() shell.
 *
 * @par Called once at the start of a DeltaSpin calculation
 * init_sc_state() fills ScState from UnitCell/INPUT data. The
 * SpinConstrain::init_sc() member then only stores solver-side pointers
 * (Hamiltonian, psi, electronic state, density matrix, PW basis).
 *
 * @par Initialization order (critical):
 * 1. Input parameters (convergence thresholds, step sizes)
 * 2. Atom/orbital/lnchi counts (needed for array sizing)
 * 3. nspin and npol (determines which code paths are used)
 * 4. target_mag, lambda, constrain (from STRU parsing)
 * 5. For nspin=2: force x,y constraint flags to 0 (collinear: only z constrained)
 * 6. Parallel orbitals info (LCAO-specific)
 * 7. Solver parameters (Hamiltonian, psi, electronic state pointers)
 *
 * @par Error conditions
 * - If UnitCell.get_atom_Counts() returns empty map, subsequent operations will
 *   fail with "atomCounts is not set" in check_atomCounts()
 * - If nspin is not 2 or 4, set_nspin() will call WARNING_QUIT
 */
#include "deltaspin_init.h"

#include "source_cell/cell_tools.h"
#include "source_cell/unitcell.h"
#include "spin_constrain.h"

namespace spinconstrain
{

void init_sc_state(const ScInitParams& params, const UnitCell& ucell, ScState& state)
{
    // Step 1: Set input parameters for lambda loop
    // - sc_thr: convergence threshold for RMS(Mi - M_target) in uB
    // - nsc: maximum inner optimization steps
    // - nsc_min: minimum steps before early exit checks
    // - alpha_trial: initial trial step size (eV/uB^2), converted to Ry/uB^2
    // - sccut: maximum lambda change per step (eV/uB), converted to Ry/uB
    // - sc_drop_thr: fraction of initial RMS for adaptive threshold
    state.set_input_parameters(params.sc_thr, params.nsc, params.nsc_min,
                               params.alpha_trial, params.sccut, params.sc_drop_thr);

    // Step 2: Get atom/orbital/lnchi counts from UnitCell for indexing
    // atomCounts: {element_type_index -> number_of_atoms_of_this_type}
    // orbitalCounts: {element_type_index -> number_of_orbitals_per_atom}
    // lnchiCounts: {element_type_index -> {angular_momentum_L -> number_of_chi_functions}}
    state.set_atomCounts(ucell.get_atom_Counts());
    state.set_orbitalCounts(ucell.get_orbital_Counts());
    state.set_lnchiCounts(ucell.get_lnchi_Counts());

    // Step 3: Set spin configuration
    // nspin=2: collinear (spin-up/down separate k-points), npol=1
    // nspin=4: non-collinear (full spinor), npol=2
    state.set_nspin(params.nspin);
    state.set_npol((params.nspin == 4) ? 2 : 1);

    // Step 4: Load target magnetic moments and initial lambda from UnitCell
    // These are parsed from the STRU file's "sc_mag" and "lambda" keywords
    state.set_target_mag(unitcell::get_target_mag(ucell.atoms, ucell.ntype, ucell.nat));
    state.lambda_ = unitcell::get_lambda(ucell.atoms, ucell.ntype, ucell.nat);
    state.constrain_ = unitcell::get_constrain(ucell.atoms, ucell.ntype, ucell.nat);

    // Step 5: CRITICAL FIX for collinear spin (nspin=2)
    // In collinear mode, spins are constrained along the z-axis only.
    // The x and y components must be set to 0 to prevent the lambda optimizer
    // from trying to constrain non-existent transverse components.
    // Without this fix, the optimizer would waste iterations trying to
    // drive Mx and My to their (usually non-zero) target values, which
    // is physically meaningless for collinear calculations.
    if (params.nspin == 2)
    {
        for (int iat = 0; iat < static_cast<int>(state.constrain_.size()); iat++)
        {
            state.constrain_[iat].x = 0;
            state.constrain_[iat].y = 0;
        }
    }

    // Step 6: Set auxiliary parameters
    state.atomLabels_ = unitcell::get_atomLabels(ucell.atoms, ucell.ntype); // "Fe_0", "Fe_1", etc.
    state.direction_only_ = params.direction_only;                          // Only optimize spin direction
    state.tpiba = ucell.tpiba;                                              // 2*pi/a lattice scaling
    state.set_decay_grad();                                                 // Initialize gradient decay thresholds
}

} // namespace spinconstrain

template <typename TK>
void spinconstrain::SpinConstrain<TK>::init_sc(double sc_thr_in,
		int nsc_in,
		int nsc_min_in,
		double alpha_trial_in,
		double sccut_in,
		double sc_drop_thr_in,
		const UnitCell& ucell,
		bool direction_only_in,
		Parallel_Orbitals* ParaV_in,
		int nspin_in,
		const K_Vectors& kv_in,
		void* p_hamilt_in,
		void* psi_in,
#ifdef __LCAO
		elecstate::DensityMatrix<TK, double>* dm_in, // mohan add 2025-11-03
#endif
		elecstate::ElecState* pelec_in,
		ModulePW::PW_Basis_K* pw_wfc_in)
{
    // Steps 1-6: solver-independent state initialization
    const spinconstrain::ScInitParams params{sc_thr_in, nsc_in, nsc_min_in,
                                             alpha_trial_in, sccut_in, sc_drop_thr_in,
                                             direction_only_in, nspin_in};
    spinconstrain::init_sc_state(params, ucell, this->state_);

    // Step 7: Solver-side pointers and parallel orbitals info
    this->pw_wfc_ = pw_wfc_in;                        // PW basis (PW mode only)
    if(ParaV_in != nullptr) this->set_ParaV(ParaV_in);
    this->set_solver_parameters(kv_in, p_hamilt_in, psi_in, pelec_in);

    // Step 8: Set density matrix pointer (LCAO mode only)
#ifdef __LCAO
    this->dm_ = dm_in; // mohan add 2025-11-03
#endif
}

// Explicit template instantiations for both spin types
template class spinconstrain::SpinConstrain<std::complex<double>>;
template class spinconstrain::SpinConstrain<double>;
