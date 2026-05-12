#include "spin_constrain.h"

/**
 * @file init_sc.cpp
 * @brief Master initialization for the SpinConstrain singleton.
 *
 * @par Called once at the start of a DeltaSpin calculation
 * This function bridges the UnitCell/InputPara data from the ESolver layer
 * to the internal SpinConstrain state. After init_sc(), the singleton is
 * fully configured and ready for the SCF lambda optimization loop.
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
    // Step 1: Set input parameters for lambda loop
    // - sc_thr: convergence threshold for RMS(Mi - M_target) in uB
    // - nsc: maximum inner optimization steps
    // - nsc_min: minimum steps before early exit checks
    // - alpha_trial: initial trial step size (eV/uB^2), converted to Ry/uB^2
    // - sccut: maximum lambda change per step (eV/uB), converted to Ry/uB
    // - sc_drop_thr: fraction of initial RMS for adaptive threshold
    this->set_input_parameters(sc_thr_in, nsc_in, nsc_min_in, alpha_trial_in, sccut_in, sc_drop_thr_in);

    // Step 2: Get atom/orbital/lnchi counts from UnitCell for indexing
    // atomCounts: {element_type_index -> number_of_atoms_of_this_type}
    // orbitalCounts: {element_type_index -> number_of_orbitals_per_atom}
    // lnchiCounts: {element_type_index -> {angular_momentum_L -> number_of_chi_functions}}
    this->set_atomCounts(ucell.get_atom_Counts());
    this->set_orbitalCounts(ucell.get_orbital_Counts());
    this->set_lnchiCounts(ucell.get_lnchi_Counts());

    // Step 3: Set spin configuration
    // nspin=2: collinear (spin-up/down separate k-points), npol=1
    // nspin=4: non-collinear (full spinor), npol=2
    this->set_nspin(nspin_in);
    this->set_npol((nspin_in == 4) ? 2 : 1);

    // Step 4: Load target magnetic moments and initial lambda from UnitCell
    // These are parsed from the STRU file's "sc_mag" and "lambda" keywords
    this->set_target_mag(ucell.get_target_mag());
    this->lambda_ = ucell.get_lambda();
    this->constrain_ = ucell.get_constrain();

    // Step 5: CRITICAL FIX for collinear spin (nspin=2)
    // In collinear mode, spins are constrained along the z-axis only.
    // The x and y components must be set to 0 to prevent the lambda optimizer
    // from trying to constrain non-existent transverse components.
    // Without this fix, the optimizer would waste iterations trying to
    // drive Mx and My to their (usually non-zero) target values, which
    // is physically meaningless for collinear calculations.
    if (nspin_in == 2)
    {
        for (int iat = 0; iat < static_cast<int>(this->constrain_.size()); iat++)
        {
            this->constrain_[iat].x = 0;
            this->constrain_[iat].y = 0;
        }
    }

    // Step 6: Set auxiliary parameters
    this->atomLabels_ = ucell.get_atomLabels();      // "Fe_0", "Fe_1", etc.
    this->direction_only_ = direction_only_in;        // Only optimize spin direction
    this->tpiba = ucell.tpiba;                        // 2*pi/a lattice scaling
    this->pw_wfc_ = pw_wfc_in;                        // PW basis (PW mode only)
    this->set_decay_grad();                           // Initialize gradient decay thresholds

    // Step 7: Set parallel orbitals info (for ScaLAPACK distributed matrices)
    if(ParaV_in != nullptr) this->set_ParaV(ParaV_in);

    // Step 8: Set solver parameters (pointers to external objects)
    this->set_solver_parameters(kv_in, p_hamilt_in, psi_in, pelec_in);

    // Step 9: Set density matrix pointer (LCAO mode only)
#ifdef __LCAO
    this->dm_ = dm_in; // mohan add 2025-11-03
#endif
}

// Explicit template instantiations for both spin types
template class spinconstrain::SpinConstrain<std::complex<double>>;
template class spinconstrain::SpinConstrain<double>;
