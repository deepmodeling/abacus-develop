#ifndef CAL_NELEC_NBAND_H
#define CAL_NELEC_NBAND_H

#include "source_cell/atom_spec.h"

namespace unitcell {

    /**
     * @brief calculate the total number of electrons in system
     *
     * @param atoms [in] atom pointer
     * @param ntype [in] number of atom types
     * @param nelec [out] total number of electrons
     */
    void cal_nelec(const Atom* atoms, const int& ntype, double& nelec, const double nelec_delta);

    /**
     * @brief Calculate the total number of local numerical atomic orbitals.
     *
     * nlocal = sum over all atom types of (atoms[it].nw * atoms[it].na).
     * For nspin == 4 (non-collinear) each basis function carries 2 polarizations,
     * so nlocal is doubled.
     *
     * Shared by cal_atoms_info() (which stores the result in PARAM.globalv.nlocal)
     * and GintInfo::init_trace_lo_(), so those two can no longer drift apart.
     * cal_wfc() still repeats the loop inline because it also needs the per-type
     * prefix sums for Atom::stapos_wf, and asserts its own total against the value
     * cal_atoms_info() produced.
     *
     * @note atoms[it].nw must already be populated, i.e. Atom::set_index() must have
     *       run for every type before calling this.
     *
     * @param atoms [in] atom pointer
     * @param ntype [in] number of atom types
     * @param nspin [in] number of spin components
     * @return total number of local basis functions
     */
    int cal_nlocal(const Atom* atoms, const int ntype, const int nspin);

    /**
     * @brief Calculate the number of bands.
     *
     * IMPORTANT: The nbands parameter must be the user-specified value from INPUT file.
     * If nbands is 0, this function will auto-calculate a default value based on nelec.
     * If nbands is non-zero (user-specified), this function will validate and use it.
     * 
     * BUG FIX NOTE: Previously, cal_atoms_info() did not pass the user-specified nbands,
     * causing result.nbands to always be 0 and triggering auto-calculation regardless
     * of user input. This led to incorrect energy calculations (deviation ~139 eV).
     *
     * @param nelec [in] total number of electrons
     * @param nlocal [in] total number of local basis
     * @param nelec_spin [in] number of electrons for each spin
     * @param nbands [in/out] number of bands - must be user-specified value on input,
     *                         will be updated if auto-calculation is triggered (nbands==0)
     * @param esolver_type [in] solver type
     * @param lspinorb [in] spin-orbit coupling flag
     * @param nspin [in] number of spin components
     * @param basis_type [in] basis type
     * @param smearing_method [in] smearing method
     */
    void cal_nbands(const int& nelec, const int& nlocal, const std::vector<double>& nelec_spin, int& nbands,
                    const std::string& esolver_type, const bool lspinorb, const int nspin,
                    const std::string& basis_type, const std::string& smearing_method);

}

#endif