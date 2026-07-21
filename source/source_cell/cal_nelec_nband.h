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