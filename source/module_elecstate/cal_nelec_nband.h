#ifndef CAL_NELEC_NBAND_H
#define CAL_NELEC_NBAND_H

#include "module_cell/atom_spec.h"

namespace elecstate {

    /**
     * @brief calculate the total number of electrons in system
     *
     * @param atoms [in] atom pointer
     * @param ntype [in] number of atom types
     * @param nelec [out] total number of electrons
     */
    void cal_nelec(const Atom* atoms, const int& ntype, double& nelec);

    /**
     * @brief Calculate the number of bands.
     *
     * @param nelec [in] total number of electrons
     * @param nlocal [in] total number of local basis
     * @param nelec_spin [in] number of electrons for each spin
     * @param nbands  [out] number of bands
     */
    void cal_nbands(const int& nelec, const int& nlocal, const std::vector<double>& nelec_spin, int& nbands);

}

#endif