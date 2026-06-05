#ifndef DIAGO_PARAMS_H
#define DIAGO_PARAMS_H

#include "source_io/module_parameter/input_parameter.h"

namespace hsolver
{

/**
 * @brief Setup diagonalization parameters for PW method
 *
 * This function sets up the diagonalization parameters for plane wave method,
 * including subspace diagonalization flag, SCF iteration number, diagonalization
 * threshold, and maximum number of diagonalization steps.
 *
 * @param istep Current ionic step
 * @param iter Current SCF iteration
 * @param ethr Diagonalization threshold
 * @param inp Input parameters
 */
template <typename T, typename Device>
void setup_diago_params_pw(const int istep,
                            const int iter,
                            const double ethr,
                            const Input_para& inp);

/**
 * @brief Setup diagonalization parameters for SDFT method
 *
 * This function sets up the diagonalization parameters for stochastic DFT method,
 * including subspace diagonalization flag, diagonalization threshold, and
 * maximum number of diagonalization steps.
 *
 * @param istep Current ionic step
 * @param iter Current SCF iteration
 * @param ethr Diagonalization threshold
 * @param inp Input parameters
 */
template <typename T, typename Device>
void setup_diago_params_sdft(const int istep,
                              const int iter,
                              const double ethr,
                              const Input_para& inp);

} // namespace hsolver

#endif // DIAGO_PARAMS_H
