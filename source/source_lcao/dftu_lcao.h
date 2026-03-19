#ifndef DFTU_LCAO_H
#define DFTU_LCAO_H

#include "source_cell/unitcell.h"
#include "source_cell/klist.h"
#include "source_io/module_parameter/input_parameter.h"

namespace ModuleESolver
{

/**
 * @brief Initialize DFT+U for LCAO method in iter_init
 *
 * This function handles the DFT+U initialization during the SCF iteration.
 * It sets the density matrix and calculates Slater integrals if needed.
 *
 * @param istep Current ionic step
 * @param iter Current SCF iteration
 * @param inp Input parameters
 * @param dftu DFT+U object
 * @param dm Density matrix
 * @param ucell Unit cell
 * @param rho Charge density
 * @param nrxx Number of real space grid points
 */
template <typename TK>
void init_dftu_lcao(const int istep,
                     const int iter,
                     const Input_para& inp,
                     void* dftu,
                     void* dm,
                     const UnitCell& ucell,
                     double** rho,
                     const int nrxx);

/**
 * @brief Finish DFT+U calculation for LCAO method in iter_finish
 *
 * This function handles the DFT+U finalization during the SCF iteration.
 * It calculates the occupation matrix and energy correction if needed.
 *
 * @param iter Current SCF iteration
 * @param conv_esolver Whether ESolver has converged
 * @param inp Input parameters
 * @param dftu DFT+U object
 * @param ucell Unit cell
 * @param dm_vec Density matrix vector
 * @param kv K-vectors
 * @param mixing_beta Mixing beta parameter
 * @param hamilt_lcao Hamiltonian LCAO object
 */
template <typename TK>
void finish_dftu_lcao(const int iter,
                       const bool conv_esolver,
                       const Input_para& inp,
                       void* dftu,
                       const UnitCell& ucell,
                       const std::vector<std::vector<TK>>& dm_vec,
                       const K_Vectors& kv,
                       const double mixing_beta,
                       void* hamilt_lcao);

} // namespace ModuleESolver

#endif // DFTU_LCAO_H
