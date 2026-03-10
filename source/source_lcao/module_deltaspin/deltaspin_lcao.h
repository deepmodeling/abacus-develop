#ifndef DELTASPIN_LCAO_H
#define DELTASPIN_LCAO_H

#include "source_cell/unitcell.h"
#include "source_cell/klist.h"
#include "source_io/module_parameter/input_parameter.h"

namespace ModuleESolver
{

/**
 * @brief Initialize DeltaSpin for LCAO method
 *
 * This function initializes the DeltaSpin calculation by setting up
 * the SpinConstrain object with input parameters.
 *
 * @param ucell Unit cell
 * @param inp Input parameters
 * @param pv Parallel orbitals
 * @param kv K-vectors
 * @param p_hamilt Pointer to Hamiltonian
 * @param psi Pointer to wave functions
 * @param dm Density matrix
 * @param pelec Pointer to electronic state
 */
template <typename TK>
void init_deltaspin_lcao(const UnitCell& ucell,
                          const Input_para& inp,
                          void* pv,
                          const K_Vectors& kv,
                          void* p_hamilt,
                          void* psi,
                          void* dm,
                          void* pelec);

/**
 * @brief Calculate magnetic moments for DeltaSpin in LCAO method
 *
 * This function calculates the magnetic moments for each atom
 * in the DeltaSpin method.
 *
 * @param iter Current iteration number
 * @param inp Input parameters
 */
template <typename TK>
void cal_mi_lcao_wrapper(const int iter, const Input_para& inp);

/**
 * @brief Run DeltaSpin lambda loop for LCAO method
 *
 * This function handles the lambda loop optimization for the DeltaSpin method
 * in LCAO calculations. It determines whether to skip the Hamiltonian solve
 * based on the convergence status of magnetic moments.
 *
 * @param iter Current iteration number
 * @param drho Charge density convergence criterion
 * @param inp Input parameters
 * @return bool Whether to skip the Hamiltonian solve
 */
template <typename TK>
bool run_deltaspin_lambda_loop_lcao(const int iter,
                                     const double drho,
                                     const Input_para& inp);

} // namespace ModuleESolver

#endif // DELTASPIN_LCAO_H
