#ifndef DELTASPIN_LCAO_H
#define DELTASPIN_LCAO_H

#include "source_cell/unitcell.h"
#include "source_io/module_parameter/input_parameter.h"

namespace ModuleESolver
{

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
