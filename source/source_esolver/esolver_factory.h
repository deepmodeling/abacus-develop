#ifndef ESOLVER_FACTORY_H
#define ESOLVER_FACTORY_H

#include <string>

struct Input_para;

namespace ModuleESolver
{

class ESolver;

/**
 * @brief Determine the ESolver type string from input parameters.
 *
 * The type is decided based on PARAM.inp.basis_type and PARAM.inp.esolver_type,
 * together with device/precision hints written to the running log.
 *
 * @return [out] std::string The type label consumed by init_esolver().
 */
std::string determine_type();

/**
 * @brief Determine and initialize an ESolver based on input information.
 *
 * This function determines the type of ESolver to create based on input
 * information and initializes the corresponding ESolver child class. It
 * supports various ESolver types including ksdft_pw, ksdft_lcao,
 * ksdft_lcao_tddft, sdft_pw, ofdft, lj_pot, and dp_pot.
 *
 * @param [in] inp  Input parameters used to select and configure the ESolver.
 * @return [out]    A pointer to the newly created ESolver object.
 */
ESolver* init_esolver(const Input_para& inp);

} // namespace ModuleESolver

#endif
