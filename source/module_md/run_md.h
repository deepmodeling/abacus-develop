#ifndef RUN_MD_H
#define RUN_MD_H

#include "module_esolver/esolver.h"
#include "module_parameter/parameter.h"

/**
 * @brief the md loop line
 *
 */
namespace Run_MD
{
/**
 * @brief the md loop line
 *
 * @param unit_in unitcell information
 * @param p_esolver energy solver
 * @param md_para input parameters used in md
 */
void md_line(UnitCell& unit_in, ModuleESolver::ESolver* p_esolver, const Parameter& param_in);
} // namespace Run_MD

#endif