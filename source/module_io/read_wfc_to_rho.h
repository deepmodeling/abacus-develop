#ifndef READ_WFC_TO_RHO_H
#define READ_WFC_TO_RHO_H

#include "module_basis/module_pw/pw_basis_k.h"
#include "module_elecstate/module_charge/charge.h"

#include <string>

namespace ModuleIO
{
/**
 * @brief read wave functions and occupation numbers to charge density
 *
 * @param pw_wfc pw basis for wave functions
 * @param symm symmetry
 * @param nkstot total number of k points
 * @param isk k index to spin index
 * @param chg charge density
 */
void read_wfc_to_rho(const ModulePW::PW_Basis_K* pw_wfc,
                     ModuleSymmetry::Symmetry& symm,
                     const int* ik2iktot,
                     const int nkstot,
                     const std::vector<int>& isk,
                     Charge& chg);

} // namespace ModuleIO

#endif
