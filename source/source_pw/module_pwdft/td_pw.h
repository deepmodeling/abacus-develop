#ifndef TD_PW_H
#define TD_PW_H

#include "source_basis/module_pw/pw_basis_k.h"
#include "source_cell/unitcell.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_pw/module_pwdft/vnl_pw.h"

namespace pw
{

/** @brief Check PW propagation restrictions after resolving the XC functional and EXX configuration. */
void check_td_input(const Input_para& input, const UnitCell& cell, const bool needs_ked, const bool cal_exx);

/** @brief Return the global unshifted momentum bound after updating the distributed basis. */
double td_momentum_bound(const ModulePW::PW_Basis_K& basis, const double tpiba);

/** @brief Extend the projector table to cover both endpoint and propagation fields in Hartree atomic units. */
void ensure_td_vnl(const UnitCell& cell, const double unshifted_bound,
                   const ModuleBase::Vector3<double>& A_right_ha,
                   const ModuleBase::Vector3<double>& A_prop_ha,
                   pseudopot_cell_vnl* projectors);

} // namespace pw
#endif
