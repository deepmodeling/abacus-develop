#ifndef POTENTIAL_IO_H
#define POTENTIAL_IO_H
#include <string>
#include "module_cell/unitcell.h"
#include "module_basis/module_pw/pw_basis.h"
#include "module_elecstate/module_charge/charge.h"

namespace ModuleIO
{

/// @brief write potential to file
/// @param[in] bz 
/// @param[in] nbz 
/// @param[in] nplane 
/// @param[in] startz_current 
/// @param[in] is 
/// @param[in] iter 
/// @param[in] fn 
/// @param[in] nx 
/// @param[in] ny 
/// @param[in] nz 
/// @param[in] v 
/// @param[in] precision 
/// @param[in] hartree 
void write_potential(
#ifdef __MPI
    const int& bz,
    const int& nbz,
    const int& nplane,
    const int& startz_current,
#endif
    const int& is,
    const int& iter,
    const std::string& fn,
    const int& nx,
    const int& ny,
    const int& nz,
    const ModuleBase::matrix& v,
    const int& precision,
    const int& hartree = 0);

/// @brief write electric static potential to file
/// @param bz 
/// @param nbz 
/// @param fn 
/// @param rho_basis 
/// @param chr 
/// @param ucell_ 
/// @param v_effective_fixed 
void write_elecstat_pot(
#ifdef __MPI
    const int& bz,
    const int& nbz,
#endif
    const std::string& fn,
    ModulePW::PW_Basis* rho_basis,
    const Charge* const chr,
    const UnitCell* ucell_,
    const double* v_effective_fixed);

} // namespace ModuleIO

#endif // POTENTIAL_IO_H