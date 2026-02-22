#ifndef POTENTIAL_IO_H
#define POTENTIAL_IO_H
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge.h"
#include "source_hamilt/module_surchem/surchem.h"

#include <string>

namespace ModuleIO
{

/// @brief write electric static potential to file
/// @param bz
/// @param nbz
/// @param fn
/// @param istep
/// @param rho_basis
/// @param chr
/// @param ucell_
/// @param v_eff_fixed
/// @param solvent: for solvation model
/// #param precision: output precision
void write_elecstat_pot(
#ifdef __MPI
    const int& bz,
    const int& nbz,
#endif
    const std::string& fn,
    const int& istep,
    ModulePW::PW_Basis* rho_basis,
    const Charge* const chr,
    const UnitCell* ucell_,
    const double* v_eff_fixed,
    const surchem& solvent,
    const int precision);

} // namespace ModuleIO

#endif // POTENTIAL_IO_H
