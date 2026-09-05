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
/// @param precision: output precision
/// @param nspin: number of spin channels (1, 2, or 4)
/// @param efield_flag: whether electric field is applied
/// @param dip_cor_flag: whether dipole correction is applied
/// @param imp_sol: whether implicit solvation model is used
/// @param two_fermi: whether two Fermi levels are used
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
    const int precision,
    const int nspin,
    const bool efield_flag,
    const bool dip_cor_flag,
    const bool imp_sol,
    const bool two_fermi);

} // namespace ModuleIO

#endif // POTENTIAL_IO_H
