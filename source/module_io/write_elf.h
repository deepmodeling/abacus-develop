#ifndef WRITE_ELF_H
#define WRITE_ELF_H
#include <string>
#include "module_cell/unitcell.h"
#include "module_basis/module_pw/pw_basis.h"
#include "module_elecstate/module_charge/charge.h"

namespace ModuleIO
{
void write_elf(
#ifdef __MPI
    const int& bz,
    const int& nbz,
#endif
    const std::string& out_dir,
    const int& istep,
    const int& nspin,
    const double* const* rho,
    const double* const* tau,
    ModulePW::PW_Basis* rho_basis,
    const Parallel_Grid& pgrid,
    const UnitCell* ucell_,
    const int& precision);
}

#endif