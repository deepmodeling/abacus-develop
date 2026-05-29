#ifndef WRITE_ELF_H
#define WRITE_ELF_H
#include <string>
#include "source_cell/unitcell.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_estate/module_charge/charge.h"

namespace ModuleIO
{
void write_elf(
    const std::string& out_dir,
    const int& istep_in,
    const int& nspin,
    const double* const* rho,
    const double* const* tau,
    ModulePW::PW_Basis* rho_basis,
    const Parallel_Grid& pgrid,
    const UnitCell* ucell_,
    const int& precision,
    const std::string& geom_block,
    const bool two_fermi);
}

#endif
