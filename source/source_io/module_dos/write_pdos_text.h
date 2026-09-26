#ifndef WRITE_PDOS_TEXT_H
#define WRITE_PDOS_TEXT_H

#include "source_base/matrix.h"
#include "source_cell/unitcell.h"

namespace ModuleIO
{

/// Write PDOS in plain text long-table format.
/// File name: pdoss{spin}g{geom}_{basis}.txt
/// Columns: energy(eV)  atom  species  l  z  m  pdos(1/eV)
/// For nspin=4, the two spinor components are summed.
void write_pdos_text(
        const UnitCell& ucell,
        const ModuleBase::matrix* pdos,
        const int nspin,
        const int nlocal,
        const int npoints,
        const double& emin,
        const double& dos_edelta_ev,
        const std::string& out_dir,
        const std::string& basis,
        const int istep);

}

#endif
