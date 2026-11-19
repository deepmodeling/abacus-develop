#ifndef WRITE_PDOS_TEXT_H
#define WRITE_PDOS_TEXT_H

#include "source_base/matrix.h"
#include "source_cell/unitcell.h"

namespace ModuleIO
{

/// Write PDOS in plain text long-table format.
/// Mode 2 adds g{geom} to pdoss{spin}.txt for numbered ionic steps.
/// Rows contain energy (eV), atom, species and one PDOS value per (l,m), summed over zeta.
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
        const int istep,
        const int out_app_flag);

}

#endif
