#ifndef NSCF_FERMI_SURF_H
#define NSCF_FERMI_SURF_H

#include "source_base/matrix.h"
#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include "source_cell/parallel_kpoints.h"

namespace ModuleIO
{
void nscf_fermi_surface(const std::string& out_band_dir,
                        const int& nband,
                        const double& ef,
                        const K_Vectors& kv,
                        const UnitCell& ucell,
                        const ModuleBase::matrix& ekb);
}
#endif

