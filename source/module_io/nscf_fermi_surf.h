#ifndef NSCF_FERMI_SURF_H
#define NSCF_FERMI_SURF_H
#include "module_base/matrix.h"
#include "module_cell/klist.h"
#include "module_cell/unitcell.h"
#include "module_cell/parallel_kpoints.h"

namespace ModuleIO
{
void nscf_fermi_surface(const std::string& out_band_dir,
                        const int& nband,
                        const double& ef,
                        const K_Vectors& kv,
                        const Parallel_Kpoints& Pkpoints,
                        const UnitCell& ucell,
                        const ModuleBase::matrix& ekb);
}
#endif

