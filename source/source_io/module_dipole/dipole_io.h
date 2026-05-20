#ifndef DIPOLE_IO_H
#define DIPOLE_IO_H

#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"

#include <fstream>
#include <string>

namespace ModuleIO
{
void write_dipole(const UnitCell& ucell,
                  const double* rho_save,
                  const ModulePW::PW_Basis* rhopw,
                  const int& istep,
                  const std::string& fn,
                  std::ofstream& ofs_running,
                  const int& precision = 11);

double prepare(const UnitCell& cell, int& dir);

} // namespace ModuleIO

#endif
