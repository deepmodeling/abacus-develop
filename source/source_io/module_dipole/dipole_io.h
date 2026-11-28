#ifndef DIPOLE_IO_H
#define DIPOLE_IO_H

#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"

#include <fstream>
#include <string>

namespace ModuleIO
{
/** @brief Write aligned electronic dipoles; precision counts significant digits.
 *  @param append Append one record instead of truncating the file. Used by mode 1
 *      (fixed filename) and by the RT-TDDFT time series; mode 0 and mode 2 truncate.
 */
void write_dipole(const UnitCell& ucell,
                  const double* rho,
                  const ModulePW::PW_Basis* rhopw,
                  const int& istep,
                  const std::string& fn,
                  std::ofstream& ofs_running,
                  const bool append,
                  const int& precision = 13);

} // namespace ModuleIO

#endif
