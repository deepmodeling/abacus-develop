#ifndef DIPOLE_IO_H
#define DIPOLE_IO_H

#include <string>

namespace ModuleIO
{
void write_dipole(const double *rho_save,
                          const int &is,
                          const int &istep,
                          const std::string &fn,
                          const int &precision = 11,
                          const bool for_plot = false);
} // namespace ModuleIO

#endif
