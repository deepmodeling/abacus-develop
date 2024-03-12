#ifndef DIPOLE_IO_H
#define DIPOLE_IO_H

#include <string>

#include "module_basis/module_pw/pw_basis.h"

namespace ModuleIO
{
void write_dipole(const double* rho_save,
                  const ModulePW::PW_Basis* rhopw,
                  const int& is,
                  const int& istep,
                  const std::string& fn,
                  const int& precision = 11,
                  const bool for_plot = false);

double prepare(const UnitCell &cell, int &dir);

} // namespace ModuleIO

#endif
