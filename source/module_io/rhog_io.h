#ifndef RHOG_IO_H
#define RHOG_IO_H

#include <string>

#include "module_basis/module_pw/pw_basis.h"

namespace ModuleIO
{

/**
 * @brief read charge density in G space from binary file
 *
 * @param filename binary filename containing charge density in G space
 * @param pw_rhod charge density FFT grids
 * @param rhog charge density in G space
 * @return true if read successfully
 * @return false if read failed
 */
bool read_rhog(const std::string& filename, const ModulePW::PW_Basis* pw_rhod, std::complex<double>** rhog);

void write_rhog();

} // namespace ModuleIO

#endif
