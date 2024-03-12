#ifndef READ_WFC_PW_H
#define READ_WFC_PW_H

#include <string>

#include "module_basis/module_pw/pw_basis_k.h"

namespace ModuleIO
{

/**
 * @brief read wave functions from binary file
 *
 * @param filename filename containing wave functions
 * @param pw_wfc wave function FFT grids
 * @param ik k index
 * @param nkstot total number of k points
 * @param wfc wave functions
 * @return true if read successfully
 * @return false if read failed
 */
bool read_wfc_pw(const std::string& filedir,
                 const ModulePW::PW_Basis_K* pw_wfc,
                 const int& ik,
                 const int& nkstot,
                 ModuleBase::ComplexMatrix& wfc);

} // namespace ModuleIO

#endif
