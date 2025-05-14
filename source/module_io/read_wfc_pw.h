#ifndef READ_WFC_PW_H
#define READ_WFC_PW_H

#include "module_basis/module_pw/pw_basis_k.h"
#include "module_elecstate/module_charge/charge.h"

#include <string>

namespace ModuleIO
{

/**
 * @brief read wave functions from binary file
 *
 * @param filename filename containing wave functions
 * @param pw_wfc wave function FFT grids
 * @param ik k index
 * @param ikstot total index of k points
 * @param nkstot total number of k points
 * @param wfc wave functions
 */
void read_wfc_pw(const std::string& filedir,
                 const ModulePW::PW_Basis_K* pw_wfc,
                 const int& ik,
                 const int& ikstot,
                 const int& nkstot,
                 ModuleBase::ComplexMatrix& wfc);
} // namespace ModuleIO

#endif
