#ifndef SPARSE_FORMAT_EXX_H
#define SPARSE_FORMAT_EXX_H

#ifdef __EXX

// --------------------------------------------------------
// Header files - minimal set for declaration only
// --------------------------------------------------------

#include <RI/global/Tensor.h>
#include <array>
#include <map>
#include <vector>

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/unitcell.h"
#include "source_lcao/LCAO_HS_arrays.hpp"

// --------------------------------------------------------
// Namespace - merged into one block
// --------------------------------------------------------

namespace sparse_format
{

/**
 * @brief Calculate the Hamiltonian matrix elements in real space using EXX data in sparse format
 * 
 * This function computes the Hamiltonian matrix elements in real space (HR) from EXX data,
 * which is stored in a sparse tensor format. The results are added to the HS_Arrays structure.
 * 
 * @tparam Tdata Data type for the matrix elements (double or std::complex<double>)
 * @param ucell Unit cell information
 * @param pv Parallel orbitals information for distributed computation
 * @param HS_Arrays Structure to store Hamiltonian and overlap matrix arrays
 * @param current_spin Current spin channel (0 or 1 for spin-polarized calculations)
 * @param sparse_thr Threshold for sparse matrix construction
 * @param nmp Periodic boundary conditions in reciprocal space
 * @param Hexxs EXX data stored as a nested map structure of tensors
 */
template <typename Tdata>
void cal_HR_exx(
    const UnitCell& ucell,
    const Parallel_Orbitals& pv,
    LCAO_HS_Arrays& HS_Arrays,
    const int& current_spin,
    const double& sparse_thr,
    const int (&nmp)[3],
    const std::vector<std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<Tdata>>>>& Hexxs);

// Explicit instantiations for double and complex<double> types
extern template void cal_HR_exx<double>(
    const UnitCell& ucell,
    const Parallel_Orbitals& pv,
    LCAO_HS_Arrays& HS_Arrays,
    const int& current_spin,
    const double& sparse_thr,
    const int (&nmp)[3],
    const std::vector<std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>>>& Hexxs);

extern template void cal_HR_exx<std::complex<double>>(
    const UnitCell& ucell,
    const Parallel_Orbitals& pv,
    LCAO_HS_Arrays& HS_Arrays,
    const int& current_spin,
    const double& sparse_thr,
    const int (&nmp)[3],
    const std::vector<std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<std::complex<double>>>>>& Hexxs);

}

#endif // __EXX
#endif // SPARSE_FORMAT_EXX_H
