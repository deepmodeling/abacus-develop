/**
 * @file boundary_fix.h
 * @brief Correct the discontinuity that occurs when crossing periodic boundary conditions
 */
#ifndef BOUNDARY_FIX_H
#define BOUNDARY_FIX_H

#include "source_cell/unitcell.h"
#include "source_base/module_container/ATen/core/tensor.h"
#include "source_psi/psi.h"

#include <cstddef>
#include <complex>
#include <vector>

class BasisFlattenedIndexing;
class BlockCyclicDistributionInformation;

namespace module_rt{

/**
*  @brief Add phases to the matrix and coefficient from the previous step to correct the boundary discontinuity.
*
* @param[in] ucell Unitcell information
* @param[in] kvec_d Direct k-point vectors
* @param[in] cart_At_last Saved hybrid-gauge Cartesian A in Hartree atomic units; zero in other gauges
* @param[in] basis canonical global LCAO basis indexing
* @param[in] layout block-cyclic distribution of Psi; gathered H/S tensors are detected from their shape
* @param[in] hk_last Hamiltonian matrix from last step 
* @param[in] sk_last Overlap matrix from last step 
* @param[in] psi_last Wavefunctions from last step
* @param[out] hk_last the fixed hk matrix
* @param[out] sk_last the fixed sk matrix
* @param[out] sk_last the fixed wavefunctions
*/
void reset_matrix_boundary(const UnitCell& ucell,
                           const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                           const ModuleBase::Vector3<double>& cart_At_last,
                           const BasisFlattenedIndexing& basis,
                           const BlockCyclicDistributionInformation& layout,
                           ct::Tensor& hk_last,
                           ct::Tensor& sk_last,
                           psi::Psi<std::complex<double>>* psi_last);

}// namespace module_rt
#endif // BOUNDARY_FIX_H
