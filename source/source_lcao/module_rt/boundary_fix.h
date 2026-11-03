/**
 * @file boundary_fix.h
 * @brief Correct the discontinuity that occurs when crossing periodic boundary conditions
 */
#ifndef BOUNDARY_FIX_H
#define BOUNDARY_FIX_H

#include "source_cell/unitcell.h"
#include <vector>
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_psi/psi.h"
#include "source_base/module_container/ATen/core/tensor.h"
namespace module_rt{

/**
*  @brief Add phases to the matrix and coefficient from the previous step to correct the boundary discontinuity.
*
* @param[in] ucell Unitcell information
* @param[in] kvec_d Direct k-point vectors
* @param[in] nks Number of local k points
* @param[in] cart_At_last Saved hybrid-gauge Cartesian A in Hartree atomic units; zero in other gauges
* @param[in] pv information of parallel
* @param[in] hk_last Hamiltonian matrix from last step 
* @param[in] sk_last Overlap matrix from last step 
* @param[in] psi_last Wavefunctions from last step
* @param[in] len_hs size of matrix element in this processor
* @param[out] hk_last the fixed hk matrix
* @param[out] sk_last the fixed sk matrix
* @param[out] sk_last the fixed wavefunctions
*/
void reset_matrix_boundary(const UnitCell& ucell,
                           const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                           const int nks,
                           const ModuleBase::Vector3<double>& cart_At_last,
                           const Parallel_Orbitals* pv,
                           ct::Tensor& hk_last,
                           ct::Tensor& sk_last,
                           psi::Psi<std::complex<double>>* psi_last,
                           const size_t len_hs);

/**
*  @brief Add extra phase to the matrix element belong to iat
*
* @param[in] phase extra phase
* @param[in] matk the matrix need to be fixed
* @param[in] pv information of parallel
* @param[in] iat atom index
* @param[out] matk the fixed matrix
*/
void boundary_shift_mat(const std::complex<double>& phase,
                        std::complex<double>* matk,
                        const Parallel_Orbitals* pv,
                        const size_t iat);
/**
*  @brief Add extra phase to the wfc coefficient belong to iat
*
* @param[in] phase extra phase
* @param[in] psi_k_last psi of last step
* @param[in] pv information of parallel
* @param[in] iat atom index
* @param[out] psi_k_last fixed psi of last step
*/
void boundary_shift_c(const std::complex<double>& phase,
                      std::complex<double>* psi_k_last,
                      const Parallel_Orbitals* pv,
                      const size_t iat);
}// namespace module_rt
#endif // BOUNDARY_FIX_H
