/**
 * @file norm_psi.h
 * @brief normalize the wave function
 *  This file originally belonged to file LCAO_evolve.cpp
 */
#ifndef NORM_PSI_H
#define NORM_PSI_H

#include "source_base/module_container/ATen/core/tensor.h" // ct::Tensor
#include "source_lcao/module_rt/kernels/cublasmp_context.h"

#include <complex>

class BlockCyclicDistribution; // Defined in source/source_base/parallel_blkcyc.h

namespace module_rt
{
#ifdef __MPI
/**
 *  @brief normalize the wave function
 *
 * @param[in] matrix_layout block-cyclic distribution of the square matrix (S)
 * @param[in] wfc_layout block-cyclic distribution of the wavefunction (nbasis x nbands)
 * @param[in] mo_repr_layout block-cyclic layout of the MO representation of S
 * @param[in] nband number of bands
 * @param[in] nlocal number of orbitals
 * @param[in] Stmp overlap matrix
 * @param[in] psi_k psi of this step
 * @param[in] print_matirx print internal matrix or not
 * @param[out] psi_k psi of this step after normalization
 */
void norm_psi(const BlockCyclicDistribution& matrix_layout,
              const BlockCyclicDistribution& wfc_layout,
              const BlockCyclicDistribution& mo_repr_layout,
              const int nband,
              const int nlocal,
              const std::complex<double>* Stmp,
              std::complex<double>* psi_k,
              std::ofstream& ofs_running,
              const int print_matrix);

void norm_psi_tensor(const BlockCyclicDistribution& matrix_layout,
                     const BlockCyclicDistribution& wfc_layout,
                     const BlockCyclicDistribution& mo_repr_layout,
                     const int nband,
                     const int nlocal,
                     const ct::Tensor& Stmp,
                     ct::Tensor& psi_k,
                     std::ofstream& ofs_running,
                     const int print_matrix,
                     CublasMpResources& cublas_res);

template <typename Device>
void norm_psi_tensor_lapack(const int nband,
                            const int nlocal,
                            const ct::Tensor& Stmp,
                            ct::Tensor& psi_k,
                            std::ofstream& ofs_running,
                            const int print_matrix);

#endif // __MPI
} // namespace module_rt

#endif
