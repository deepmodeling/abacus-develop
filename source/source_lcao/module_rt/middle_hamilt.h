/**
 * @file middle_hamilt.h
 * @brief  compute H(t+dt/2)
 *  This file originally belonged to file LCAO_evolve.cpp
 */
#ifndef MIDDLE_HAMILT_H
#define MIDDLE_HAMILT_H

#include "source_base/module_container/ATen/core/tensor.h" // ct::Tensor
#include "source_lcao/module_rt/kernels/cublasmp_context.h"

#include <complex>
#include <fstream>

class Parallel_Blkcyclic; // Defined in source/source_base/parallel_blkcyc.h

template <typename T>
class Matrix2dView;

namespace module_rt
{
#ifdef __MPI
/**
 *  @brief compute H(t+dt/2)
 *
 * The four matrix views must use the same host or distributed layout.
 *
 * @param[in] previous_hamiltonian H(t).
 * @param[in] previous_overlap S(t).
 * @param[in,out] hamiltonian H(t+dt) on input, H(t+dt/2) on return.
 * @param[in,out] overlap S(t+dt) on input, S(t+dt/2) on return.
 * @param[in,out] ofs_running Running log.
 * @param[in] print_matrix Whether to print intermediate matrices.
 */
void half_Hmatrix(
                  Matrix2dView<const std::complex<double>> previous_hamiltonian,
                  Matrix2dView<const std::complex<double>> previous_overlap,
                  Matrix2dView<std::complex<double>> hamiltonian,
                  Matrix2dView<std::complex<double>> overlap,
                  std::ofstream& ofs_running,
                  const int print_matrix);

void half_Hmatrix_tensor(const Parallel_Blkcyclic& layout,
                         const int nband,
                         const int nlocal,
                         ct::Tensor& Htmp,
                         ct::Tensor& Stmp,
                         const ct::Tensor& H_laststep,
                         const ct::Tensor& S_laststep,
                         std::ofstream& ofs_running,
                         const int print_matrix,
                         CublasMpResources& cublas_res);

template <typename Device>
void half_Hmatrix_tensor_lapack(const int nband,
                                const int nlocal,
                                ct::Tensor& Htmp,
                                ct::Tensor& Stmp,
                                const ct::Tensor& H_laststep,
                                const ct::Tensor& S_laststep,
                                std::ofstream& ofs_running,
                                const int print_matrix);

#endif // __MPI
} // namespace module_rt

#endif
