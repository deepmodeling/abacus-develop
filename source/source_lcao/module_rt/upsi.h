/**
 * @file upsi.h
 * @brief apply U_operator to the wave function of the previous step for new wave function
 *  This file originally belonged to file LCAO_evolve.cpp
 */

#ifndef UPSI_H
#define UPSI_H

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
 *  @brief apply U_operator to the wave function of the previous step for new wave function
 *
 * All views must use compatible process grids and block sizes.
 *
 * @param[in] propagator Time-evolution operator.
 * @param[in] previous_psi Wave function at the previous step.
 * @param[out] psi Wave function at this step.
 * @param[in,out] ofs_running Running log.
 * @param[in] print_matrix Whether to print intermediate matrices.
 */
void upsi(Matrix2dView<const std::complex<double>> propagator,
          Matrix2dView<const std::complex<double>> previous_psi,
          Matrix2dView<std::complex<double>> psi,
          std::ofstream& ofs_running,
          const int print_matrix);

void upsi_tensor(const Parallel_Blkcyclic& matrix_layout,
                 const Parallel_Blkcyclic& wfc_layout,
                 const int nband,
                 const int nlocal,
                 const ct::Tensor& U_operator,
                 const ct::Tensor& psi_k_laststep,
                 ct::Tensor& psi_k,
                 std::ofstream& ofs_running,
                 const int print_matrix,
                 CublasMpResources& cublas_res);

template <typename Device>
void upsi_tensor_lapack(const int nband,
                        const int nlocal,
                        const ct::Tensor& U_operator,
                        const ct::Tensor& psi_k_laststep,
                        ct::Tensor& psi_k,
                        std::ofstream& ofs_running,
                        const int print_matrix);

#endif // __MPI
} // namespace module_rt

#endif
