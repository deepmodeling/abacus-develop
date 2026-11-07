/**
 * @file band_energy.h
 * @brief compute band energy ekb
 *  This file originally belonged to file LCAO_evolve.cpp
 */
#ifndef BANDENERGY_H
#define BANDENERGY_H

#include "source_base/module_container/ATen/core/tensor.h" // ct::Tensor
#include "source_io/module_energy/td_eij_io.h"
#include "source_lcao/module_rt/kernels/cublasmp_context.h"

#include <complex>

class BlockCyclicDistribution;

namespace module_rt
{
using ModuleIO::EijOutputInfo;

#ifdef __MPI
/**
 *  @brief compute band energy ekb <psi_i|H|psi_i>
 *
 * @param[in] matrix_layout block-cyclic distribution of the square matrix (H)
 * @param[in] wfc_layout block-cyclic distribution of the wavefunction (nbasis x nbands)
 * @param[in] mo_repr_layout block-cyclic layout of the MO representation of H
 * @param[in] nband number of bands
 * @param[in] nlocal number of orbitals
 * @param[in] Htmp Hamiltonian
 * @param[in] psi_k psi of this step
 * @param[out] ekb band energy
 */
void compute_ekb(const BlockCyclicDistribution& matrix_layout,
                 const BlockCyclicDistribution& wfc_layout,
                 const BlockCyclicDistribution& mo_repr_layout,
                 const int nband,
                 const int nlocal,
                 const std::complex<double>* Htmp,
                 const std::complex<double>* psi_k,
                 double* ekb,
                 std::ofstream& ofs_running,
                 const EijOutputInfo& output_info);

void compute_ekb_tensor(const BlockCyclicDistribution& matrix_layout,
                        const BlockCyclicDistribution& wfc_layout,
                        const BlockCyclicDistribution& mo_repr_layout,
                        const int nband,
                        const int nlocal,
                        const ct::Tensor& Htmp,
                        const ct::Tensor& psi_k,
                        ct::Tensor& ekb,
                        std::ofstream& ofs_running,
                        CublasMpResources& cublas_res);

template <typename Device>
void compute_ekb_tensor_lapack(const int nband,
                               const int nlocal,
                               const ct::Tensor& Htmp,
                               const ct::Tensor& psi_k,
                               ct::Tensor& ekb,
                               std::ofstream& ofs_running,
                               const EijOutputInfo& output_info);
#endif // __MPI
} // namespace module_rt
#endif
