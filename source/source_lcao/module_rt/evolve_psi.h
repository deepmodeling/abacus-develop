/**
 * @file evolve_psi.h
 * @brief evolve the wave function
 *  This file originally belonged to file LCAO_evolve.cpp
 */
#ifndef ELEC_PSI_H
#define ELEC_PSI_H

#include "source_base/module_container/ATen/core/tensor.h"     // ct::Tensor
#include "source_base/module_container/ATen/core/tensor_map.h" // TensorMap
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_lcao/hamilt_lcao.h"

namespace module_rt
{
void evolve_psi(const int nband,
                const int nlocal,
                const Parallel_Orbitals* pv,
                hamilt::Hamilt<std::complex<double>>* p_hamilt,
                std::complex<double>* psi_k,
                std::complex<double>* psi_k_laststep,
                std::complex<double>* H_laststep,
                std::complex<double>* S_laststep,
                double* ekb,
                int htype,
                int propagator,
                std::ofstream& ofs_running,
                const int print_matrix);

template <typename Device>
void evolve_psi_tensor(const int nband,
                       const int nlocal,
                       const Parallel_Orbitals* pv,
                       hamilt::Hamilt<std::complex<double>>* p_hamilt,
                       ct::Tensor& psi_k,
                       ct::Tensor& psi_k_laststep,
                       ct::Tensor& H_laststep,
                       ct::Tensor& S_laststep,
                       ct::Tensor& ekb,
                       int htype,
                       int propagator,
                       std::ofstream& ofs_running,
                       const int print_matrix,
                       const bool use_lapack);
} // namespace module_rt

#endif