#ifndef DFTU_LCAO_OCC_H
#define DFTU_LCAO_OCC_H

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_base/matrix.h"
#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include "source_hamilt/hamilt.h"

#include <complex>
#include <string>
#include <vector>

class Plus_U;

#ifdef __LCAO
namespace DFTU_LCAO {

/// @brief Compute the occupation matrix and delegate to Plus_U member.
/// Dispatches to DFTU_LCAO::cal_occ_mat_gamma (gamma-only, double) or
/// DFTU_LCAO::cal_occ_mat_k (multi-k, std::complex<double>) via template.
/// @param pv parallel-orbitals descriptor that owns BLACS context and
///        global<->local index maps; sourced by the caller from the same
///        Parallel_Orbitals used to build the Hamiltonian and density matrix.
template <typename T>
void cal_occ_mat(const Parallel_Orbitals* pv,
                 const int iter,
                 const UnitCell& ucell,
                 const std::vector<std::vector<T>>& dm,
                 const K_Vectors& kv,
                 const double& mixing_beta,
                 hamilt::Hamilt<T>* p_ham,
                 Plus_U& dftu,
                 const bool gamma_only_local,
                 const int nspin);

// calculate the local occupation number matrix (k-point version)
void cal_occ_mat_k(const Parallel_Orbitals* pv,
                   const int iter,
                   const UnitCell& ucell,
                   const std::vector<std::vector<std::complex<double>>>& dm_k,
                   const K_Vectors& kv,
                   const double& mixing_beta,
                   hamilt::Hamilt<std::complex<double>>* p_ham,
                   const bool gamma_only_local,
                   Plus_U& dftu);

// calculate the local occupation number matrix (gamma-point version)
void cal_occ_mat_gamma(const Parallel_Orbitals* pv,
                       const int iter,
                       const UnitCell& ucell,
                       const std::vector<std::vector<double>>& dm_gamma,
                       const double& mixing_beta,
                       hamilt::Hamilt<double>* p_ham,
                       Plus_U& dftu);

} // namespace DFTU_LCAO
#endif

#endif
