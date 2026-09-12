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

class Plus_U_Base;
class OccupationMatrix;

namespace DFTU_LCAO {

/// @brief Compute the occupation matrix and delegate to Plus_U member.
/// Dispatches to DFTU_LCAO::cal_occ_mat_gamma (gamma-only, double) or
/// DFTU_LCAO::cal_occ_mat_k (multi-k, std::complex<double>) via template.
/// @param pv parallel-orbitals descriptor that owns BLACS context and
///        global<->local index maps; sourced by the caller from the same
///        Parallel_Orbitals used to build the Hamiltonian and density matrix.
/// @param ks_solver KS solver name (e.g. "scalapack"); forwarded to the
///        multi-k path for folding-matrix selection.
template <typename T>
void cal_occ_mat(const Parallel_Orbitals* pv,
                 const UnitCell& ucell,
                 const std::vector<std::vector<T>>& dm,
                 const K_Vectors& kv,
                 const double& mixing_beta,
                 hamilt::Hamilt<T>* p_ham,
                 Plus_U_Base& dftu,
                 const bool gamma_only_local,
                 const int nspin,
                 const std::string& ks_solver);

/// @brief Accumulate one (iat, l, n, spin) channel of the occupation matrix
///        from the complex S*DM product srho for the multi-k case. Reads npol
///        and the iatlnmipol2iwt lookup directly from occmat so callers do
///        not need to thread those scalars through.
void accumulate_occ_channel_k(OccupationMatrix& occmat,
                              const Parallel_Orbitals& pv,
                              const std::complex<double>* srho,
                              int iat,
                              int l,
                              int n,
                              int spin);

/// @brief Accumulate one (iat, l, n, spin) channel of the occupation matrix
///        from the real S*DM product srho for the gamma-only case. Reads npol
///        and the iatlnmipol2iwt lookup directly from occmat so callers do
///        not need to thread those scalars through.
void accumulate_occ_channel_gamma(OccupationMatrix& occmat,
                                  const Parallel_Orbitals& pv,
                                  const double* srho,
                                  int iat,
                                  int l,
                                  int n,
                                  int spin);

/// @brief MPI Allreduce each (iat, l, n=0) channel of occmat across all ranks
///        and symmetrize it (Hermitian average) per the nspin convention:
///        nspin=1 mirrors spin-0 into spin-1; nspin=2 symmetrizes each spin;
///        nspin=4 symmetrizes the single Pauli block. Reads nspin and npol
///        from occmat so callers do not thread them through.
void reduce_and_symmetrize_occ_k(OccupationMatrix& occmat,
                                 const UnitCell& ucell,
                                 const std::vector<int>& l_channel);

/// @brief Walk the (it, ia, l, n=0) atom mesh for one k-point and accumulate
///        each qualifying channel of occmat from the complex S*DM product
///        srho. Reads npol and the iatlnmipol2iwt lookup from occmat so
///        callers do not thread them through.
void accumulate_occ_k_for_ik(OccupationMatrix& occmat,
                             const UnitCell& ucell,
                             const Parallel_Orbitals& pv,
                             const std::complex<double>* srho,
                             int spin,
                             const std::vector<int>& l_channel);

/// @brief Process one (it, ia, l, n=0, spin) block of the gamma-only
///        occupation matrix: accumulate from the real S*DM product srho,
///        MPI-Allreduce across ranks, then symmetrize per the nspin
///        convention. Reads nspin and npol from occmat so callers do not
///        thread them through.
void process_occ_channel_gamma(OccupationMatrix& occmat,
                               const UnitCell& ucell,
                               const Parallel_Orbitals& pv,
                               const double* srho,
                               int spin,
                               const std::vector<int>& l_channel);

// calculate the local occupation number matrix (k-point version)
void cal_occ_mat_k(const Parallel_Orbitals* pv,
                   const UnitCell& ucell,
                   const std::vector<std::vector<std::complex<double>>>& dm_k,
                   const K_Vectors& kv,
                   const double& mixing_beta,
                   hamilt::Hamilt<std::complex<double>>* p_ham,
                   const bool gamma_only_local,
                   Plus_U_Base& dftu,
                   const std::string& ks_solver);

// calculate the local occupation number matrix (gamma-point version)
void cal_occ_mat_gamma(const Parallel_Orbitals* pv,
                       const UnitCell& ucell,
                       const std::vector<std::vector<double>>& dm_gamma,
                       const double& mixing_beta,
                       hamilt::Hamilt<double>* p_ham,
                       Plus_U_Base& dftu);

} // namespace DFTU_LCAO

#endif
