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

/// @brief Compute the occupation matrix
///        occ(m,m') = sum_R DMR(I,J,R) * <phi_0|chi_m(I)> * <chi_m'(J)|phi_R>
///        and delegate to the Plus_U member.
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

// calculate the local occupation number matrix (k-point version)
//
// @note Hard to unit-test: requires a full wavefunction (Psi), two-center
// integrator, and PSI-to-2D distribution (p2s_dist) to build srho.
// Consider extracting the srho computation into an injectable interface
// if unit-test coverage is needed.
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
//
// @note Hard to unit-test: requires a full wavefunction (Psi) and
// two-center integrator to build srho. Consider extracting the srho
// computation into an injectable interface if unit-test coverage is needed.
void cal_occ_mat_gamma(const Parallel_Orbitals* pv,
                       const UnitCell& ucell,
                       const std::vector<std::vector<double>>& dm_gamma,
                       const double& mixing_beta,
                       hamilt::Hamilt<double>* p_ham,
                       Plus_U_Base& dftu);

} // namespace DFTU_LCAO

#endif
