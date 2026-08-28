#ifndef DFTU_OCCUP_H
#define DFTU_OCCUP_H

#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include "source_hamilt/hamilt.h"

#include <vector>

class Plus_U;

#ifdef __LCAO
namespace DFTU_LCAO {

/// @brief Compute the occupation matrix and delegate to Plus_U member.
/// Dispatches to Plus_U::cal_occ_mat_gamma (gamma-only, double) or
/// Plus_U::cal_occ_mat_k (multi-k, std::complex<double>) via template.
template <typename T>
void cal_occ_mat(const int iter,
                 const UnitCell& ucell,
                 const std::vector<std::vector<T>>& dm,
                 const K_Vectors& kv,
                 const double& mixing_beta,
                 hamilt::Hamilt<T>* p_ham,
                 Plus_U& dftu);

} // namespace DFTU_LCAO
#endif

#endif
