//=====================
// REFACTORING AUTHOR : Xiaoyang Zhang
// DATE : 2024-6-24
//=====================

// This is fully refactored according to original diago_scalapack

#ifndef DIAGOLAPACK_H
#define DIAGOLAPACK_H

#include "module_base/macros.h"   // GetRealType
#include "module_hamilt_general/hamilt.h"
#include "module_base/matrix.h"
#include "module_basis/module_ao/parallel_orbitals.h"

#include <complex>
#include <utility>
#include <vector>

namespace hsolver
{
template <typename T>
class DiagoLapack
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    void diag(hamilt::Hamilt<T>* phm_in, psi::Psi<T>& psi, Real* eigenvalue_in);

    void dsygvx_diag(const int ncol,
                     const int nrow,
                     const double* const h_mat,
                     const double* const s_mat,
                     double* const ekb,
                     psi::Psi<double>& wfc_2d);
    void zhegvx_diag(const int ncol,
                     const int nrow,
                     const std::complex<double>* const h_mat,
                     const std::complex<double>* const s_mat,
                     double* const ekb,
                     psi::Psi<std::complex<double>>& wfc_2d);

    int dsygvx_once(const int ncol,
                    const int nrow,
                    const double* const h_mat,
                    const double* const s_mat,
                    double* const ekb,
                    psi::Psi<double>& wfc_2d) const;
    int zhegvx_once(const int ncol,
                    const int nrow,
                    const std::complex<double>* const h_mat,
                    const std::complex<double>* const s_mat,
                    double* const ekb,
                    psi::Psi<std::complex<double>>& wfc_2d) const;

    int degeneracy_max = 12; // For reorthogonalized memory. 12 followes siesta.

    void post_processing(const int info, const std::vector<int>& vec);
};

} // namespace hsolver

#endif