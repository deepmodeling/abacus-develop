//=====================
// AUTHOR : Peize Lin
// DATE : 2021-11-02
// REFACTORING AUTHOR : Daye Zheng
// DATE : 2022-04-14
//=====================

#ifndef DIAGO_SCALAPACK_H
#define DIAGO_SCALAPACK_H

#include <complex>
#include <utility>
#include <vector>

#include "diagh.h"
#include "module_base/complexmatrix.h"
#include "module_base/matrix.h"
#include "module_basis/module_ao/parallel_orbitals.h"

namespace hsolver
{
    template<typename T>
    class DiagoScalapack
{
private:
    using Real = typename GetTypeReal<T>::type;
  public:
    void diag(hamilt::Hamilt<T>* phm_in, psi::Psi<T>& psi, Real* eigenvalue_in);
#ifdef __MPI
    // diagnolization used in parallel-k case
    void diag_pool(hamilt::MatrixBlock<T>& h_mat, hamilt::MatrixBlock<T>& s_mat, psi::Psi<T>& psi, Real* eigenvalue_in, MPI_Comm& comm);
#endif

  private:
    void pdsygvx_diag(const int *const desc,
                      const int ncol,
                      const int nrow,
                      const double *const h_mat,
                      const double *const s_mat,
                      double *const ekb,
                      psi::Psi<double> &wfc_2d);
    void pzhegvx_diag(const int *const desc,
                      const int ncol,
                      const int nrow,
                      const std::complex<double> *const h_mat,
                      const std::complex<double> *const s_mat,
                      double *const ekb,
                      psi::Psi<std::complex<double>> &wfc_2d);

    std::pair<int, std::vector<int>> pdsygvx_once(const int *const desc,
                                                  const int ncol,
                                                  const int nrow,
                                                  const double *const h_mat,
                                                  const double *const s_mat,
                                                  double *const ekb,
                                                  psi::Psi<double> &wfc_2d) const;
    std::pair<int, std::vector<int>> pzhegvx_once(const int *const desc,
                                                  const int ncol,
                                                  const int nrow,
                                                  const std::complex<double> *const h_mat,
                                                  const std::complex<double> *const s_mat,
                                                  double *const ekb,
                                                  psi::Psi<std::complex<double>> &wfc_2d) const;

    int degeneracy_max = 12; // For reorthogonalized memory. 12 followes siesta.

    void post_processing(const int info, const std::vector<int> &vec);
};

} // namespace hsolver

#endif