#ifndef HSOLVER_DIAGO_PXXXGVX_H
#define HSOLVER_DIAGO_PXXXGVX_H
#include <complex>
#include "module_base/macros.h"

namespace hsolver
{

#ifdef __MPI
/**
 * @brief Wrapper function for Scalapack's generalized eigensolver routines: pdsygvx_, pzhegvx_, pssygvx_, pchegvx_.
 * 
 * @param desc the descriptor of scalapack descriptor
 * @param ncol the number of columns of the H/S matrix in current processor
 * @param nrow the number of rows of the H/S matrix in current processor
 * @param nbands the number of bands to be solved
 * @param h_mat the Hamiltonian matrix
 * @param s_mat the overlap matrix
 * @param ekb the eigenvalues
 * @param wfc_2d the eigenvectors in 2D block cyclic distribution
 * 
 */

template <typename T>
void pxxxgvx_diag(const int* const desc,
                  const int ncol,
                  const int nrow,
                  const int nbands,
                  const T* const h_mat,
                  const T* const s_mat,
                  typename GetTypeReal<T>::type* const ekb,
                  T* const wfc_2d);
#endif 

} // namespace hsolver

#endif // HSOLVER_DIAGO_PXXXGVX_H