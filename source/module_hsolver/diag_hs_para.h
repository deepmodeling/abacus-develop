#include "module_base/parallel_2d.h"
#include "module_base/macros.h"

#ifdef __MPI
#include <mpi.h>
#endif

namespace hsolver
{


#ifdef __MPI

/**
 * @brief Parallel do the generalized eigenvalue problem
 * 
 * @tparam T double or complex<double> or float or complex<float>
 * @param H the hermitian matrix H.
 * @param S the overlap matrix S.
 * @param lda the leading dimension of H and S
 * @param nband the number of bands to be calculated
 * @param ekb to store the eigenvalues.
 * @param wfc to store the eigenvectors
 * @param comm the communicator
 * @param diag_subspace the method to solve the generalized eigenvalue problem
 * @param block_size the block size in 2d block cyclic distribution if use elpa or scalapack. 
 * 
 * @note 1. h and s should be full matrix in rank 0 of the communicator, and the other ranks is not concerned.
 * @note 2. wfc is complete in rank 0, and not store in other ranks.
 * @note 3. diag_subspace should be 1: by elpa, 2: by scalapack
 * @note 4. block_size should be 0 or a positive integer. If it is 0, then will use a value as large as possible that is allowed
 */
template <typename T>
void diago_hs_para(T* h,
                   T* s,
                   const int lda,
                   const int nband,
                   typename GetTypeReal<T>::type* const ekb,
                   T* const wfc,
                   const MPI_Comm& comm,
                   const int diag_subspace,
                   const int block_size = 0);
#endif

} // namespace hsolver                 
              