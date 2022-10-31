// TODO: This is a temperary location for these functions.
// And will be moved to a global module(module base) later.
#ifndef MODULE_HSOLVER_DNGVD_H
#define MODULE_HSOLVER_DNGVD_H

#include "module_base/lapack_connector.h"
#include "module_psi/include/memory.h"


namespace hsolver
{

    

template <typename FPTYPE, typename Device>
struct dngvd_op
{
    /// @brief dngvd computes selected eigenvalues, and optionally, eigenvectors of
    /// a complex generalized Hermitian-definite eigenproblem
    ///
    /// Input Parameters
    ///     @param d : the type of device
    ///     @param row : the number of rows of the matrix 
    ///     @param col : the number of cols of the matrix  
    ///     @param A : the hermitian matrix A in A x=lambda B x (row major) 
    ///     @param B : the SPD matrix B in A x=lambda B x (row major) 
    ///     @param m : the number of the first m eigenvalues to calculate
    /// Output Parameter
    ///     @param W : calculated eigenvalues
    ///     @param V : calculated eigenvectors (row major)
    void operator()(
        const Device* d,
        const int row,
        const int col,
        const std::complex<FPTYPE>* A,
        const std::complex<FPTYPE>* B,
        const int m,
        double* W, 
        std::complex<FPTYPE>* V);
};



} // namespace hsolver

#endif // !MODULE_HSOLVER_DNGVD_H