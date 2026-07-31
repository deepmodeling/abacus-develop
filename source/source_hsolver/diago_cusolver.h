#ifndef DIAGOCUSOLVER_H
#define DIAGOCUSOLVER_H

#include "source_base/macros.h"   // GetRealType
#include "source_hamilt/hamilt.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_hsolver/kernels/cuda/diag_cusolver.cuh"

namespace hsolver
{

// DiagoCusolver class for diagonalization using CUSOLVER
template <typename T>
class DiagoCusolver
{
  private:
    // Real is the real part of the complex type T
    using Real = typename GetTypeReal<T>::type;

  public:

    /// @param nlocal_in global dimension of the NAO Hamiltonian
    /// @param nbands_in number of lowest eigenpairs to compute
    DiagoCusolver(const int nlocal_in, const int nbands_in);
    ~DiagoCusolver();

    // Override the diag function for CUSOLVER diagonalization
    void diag(
      hamilt::MatrixBlock<T>& h_mat,
      hamilt::MatrixBlock<T>& s_mat,
      psi::Psi<T>& psi,
      Real* eigenvalue_in);

    // Static variable to keep track of the decomposition state
    static int DecomposedState;

    // Diag_Cusolver_gvd object for CUSOLVER operations
    Diag_Cusolver_gvd dc;

  private:
#ifdef __MPI
    // Function to check if ELPA handle needs to be created or reused in MPI settings
    bool ifElpaHandle(const bool& newIteration, const bool& ifNSCF) const;
#endif

    const int nlocal;
    const int nbands;
};

} // namespace hsolver

#endif
