#ifndef DIAGOCUSOLVER_H
#define DIAGOCUSOLVER_H

#include "diagh.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_hsolver/kernels/cuda/diag_cusolver.cuh"
// #include "module_hsolver/kernels/cuda/dngvd_op.cu"

namespace hsolver
{

// DiagoCusolver class, derived from DiagH, for diagonalization using CUSOLVER
template <typename T>
class DiagoCusolver : public DiagH<T>
{
  private:
    // Real is the real part of the complex type T
    using Real = typename GetTypeReal<T>::type;
    Parallel_Orbitals const * ParaV;

  public:

    DiagoCusolver(const Parallel_Orbitals* ParaV = nullptr);
    ~DiagoCusolver();
    
    // Override the diag function for CUSOLVER diagonalization
    void diag(hamilt::Hamilt<T>* phm_in, psi::Psi<T>& psi, Real* eigenvalue_in) override;

    // Static variable to keep track of the decomposition state
    static int DecomposedState;

    // Diag_Cusolver_gvd object for CUSOLVER operations
    Diag_Cusolver_gvd dc;

  private:
#ifdef __MPI
    // Function to check if ELPA handle needs to be created or reused in MPI settings
    bool ifElpaHandle(const bool& newIteration, const bool& ifNSCF);
#endif
};

} // namespace hsolver

#endif
