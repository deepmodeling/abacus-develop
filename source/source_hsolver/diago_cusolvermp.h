#ifndef DIAGO_CUSOLVERMPH
#define DIAGO_CUSOLVERMPH

#ifdef __CUSOLVERMP
#include "source_hamilt/hamilt.h"
#include "source_base/macros.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_hsolver/kernels/cuda/diag_cusolvermp.cuh"
namespace hsolver
{
// DiagoCusolverMP class, for diagonalization using CUSOLVERMP
template <typename T>
class DiagoCusolverMP
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    /// @param nlocal_in global dimension of the NAO Hamiltonian
    /// @param nbands_in number of lowest eigenpairs to compute
    DiagoCusolverMP(const int nlocal_in, const int nbands_in) : nlocal(nlocal_in), nbands(nbands_in)
    {
    }
    // the diag function for CUSOLVERMP diagonalization
    void diag(hamilt::Hamilt<T>* phm_in, psi::Psi<T>& psi, Real* eigenvalue_in);

  private:
    const int nlocal;
    const int nbands;
};
} // namespace hsolver
#endif // __CUSOLVERMP
#endif // DIAGO_CUSOLVERMPH
