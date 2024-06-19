#ifndef DIAGO_CUSOLVERMPH
#define DIAGO_CUSOLVERMPH

#ifdef __CUSOLVERMP
#include "diagh.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_hsolver/kernels/cuda/diag_cusolvermp.cuh"
#include "module_base/macros.h"
namespace hsolver
{
// DiagoCusolver class, derived from DiagH, for diagonalization using CUSOLVER
template <typename T>
class DiagoCusolverMP : public DiagH<T>
{
  private:
    using Real = typename GetTypeReal<T>::type;
  public:
    DiagoCusolverMP()
    {

    }
    // Override the diag function for CUSOLVER diagonalization
    void diag(hamilt::Hamilt<T>* phm_in, psi::Psi<T>& psi, Real* eigenvalue_in) override;
};

} // namespace hsolver
#endif // __CUSOLVERMP
#endif // DIAGO_CUSOLVERMPH
