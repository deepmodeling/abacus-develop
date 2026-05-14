#ifndef DIAGO_CUSOLVERMPH
#define DIAGO_CUSOLVERMPH

#ifdef __CUSOLVERMP
#include "source_hamilt/hamilt.h"
#include "source_base/macros.h"
namespace hsolver
{
// DiagoCusolverMP class, for diagonalization using CUSOLVERMP
template <typename T>
class DiagoCusolverMP
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    DiagoCusolverMP()
    {
    }
    // the diag function for CUSOLVERMP diagonalization
    void diag(hamilt::Hamilt<T>* phm_in, psi::Psi<T>& psi, Real* eigenvalue_in);
};
} // namespace hsolver
#endif // __CUSOLVERMP
#endif // DIAGO_CUSOLVERMPH
