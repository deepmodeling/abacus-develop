#ifndef DIAGOELPA_H
#define DIAGOELPA_H

#include "diagh.h"
#include "module_basis/module_ao/parallel_orbitals.h"

namespace hsolver
{

class DiagoElpa : public DiagH<std::complex<double>>
{

  public:
    void diag(hamilt::Hamilt<std::complex<double>>* phm_in, psi::Psi<double>& psi, double* eigenvalue_in) override;

    void diag(hamilt::Hamilt<std::complex<double>>* phm_in, psi::Psi<std::complex<double>>& psi, double* eigenvalue_in) override;
    
    static int DecomposedState;

  private:
#ifdef __MPI
    bool ifElpaHandle(const bool& newIteration, const bool& ifNSCF);
#endif
};

} // namespace hsolver

#endif
