#ifndef DIAGOLAPACK_H
#define DIAGOLAPACK_H

#include "diagh.h"
#include "module_basis/module_ao/parallel_orbitals.h"

namespace hsolver
{

class DiagoLapack : public DiagH<std::complex<double>>
{

  public:
    void diag(hamilt::Hamilt<std::complex<double>>* phm_in, psi::Psi<double>& psi, double* eigenvalue_in) override;

    void diag(hamilt::Hamilt<std::complex<double>>* phm_in, psi::Psi<std::complex<double>>& psi, double* eigenvalue_in) override;
    
  private:
};

} // namespace hsolver

#endif