#ifndef DIAGOLAPACK_H
#define DIAGOLAPACK_H

#include "diagh.h"
#include "module_orbital/parallel_orbitals.h"

namespace hsolver
{

class DiagoLapack : public DiagH
{

  public:
    void diag(hamilt::Hamilt* phm_in, psi::Psi<double>& psi, double* eigenvalue_in) override;

    void diag(hamilt::Hamilt* phm_in, psi::Psi<std::complex<double>>& psi, double* eigenvalue_in) override;
    
  private:
};

} // namespace hsolver

#endif