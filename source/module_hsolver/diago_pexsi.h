#ifndef DIGAOPEXSI_H
#define DIGAOPEXSI_H

#ifdef  __PEXSI

#include "diagh.h"

namespace hsolver
{

class DiagoPexsi : public DiagH<double>
{
  public:
    void diag(hamilt::Hamilt<double>* phm_in, psi::Psi<double>& psi, double* eigenvalue_in) override;
    void diag(hamilt::Hamilt<double>* phm_in, psi::Psi<std::complex<double>>& psi, double* eigenvalue_in) override;
};

}

#endif
#endif
