#ifndef DIAGH_H
#define DIAGH_H

#include "module_base/complexmatrix.h"
#include "module_hamilt/hamilt.h"
#include "module_psi/psi.h"
#include "string"

namespace hsolver
{

class DiagH
{
  public:
    virtual ~DiagH(){};
    // virtual void init()=0;
    std::string method = "none";

    virtual void diag(hamilt::Hamilt *phm_in, psi::Psi<std::complex<double>> &psi, double *eigenvalue_in) = 0;

    virtual void diag(hamilt::Hamilt *phm_in, psi::Psi<double> &psi, double *eigenvalue_in)
    {
        return;
    }
};

} // namespace hsolver

#endif