#ifndef DIAGH_H
#define DIAGH_H

#include "module_base/complexmatrix.h"
#include "module_hamilt/hamilt.h"
#include "module_psi/psi.h"
#include "string"

namespace hsolver
{

template <typename FPTYPE, typename Device = psi::DEVICE_CPU>
class DiagH
{
  public:
    virtual ~DiagH(){};
    // virtual void init()=0;
    std::string method = "none";

    virtual void diag(hamilt::Hamilt<FPTYPE, Device> *phm_in, psi::Psi<std::complex<FPTYPE>, Device> &psi, FPTYPE *eigenvalue_in) = 0;

    virtual void diag(hamilt::Hamilt<FPTYPE, Device> *phm_in, psi::Psi<FPTYPE, Device> &psi, FPTYPE *eigenvalue_in) {
        return;
    }
};

} // namespace hsolver

#endif