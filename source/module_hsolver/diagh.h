#ifndef DIAGH_H
#define DIAGH_H

#include <string>

#include "module_base/complexmatrix.h"
#include "module_base/macros.h"
#include "module_hamilt_general/hamilt.h"
#include "module_psi/psi.h"

namespace hsolver
{

template <typename T, typename Device = psi::DEVICE_CPU>
class DiagH
{
  private:
    using Real = typename GetTypeReal<T>::type;
  public:
    virtual ~DiagH(){};
    // virtual void init()=0;
    std::string method = "none";

    virtual void diag(hamilt::Hamilt<T, Device> *phm_in, psi::Psi<T, Device> &psi, Real *eigenvalue_in) = 0;

    virtual void diag(hamilt::Hamilt<T, Device> *phm_in, psi::Psi<Real, Device> &psi, Real *eigenvalue_in) {
        return;
    }
};

} // namespace hsolver

#endif