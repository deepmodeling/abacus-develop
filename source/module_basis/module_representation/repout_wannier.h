#ifndef REPIN_WANNIER_H
#define REPIN_WANNIER_H

#include "module_basis/module_representation/repout.h"

template<typename T, typename Device>
class RepOut_Wannier : public RepOut<T, Device>
{
    public:
        RepOut_Wannier();
        ~RepOut_Wannier();
        void project(const psi::Psi<T, Device>& psi_in,
                     const psi::Psi<T, Device>* psig,
                     psi::Psi<T, Device>& psi_out) override;
};

#endif // REPIN_WANNIER_H