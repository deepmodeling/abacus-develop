#ifndef REPOUT_QO_H
#define REPOUT_QO_H

#include "module_basis/module_representation/repout.h"

template<typename T, typename Device>
class RepOut_QO : public RepOut<T, Device>
{
    public:
        RepOut_QO();
        ~RepOut_QO();
        void project(const psi::Psi<T, Device>& psi_in,
                     const psi::Psi<T, Device>* psig,
                     psi::Psi<T, Device>& psi_out) override;
};

#endif // REPOUT_QO_H