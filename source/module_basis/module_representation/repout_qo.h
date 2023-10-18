#ifndef REPOUT_QO_H
#define REPOUT_QO_H

#include "module_basis/module_representation/repout.h"
/*
    RepOut_QO: output wavefunction in QO representation (from pw representation)
    
*/
template<typename T, typename Device>
class RepOut_QO : public RepOut<T, Device>
{
    private:
        using Real = typename GetTypeReal<T>::type;
    public:
        RepOut_QO();
        ~RepOut_QO();
        void project(psi::Psi<T, Device>* psi_in,
                     psi::Psi<T, Device>* psig,
                     psi::Psi<T, Device>* psi_out) override;
};

#endif // REPOUT_QO_H