#ifndef REPOUT_PW_H
#define REPOUT_PW_H

#include "module_basis/module_representation/repout.h"
#include "module_base/macros.h"
/*
    RepOut_PW: output representation of plane wave basis
    (note that in class Representation, the intermediate representation is pw, so 
    this class is the simplest case that only need to directly output, copy psig to psi_out,
    but needed to multiply by psi_in.

    For psi initialization, the psi_in would be identity matrix,
    For lcao to transform to pw, the psi_in would be the one in lcao basis,
    )
*/
template<typename T, typename Device>
class RepOut_PW : public RepOut<T, Device>
{
    private:
        using Real = GetTypeReal<T>::type;
        constexpr static const Device * ctx = {};
    public:
        RepOut_PW();
        ~RepOut_PW();
        void project(const psi::Psi<T, Device>& psi_in,
                     const psi::Psi<T, Device>* psig,
                     psi::Psi<T, Device>& psi_out) override;
};

#endif // REPOUT_PW_H