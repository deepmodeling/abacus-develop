#ifndef REPOUT_PW_H
#define REPOUT_PW_H

#include "module_basis/module_representation/repout.h"
#include "module_base/macros.h"
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