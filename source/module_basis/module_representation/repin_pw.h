#ifndef REPIN_PW_H
#define REPIN_PW_H

#include "module_basis/module_representation/repin.h"

template<typename T, typename Device>
class RepIn_PW : public RepIn<T, Device>
{
    private:
        using Real = typename GetTypeReal<T>::type;
    public:
        RepIn_PW();
        ~RepIn_PW();
        void cal_psig(psi::Psi<T, Device>* psig) override;
};

#endif // REPIN_PW_H