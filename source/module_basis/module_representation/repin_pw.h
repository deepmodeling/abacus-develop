#ifndef REPIN_PW_H
#define REPIN_PW_H

#include "module_basis/module_representation/repin.h"

template<typename T, typename Device>
class RepIn_PW : public RepIn<T, Device>
{
    public:
        RepIn_PW();
        ~RepIn_PW();
        void project_in(const psi::Psi<T, Device>& psi_in, 
                                psi::Psi<T, Device>& psi_out) override;
};

#endif // REPIN_PW_H