#ifndef REPOUT_H
#define REPOUT_H

// data structure support
#include <string>
#include "module_psi/psi.h"

template<typename T, typename Device>
class RepOut
{
    public:
        RepOut();
        ~RepOut();
        virtual void project(const psi::Psi<T, Device>& psi_in,
                             const psi::Psi<T, Device>* psig,
                                   psi::Psi<T, Device>& psi_out) = 0;
        void set_rep_from(const std::string rep_from_in) { rep_from = rep_from_in; }
    protected:
        std::string rep_from = "";
        T one = static_cast<T>(1.0);
        T zero = static_cast<T>(0.0);
};

#endif // REPOUT_H