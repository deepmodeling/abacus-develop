#ifndef REPOUT_H
#define REPOUT_H

// data structure support
#include <string>
#include "module_psi/psi.h"
/*
    RepOut: output wavefunction from one representation via pw to another
    (Interface to various postprocessing software which supports pw basis)

    project: project psi_in to psi_out
             this function should do two things:
             1. calculate pw basis in output representation, get transformation matrix pw->output
             2. multiply psi_in with psig, get bands in pw representation, state->pw
             3. multiply bands in pw representation with pw basis in output representation, get psi_out,
                state->output = state->pw * pw->output
*/
template<typename T, typename Device>
class RepOut
{
    public:
        RepOut();
        ~RepOut();
        virtual void project(const psi::Psi<T, Device>& psi_in,
                             const psi::Psi<T, Device>* psig,
                                   psi::Psi<T, Device>& psi_out) = 0;
        void set_kpoint(int ik_in) { this->ik = ik_in; }
    protected:
        int ik = 0;
        T one = static_cast<T>(1.0);
        T zero = static_cast<T>(0.0);
};

#endif // REPOUT_H