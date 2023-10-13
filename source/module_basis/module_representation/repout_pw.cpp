#include "repout_pw.h"
#include "module_hsolver/kernels/math_kernel_op.h"

template<typename T, typename Device>
RepOut_PW<T, Device>::RepOut_PW()
{
}

template<typename T, typename Device>
RepOut_PW<T, Device>::~RepOut_PW()
{
}

template<typename T, typename Device>
void RepOut_PW<T, Device>::project(const psi::Psi<T, Device>& psi_in,
                                   const psi::Psi<T, Device>* psig,
                                         psi::Psi<T, Device>& psi_out)
{
    // default the psi_in to be identity matrix
}