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
    if(this->rep_from == "pw")
    {
        // it is from pw to pw, do nothing... check if it is needed to truncate psi
    }
    else if(this->rep_from == "nao")
    {
        // transform nao to pw, can be used in COHP post-processing, multiply psi_in by psig
        // CAREFULLY CHECK!
        hsolver::gemm_op<Real, Device>()(ctx, 'N', 'N', 
                                         psi_in.get_nbands(), psi_out.get_nbands(), psig->get_nbasis(),
                                         &(this->one), psi_in, psi_in.get_nbands(),
                                         psig->get_pointer(), psig->get_nbasis(),
                                         &(this->zero), psi_out, psi_out.get_nbands());
    }
}