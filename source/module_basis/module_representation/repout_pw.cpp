#include "repout_pw.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_base/timer.h"
template<typename T, typename Device>
RepOut_PW<T, Device>::RepOut_PW()
{
}

template<typename T, typename Device>
RepOut_PW<T, Device>::~RepOut_PW()
{
}

template<typename T, typename Device>
void RepOut_PW<T, Device>::matmul(
            T* mat1, T* mat2, T* mat3,
            const int nrow1, const int ncol1, 
            const int nrow2, const int ncol2,
            const int nrow3, const int ncol3
        )
{
    hsolver::gemm_op<T, Device>()(
        this->ctx,
        'N', 'N',
        ncol3, nrow3, ncol1,
        &this->one,
        mat1, ncol1,
        mat2, ncol2,
        &this->zero,
        mat3, ncol3
    );
}

template<typename T, typename Device>
void RepOut_PW<T, Device>::project(psi::Psi<T, Device>* psi_in,
                                   psi::Psi<T, Device>* psig,
                                   psi::Psi<T, Device>* psi_out)
{
    ModuleBase::timer::tick("RepOut_PW", "project");
    this->matmul(
        psi_in->get_pointer(), psig->get_pointer(), psi_out->get_pointer(),
        psi_in->get_nbands(), psi_in->get_nbasis(),
        psig->get_nbands(), psig->get_nbasis(),
        psi_out->get_nbasis(), psi_out->get_nbands()
    );
    ModuleBase::timer::tick("RepOut_PW", "project");
}

// Explicit instantiation
template class RepOut_PW<std::complex<double>, psi::DEVICE_CPU>;
template class RepOut_PW<std::complex<float>, psi::DEVICE_CPU>;
// gamma point support
template class RepOut_PW<double, psi::DEVICE_CPU>;
template class RepOut_PW<float, psi::DEVICE_CPU>;
// gpu support
#if ((defined __CUDA) || (defined __ROCM))
template class RepOut_PW<std::complex<double>, psi::DEVICE_GPU>;
template class RepOut_PW<std::complex<float>, psi::DEVICE_GPU>;
// gamma point support
template class RepOut_PW<double, psi::DEVICE_GPU>;
template class RepOut_PW<float, psi::DEVICE_GPU>;
#endif