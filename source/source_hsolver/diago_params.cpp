#include "diago_params.h"
#include "diago_iter_assist.h"

namespace hsolver
{

template <typename T, typename Device>
void setup_diago_params_pw(const int istep,
                            const int iter,
                            const double ethr,
                            const Input_para& inp)
{
    /// choose if psi should be diag in subspace
    /// be careful that istep start from 0 and iter start from 1
    DiagoIterAssist<T, Device>::need_subspace = ((istep == 0 || istep == 1) && iter == 1) ? false : true;
    DiagoIterAssist<T, Device>::SCF_ITER = iter;
    DiagoIterAssist<T, Device>::PW_DIAG_THR = ethr;

    if (inp.calculation != "nscf")
    {
        DiagoIterAssist<T, Device>::PW_DIAG_NMAX = inp.pw_diag_nmax;
    }
}

template <typename T, typename Device>
void setup_diago_params_sdft(const int istep,
                              const int iter,
                              const double ethr,
                              const Input_para& inp)
{
    /// choose if psi should be diag in subspace
    /// be careful that istep start from 0 and iter start from 1
    if (istep == 0 && iter == 1 || inp.calculation == "nscf")
    {
        DiagoIterAssist<T, Device>::need_subspace = false;
    }
    else
    {
        DiagoIterAssist<T, Device>::need_subspace = true;
    }

    DiagoIterAssist<T, Device>::PW_DIAG_THR = ethr;
    DiagoIterAssist<T, Device>::PW_DIAG_NMAX = inp.pw_diag_nmax;
}

/// Template instantiation for CPU
template void setup_diago_params_pw<std::complex<float>, base_device::DEVICE_CPU>(const int istep,
                                                                                    const int iter,
                                                                                    const double ethr,
                                                                                    const Input_para& inp);
template void setup_diago_params_pw<std::complex<double>, base_device::DEVICE_CPU>(const int istep,
                                                                                     const int iter,
                                                                                     const double ethr,
                                                                                     const Input_para& inp);
template void setup_diago_params_pw<double, base_device::DEVICE_CPU>(const int istep,
                                                                       const int iter,
                                                                       const double ethr,
                                                                       const Input_para& inp);

/// Template instantiation for GPU
#if ((defined __CUDA) || (defined __ROCM))
template void setup_diago_params_pw<std::complex<float>, base_device::DEVICE_GPU>(const int istep,
                                                                                    const int iter,
                                                                                    const double ethr,
                                                                                    const Input_para& inp);
template void setup_diago_params_pw<std::complex<double>, base_device::DEVICE_GPU>(const int istep,
                                                                                     const int iter,
                                                                                     const double ethr,
                                                                                     const Input_para& inp);
template void setup_diago_params_pw<double, base_device::DEVICE_GPU>(const int istep,
                                                                       const int iter,
                                                                       const double ethr,
                                                                       const Input_para& inp);
#endif

/// Template instantiation for SDFT CPU
template void setup_diago_params_sdft<std::complex<float>, base_device::DEVICE_CPU>(const int istep,
                                                                                      const int iter,
                                                                                      const double ethr,
                                                                                      const Input_para& inp);
template void setup_diago_params_sdft<std::complex<double>, base_device::DEVICE_CPU>(const int istep,
                                                                                       const int iter,
                                                                                       const double ethr,
                                                                                       const Input_para& inp);
template void setup_diago_params_sdft<double, base_device::DEVICE_CPU>(const int istep,
                                                                         const int iter,
                                                                         const double ethr,
                                                                         const Input_para& inp);

/// Template instantiation for SDFT GPU
#if ((defined __CUDA) || (defined __ROCM))
template void setup_diago_params_sdft<std::complex<float>, base_device::DEVICE_GPU>(const int istep,
                                                                                      const int iter,
                                                                                      const double ethr,
                                                                                      const Input_para& inp);
template void setup_diago_params_sdft<std::complex<double>, base_device::DEVICE_GPU>(const int istep,
                                                                                       const int iter,
                                                                                       const double ethr,
                                                                                       const Input_para& inp);
template void setup_diago_params_sdft<double, base_device::DEVICE_GPU>(const int istep,
                                                                         const int iter,
                                                                         const double ethr,
                                                                         const Input_para& inp);
#endif

} // namespace hsolver
