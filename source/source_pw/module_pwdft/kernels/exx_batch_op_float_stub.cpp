#include "source_pw/module_pwdft/kernels/exx_batch_op.h"

#include "source_base/tool_quit.h"

// Stubs for the single-precision batched FFTW plans, compiled when
// ENABLE_FLOAT_FFTW is off (see module_pwdft/CMakeLists.txt); the real
// implementation lives in exx_batch_op_float.cpp. Float PW on CPU needs the
// single-precision FFTW library anyway, so reaching these means the build
// lacks it: fail with a clear message instead of a null-symbol crash.

namespace hamilt
{

namespace
{
[[noreturn]] void exx_fftw_float_missing()
{
    ModuleBase::WARNING_QUIT("exx_batch_fft",
                             "single-precision FFTW is required for the float EXX FFT on CPU; "
                             "rebuild with ENABLE_FLOAT_FFTW=ON");
}
} // namespace

template <>
void exx_batch_fft_plan_create<std::complex<float>, base_device::DEVICE_CPU>(void**,
                                                                             const int,
                                                                             const int,
                                                                             const int,
                                                                             const int)
{
    exx_fftw_float_missing();
}

template <>
void exx_batch_fft_exec<std::complex<float>, base_device::DEVICE_CPU>(void*, std::complex<float>*, const bool)
{
    exx_fftw_float_missing();
}

template <>
void exx_batch_fft_plan_destroy<std::complex<float>, base_device::DEVICE_CPU>(void** plan)
{
    // no plan could have been created, but stay harmless at process teardown
    *plan = nullptr;
}

#if !defined(__CUDA)
// link fodder for the DEVICE_GPU instantiations, never executed
template <>
void exx_batch_fft_plan_create<std::complex<float>, base_device::DEVICE_GPU>(void**,
                                                                             const int,
                                                                             const int,
                                                                             const int,
                                                                             const int)
{
    exx_fftw_float_missing();
}

template <>
void exx_batch_fft_exec<std::complex<float>, base_device::DEVICE_GPU>(void*, std::complex<float>*, const bool)
{
    exx_fftw_float_missing();
}

template <>
void exx_batch_fft_plan_destroy<std::complex<float>, base_device::DEVICE_GPU>(void** plan)
{
    *plan = nullptr;
}
#endif

} // namespace hamilt
