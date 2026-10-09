#include "source_pw/module_pwdft/kernels/exx_batch_op.h"

#include <fftw3.h>

// Single-precision batched FFTW plans for the EXX batched path, compiled only
// when ENABLE_FLOAT_FFTW is on (see module_pwdft/CMakeLists.txt); otherwise
// exx_batch_op_float_stub.cpp provides the same symbols with a runtime error.
// Mirrors the double-precision implementation in exx_batch_op.cpp.

namespace hamilt
{

namespace
{
struct ExxFftwPlanF
{
    fftwf_plan fwd = nullptr;
    fftwf_plan bac = nullptr;
};

void exx_fftw_plan_create_f(void** plan, const int nx, const int ny, const int nz, const int batch)
{
    const int n[3] = {nx, ny, nz};
    const int nxyz = nx * ny * nz;
    fftwf_complex* tmp = reinterpret_cast<fftwf_complex*>(fftwf_malloc(sizeof(fftwf_complex) * nxyz * batch));
    auto* p = new ExxFftwPlanF;
    const unsigned flag = FFTW_ESTIMATE | FFTW_UNALIGNED;
    p->fwd = fftwf_plan_many_dft(3, n, batch, tmp, nullptr, 1, nxyz, tmp, nullptr, 1, nxyz, FFTW_FORWARD, flag);
    p->bac = fftwf_plan_many_dft(3, n, batch, tmp, nullptr, 1, nxyz, tmp, nullptr, 1, nxyz, FFTW_BACKWARD, flag);
    fftwf_free(tmp);
    *plan = p;
}

void exx_fftw_exec_f(void* plan, std::complex<float>* data, const bool forward)
{
    auto* p = reinterpret_cast<ExxFftwPlanF*>(plan);
    fftwf_execute_dft(forward ? p->fwd : p->bac,
                      reinterpret_cast<fftwf_complex*>(data),
                      reinterpret_cast<fftwf_complex*>(data));
}

void exx_fftw_plan_destroy_f(void** plan)
{
    if (*plan != nullptr)
    {
        auto* p = reinterpret_cast<ExxFftwPlanF*>(*plan);
        fftwf_destroy_plan(p->fwd);
        fftwf_destroy_plan(p->bac);
        delete p;
        *plan = nullptr;
    }
}
} // namespace

template <>
void exx_batch_fft_plan_create<std::complex<float>, base_device::DEVICE_CPU>(void** plan,
                                                                             const int nx,
                                                                             const int ny,
                                                                             const int nz,
                                                                             const int batch)
{
    exx_fftw_plan_create_f(plan, nx, ny, nz, batch);
}

template <>
void exx_batch_fft_exec<std::complex<float>, base_device::DEVICE_CPU>(void* plan,
                                                                      std::complex<float>* data,
                                                                      const bool forward)
{
    exx_fftw_exec_f(plan, data, forward);
}

template <>
void exx_batch_fft_plan_destroy<std::complex<float>, base_device::DEVICE_CPU>(void** plan)
{
    exx_fftw_plan_destroy_f(plan);
}

#if !defined(__CUDA)
// link fodder for the DEVICE_GPU instantiations, never executed (ROCm keeps
// the batched path disabled)
template <>
void exx_batch_fft_plan_create<std::complex<float>, base_device::DEVICE_GPU>(void** plan,
                                                                             const int nx,
                                                                             const int ny,
                                                                             const int nz,
                                                                             const int batch)
{
    exx_fftw_plan_create_f(plan, nx, ny, nz, batch);
}

template <>
void exx_batch_fft_exec<std::complex<float>, base_device::DEVICE_GPU>(void* plan,
                                                                      std::complex<float>* data,
                                                                      const bool forward)
{
    exx_fftw_exec_f(plan, data, forward);
}

template <>
void exx_batch_fft_plan_destroy<std::complex<float>, base_device::DEVICE_GPU>(void** plan)
{
    exx_fftw_plan_destroy_f(plan);
}
#endif

} // namespace hamilt
