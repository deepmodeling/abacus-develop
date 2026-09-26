#include "source_pw/module_pwdft/kernels/exx_batch_op.h"

#include <fftw3.h>

namespace hamilt
{

// Generic host implementation of the EXX batched kernels: plain loops plus
// FFTW plans. The CUDA explicit specializations in kernels/cuda/exx_batch_op.cu
// provide the DEVICE_GPU versions; on non-CUDA builds the DEVICE_GPU
// instantiations below exist only to satisfy the linker (the ROCm GPU path
// never activates the batched path, so they are never executed).

template <typename T, typename Device>
void exx_batch_density_real(const int nbands,
                            const int nrxx,
                            const T* nk_real_all,
                            const T* mq_real,
                            const double omega,
                            T* out)
{
    using Real = typename T::value_type;
    const Real omega_inv = static_cast<Real>(1.0 / omega);
    for (int n = 0; n < nbands; n++)
    {
        const T* nk_real = nk_real_all + n * nrxx;
        T* out_n = out + n * nrxx;
        for (int i = 0; i < nrxx; i++)
        {
            out_n[i] = nk_real[i] * std::conj(mq_real[i]) * omega_inv;
        }
    }
}

template <typename T, typename Device>
void exx_batch_gather_pw(const int nbands,
                         const int npw,
                         const int nxyz,
                         const int* box_map,
                         const T* box,
                         T* pw)
{
    using Real = typename T::value_type;
    const Real nxyz_inv = static_cast<Real>(1.0) / static_cast<Real>(nxyz);
    for (int n = 0; n < nbands; n++)
    {
        const T* box_n = box + n * nxyz;
        T* pw_n = pw + n * npw;
        for (int ig = 0; ig < npw; ig++)
        {
            pw_n[ig] = box_n[box_map[ig]] * nxyz_inv;
        }
    }
}

template <typename T, typename Real, typename Device>
void exx_batch_mul_pot(const int nbands, const int npw, const Real* pot, T* pw)
{
    for (int n = 0; n < nbands; n++)
    {
        T* pw_n = pw + n * npw;
        for (int ig = 0; ig < npw; ig++)
        {
            pw_n[ig] *= pot[ig];
        }
    }
}

template <typename T, typename Device>
void exx_batch_scatter_pw(const int nbands,
                          const int npw,
                          const int nxyz,
                          const int* box_map,
                          const T* pw,
                          T* box)
{
    for (int n = 0; n < nbands; n++)
    {
        const T* pw_n = pw + n * npw;
        T* box_n = box + n * nxyz;
        for (int ig = 0; ig < npw; ig++)
        {
            box_n[box_map[ig]] = pw_n[ig];
        }
    }
}

template <typename T, typename Device>
void exx_batch_mul_real(const int nbands, const int nrxx, const T* mq_real, T* data)
{
    for (int n = 0; n < nbands; n++)
    {
        T* data_n = data + n * nrxx;
        for (int i = 0; i < nrxx; i++)
        {
            data_n[i] *= mq_real[i];
        }
    }
}

template <typename T, typename Real, typename Device>
void exx_batch_gather_accum(const int nbands,
                            const int npwk,
                            const int nxyz,
                            const int* box_map,
                            const T* box,
                            const Real factor,
                            T* out,
                            const int out_stride)
{
    const Real factor_scaled = factor / static_cast<Real>(nxyz);
    for (int n = 0; n < nbands; n++)
    {
        const T* box_n = box + n * nxyz;
        T* out_n = out + n * out_stride;
        for (int ig = 0; ig < npwk; ig++)
        {
            out_n[ig] += factor_scaled * box_n[box_map[ig]];
        }
    }
}

template <typename T, typename Device>
void exx_batch_scatter_wfc(const int nbands,
                           const int npwk,
                           const int nxyz,
                           const int* map,
                           const T* psi,
                           const int psi_stride,
                           T* box)
{
    for (int n = 0; n < nbands; n++)
    {
        const T* psi_n = psi + n * psi_stride;
        T* box_n = box + n * nxyz;
        for (int ig = 0; ig < npwk; ig++)
        {
            box_n[map[ig]] = psi_n[ig];
        }
    }
}

#define EXX_BATCH_INSTANTIATE(T, Real, DEV)                                                                 \
    template void exx_batch_density_real<T, DEV>(int, int, const T*, const T*, double, T*);                 \
    template void exx_batch_gather_pw<T, DEV>(int, int, int, const int*, const T*, T*);                     \
    template void exx_batch_mul_pot<T, Real, DEV>(int, int, const Real*, T*);                               \
    template void exx_batch_scatter_pw<T, DEV>(int, int, int, const int*, const T*, T*);                    \
    template void exx_batch_mul_real<T, DEV>(int, int, const T*, T*);                                       \
    template void exx_batch_gather_accum<T, Real, DEV>(int, int, int, const int*, const T*, Real, T*, int); \
    template void exx_batch_scatter_wfc<T, DEV>(int, int, int, const int*, const T*, int, T*);

EXX_BATCH_INSTANTIATE(std::complex<double>, double, base_device::DEVICE_CPU)
EXX_BATCH_INSTANTIATE(std::complex<float>, float, base_device::DEVICE_CPU)
#if !defined(__CUDA)
// Non-CUDA builds (including ROCm): the DEVICE_GPU instantiations of
// OperatorEXXPW reference these symbols, but the batched path is never
// activated there, so the host implementations are only link fodder.
EXX_BATCH_INSTANTIATE(std::complex<double>, double, base_device::DEVICE_GPU)
EXX_BATCH_INSTANTIATE(std::complex<float>, float, base_device::DEVICE_GPU)
#endif

// Batched FFTW c2c plans (double precision; the single-precision ones live in
// exx_batch_op_float.cpp / exx_batch_op_float_stub.cpp). FFTW encodes the
// transform direction in the plan (unlike cuFFT), so the handle holds one
// plan per direction. The plans are created on a scratch buffer with
// FFTW_UNALIGNED and run through fftw_execute_dft, so any buffer of the right
// size can be transformed.
namespace
{
struct ExxFftwPlanD
{
    fftw_plan fwd = nullptr;
    fftw_plan bac = nullptr;
};

void exx_fftw_plan_create_d(void** plan, const int nx, const int ny, const int nz, const int batch)
{
    const int n[3] = {nx, ny, nz};
    const int nxyz = nx * ny * nz;
    fftw_complex* tmp = reinterpret_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * nxyz * batch));
    auto* p = new ExxFftwPlanD;
    const unsigned flag = FFTW_ESTIMATE | FFTW_UNALIGNED;
    p->fwd = fftw_plan_many_dft(3, n, batch, tmp, nullptr, 1, nxyz, tmp, nullptr, 1, nxyz, FFTW_FORWARD, flag);
    p->bac = fftw_plan_many_dft(3, n, batch, tmp, nullptr, 1, nxyz, tmp, nullptr, 1, nxyz, FFTW_BACKWARD, flag);
    fftw_free(tmp);
    *plan = p;
}

void exx_fftw_exec_d(void* plan, std::complex<double>* data, const bool forward)
{
    auto* p = reinterpret_cast<ExxFftwPlanD*>(plan);
    fftw_execute_dft(forward ? p->fwd : p->bac,
                     reinterpret_cast<fftw_complex*>(data),
                     reinterpret_cast<fftw_complex*>(data));
}

void exx_fftw_plan_destroy_d(void** plan)
{
    if (*plan != nullptr)
    {
        auto* p = reinterpret_cast<ExxFftwPlanD*>(*plan);
        fftw_destroy_plan(p->fwd);
        fftw_destroy_plan(p->bac);
        delete p;
        *plan = nullptr;
    }
}
} // namespace

template <>
void exx_batch_fft_plan_create<std::complex<double>, base_device::DEVICE_CPU>(void** plan,
                                                                              const int nx,
                                                                              const int ny,
                                                                              const int nz,
                                                                              const int batch)
{
    exx_fftw_plan_create_d(plan, nx, ny, nz, batch);
}

template <>
void exx_batch_fft_exec<std::complex<double>, base_device::DEVICE_CPU>(void* plan,
                                                                       std::complex<double>* data,
                                                                       const bool forward)
{
    exx_fftw_exec_d(plan, data, forward);
}

template <>
void exx_batch_fft_plan_destroy<std::complex<double>, base_device::DEVICE_CPU>(void** plan)
{
    exx_fftw_plan_destroy_d(plan);
}

#if !defined(__CUDA)
// link fodder for the DEVICE_GPU instantiations, never executed (see above)
template <>
void exx_batch_fft_plan_create<std::complex<double>, base_device::DEVICE_GPU>(void** plan,
                                                                              const int nx,
                                                                              const int ny,
                                                                              const int nz,
                                                                              const int batch)
{
    exx_fftw_plan_create_d(plan, nx, ny, nz, batch);
}

template <>
void exx_batch_fft_exec<std::complex<double>, base_device::DEVICE_GPU>(void* plan,
                                                                       std::complex<double>* data,
                                                                       const bool forward)
{
    exx_fftw_exec_d(plan, data, forward);
}

template <>
void exx_batch_fft_plan_destroy<std::complex<double>, base_device::DEVICE_GPU>(void** plan)
{
    exx_fftw_plan_destroy_d(plan);
}
#endif

} // namespace hamilt
