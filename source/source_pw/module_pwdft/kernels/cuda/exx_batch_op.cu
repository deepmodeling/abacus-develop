#include "source_pw/module_pwdft/kernels/exx_batch_op.h"

#include <thrust/complex.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <base/macros/macros.h>
#include <source_base/module_device/device_check.h>

namespace hamilt
{

#define THREADS_PER_BLOCK 256

template <typename T>
struct CufftTraits;

template <>
struct CufftTraits<std::complex<double>>
{
    using cufft_t = cufftDoubleComplex;
    static constexpr cufftType type = CUFFT_Z2Z;
    static cufftResult exec(cufftHandle plan, cufft_t* data, int dir) { return cufftExecZ2Z(plan, data, data, dir); }
};

template <>
struct CufftTraits<std::complex<float>>
{
    using cufft_t = cufftComplex;
    static constexpr cufftType type = CUFFT_C2C;
    static cufftResult exec(cufftHandle plan, cufft_t* data, int dir) { return cufftExecC2C(plan, data, data, dir); }
};

template <typename T>
static void exx_batch_fft_plan_create_gpu(void** plan, const int nx, const int ny, const int nz, const int batch)
{
    cufftHandle handle;
    int n[3] = {nx, ny, nz};
    CHECK_CUFFT(cufftPlanMany(&handle,
                              3,
                              n,
                              nullptr, 1, 0, // input: contiguous, distance = nx*ny*nz
                              nullptr, 1, 0, // output: same (in-place)
                              CufftTraits<T>::type,
                              batch));
    *plan = new cufftHandle(handle);
}

template <typename T>
static void exx_batch_fft_exec_gpu(void* plan, T* data, const bool forward)
{
    cufftHandle handle = *reinterpret_cast<cufftHandle*>(plan);
    CHECK_CUFFT(CufftTraits<T>::exec(handle,
                                     reinterpret_cast<typename CufftTraits<T>::cufft_t*>(data),
                                     forward ? CUFFT_FORWARD : CUFFT_INVERSE));
}

template <typename T>
static void exx_batch_fft_plan_destroy_gpu(void** plan)
{
    if (*plan != nullptr)
    {
        cufftHandle handle = *reinterpret_cast<cufftHandle*>(plan);
        // no CHECK_CUFFT here: at process teardown the CUDA context may already
        // be gone, in which case cufftDestroy reports CUFFT_INVALID_PLAN; that
        // is harmless and must not abort the run
        cufftDestroy(handle);
        delete reinterpret_cast<cufftHandle*>(*plan);
        *plan = nullptr;
    }
}

template <class FPTYPE>
__global__ void batch_density_real(const int nbands,
                                   const int nrxx,
                                   const thrust::complex<FPTYPE>* nk_real_all,
                                   const thrust::complex<FPTYPE>* mq_real,
                                   const FPTYPE omega_inv,
                                   thrust::complex<FPTYPE>* out)
{
    const long long idx = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx < static_cast<long long>(nbands) * nrxx)
    {
        const int i = idx % nrxx;
        out[idx] = nk_real_all[idx] * thrust::conj(mq_real[i]) * omega_inv;
    }
}

template <class FPTYPE>
__global__ void batch_gather_pw(const int nbands,
                                const int npw,
                                const int nxyz,
                                const int* box_map,
                                const thrust::complex<FPTYPE>* box,
                                thrust::complex<FPTYPE>* pw)
{
    const long long idx = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx < static_cast<long long>(nbands) * npw)
    {
        const long long n = idx / npw;
        const int ig = idx % npw;
        pw[idx] = box[n * nxyz + box_map[ig]] / static_cast<FPTYPE>(nxyz);
    }
}

template <class FPTYPE>
__global__ void batch_mul_pot(const int nbands,
                              const int npw,
                              const FPTYPE* pot,
                              thrust::complex<FPTYPE>* pw)
{
    const long long idx = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx < static_cast<long long>(nbands) * npw)
    {
        const int ig = idx % npw;
        pw[idx] *= pot[ig];
    }
}

template <class FPTYPE>
__global__ void batch_scatter_pw(const int nbands,
                                 const int npw,
                                 const int nxyz,
                                 const int* box_map,
                                 const thrust::complex<FPTYPE>* pw,
                                 thrust::complex<FPTYPE>* box)
{
    const long long idx = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx < static_cast<long long>(nbands) * npw)
    {
        const long long n = idx / npw;
        const int ig = idx % npw;
        box[n * nxyz + box_map[ig]] = pw[idx];
    }
}

template <class FPTYPE>
__global__ void batch_mul_real(const int nbands,
                               const int nrxx,
                               const thrust::complex<FPTYPE>* mq_real,
                               thrust::complex<FPTYPE>* data)
{
    const long long idx = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx < static_cast<long long>(nbands) * nrxx)
    {
        const int i = idx % nrxx;
        data[idx] *= mq_real[i];
    }
}

template <class FPTYPE>
__global__ void batch_gather_accum(const int nbands,
                                   const int npwk,
                                   const int nxyz,
                                   const int* box_map,
                                   const thrust::complex<FPTYPE>* box,
                                   const FPTYPE factor,
                                   thrust::complex<FPTYPE>* out,
                                   const int out_stride)
{
    const long long idx = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx < static_cast<long long>(nbands) * npwk)
    {
        const long long n = idx / npwk;
        const int ig = idx % npwk;
        out[n * out_stride + ig] += factor / static_cast<FPTYPE>(nxyz) * box[n * nxyz + box_map[ig]];
    }
}

template <class FPTYPE>
__global__ void batch_scatter_wfc(const int nbands,
                                  const int npwk,
                                  const int nxyz,
                                  const int* map,
                                  const thrust::complex<FPTYPE>* psi,
                                  const int psi_stride,
                                  thrust::complex<FPTYPE>* box)
{
    const long long idx = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx < static_cast<long long>(nbands) * npwk)
    {
        const long long n = idx / npwk;
        const int ig = idx % npwk;
        box[n * nxyz + map[ig]] = psi[n * psi_stride + ig];
    }
}

static inline int nblocks(const long long n)
{
    return static_cast<int>((n + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
}

template <typename T>
static void exx_batch_density_real_gpu(const int nbands,
                            const int nrxx,
                            const T* nk_real_all,
                            const T* mq_real,
                            const double omega,
                            T* out)
{
    using FPTYPE = typename T::value_type;
    batch_density_real<FPTYPE><<<nblocks(static_cast<long long>(nbands) * nrxx), THREADS_PER_BLOCK>>>(
        nbands,
        nrxx,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(nk_real_all),
        reinterpret_cast<const thrust::complex<FPTYPE>*>(mq_real),
        static_cast<FPTYPE>(1.0 / omega),
        reinterpret_cast<thrust::complex<FPTYPE>*>(out));
    CHECK_CUDA_SYNC();
}

template <typename T>
static void exx_batch_gather_pw_gpu(const int nbands,
                         const int npw,
                         const int nxyz,
                         const int* box_map,
                         const T* box,
                         T* pw)
{
    using FPTYPE = typename T::value_type;
    batch_gather_pw<FPTYPE><<<nblocks(static_cast<long long>(nbands) * npw), THREADS_PER_BLOCK>>>(
        nbands,
        npw,
        nxyz,
        box_map,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(box),
        reinterpret_cast<thrust::complex<FPTYPE>*>(pw));
    CHECK_CUDA_SYNC();
}

template <typename T, typename Real>
static void exx_batch_mul_pot_gpu(const int nbands, const int npw, const Real* pot, T* pw)
{
    using FPTYPE = typename T::value_type;
    batch_mul_pot<FPTYPE><<<nblocks(static_cast<long long>(nbands) * npw), THREADS_PER_BLOCK>>>(
        nbands,
        npw,
        pot,
        reinterpret_cast<thrust::complex<FPTYPE>*>(pw));
    CHECK_CUDA_SYNC();
}

template <typename T>
static void exx_batch_scatter_pw_gpu(const int nbands,
                          const int npw,
                          const int nxyz,
                          const int* box_map,
                          const T* pw,
                          T* box)
{
    using FPTYPE = typename T::value_type;
    batch_scatter_pw<FPTYPE><<<nblocks(static_cast<long long>(nbands) * npw), THREADS_PER_BLOCK>>>(
        nbands,
        npw,
        nxyz,
        box_map,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(pw),
        reinterpret_cast<thrust::complex<FPTYPE>*>(box));
    CHECK_CUDA_SYNC();
}

template <typename T>
static void exx_batch_mul_real_gpu(const int nbands, const int nrxx, const T* mq_real, T* data)
{
    using FPTYPE = typename T::value_type;
    batch_mul_real<FPTYPE><<<nblocks(static_cast<long long>(nbands) * nrxx), THREADS_PER_BLOCK>>>(
        nbands,
        nrxx,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(mq_real),
        reinterpret_cast<thrust::complex<FPTYPE>*>(data));
    CHECK_CUDA_SYNC();
}

template <typename T>
static void exx_batch_scatter_wfc_gpu(const int nbands,
                           const int npwk,
                           const int nxyz,
                           const int* map,
                           const T* psi,
                           const int psi_stride,
                           T* box)
{
    using FPTYPE = typename T::value_type;
    batch_scatter_wfc<FPTYPE><<<nblocks(static_cast<long long>(nbands) * npwk), THREADS_PER_BLOCK>>>(
        nbands,
        npwk,
        nxyz,
        map,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(psi),
        psi_stride,
        reinterpret_cast<thrust::complex<FPTYPE>*>(box));
    CHECK_CUDA_SYNC();
}

template <typename T, typename Real>
static void exx_batch_gather_accum_gpu(const int nbands,
                            const int npwk,
                            const int nxyz,
                            const int* box_map,
                            const T* box,
                            const Real factor,
                            T* out,
                            const int out_stride)
{
    using FPTYPE = typename T::value_type;
    batch_gather_accum<FPTYPE><<<nblocks(static_cast<long long>(nbands) * npwk), THREADS_PER_BLOCK>>>(
        nbands,
        npwk,
        nxyz,
        box_map,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(box),
        static_cast<FPTYPE>(factor),
        reinterpret_cast<thrust::complex<FPTYPE>*>(out),
        out_stride);
    CHECK_CUDA_SYNC();
}

// DEVICE_GPU explicit specializations of the API declared in exx_batch_op.h,
// forwarding to the CUDA implementations above.
#define EXX_BATCH_GPU_SPEC(T, Real)                                                                            \
    template <>                                                                                                \
    void exx_batch_density_real<T, base_device::DEVICE_GPU>(const int nbands,                                  \
                                                            const int nrxx,                                    \
                                                            const T* nk_real_all,                              \
                                                            const T* mq_real,                                  \
                                                            const double omega,                                \
                                                            T* out)                                            \
    {                                                                                                          \
        exx_batch_density_real_gpu<T>(nbands, nrxx, nk_real_all, mq_real, omega, out);                         \
    }                                                                                                          \
    template <>                                                                                                \
    void exx_batch_gather_pw<T, base_device::DEVICE_GPU>(const int nbands,                                     \
                                                         const int npw,                                        \
                                                         const int nxyz,                                       \
                                                         const int* box_map,                                   \
                                                         const T* box,                                         \
                                                         T* pw)                                                \
    {                                                                                                          \
        exx_batch_gather_pw_gpu<T>(nbands, npw, nxyz, box_map, box, pw);                                       \
    }                                                                                                          \
    template <>                                                                                                \
    void exx_batch_mul_pot<T, Real, base_device::DEVICE_GPU>(const int nbands,                                 \
                                                             const int npw,                                    \
                                                             const Real* pot,                                  \
                                                             T* pw)                                            \
    {                                                                                                          \
        exx_batch_mul_pot_gpu<T, Real>(nbands, npw, pot, pw);                                                  \
    }                                                                                                          \
    template <>                                                                                                \
    void exx_batch_scatter_pw<T, base_device::DEVICE_GPU>(const int nbands,                                    \
                                                          const int npw,                                       \
                                                          const int nxyz,                                      \
                                                          const int* box_map,                                  \
                                                          const T* pw,                                         \
                                                          T* box)                                              \
    {                                                                                                          \
        exx_batch_scatter_pw_gpu<T>(nbands, npw, nxyz, box_map, pw, box);                                      \
    }                                                                                                          \
    template <>                                                                                                \
    void exx_batch_mul_real<T, base_device::DEVICE_GPU>(const int nbands,                                      \
                                                        const int nrxx,                                        \
                                                        const T* mq_real,                                      \
                                                        T* data)                                               \
    {                                                                                                          \
        exx_batch_mul_real_gpu<T>(nbands, nrxx, mq_real, data);                                                \
    }                                                                                                          \
    template <>                                                                                                \
    void exx_batch_gather_accum<T, Real, base_device::DEVICE_GPU>(const int nbands,                            \
                                                                  const int npwk,                              \
                                                                  const int nxyz,                              \
                                                                  const int* box_map,                          \
                                                                  const T* box,                                \
                                                                  const Real factor,                           \
                                                                  T* out,                                      \
                                                                  const int out_stride)                        \
    {                                                                                                          \
        exx_batch_gather_accum_gpu<T, Real>(nbands, npwk, nxyz, box_map, box, factor, out, out_stride);        \
    }                                                                                                          \
    template <>                                                                                                \
    void exx_batch_scatter_wfc<T, base_device::DEVICE_GPU>(const int nbands,                                   \
                                                           const int npwk,                                     \
                                                           const int nxyz,                                     \
                                                           const int* map,                                     \
                                                           const T* psi,                                       \
                                                           const int psi_stride,                               \
                                                           T* box)                                             \
    {                                                                                                          \
        exx_batch_scatter_wfc_gpu<T>(nbands, npwk, nxyz, map, psi, psi_stride, box);                           \
    }                                                                                                          \
    template <>                                                                                                \
    void exx_batch_fft_plan_create<T, base_device::DEVICE_GPU>(void** plan,                                    \
                                                               const int nx,                                   \
                                                               const int ny,                                   \
                                                               const int nz,                                   \
                                                               const int batch)                                \
    {                                                                                                          \
        exx_batch_fft_plan_create_gpu<T>(plan, nx, ny, nz, batch);                                             \
    }                                                                                                          \
    template <>                                                                                                \
    void exx_batch_fft_exec<T, base_device::DEVICE_GPU>(void* plan, T* data, const bool forward)               \
    {                                                                                                          \
        exx_batch_fft_exec_gpu<T>(plan, data, forward);                                                        \
    }                                                                                                          \
    template <>                                                                                                \
    void exx_batch_fft_plan_destroy<T, base_device::DEVICE_GPU>(void** plan)                                   \
    {                                                                                                          \
        exx_batch_fft_plan_destroy_gpu<T>(plan);                                                               \
    }

EXX_BATCH_GPU_SPEC(std::complex<double>, double)
EXX_BATCH_GPU_SPEC(std::complex<float>, float)

} // namespace hamilt
