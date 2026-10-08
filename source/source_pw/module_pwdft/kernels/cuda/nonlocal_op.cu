#include "source_pw/module_pwdft/kernels/nonlocal_op.h"

#include <complex>
#include <cstdint>

#include <cuda_runtime.h>
#include <thrust/complex.h>
#include <base/macros/macros.h>

namespace hamilt {

#define THREADS_PER_BLOCK 256

template <typename FPTYPE>
__global__ void nonlocal_pw(
    const int l1,
    const int l2,
    const int l3,
    const int sum,
    const int iat,
    const int spin,
    const int nkb,
    const int deeq_x,
    const int deeq_y,
    const int deeq_z,
    const FPTYPE* deeq,
    thrust::complex<FPTYPE>* ps,
    const thrust::complex<FPTYPE>* becp)
{
  const int ii = blockIdx.x / l2;
  const int jj = blockIdx.x % l2;
  for (int kk = threadIdx.x; kk < l3; kk += blockDim.x) {
    thrust::complex<FPTYPE> res(0.0, 0.0);
    const std::int64_t projector_offset = sum + static_cast<std::int64_t>(ii) * l3;
    for (int xx = 0; xx < l3; xx++) {
      const std::int64_t deeq_index = ((static_cast<std::int64_t>(spin) * deeq_x + iat + ii) * deeq_y + xx) * deeq_z + kk;
      const std::int64_t input_index = static_cast<std::int64_t>(jj) * nkb + projector_offset + xx;
      res
        += deeq[deeq_index]
        *  becp[input_index];
    }
    const std::int64_t output_index = (projector_offset + kk) * l2 + jj;
    ps[output_index] += res;
  }
}

template <typename FPTYPE>
__global__ void nonlocal_pw(
    const int l1,
    const int l2,
    const int l3,
    const int sum,
    const int iat,
    const int nkb,
    const int deeq_x,
    const int deeq_y,
    const int deeq_z,
    const thrust::complex<FPTYPE>* deeq_nc,
    thrust::complex<FPTYPE>* ps,
    const thrust::complex<FPTYPE>* becp)
{
  const std::int64_t spinor_index = static_cast<std::int64_t>(blockIdx.x) * 2;
  const int ii = spinor_index / l2;
  const int jj = spinor_index % l2;
  for (int kk = threadIdx.x; kk < l3; kk += blockDim.x) {
    thrust::complex<FPTYPE> res1(0.0, 0.0);
    thrust::complex<FPTYPE> res2(0.0, 0.0);
    const std::int64_t projector_offset = sum + static_cast<std::int64_t>(ii) * l3;
    const std::int64_t psind = (projector_offset + kk) * l2 + jj;
    const std::int64_t spin_stride = static_cast<std::int64_t>(deeq_x) * deeq_y * deeq_z;
    for (int xx = 0; xx < l3; xx++) {
      const std::int64_t becpind = static_cast<std::int64_t>(jj) * nkb + projector_offset + xx;
      const std::int64_t deeq_index = ((static_cast<std::int64_t>(iat) + ii) * deeq_y + kk) * deeq_z + xx;
      thrust::complex<FPTYPE> becp1 = becp[becpind];
      thrust::complex<FPTYPE> becp2 = becp[becpind + nkb];
      res1 += deeq_nc[deeq_index] * becp1
                   + deeq_nc[spin_stride + deeq_index] * becp2;
      res2 += deeq_nc[2 * spin_stride + deeq_index] * becp1
                       + deeq_nc[3 * spin_stride + deeq_index] * becp2;
    }
    ps[psind] += res1;
    ps[psind + 1] += res2;
  }
}

template <typename FPTYPE>
void hamilt::nonlocal_pw_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* dev,
                                                                         const int& l1,
                                                                         const int& l2,
                                                                         const int& l3,
                                                                         int& sum,
                                                                         int& iat,
                                                                         const int& spin,
                                                                         const int& nkb,
                                                                         const int& deeq_x,
                                                                         const int& deeq_y,
                                                                         const int& deeq_z,
                                                                         const FPTYPE* deeq,
                                                                         std::complex<FPTYPE>* ps,
                                                                         const std::complex<FPTYPE>* becp)
{
  // denghui implement 20221019
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  nonlocal_pw<FPTYPE><<<l1 * l2, THREADS_PER_BLOCK>>>(
    l1, l2, l3, // loop size
    sum, iat, spin, nkb,   // control params
    deeq_x, deeq_y, deeq_z, deeq,  // deeq realArray operator()
    reinterpret_cast<thrust::complex<FPTYPE>*>(ps), // array of data
    reinterpret_cast<const thrust::complex<FPTYPE>*>(becp)); // array of data

  CHECK_CUDA_SYNC();
  iat += l1;
  sum += l1 * l3;
}

template <typename FPTYPE>
void hamilt::nonlocal_pw_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* dev,
                                                                         const int& l1,
                                                                         const int& l2,
                                                                         const int& l3,
                                                                         int& sum,
                                                                         int& iat,
                                                                         const int& nkb,
                                                                         const int& deeq_x,
                                                                         const int& deeq_y,
                                                                         const int& deeq_z,
                                                                         const std::complex<FPTYPE>* deeq_nc,
                                                                         std::complex<FPTYPE>* ps,
                                                                         const std::complex<FPTYPE>* becp)
{
  // denghui implement 20221109
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  const std::int64_t block_count = static_cast<std::int64_t>(l1) * l2 / 2;
  nonlocal_pw<FPTYPE><<<block_count, THREADS_PER_BLOCK>>>(
    l1, l2, l3, // loop size
    sum, iat, nkb,   // control params
    deeq_x, deeq_y, deeq_z,
    reinterpret_cast<const thrust::complex<FPTYPE>*>(deeq_nc),  // deeq realArray operator()
    reinterpret_cast<thrust::complex<FPTYPE>*>(ps), // array of data
    reinterpret_cast<const thrust::complex<FPTYPE>*>(becp)); // array of data

  CHECK_CUDA_SYNC();
  iat += l1;
  sum += l1 * l3;
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
}

template struct nonlocal_pw_op<float, base_device::DEVICE_GPU>;
template struct nonlocal_pw_op<double, base_device::DEVICE_GPU>;

}  // namespace hamilt
