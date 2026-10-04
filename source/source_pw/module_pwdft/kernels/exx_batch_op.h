#ifndef EXX_BATCH_OP_H
#define EXX_BATCH_OP_H

#include "source_base/module_device/types.h"

#include <complex>

namespace hamilt
{

// Batched elementwise / gather-scatter kernels and batched FFT helpers for
// the n-band inner loop of OperatorEXXPW. All bands of one (iq, m_iband)
// pair are processed with a single kernel launch / batched FFT instead of
// one launch per band. The generic implementation (exx_batch_op.cpp) runs on
// the host with plain loops and FFTW plans; the CUDA explicit
// specializations (kernels/cuda/exx_batch_op.cu) run the same operations on
// the GPU with kernels and cuFFT. Instantiated for std::complex<float> and
// std::complex<double>.

// out[n*nrxx + i] = nk_real_all[n*nrxx + i] * conj(mq_real[i]) / omega
template <typename T, typename Device>
void exx_batch_density_real(const int nbands,
                            const int nrxx,
                            const T* nk_real_all,
                            const T* mq_real,
                            const double omega,
                            T* out);

// pw[n*npw + ig] = box[n*nxyz + box_map[ig]] / nxyz   (real -> recip, gather)
template <typename T, typename Device>
void exx_batch_gather_pw(const int nbands,
                         const int npw,
                         const int nxyz,
                         const int* box_map,
                         const T* box,
                         T* pw);

// pw[n*npw + ig] *= pot[ig]
template <typename T, typename Real, typename Device>
void exx_batch_mul_pot(const int nbands, const int npw, const Real* pot, T* pw);

// box[n*nxyz + box_map[ig]] = pw[n*npw + ig]   (recip -> real, scatter; box pre-zeroed)
template <typename T, typename Device>
void exx_batch_scatter_pw(const int nbands,
                          const int npw,
                          const int nxyz,
                          const int* box_map,
                          const T* pw,
                          T* box);

// data[n*nrxx + i] *= mq_real[i]
template <typename T, typename Device>
void exx_batch_mul_real(const int nbands, const int nrxx, const T* mq_real, T* data);

// out[n*out_stride + ig] += factor / nxyz * box[n*nxyz + box_map[ig]]   (real -> recip, gather + accumulate)
template <typename T, typename Real, typename Device>
void exx_batch_gather_accum(const int nbands,
                            const int npwk,
                            const int nxyz,
                            const int* box_map,
                            const T* box,
                            const Real factor,
                            T* out,
                            const int out_stride);

// box[n*nxyz + map[ig]] = psi[n*psi_stride + ig]   (scatter PW coefficients into the box; box pre-zeroed)
template <typename T, typename Device>
void exx_batch_scatter_wfc(const int nbands,
                           const int npwk,
                           const int nxyz,
                           const int* map,
                           const T* psi,
                           const int psi_stride,
                           T* box);

// Batched in-place 3D complex-to-complex FFT over the (nx, ny, nz) grid,
// batch = nbands transforms with contiguous distance nx*ny*nz, matching the
// layout of cufftPlan3d used in FFT_CUDA.
template <typename T, typename Device>
void exx_batch_fft_plan_create(void** plan, const int nx, const int ny, const int nz, const int batch);

// forward = FFT_FORWARD, otherwise FFT_BACKWARD (unnormalized, as in FFT_CUDA)
template <typename T, typename Device>
void exx_batch_fft_exec(void* plan, T* data, const bool forward);

template <typename T, typename Device>
void exx_batch_fft_plan_destroy(void** plan);

} // namespace hamilt

#endif // EXX_BATCH_OP_H
