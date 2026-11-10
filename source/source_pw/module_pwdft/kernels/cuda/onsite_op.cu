#include "source_pw/module_pwdft/kernels/onsite_op.h"

#include <base/macros/macros.h>
#include <complex>
#include <cuda_runtime.h>
#include <thrust/complex.h>

namespace hamilt
{

#define THREADS_PER_BLOCK 256

namespace
{

template <bool Spinor>
__global__ void apply_spin_constraint(int bands, int nkb, int sign, const double* lambda,
                                      const thrust::complex<double>* becp, thrust::complex<double>* ps)
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < bands * nkb)
    {
        const int ip = index % nkb;
        const int band = index / nkb;
        const int row = band * (Spinor ? 2 : 1) * nkb + ip;
        const double z = lambda[3 * ip + 2];
        if (Spinor)
        {
            const thrust::complex<double> diagonal(z, 0.0);
            const thrust::complex<double> plus(lambda[3 * ip], lambda[3 * ip + 1]);
            const thrust::complex<double> minus(lambda[3 * ip], -lambda[3 * ip + 1]);
            const auto up = becp[row], down = becp[row + nkb];
            ps[row] = diagonal * up + minus * down;
            ps[row + nkb] = plus * up - diagonal * down;
        }
        else
        {
            ps[row] = (z * sign) * becp[row];
        }
    }
}

} // namespace

namespace
{

template <bool Spinor>
__global__ void accumulate_moments(int bands, int nkb, int sign, const int* offsets,
                                   const double* weights, const thrust::complex<double>* becp,
                                   double* mi_xyz)
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= bands * gridDim.y)
    {
        return;
    }
    const int atom = blockIdx.y;
    const int band = index / gridDim.y;
    const double weight = weights[band];
    if (weight == 0.0)
    {
        return;
    }
    double up_norm = 0.0, down_norm = 0.0, cross_real = 0.0, cross_imag = 0.0;
    for (int projector = offsets[atom]; projector < offsets[atom + 1]; ++projector)
    {
        const int row = band * (Spinor ? 2 : 1) * nkb + projector;
        const auto up = becp[row];
        up_norm += up.real() * up.real() + up.imag() * up.imag();
        if (Spinor)
        {
            const auto down = becp[row + nkb];
            down_norm += down.real() * down.real() + down.imag() * down.imag();
            cross_real += up.real() * down.real() + up.imag() * down.imag();
            cross_imag += up.real() * down.imag() - up.imag() * down.real();
        }
    }
    if (Spinor)
    {
        atomicAdd(&mi_xyz[3 * atom], weight * 2.0 * cross_real);
        atomicAdd(&mi_xyz[3 * atom + 1], weight * 2.0 * cross_imag);
        atomicAdd(&mi_xyz[3 * atom + 2], weight * (up_norm - down_norm));
    }
    else
    {
        atomicAdd(&mi_xyz[3 * atom + 2], weight * sign * up_norm);
    }
}

} // namespace

void accumulate_projector_moments_gpu(int bands, int npol, int atoms, int nkb, int sign,
                                      const int* offsets, const double* weights,
                                      const std::complex<double>* projections, double* mi_xyz)
{
    if (atoms <= 0 || bands <= 0)
    {
        return;
    }
    const int threads = 128;
    const dim3 grid((bands + threads - 1) / threads, atoms);
    const auto* input = reinterpret_cast<const thrust::complex<double>*>(projections);
    if (npol == 2)
    {
        accumulate_moments<true><<<grid, threads>>>(bands, nkb, sign, offsets, weights, input, mi_xyz);
    }
    else
    {
        accumulate_moments<false><<<grid, threads>>>(bands, nkb, sign, offsets, weights, input, mi_xyz);
    }
    CHECK_CUDA_SYNC();
}

void spin_constraint_projectors_gpu(int bands, int npol, int nkb, int sign,
                                    const double* lambda, const std::complex<double>* becp,
                                    std::complex<double>* ps)
{
    const int blocks = (bands * nkb + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    const auto* input = reinterpret_cast<const thrust::complex<double>*>(becp);
    auto* output = reinterpret_cast<thrust::complex<double>*>(ps);
    if (npol == 2)
    {
        apply_spin_constraint<true><<<blocks, THREADS_PER_BLOCK>>>(bands, nkb, sign, lambda, input, output);
    }
    else
    {
        apply_spin_constraint<false><<<blocks, THREADS_PER_BLOCK>>>(bands, nkb, sign, lambda, input, output);
    }
    CHECK_CUDA_SYNC();
}

template <typename FPTYPE>
__global__ void onsite_op(const int npm,
                          const int npol,
                          const int* ip_iat,
                          const int tnp,
                          const thrust::complex<FPTYPE>* lambda_coeff,
                          thrust::complex<FPTYPE>* ps,
                          const thrust::complex<FPTYPE>* becp)
{
    const int ip = blockIdx.x;
    const int nbands = npm / npol;
    int iat = ip_iat[ip];
    if (npol == 2)
    {
        for (int ib = threadIdx.x; ib < nbands; ib += blockDim.x)
        {
            int ib2 = ib * npol;
            const int psind = ip * npm + ib2;
            const int becpind = ib2 * tnp + ip;
            ps[psind] += lambda_coeff[iat * 4] * becp[becpind] + lambda_coeff[iat * 4 + 2] * becp[becpind + tnp];
            ps[psind + 1] += lambda_coeff[iat * 4 + 1] * becp[becpind] + lambda_coeff[iat * 4 + 3] * becp[becpind + tnp];
        }
    }
    else
    {
        for (int ib = threadIdx.x; ib < nbands; ib += blockDim.x)
        {
            const int psind = ip * npm + ib;
            const int becpind = ib * tnp + ip;
            ps[psind] += lambda_coeff[iat] * becp[becpind];
        }
    }
}

template <typename FPTYPE>
__global__ void onsite_op(const int npm,
                          const int npol,
                          const int* orb_l_iat,
                          const int* ip_iat,
                          const int* ip_m,
                          const int* pot_onsite_begin_iat,
                          const int tnp,
                          const thrust::complex<FPTYPE>* pot_onsite,
                          thrust::complex<FPTYPE>* ps,
                          const thrust::complex<FPTYPE>* becp)
{
    const int ip = blockIdx.x;
    int m1 = ip_m[ip];
    if (m1 >= 0)
    {
        const int nbands = npm / npol;
        int iat = ip_iat[ip];
        const thrust::complex<FPTYPE>* pot_onsite_iat = pot_onsite + pot_onsite_begin_iat[iat];
        int orb_l = orb_l_iat[iat];
        int tlp1 = 2 * orb_l + 1;
        int tlp1_2 = tlp1 * tlp1;
        int ip2_begin = ip - m1;
        int ip2_end = ip - m1 + tlp1;

        if (npol == 2)
        {
            for (int ib = threadIdx.x; ib < nbands; ib += blockDim.x)
            {
                int ib2 = ib * npol;
                const int psind = ip * npm + ib2;
                for (int ip2 = ip2_begin; ip2 < ip2_end; ip2++)
                {
                    const int becpind = ib2 * tnp + ip2;
                    int m2 = ip_m[ip2];
                    const int index_mm = m1 * tlp1 + m2;
                    ps[psind] += pot_onsite_iat[index_mm] * becp[becpind] + pot_onsite_iat[index_mm + tlp1_2 * 2] * becp[becpind + tnp];
                    ps[psind + 1] += pot_onsite_iat[index_mm + tlp1_2 * 1] * becp[becpind]
                                     + pot_onsite_iat[index_mm + tlp1_2 * 3] * becp[becpind + tnp];
                }
            }
        }
        else // npol == 1
        {
            for (int ib = threadIdx.x; ib < nbands; ib += blockDim.x)
            {
                const int psind = ip * npm + ib;
                for (int ip2 = ip2_begin; ip2 < ip2_end; ip2++)
                {
                    const int becpind = ib * tnp + ip2;
                    int m2 = ip_m[ip2];
                    const int index_mm = m1 * tlp1 + m2;
                    ps[psind] += pot_onsite_iat[index_mm] * becp[becpind];
                }
            }
        }
    }
}

template <typename FPTYPE>
void hamilt::onsite_ps_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* dev,
                                                                       const int& npm,
                                                                       const int npol,
                                                                       const int* ip_iat,
                                                                       const int& tnp,
                                                                       const std::complex<FPTYPE>* lambda_coeff,
                                                                       std::complex<FPTYPE>* ps,
                                                                       const std::complex<FPTYPE>* becp)
{
    // denghui implement 20221019
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    onsite_op<FPTYPE>
        <<<tnp, THREADS_PER_BLOCK>>>(npm,
                                     npol,
                                     ip_iat,
                                     tnp,
                                     reinterpret_cast<const thrust::complex<FPTYPE>*>(lambda_coeff),
                                     reinterpret_cast<thrust::complex<FPTYPE>*>(ps),          // array of data
                                     reinterpret_cast<const thrust::complex<FPTYPE>*>(becp)); // array of data

    CHECK_CUDA_SYNC();
}

template <typename FPTYPE>
void hamilt::onsite_ps_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* dev,
                                                                       const int& npm,
                                                                       const int npol,
                                                                       const int* orb_l_iat,
                                                                       const int* ip_iat,
                                                                       const int* ip_m,
                                                                       const int* pot_onsite_begin_iat,
                                                                       const int& tnp,
                                                                       const std::complex<FPTYPE>* pot_onsite,
                                                                       std::complex<FPTYPE>* ps,
                                                                       const std::complex<FPTYPE>* becp)
{
    // denghui implement 20221109
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    onsite_op<FPTYPE>
        <<<tnp, THREADS_PER_BLOCK>>>(npm,
                                     npol,
                                     orb_l_iat,
                                     ip_iat,
                                     ip_m,
                                     pot_onsite_begin_iat,
                                     tnp,
                                     reinterpret_cast<const thrust::complex<FPTYPE>*>(pot_onsite),
                                     reinterpret_cast<thrust::complex<FPTYPE>*>(ps),          // array of data
                                     reinterpret_cast<const thrust::complex<FPTYPE>*>(becp)); // array of data

    CHECK_CUDA_SYNC();
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
}

template struct onsite_ps_op<float, base_device::DEVICE_GPU>;
template struct onsite_ps_op<double, base_device::DEVICE_GPU>;

} // namespace hamilt
