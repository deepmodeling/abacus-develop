#include "module_hamilt_pw/hamilt_pwdft/kernels/onsite_op.h"

#include <base/macros/macros.h>
#include <complex>
#include <cuda_runtime.h>
#include <thrust/complex.h>

namespace hamilt
{

#define THREADS_PER_BLOCK 256

namespace
{

__global__ void projector_band_moments(const double2* projections, const int* offsets,
                                      int atoms, int bands, int npol, int nkb, double* output)
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= atoms * bands) { return; }
    const int atom = index % atoms, band = index / atoms;
    double up_norm = 0, down_norm = 0, cross_real = 0, cross_imag = 0;
    // Separate rounding preserves the host complex-product arithmetic and
    // projector accumulation order, including cancellation between spin terms.
    for (int projector = offsets[atom]; projector < offsets[atom + 1]; ++projector)
    {
        const int row = band * npol * nkb + projector;
        const auto up = projections[row];
        up_norm = __dadd_rn(up_norm, __dadd_rn(__dmul_rn(up.x, up.x), __dmul_rn(up.y, up.y)));
        if (npol == 2)
        {
            const auto down = projections[row + nkb];
            down_norm = __dadd_rn(down_norm, __dadd_rn(__dmul_rn(down.x, down.x), __dmul_rn(down.y, down.y)));
            cross_real = __dadd_rn(cross_real, __dadd_rn(__dmul_rn(up.x, down.x), __dmul_rn(up.y, down.y)));
            cross_imag = __dadd_rn(cross_imag, __dsub_rn(__dmul_rn(up.x, down.y), __dmul_rn(up.y, down.x)));
        }
    }
    output[3 * index] = __dadd_rn(cross_real, cross_real);
    output[3 * index + 1] = __dadd_rn(cross_imag, cross_imag);
    output[3 * index + 2] = __dsub_rn(up_norm, down_norm);
}

__global__ void sum_projector_moments(const double* moments, const double* weights, const int* spins,
                                     int points, int bands, int atoms, int npol, double* output)
{
    const int component = blockIdx.x * blockDim.x + threadIdx.x;
    if (component >= atoms * 3) { return; }
    double sum = 0;
    for (int point = 0; point < points; ++point)
    {
        for (int band = 0; band < bands; ++band)
        {
            const int state = point * bands + band;
            double contribution = __dmul_rn(weights[state], moments[state * atoms * 3 + component]);
            if (npol == 1) { contribution = __dmul_rn(contribution, spins[point] == 0 ? 1.0 : -1.0); }
            sum = __dadd_rn(sum, contribution);
        }
    }
    output[component] = sum;
}

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
        else { ps[row] = (z * sign) * becp[row]; }
    }
}

} // namespace

void projector_band_moments_gpu(int bands, int npol, int atoms, int nkb, const int* offsets,
                                const std::complex<double>* projections, double* moments)
{
    projector_band_moments<<<(atoms * bands + 127) / 128, 128>>>(
        reinterpret_cast<const double2*>(projections), offsets, atoms, bands, npol, nkb, moments);
    cudaErrcheck(cudaGetLastError());
    cudaCheckOnDebug();
}

void sum_projector_moments_gpu(int points, int bands, int atoms, int npol,
                              const double* moments, const double* weights, const int* spins, double* output)
{
    sum_projector_moments<<<(atoms * 3 + 127) / 128, 128>>>(
        moments, weights, spins, points, bands, atoms, npol, output);
    cudaErrcheck(cudaGetLastError());
    cudaCheckOnDebug();
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
    cudaErrcheck(cudaGetLastError());
    cudaCheckOnDebug();
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
    for (int ib = threadIdx.x; ib < nbands; ib += blockDim.x)
    {
        int ib2 = ib * npol;
        int iat = ip_iat[ip];
        const int psind = ip * npm + ib2;
        const int becpind = ib2 * tnp + ip;
        ps[psind] += lambda_coeff[iat * 4] * becp[becpind] + lambda_coeff[iat * 4 + 2] * becp[becpind + tnp];
        ps[psind + 1] += lambda_coeff[iat * 4 + 1] * becp[becpind] + lambda_coeff[iat * 4 + 3] * becp[becpind + tnp];
    }
}

template <typename FPTYPE>
__global__ void onsite_op(const int npm,
                          const int npol,
                          const int* orb_l_iat,
                          const int* ip_iat,
                          const int* ip_m,
                          const int* vu_begin_iat,
                          const int tnp,
                          const thrust::complex<FPTYPE>* vu,
                          thrust::complex<FPTYPE>* ps,
                          const thrust::complex<FPTYPE>* becp)
{
    const int ip = blockIdx.x;
    int m1 = ip_m[ip];
    if (m1 >= 0)
    {
        const int nbands = npm / npol;
        for (int ib = threadIdx.x; ib < nbands; ib += blockDim.x)
        {
            int ib2 = ib * npol;
            int iat = ip_iat[ip];
            const thrust::complex<FPTYPE>* vu_iat = vu + vu_begin_iat[iat];
            int orb_l = orb_l_iat[iat];
            int tlp1 = 2 * orb_l + 1;
            int tlp1_2 = tlp1 * tlp1;
            int ip2_begin = ip - m1;
            int ip2_end = ip - m1 + tlp1;
            const int psind = ip * npm + ib2;
            for (int ip2 = ip2_begin; ip2 < ip2_end; ip2++)
            {
                const int becpind = ib2 * tnp + ip2;
                int m2 = ip_m[ip2];
                const int index_mm = m1 * tlp1 + m2;
                ps[psind] += vu_iat[index_mm] * becp[becpind] + vu_iat[index_mm + tlp1_2 * 2] * becp[becpind + tnp];
                ps[psind + 1] += vu_iat[index_mm + tlp1_2 * 1] * becp[becpind]
                                 + vu_iat[index_mm + tlp1_2 * 3] * becp[becpind + tnp];
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

    cudaCheckOnDebug();
}

template <typename FPTYPE>
void hamilt::onsite_ps_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* dev,
                                                                       const int& npm,
                                                                       const int npol,
                                                                       const int* orb_l_iat,
                                                                       const int* ip_iat,
                                                                       const int* ip_m,
                                                                       const int* vu_begin_iat,
                                                                       const int& tnp,
                                                                       const std::complex<FPTYPE>* vu,
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
                                     vu_begin_iat,
                                     tnp,
                                     reinterpret_cast<const thrust::complex<FPTYPE>*>(vu),
                                     reinterpret_cast<thrust::complex<FPTYPE>*>(ps),          // array of data
                                     reinterpret_cast<const thrust::complex<FPTYPE>*>(becp)); // array of data

    cudaCheckOnDebug();
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
}

template struct onsite_ps_op<float, base_device::DEVICE_GPU>;
template struct onsite_ps_op<double, base_device::DEVICE_GPU>;

} // namespace hamilt
