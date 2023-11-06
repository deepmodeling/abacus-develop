#include "module_hamilt_pw/hamilt_pwdft/kernels/wf_op.h"

#include <complex>

#include <thrust/complex.h>
#include <cuda_runtime.h>

#define THREADS_PER_BLOCK 256

namespace hamilt {

template<typename FPTYPE>
__global__ void cal_sk(
    const int ik,
    const int ntype,
    const int nx,
    const int ny,
    const int nz,
    const int rho_nx,
    const int rho_ny,
    const int rho_nz,
    const int npw,
    const int npwx,
    const int fftny,
    const int eigts1_nc,
    const int eigts2_nc,
    const int eigts3_nc,
    const int * atom_na,
    const int * igl2isz,
    const int * is2fftixy,
    const FPTYPE TWO_PI,
    const FPTYPE *kvec_c,
    const FPTYPE *atom_tau,
    thrust::complex<FPTYPE> *eigts1,
    thrust::complex<FPTYPE> *eigts2,
    thrust::complex<FPTYPE> *eigts3,
    thrust::complex<FPTYPE> *sk)
{
    int iat = 0, igl = blockIdx.x * blockDim.x + threadIdx.x;
    if (igl >= npw) {return;}
    for (int it = 0; it < ntype; it++) {
        for (int ia = 0; ia < atom_na[it]; ia++) {
            FPTYPE arg = 0.0;
            for (int ii = 0; ii < 3; ii++) {
                arg += kvec_c[ik * 3 + ii] * atom_tau[iat * 3 + ii];
            }
            arg *= TWO_PI;
            const thrust::complex<FPTYPE> kphase = thrust::complex<FPTYPE>(cos(arg), -sin(arg));
            const int isz = igl2isz[ik * npwx + igl];
            int iz = isz % nz;
            const int is = isz / nz;
            const int ixy = is2fftixy[is];
            int ix = ixy / fftny;
            int iy = ixy % fftny;
            if (ix >= int(nx / 2) + 1)
                ix -= nx;
            if (iy >= int(ny / 2) + 1)
                iy -= ny;
            if (iz >= int(nz / 2) + 1)
                iz -= nz;
            ix += rho_nx;
            iy += rho_ny;
            iz += rho_nz;
            sk[iat * npw + igl] = kphase * eigts1[iat * eigts1_nc + ix] * eigts2[iat * eigts2_nc + iy]
                                  * eigts3[iat * eigts3_nc + iz];
            iat++;
        }
    }
}

template <typename FPTYPE>
void cal_sk_op<FPTYPE, psi::DEVICE_GPU>::operator() (
    const psi::DEVICE_GPU *ctx,
    const int &ik,
    const int &ntype,
    const int &nx,
    const int &ny,
    const int &nz,
    const int& rho_nx,
    const int& rho_ny,
    const int& rho_nz,
    const int &npw,
    const int &npwx,
    const int &fftny,
    const int &eigts1_nc,
    const int &eigts2_nc,
    const int &eigts3_nc,
    const int * atom_na,
    const int * igl2isz,
    const int * is2fftixy,
    const FPTYPE &TWO_PI,
    const FPTYPE *kvec_c,
    const FPTYPE *atom_tau,
    std::complex<FPTYPE> *eigts1,
    std::complex<FPTYPE> *eigts2,
    std::complex<FPTYPE> *eigts3,
    std::complex<FPTYPE> *sk)
{
    int block = (npw + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    cal_sk<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
         ik, ntype,
         nx, ny, nz,
         rho_nx, rho_ny, rho_nz,
         npw, npwx,
         fftny,
         eigts1_nc, eigts2_nc, eigts3_nc,
         atom_na, igl2isz, is2fftixy,
         TWO_PI,
         kvec_c,
         atom_tau,
         reinterpret_cast<thrust::complex<FPTYPE>*>(eigts1),
         reinterpret_cast<thrust::complex<FPTYPE>*>(eigts2),
         reinterpret_cast<thrust::complex<FPTYPE>*>(eigts3),
         reinterpret_cast<thrust::complex<FPTYPE>*>(sk));
}

template struct cal_sk_op<float, psi::DEVICE_GPU>;
template struct cal_sk_op<double, psi::DEVICE_GPU>;

}  // namespace hamilt