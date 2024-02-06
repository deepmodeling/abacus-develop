#include "module_hamilt_pw/hamilt_pwdft/kernels/vnl_op.h"

#include <complex>

#include <thrust/complex.h>
#include <hip/hip_runtime.h>
#include <base/macros/macros.h>

#define THREADS_PER_BLOCK 256

namespace hamilt {

template <typename FPTYPE>
__device__ FPTYPE _polynomial_interpolation(
        const FPTYPE *table,
        const int &dim1,
        const int &dim2,
        const int &tab_2,
        const int &tab_3,
        const int &table_length,
        const FPTYPE &table_interval,
        const FPTYPE &x)
{
    const FPTYPE position = x / table_interval;
    const int iq = static_cast<int>(position);

    const FPTYPE x0 = position - static_cast<FPTYPE>(iq);
    const FPTYPE x1 = 1.0 - x0;
    const FPTYPE x2 = 2.0 - x0;
    const FPTYPE x3 = 3.0 - x0;
    const FPTYPE y =
            table[(dim1 * tab_2 + dim2) * tab_3 + iq + 0] * x1 * x2 * x3 / 6.0 +
            table[(dim1 * tab_2 + dim2) * tab_3 + iq + 0 + 1] * x0 * x2 * x3 / 2.0 -
            table[(dim1 * tab_2 + dim2) * tab_3 + iq + 0 + 2] * x1 * x0 * x3 / 2.0 +
            table[(dim1 * tab_2 + dim2) * tab_3 + iq + 0 + 3] * x1 * x2 * x0 / 6.0 ;

    return y;
}

template<typename FPTYPE>
__global__ void cal_vnl(
    const int ntype,
    const int npw,
    const int npwx,
    const int nhm,
    const int NQX,
    const int tab_2,
    const int tab_3,
    const int * atom_na,
    const int * atom_nb,
    const int * atom_nh,
    const FPTYPE DQ,
    const FPTYPE tpiba,
    const thrust::complex<FPTYPE> NEG_IMAG_UNIT,
    const FPTYPE *gk,
    const FPTYPE *ylm,
    const FPTYPE *indv,
    const FPTYPE *nhtol,
    const FPTYPE *nhtolm,
    const FPTYPE *tab,
    FPTYPE *vkb1,
    const thrust::complex<FPTYPE> *sk,
    thrust::complex<FPTYPE> *vkb_in)
{
    FPTYPE vq = 0.0;
    int iat = 0, jkb = 0, ig = blockIdx.x * blockDim.x + threadIdx.x;
    if (ig >= npw) {return;}
    for (int it = 0; it < ntype; it++) {
        // calculate beta in G-space using an interpolation table
        const int nh = atom_nh[it];
        const int nbeta = atom_nb[it];

        for (int nb = 0; nb < nbeta; nb++) {
            const FPTYPE gnorm = sqrt(gk[ig * 3 + 0] * gk[ig * 3 + 0] + gk[ig * 3 + 1] * gk[ig * 3 + 1] +
                                      gk[ig * 3 + 2] * gk[ig * 3 + 2]) * tpiba;

            vq = _polynomial_interpolation(
                    tab, it, nb, tab_2, tab_3, NQX, DQ, gnorm);

            // add spherical harmonic part
            for (int ih = 0; ih < nh; ih++) {
                if (nb == indv[it * nhm + ih]) {
                    const int lm = static_cast<int>(nhtolm[it * nhm + ih]);
                    vkb1[ih * npw + ig] = ylm[lm * npw + ig] * vq;
                }
            } // end ih
        } // end nbeta

        // vkb1 contains all betas including angular part for type nt
        // now add the structure factor and factor (-i)^l
        for (int ia = 0; ia < atom_na[it]; ia++) {
            for (int ih = 0; ih < nh; ih++) {
                thrust::complex<FPTYPE> pref = pow(NEG_IMAG_UNIT, nhtol[it * nhm + ih]);    //?
                thrust::complex<FPTYPE> *pvkb = vkb_in + jkb * npwx;
                pvkb[ig] = vkb1[ih * npw + ig] * sk[iat * npw + ig] * pref;
                ++jkb;
            } // end ih
            iat++;
        } // end ia
    } // enddo
}

template <typename FPTYPE>
void cal_vnl_op<FPTYPE, psi::DEVICE_GPU>::operator() (
    const psi::DEVICE_GPU *ctx,
    const int &ntype,
    const int &npw,
    const int &npwx,
    const int &nhm,
    const int &NQX,
    const int &tab_2,
    const int &tab_3,
    const int * atom_na,
    const int * atom_nb,
    const int * atom_nh,
    const FPTYPE &DQ,
    const FPTYPE &tpiba,
    const std::complex<FPTYPE> &NEG_IMAG_UNIT,
    const FPTYPE *gk,
    const FPTYPE *ylm,
    const FPTYPE *indv,
    const FPTYPE *nhtol,
    const FPTYPE *nhtolm,
    const FPTYPE *tab,
    FPTYPE *vkb1,
    const std::complex<FPTYPE> *sk,
    std::complex<FPTYPE> *vkb_in)
{
    int block = (npw + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    hipLaunchKernelGGL(HIP_KERNEL_NAME(cal_vnl<FPTYPE>), dim3(block), dim3(THREADS_PER_BLOCK), 0, 0, 
            ntype, npw, npwx, nhm, NQX, tab_2, tab_3,
            atom_na, atom_nb, atom_nh,
            DQ, tpiba,
            static_cast<thrust::complex<FPTYPE>>(NEG_IMAG_UNIT),
            gk, ylm, indv, nhtol, nhtolm, tab, vkb1,
            reinterpret_cast<const thrust::complex<FPTYPE>*>(sk),
            reinterpret_cast<thrust::complex<FPTYPE>*>(vkb_in));
    
    hipErrcheck(hipGetLastError());
    hipErrcheck(hipDeviceSynchronize());
}

template struct cal_vnl_op<float, psi::DEVICE_GPU>;
template struct cal_vnl_op<double, psi::DEVICE_GPU>;

}  // namespace hamilt