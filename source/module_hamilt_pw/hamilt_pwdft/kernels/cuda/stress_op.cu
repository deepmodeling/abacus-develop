#include "module_hamilt_pw/hamilt_pwdft/kernels/stress_op.h"
#include "module_psi/kernels/device.h"

#include <complex>
#include <thrust/complex.h>

#include <cuda_runtime.h>

#define THREADS_PER_BLOCK 256
#define FULL_MASK 0xffffffff
#define WARP_SIZE 32

namespace hamilt{

template <typename FPTYPE>
__forceinline__
__device__
void warp_reduce(FPTYPE & val) {
    for (int offset = 16; offset > 0; offset >>= 1) {
        val += __shfl_down_sync(FULL_MASK, val, offset);
    }
}

template <typename FPTYPE>
__global__ void cal_dbecp_noevc_nl(
        const int ipol,
        const int jpol,
        const int npw,
        const int npwx,
        const int ik,
        const FPTYPE tpiba,
        const FPTYPE *gcar,
        const FPTYPE *kvec_c,
        thrust::complex<FPTYPE> *vkbi,
        thrust::complex<FPTYPE> *vkbj,
        thrust::complex<FPTYPE> *vkb,
        thrust::complex<FPTYPE> *vkb1,
        thrust::complex<FPTYPE> *vkb2,
        thrust::complex<FPTYPE> *dbecp_noevc)
{
    int i = blockIdx.x;
    const thrust::complex<FPTYPE>* pvkb0i = vkbi + i * npwx;
    const thrust::complex<FPTYPE>* pvkb0j = vkbj + i * npwx;
    thrust::complex<FPTYPE>* pvkb = nullptr;
    thrust::complex<FPTYPE>* pdbecp_noevc = dbecp_noevc + i * npwx;
    // third term of dbecp_noevc
    //std::complex<FPTYPE>* pvkb = &vkb2(i,0);
    //std::complex<FPTYPE>* pdbecp_noevc = &dbecp_noevc(i, 0);
    FPTYPE qvec[3] = {0, 0, 0};
    for (int ig = threadIdx.x; ig < npw; ig += blockDim.x)
    {
        pvkb = vkb1 + i * npwx;
        qvec[ipol] = gcar[(ik * npwx + ig) * 3 + ipol] + kvec_c[ik * 3 + ipol];
        qvec[jpol] = gcar[(ik * npwx + ig) * 3 + jpol] + kvec_c[ik * 3 + jpol];
        pvkb[ig] += 0.5 * qvec[ipol] * pvkb0j[ig] +
                    0.5 * qvec[jpol] * pvkb0i[ig];
        pdbecp_noevc[ig] -= 2.0 * pvkb[ig];
        if (ipol == jpol) {
            pvkb = vkb + i * npwx;
            pdbecp_noevc[ig] -= pvkb[ig];
        }
        pvkb = vkb2 + i * npwx;
        for (int ii = 0; ii < 3; ii++) {
            qvec[ii] = gcar[(ik * npwx + ig) * 3 + ii] + kvec_c[ik * 3 + ii];
        }
        FPTYPE qvec_norm2 = qvec[0] * qvec[0] + qvec[1] * qvec[1] + qvec[2] * qvec[2];
        FPTYPE qm1 = qvec_norm2 > 1e-16 ? 1.0 / sqrt(qvec_norm2) : 0;
        pdbecp_noevc[ig] -= 2.0 * pvkb[ig] * qvec[ipol] *
                            qvec[jpol] * qm1 *	tpiba;
    } // end ig
}

template <typename FPTYPE>
__global__ void cal_stress_nl(
        const bool multi_proj,
        const int ipol,
        const int jpol,
        const int nkb,
        const int ntype,
        const int spin,
        const int wg_nc,
        const int ik,
        const int deeq_2,
        const int deeq_3,
        const int deeq_4,
        const int *atom_nh,
        const int *atom_na,
        const FPTYPE *d_wg,
        const FPTYPE *deeq,
        const thrust::complex<FPTYPE> *becp,
        const thrust::complex<FPTYPE> *dbecp,
        FPTYPE *stress)
{
    int ib = blockIdx.x / ntype;
    int it = blockIdx.x % ntype;

    int iat = 0, sum = 0;
    for (int ii = 0; ii < it; ii++) {
        iat += atom_na[ii];
        sum += atom_na[ii] * atom_nh[ii];
    }

    FPTYPE stress_var = 0, fac = d_wg[ik * wg_nc + ib] * 1.0;
    const int Nprojs = atom_nh[it];
    for (int ia = 0; ia < atom_na[it]; ia++)
    {
        for (int ii = threadIdx.x; ii < Nprojs * Nprojs; ii += blockDim.x) {
            int ip1 = ii / Nprojs, ip2 = ii % Nprojs;
            if(!multi_proj && ip1 != ip2) {
                continue;
            }
            FPTYPE ps = deeq[((spin * deeq_2 + iat) * deeq_3 + ip1) * deeq_4 + ip2];
            const int inkb1 = sum + ip1;
            const int inkb2 = sum + ip2;
            //out<<"\n ps = "<<ps;
            const FPTYPE dbb = ( conj( dbecp[ ib * nkb + inkb1] ) * becp[ ib * nkb + inkb2] ).real();
            stress_var -= ps * fac * dbb;
        }
        ++iat;
        sum+=Nprojs;
    }//ia
    __syncwarp();
    warp_reduce(stress_var);
    if (threadIdx.x % WARP_SIZE == 0) {
        atomicAdd(stress + ipol * 3 + jpol, stress_var);
    }
}

template <typename FPTYPE>
void cal_dbecp_noevc_nl_op<FPTYPE, psi::DEVICE_GPU>::operator() (
        const psi::DEVICE_GPU *ctx,
        const int &ipol,
        const int &jpol,
        const int &nkb,
        const int &npw,
        const int &npwx,
        const int &ik,
        const FPTYPE &tpiba,
        const FPTYPE *gcar,
        const FPTYPE *kvec_c,
        std::complex<FPTYPE> *vkbi,
        std::complex<FPTYPE> *vkbj,
        std::complex<FPTYPE> *vkb,
        std::complex<FPTYPE> *vkb1,
        std::complex<FPTYPE> *vkb2,
        std::complex<FPTYPE> *dbecp_noevc)
{
    cal_dbecp_noevc_nl<FPTYPE><<<nkb, THREADS_PER_BLOCK>>>(
            ipol,
            jpol,
            npw,
            npwx,
            ik,
            tpiba,
            gcar,
            kvec_c,
            reinterpret_cast<thrust::complex<FPTYPE>*>(vkbi),
            reinterpret_cast<thrust::complex<FPTYPE>*>(vkbj),
            reinterpret_cast<thrust::complex<FPTYPE>*>(vkb),
            reinterpret_cast<thrust::complex<FPTYPE>*>(vkb1),
            reinterpret_cast<thrust::complex<FPTYPE>*>(vkb2),
            reinterpret_cast<thrust::complex<FPTYPE>*>(dbecp_noevc));
}

template <typename FPTYPE>
void cal_stress_nl_op<FPTYPE, psi::DEVICE_GPU>::operator() (
        const psi::DEVICE_GPU *ctx,
        const bool &multi_proj,
        const int &ipol,
        const int &jpol,
        const int &nkb,
        const int &nbands_occ,
        const int &ntype,
        const int &spin,
        const int &wg_nc,
        const int &ik,
        const int &deeq_2,
        const int &deeq_3,
        const int &deeq_4,
        const int *atom_nh,
        const int *atom_na,
        const FPTYPE *d_wg,
        const FPTYPE *deeq,
        const std::complex<FPTYPE> *becp,
        const std::complex<FPTYPE> *dbecp,
        FPTYPE *stress)
{
     cal_stress_nl<FPTYPE><<<nbands_occ * ntype, THREADS_PER_BLOCK>>>(
             multi_proj,
             ipol,
             jpol,
             nkb,
             ntype,
             spin,
             wg_nc,
             ik,
             deeq_2,
             deeq_3,
             deeq_4,
             atom_nh,
             atom_na,
             d_wg,
             deeq,
             reinterpret_cast<const thrust::complex<FPTYPE>*>(becp),
             reinterpret_cast<const thrust::complex<FPTYPE>*>(dbecp),
             stress);// array of data
}

template struct cal_dbecp_noevc_nl_op<float, psi::DEVICE_GPU>;
template struct cal_stress_nl_op<float, psi::DEVICE_GPU>;

template struct cal_dbecp_noevc_nl_op<double, psi::DEVICE_GPU>;
template struct cal_stress_nl_op<double, psi::DEVICE_GPU>;

}  // namespace hamilt