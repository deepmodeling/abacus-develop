#include "module_hamilt_pw/hamilt_pwdft/kernels/stress_op.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/cuda/vnl_tools_cu.hpp"

#include <complex>

#include <thrust/complex.h>
#include <hip/hip_runtime.h>
#include <base/macros/macros.h>

#define THREADS_PER_BLOCK 256
#define FULL_MASK 0xffffffff
#define WARP_SIZE 64

namespace hamilt {

template <typename FPTYPE>
__forceinline__
__device__
void warp_reduce(FPTYPE & val) {
    for (int offset = 32; offset > 0; offset >>= 1) {
        val += __shfl_down( val, offset);//########????
    }
}


template <typename T>
__global__ void cal_stress_mgga(
    const int spin,
    const int nrxx,
    const T w1,
    const thrust::complex<T> * gradwfc,
    T * crosstaus)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= nrxx) { return; }
    int ipol = 0;
    for (int ix = 0; ix < 3; ix++) {
        for (int iy = 0; iy < ix + 1; iy++) {
            crosstaus[spin * nrxx * 6 + ipol * nrxx + idx]
                += 2.0 * w1
                * (gradwfc[ix * nrxx + idx].real() * gradwfc[iy*nrxx + idx].real()
                +  gradwfc[ix * nrxx + idx].imag() * gradwfc[iy*nrxx + idx].imag());
            ipol += 1;
        }
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
        const bool nondiagonal,
        const int ipol,
        const int jpol,
        const int nkb,
        const int ntype,
        const int spin,
        const int deeq_2,
        const int deeq_3,
        const int deeq_4,
        const int *atom_nh,
        const int *atom_na,
        const FPTYPE *d_wg,
        const bool occ,
        const FPTYPE* d_ekb,
        const FPTYPE* qq_nt,
        const FPTYPE *deeq,
        const thrust::complex<FPTYPE> *becp,
        const thrust::complex<FPTYPE> *dbecp,
        FPTYPE *stress)
{
    int ib = blockIdx.x / ntype;
    int it = blockIdx.x % ntype;

    int iat = 0;
    int sum = 0;
    for (int ii = 0; ii < it; ii++) {
        iat += atom_na[ii];
        sum += atom_na[ii] * atom_nh[ii];
    }

    FPTYPE stress_var = 0;
    FPTYPE fac;
    if (occ)
    {
        fac = d_wg[ib];
    }
    else
    {
        fac = d_wg[0];
    }
    FPTYPE ekb_now = 0.0;
    if (d_ekb != nullptr)
    {
        ekb_now = d_ekb[ib];
    }
    const int nproj = atom_nh[it];
    for (int ia = 0; ia < atom_na[it]; ia++)
    {
        for (int ii = threadIdx.x; ii < nproj * nproj; ii += blockDim.x) {
            int ip1 = ii / nproj, ip2 = ii % nproj;
            if(!nondiagonal && ip1 != ip2) {
                continue;
            }
            FPTYPE ps_qq = 0;
            if (ekb_now != 0)
            {
                ps_qq = -ekb_now * qq_nt[it * deeq_3 * deeq_4 + ip1 * deeq_4 + ip2];
            }
            FPTYPE ps = deeq[((spin * deeq_2 + iat) * deeq_3 + ip1) * deeq_4 + ip2] + ps_qq;
            const int inkb1 = sum + ip1;
            const int inkb2 = sum + ip2;
            //out<<"\n ps = "<<ps;
            const FPTYPE dbb = ( conj( dbecp[ ib * nkb + inkb1] ) * becp[ ib * nkb + inkb2] ).real();
            stress_var -= ps * fac * dbb;
        }
        ++iat;
        sum+=nproj;
    }//ia
    warp_reduce(stress_var);
    if (threadIdx.x % WARP_SIZE == 0) {
        atomicAdd(stress + ipol * 3 + jpol, stress_var);
    }
}

template <typename FPTYPE>
void cal_dbecp_noevc_nl_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* ctx,
                                                                        const int& ipol,
                                                                        const int& jpol,
                                                                        const int& nkb,
                                                                        const int& npw,
                                                                        const int& npwx,
                                                                        const int& ik,
                                                                        const FPTYPE& tpiba,
                                                                        const FPTYPE* gcar,
                                                                        const FPTYPE* kvec_c,
                                                                        std::complex<FPTYPE>* vkbi,
                                                                        std::complex<FPTYPE>* vkbj,
                                                                        std::complex<FPTYPE>* vkb,
                                                                        std::complex<FPTYPE>* vkb1,
                                                                        std::complex<FPTYPE>* vkb2,
                                                                        std::complex<FPTYPE>* dbecp_noevc)
{
    hipLaunchKernelGGL(HIP_KERNEL_NAME(cal_dbecp_noevc_nl<FPTYPE>), dim3(nkb), dim3(THREADS_PER_BLOCK), 0, 0,
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

    hipCheckOnDebug();
}

template <typename FPTYPE>
void cal_stress_nl_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* ctx,
                                                                   const bool& nondiagonal,
                                                                   const int& ipol,
                                                                   const int& jpol,
                                                                   const int& nkb,
                                                                   const int& nbands_occ,
                                                                   const int& ntype,
                                                                   const int& spin,
                                                                   const int& deeq_2,
                                                                   const int& deeq_3,
                                                                   const int& deeq_4,
                                                                   const int* atom_nh,
                                                                   const int* atom_na,
                                                                   const FPTYPE* d_wg,
                                                                   const bool& occ,
                                                                   const FPTYPE* d_ekb,
                                                                   const FPTYPE* qq_nt,
                                                                   const FPTYPE* deeq,
                                                                   const std::complex<FPTYPE>* becp,
                                                                   const std::complex<FPTYPE>* dbecp,
                                                                   FPTYPE* stress)
{
     hipLaunchKernelGGL(HIP_KERNEL_NAME(cal_stress_nl<FPTYPE>), dim3(nbands_occ * ntype), dim3(THREADS_PER_BLOCK), 0, 0,
             nondiagonal,
             ipol,
             jpol,
             nkb,
             ntype,
             spin,
             deeq_2,
             deeq_3,
             deeq_4,
             atom_nh,
             atom_na,
             d_wg,
             occ,
             d_ekb,
             qq_nt,
             deeq,
             reinterpret_cast<const thrust::complex<FPTYPE>*>(becp),
             reinterpret_cast<const thrust::complex<FPTYPE>*>(dbecp),
             stress);// array of data

    hipCheckOnDebug();
}

template <typename FPTYPE>
__global__ void cal_stress_nl(
        const int ipol,
        const int jpol,
        const int nkb,
        const int ntype,
        const int deeq_2,
        const int deeq_3,
        const int deeq_4,
        const int *atom_nh,
        const int *atom_na,
        const FPTYPE *d_wg,
        const bool occ,
        const FPTYPE* d_ekb,
        const FPTYPE* qq_nt,
        const thrust::complex<FPTYPE> *deeq_nc,
        const thrust::complex<FPTYPE> *becp,
        const thrust::complex<FPTYPE> *dbecp,
        FPTYPE *stress)
{
    int ib = blockIdx.x / ntype;
    const int ib2  = ib * 2;
    int it = blockIdx.x % ntype;

    int iat = 0;
    int sum = 0;
    for (int ii = 0; ii < it; ii++) {
        iat += atom_na[ii];
        sum += atom_na[ii] * atom_nh[ii];
    }

    FPTYPE stress_var = 0;
    FPTYPE fac;
    if (occ)
    {
        fac = d_wg[ib];
    }
    else
    {
        fac = d_wg[0];
    }
    FPTYPE ekb_now = 0.0;
    if (d_ekb != nullptr)
    {
        ekb_now = d_ekb[ib];
    }
    const int nproj = atom_nh[it];
    for (int ia = 0; ia < atom_na[it]; ia++)
    {
        for (int ii = threadIdx.x; ii < nproj * nproj; ii += blockDim.x) {
            const int ip1 = ii / nproj;
            const int ip2 = ii % nproj;
            thrust::complex<FPTYPE> ps_qq = 0;
            if (ekb_now != 0)
            {
                ps_qq = thrust::complex<FPTYPE>(- ekb_now * qq_nt[it * deeq_3 * deeq_4 + ip1 * deeq_4 + ip2], 0.0);
            }
            const thrust::complex<FPTYPE> ps0 = deeq_nc[((iat + ia) * deeq_3 + ip1) * deeq_4 + ip2] + ps_qq;
            const thrust::complex<FPTYPE> ps1 = deeq_nc[((1 * deeq_2 + iat + ia) * deeq_3 + ip1) * deeq_4 + ip2];
            const thrust::complex<FPTYPE> ps2 = deeq_nc[((2 * deeq_2 + iat + ia) * deeq_3 + ip1) * deeq_4 + ip2];
            const thrust::complex<FPTYPE> ps3 = deeq_nc[((3 * deeq_2 + iat + ia) * deeq_3 + ip1) * deeq_4 + ip2] + ps_qq;
            const int inkb1 = sum + ip1;
            const int inkb2 = sum + ip2;
            //out<<"\n ps = "<<ps;
            const thrust::complex<FPTYPE> dbb0 = conj(dbecp[ib2 * nkb + inkb1]) * becp[ib2 * nkb + inkb2];
            const thrust::complex<FPTYPE> dbb1 = conj(dbecp[ib2 * nkb + inkb1]) * becp[(ib2+1) * nkb + inkb2];
            const thrust::complex<FPTYPE> dbb2 = conj(dbecp[(ib2+1) * nkb + inkb1]) * becp[ib2 * nkb + inkb2];
            const thrust::complex<FPTYPE> dbb3 = conj(dbecp[(ib2+1) * nkb + inkb1]) * becp[(ib2+1) * nkb + inkb2];
            stress_var -= fac * (ps0 * dbb0 + ps1 * dbb1 + ps2 * dbb2 + ps3 * dbb3).real();
        }
        ++iat;
        sum+=nproj;
    }//ia
    warp_reduce(stress_var);
    if (threadIdx.x % WARP_SIZE == 0) {
        atomicAdd(stress + ipol * 3 + jpol, stress_var);
    }
}

template <typename FPTYPE>
__global__ void cal_multi_dot(const int npw,
                              const FPTYPE fac,
                              const FPTYPE* gk1,
                              const FPTYPE* gk2,
                              const FPTYPE* d_kfac,
                              const thrust::complex<FPTYPE>* psi,
                              FPTYPE* sum)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < npw) {
        atomicAdd(sum, fac * gk1[idx] * gk2[idx] * d_kfac[idx] * thrust::norm(psi[idx]));
    }
}

template <typename FPTYPE>
void cal_stress_nl_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* ctx,
                                                                   const int& ipol,
                                                                   const int& jpol,
                                                                   const int& nkb,
                                                                   const int& nbands_occ,
                                                                   const int& ntype,
                                                                   const int& deeq_2,
                                                                   const int& deeq_3,
                                                                   const int& deeq_4,
                                                                   const int* atom_nh,
                                                                   const int* atom_na,
                                                                   const FPTYPE* d_wg,
                                                                   const bool& occ,
                                                                   const FPTYPE* d_ekb,
                                                                   const FPTYPE* qq_nt,
                                                                   const std::complex<FPTYPE>* deeq_nc,
                                                                   const std::complex<FPTYPE>* becp,
                                                                   const std::complex<FPTYPE>* dbecp,
                                                                   FPTYPE* stress)
{
     hipLaunchKernelGGL(HIP_KERNEL_NAME(cal_stress_nl<FPTYPE>), dim3(nbands_occ * ntype), dim3(THREADS_PER_BLOCK), 0, 0,
             ipol,
             jpol,
             nkb,
             ntype,
             deeq_2,
             deeq_3,
             deeq_4,
             atom_nh,
             atom_na,
             d_wg,
             occ,
             d_ekb,
             qq_nt,
             reinterpret_cast<const thrust::complex<FPTYPE>*>(deeq_nc),
             reinterpret_cast<const thrust::complex<FPTYPE>*>(becp),
             reinterpret_cast<const thrust::complex<FPTYPE>*>(dbecp),
             stress);// array of data

    hipCheckOnDebug();
}

template <typename T, typename Device>
void cal_stress_mgga_op<T, Device>::operator()(
    const int& spin,
    const int& nrxx,
    const Real& w1,
    const T * gradwfc,
    Real * crosstaus)
{
    auto gradwfc_ = reinterpret_cast<const thrust::complex<Real>*>(gradwfc);
    const int block = (nrxx + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    cal_stress_mgga<Real><<<block, THREADS_PER_BLOCK>>>(
        spin, nrxx, w1, gradwfc_, crosstaus);

    hipCheckOnDebug();
}




template <typename FPTYPE>
__global__ void cal_vkb(
    const int npw,
    const int* indexes,
    const FPTYPE* vqs_in,
    const FPTYPE* ylms_in,
    const thrust::complex<FPTYPE>* sk_in,
    const thrust::complex<FPTYPE>* pref_in,
    thrust::complex<FPTYPE>* vkbs_out
){
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int ih =  blockIdx.y;

    thrust::complex<FPTYPE>* vkb_ptr = vkbs_out + ih * npw;
    const FPTYPE* ylm_ptr = ylms_in + indexes[ih*4] * npw;
    const FPTYPE* vq_ptr = vqs_in + indexes[ih*4+1] * npw;
    if(idx<npw) vkb_ptr[idx] = ylm_ptr[idx] * vq_ptr[idx] * sk_in[idx] * pref_in[ih];              
    
}

template <typename FPTYPE>
__global__ void cal_vkb_deri(
        const int npw,
        const int ipol,
        const int jpol,
        const int* indexes,
        const FPTYPE* vqs_in, const FPTYPE* vqs_deri_in,
        const FPTYPE* ylms_in, const FPTYPE* ylms_deri_in,
        const thrust::complex<FPTYPE>* sk_in,
        const thrust::complex<FPTYPE>* pref_in,
        const FPTYPE* gk_in,
        thrust::complex<FPTYPE>* vkbs_out
){
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int ih =  blockIdx.y;

    thrust::complex<FPTYPE>* vkb_ptr = vkbs_out + ih * npw;
    const FPTYPE* ylm_ptr = ylms_in + indexes[ih*4] * npw;
    const FPTYPE* vq_ptr = vqs_in + indexes[ih*4 + 1] * npw;

    const FPTYPE* ylm_deri_ptr1 = ylms_deri_in + indexes[ih*4+2] * npw;
    const FPTYPE* ylm_deri_ptr2 = ylms_deri_in + indexes[ih*4+3] * npw;
    const FPTYPE* vq_deri_ptr = vqs_deri_in + indexes[ih*4+1] * npw;
    const FPTYPE* gkn = &gk_in[4 * npw];
    const FPTYPE* gk = &gk_in[idx * 3];

    if(idx<npw) {
        vkb_ptr[idx] = thrust::complex<FPTYPE>(0.0, 0.0);
        if(ipol == jpol)
        {
            vkb_ptr[idx] -= ylm_ptr[idx] * vq_ptr[idx] * sk_in[idx] * pref_in[ih];
        }
        vkb_ptr[idx] -= (gk[ipol] * ylm_deri_ptr2[idx] 
                        + gk[jpol] * ylm_deri_ptr1[idx]) 
                        * vq_ptr[idx] * sk_in[idx] * pref_in[ih];

        vkb_ptr[idx] -= 2.0 * ylm_ptr[idx] * vq_deri_ptr[idx] * sk_in[idx] * pref_in[ih]
                    * gk[ipol] * gk[jpol] * gkn[idx];  
    }
}


template <typename FPTYPE>
__global__ void cal_vq(
        const FPTYPE* tab,
        int it, const FPTYPE* gk, int npw,
        const int tab_2,const int tab_3,  const FPTYPE table_interval, 
        const int nbeta, FPTYPE* vq
){
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int ib =  blockIdx.y;

    FPTYPE* vq_ptr = &vq[ib * npw];
    const FPTYPE* gnorm = &gk[3 * npw];
    if(idx<npw) vq_ptr[idx] = _polynomial_interpolation(
        tab, it, ib, tab_2, tab_3, table_interval, gnorm[idx]);
}

template <typename FPTYPE>
__global__ void cal_vq_deri(
        const FPTYPE* tab,
        int it, const FPTYPE* gk, int npw,
        const int tab_2,const int tab_3,  const FPTYPE table_interval, 
        const int nbeta, FPTYPE* vq
){
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int ib =  blockIdx.y;

    FPTYPE* vq_ptr = &vq[ib * npw];
    const FPTYPE* gnorm = &gk[3 * npw];
    if(idx<npw) vq_ptr[idx] = _polynomial_interpolation_nl(
        tab, it, ib, tab_2, tab_3, table_interval, gnorm[idx]);
}


template <typename FPTYPE>
__global__ void cal_stress_drhoc_aux0(
        const FPTYPE* r, const FPTYPE* rhoc, 
        const FPTYPE *gx_arr, const FPTYPE *rab, FPTYPE *drhocg, 
        const int mesh, const int igl0, const int ngg, const double omega
){
    const double FOUR_PI =  4.0 * 3.14159265358979323846;

    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx >= ngg) {return;}

    FPTYPE rhocg1=0.0;
    FPTYPE gx = gx_arr[idx];

    auto aux = [](FPTYPE r, FPTYPE rhoc, FPTYPE gx, FPTYPE rab) -> FPTYPE{
        return r * rhoc * (r * cos (gx * r) / gx - sin (gx * r) / (gx * gx)) * rab;
    };

    FPTYPE f_0 = aux(r[0],rhoc[0], gx, rab[0]);
    for( int ir = 1 ; ir< mesh - 2; ir+=2)
    {
        rhocg1 += 2 * aux(r[ir],rhoc[ir], gx, rab[ir]) + aux(r[ir+1],rhoc[ir+1], gx, rab[ir+1]);
    }//ir
    FPTYPE f_2 = aux(r[mesh - 2],rhoc[mesh - 2], gx, rab[mesh - 2]);
    FPTYPE f_1 = aux(r[mesh - 1],rhoc[mesh - 1], gx, rab[mesh - 1]);

    rhocg1 += f_2+f_2;
    rhocg1 += rhocg1;
    rhocg1 += f_0 + f_1;
    rhocg1/=3.0;

    drhocg [idx] = FOUR_PI / omega * rhocg1;
}

template <typename FPTYPE>
__global__ void cal_stress_drhoc_aux1(
        const FPTYPE* r, const FPTYPE* rhoc, 
        const FPTYPE *gx_arr, const FPTYPE *rab, FPTYPE *drhocg, 
        const int mesh, const int igl0, const int ngg, const double omega
){
    const double FOUR_PI =  4.0 * 3.14159265358979323846;

    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx >= ngg) {return;}

    FPTYPE rhocg1=0.0;
    FPTYPE gx = gx_arr[idx];

    auto aux = [](FPTYPE r, FPTYPE rhoc, FPTYPE gx, FPTYPE rab) -> FPTYPE{
        return sin (gx * r) / (gx * r) * r * r * rhoc * rab;
    };

    FPTYPE f_0 = r[0] * r[0] * rhoc[0] * rab[0];
    for( int ir = 1 ; ir< mesh - 2; ir+=2)
    {
        rhocg1 += 2 * aux(r[ir],rhoc[ir], gx, rab[ir]) + aux(r[ir+1],rhoc[ir+1], gx, rab[ir+1]);
    }//ir
    
    FPTYPE f_2 = aux(r[mesh - 2],rhoc[mesh - 2], gx, rab[mesh - 2]);
    FPTYPE f_1 = aux(r[mesh - 1],rhoc[mesh - 1], gx, rab[mesh - 1]);

    rhocg1 += f_2+f_2;
    rhocg1 += rhocg1;
    rhocg1 += f_0 + f_1;
    rhocg1/=3.0;

    drhocg [idx] = FOUR_PI * rhocg1 / omega;
}


template <typename FPTYPE>
__global__ void cal_stress_drhoc_aux2(
        const FPTYPE* r, const FPTYPE* rhoc, 
        const FPTYPE *gx_arr, const FPTYPE *rab, FPTYPE *drhocg, 
        const int mesh, const int igl0, const int ngg, const double omega
){


    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx >= ngg) {return;}

    FPTYPE rhocg1=0.0;
    FPTYPE gx = gx_arr[idx];    

    auto aux = [](FPTYPE r, FPTYPE rhoc, FPTYPE gx, FPTYPE rab) -> FPTYPE{
        return r < 1.0e-8 ? rab * rhoc : rab * rhoc * sin(gx * r) / (gx * r);
    };


    FPTYPE f_0 = r[0] * r[0] * rhoc[0] * rab[0];
    for( int ir = 1 ; ir< mesh - 2; ir+=2)
    {
        rhocg1 += 2 * aux(r[ir],rhoc[ir], gx, rab[ir]) + aux(r[ir+1],rhoc[ir+1], gx, rab[ir+1]);
    }//ir
    FPTYPE f_2 = aux(r[mesh - 2],rhoc[mesh - 2], gx, rab[mesh - 2]);
    FPTYPE f_1 = aux(r[mesh - 1],rhoc[mesh - 1], gx, rab[mesh - 1]);
    
    rhocg1 += f_2+f_2;
    rhocg1 += rhocg1;
    rhocg1 += f_0 + f_1;
    rhocg1/=3.0;

    drhocg [idx] = rhocg1;
}


template <typename FPTYPE>
__global__ void cal_stress_drhoc_aux3(
        const FPTYPE* r, const FPTYPE* rhoc, 
        const FPTYPE *gx_arr, const FPTYPE *rab, FPTYPE *drhocg, 
        const int mesh, const int igl0, const int ngg, const double omega
){
    const double FOUR_PI =  4.0 * 3.14159265358979323846;

    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx >= ngg) {return;}

    FPTYPE rhocg1=0.0;
    FPTYPE gx = gx_arr[idx];    
    const FPTYPE pow_gx = gx * gx;

    auto aux = [](FPTYPE r, FPTYPE rhoc, FPTYPE gx, FPTYPE rab) -> FPTYPE{
        return rab * rhoc * (r * cos(gx * r)/gx - sin(gx * r)/(gx * gx));
    };

    FPTYPE f_0 = r[0] * r[0] * rhoc[0] * rab[0];
    for( int ir = 1 ; ir< mesh - 2; ir+=2)
    {
        rhocg1 += 2 * aux(r[ir],rhoc[ir], gx, rab[ir]) + aux(r[ir+1],rhoc[ir+1], gx, rab[ir+1]);
    }//ir
    FPTYPE f_2 = aux(r[mesh - 2],rhoc[mesh - 2], gx, rab[mesh - 2]);
    FPTYPE f_1 = aux(r[mesh - 1],rhoc[mesh - 1], gx, rab[mesh - 1]);
    
    rhocg1 += f_2+f_2;
    rhocg1 += rhocg1;
    rhocg1 += f_0 + f_1;
    rhocg1/=3.0;

    // calculations after Simpson Integral
    const double g2a = pow_gx / 4.0;
    rhocg1 *= FOUR_PI / omega / 2.0 / gx;
    rhocg1 += FOUR_PI / omega * gx_arr[ngg] * exp(-g2a) * (g2a + 1) / (pow_gx*pow_gx);
    drhocg [idx] = rhocg1;
}


template <typename FPTYPE>
__global__ void cal_force_npw(
        const thrust::complex<FPTYPE> *psiv,
        const FPTYPE* gv_x, const FPTYPE* gv_y, const FPTYPE* gv_z,
        const FPTYPE* rhocgigg_vec,
        FPTYPE* force,
        const FPTYPE pos_x, const FPTYPE pos_y, const FPTYPE pos_z,
        const int npw,
        const FPTYPE omega, const FPTYPE tpiba
){
    const double TWO_PI = 2.0 * 3.14159265358979323846;
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int begin_idx = tid * 1024;
    if(begin_idx > npw) return;

    FPTYPE t_force0 = 0;
    FPTYPE t_force1 = 0;
    FPTYPE t_force2 = 0;
    for(int ig = begin_idx; ig<begin_idx+1024 && ig<npw;ig++) {
        const thrust::complex<FPTYPE> psiv_conj = conj(psiv[ig]);

        const FPTYPE arg = TWO_PI * (gv_x[ig] * pos_x + gv_y[ig] * pos_y + gv_z[ig] * pos_z);
        const FPTYPE sinp = sin(arg);
        const FPTYPE cosp = cos(arg);
        const thrust::complex<FPTYPE> expiarg = thrust::complex<FPTYPE>(sinp, cosp);

        const thrust::complex<FPTYPE> tmp_var = psiv_conj * expiarg * tpiba * omega * rhocgigg_vec[ig];

        const thrust::complex<FPTYPE> ipol0 = tmp_var * gv_x[ig];
        t_force0 += ipol0.real();

        const thrust::complex<FPTYPE> ipol1 = tmp_var * gv_y[ig];
        t_force1 += ipol1.real();

        const thrust::complex<FPTYPE> ipol2 = tmp_var * gv_z[ig];
        t_force2 += ipol2.real();
    }
    atomicAdd(&force[0], t_force0);
    atomicAdd(&force[1], t_force1);
    atomicAdd(&force[2], t_force2);
}


template <typename FPTYPE>
void cal_vkb_op<FPTYPE, base_device::DEVICE_GPU>::operator()(
        const base_device::DEVICE_GPU* ctx,
        const int nh,
        const int npw,
        const int* indexes,
        const FPTYPE* vqs_in,
        const FPTYPE* ylms_in,
        const std::complex<FPTYPE>* sk_in,
        const std::complex<FPTYPE>* pref_in,
        std::complex<FPTYPE>* vkbs_out
    )
{
    const int block = (npw + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    dim3 gridsize(block,nh);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(cal_vkb<FPTYPE>),gridsize,THREADS_PER_BLOCK,0,0,
        npw, indexes, vqs_in, ylms_in,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(sk_in), 
        reinterpret_cast<const thrust::complex<FPTYPE>*>(pref_in), 
        reinterpret_cast<thrust::complex<FPTYPE>*>(vkbs_out)
        
    );

}

template <typename FPTYPE>
void cal_vkb_deri_op<FPTYPE, base_device::DEVICE_GPU>::operator()(
        const base_device::DEVICE_GPU* ctx,
        const int nh,
        const int npw,
        const int ipol,
        const int jpol,
        const int* indexes,
        const FPTYPE* vqs_in,
        const FPTYPE* vqs_deri_in,
        const FPTYPE* ylms_in,
        const FPTYPE* ylms_deri_in,
        const std::complex<FPTYPE>* sk_in,
        const std::complex<FPTYPE>* pref_in,
        const FPTYPE* gk_in,
        std::complex<FPTYPE>* vkbs_out
    )
{
    const int block = (npw + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    dim3 gridsize(block,nh);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(cal_vkb_deri<FPTYPE>),gridsize,THREADS_PER_BLOCK,0,0,
        npw, ipol, jpol, indexes,
        vqs_in, vqs_deri_in, ylms_in, ylms_deri_in,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(sk_in), 
        reinterpret_cast<const thrust::complex<FPTYPE>*>(pref_in),       
        gk_in,
        reinterpret_cast<thrust::complex<FPTYPE>*>(vkbs_out)
    );
}

template <typename FPTYPE>
void cal_vq_op<FPTYPE, base_device::DEVICE_GPU>::operator()(
        const base_device::DEVICE_GPU *ctx,
        const FPTYPE* tab,
        int it, const FPTYPE* gk, int npw,
        const int tab_2, const int tab_3, const FPTYPE table_interval, 
        const int nbeta, FPTYPE* vq
    )
{
    const int block = (npw + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    dim3 gridsize(block,nbeta);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(cal_vq<FPTYPE>),gridsize,THREADS_PER_BLOCK,0,0,
        tab, it, gk, npw, tab_2, tab_3,
        table_interval, nbeta, vq
    );
}


template <typename FPTYPE>
void cal_vq_deri_op<FPTYPE, base_device::DEVICE_GPU>::operator()(
        const base_device::DEVICE_GPU *ctx,
        const FPTYPE* tab,
        int it, const FPTYPE* gk, int npw,
        const int tab_2, const int tab_3, const FPTYPE table_interval, 
        const int nbeta, FPTYPE* vq
    )
{
    const int block = (npw + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    dim3 gridsize(block,nbeta);
    hipLaunchKernelGGL(HIP_KERNEL_NAME(cal_vq_deri<FPTYPE>),gridsize,THREADS_PER_BLOCK,0,0,
        tab, it, gk, npw, tab_2, tab_3,
        table_interval, nbeta, vq
    );

    return ;
}

template <typename FPTYPE>
void cal_stress_drhoc_aux_op<FPTYPE, base_device::DEVICE_GPU>::operator()(
        const FPTYPE* r, const FPTYPE* rhoc,  
        const FPTYPE *gx_arr, const FPTYPE *rab, FPTYPE *drhocg, 
        const int mesh, const int igl0, const int ngg, const double omega,
        int type
    )
{
    const int block = (ngg + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;



    if(type == 0) {
        hipLaunchKernelGGL(HIP_KERNEL_NAME(cal_stress_drhoc_aux0<FPTYPE>),block,THREADS_PER_BLOCK,0,0,
            r,rhoc,gx_arr,rab,drhocg,mesh,igl0,ngg,omega
        );
    } else if(type == 1 ){
        hipLaunchKernelGGL(HIP_KERNEL_NAME(cal_stress_drhoc_aux1<FPTYPE>),block,THREADS_PER_BLOCK,0,0,
            r,rhoc,gx_arr,rab,drhocg,mesh,igl0,ngg,omega
        );     
    } else if(type == 2 ){
        hipLaunchKernelGGL(HIP_KERNEL_NAME(cal_stress_drhoc_aux2<FPTYPE>),block,THREADS_PER_BLOCK,0,0,
            r,rhoc,gx_arr,rab,drhocg,mesh,igl0,ngg,omega
        );     
    } else if(type == 3 ){
        hipLaunchKernelGGL(HIP_KERNEL_NAME(cal_stress_drhoc_aux3<FPTYPE>),block,THREADS_PER_BLOCK,0,0,
            r,rhoc,gx_arr,rab,drhocg,mesh,igl0,ngg,omega
        );     
    }

    return ;
}
template <typename FPTYPE>
void cal_force_npw_op<FPTYPE, base_device::DEVICE_GPU>::operator()(
        const std::complex<FPTYPE> *psiv,
        const FPTYPE* gv_x, const FPTYPE* gv_y, const FPTYPE* gv_z,
        const FPTYPE* rhocgigg_vec,
        FPTYPE* force,
        const FPTYPE pos_x, const FPTYPE pos_y, const FPTYPE pos_z,
        const int npw,
        const FPTYPE omega, const FPTYPE tpiba
    )
{
    int t_size = 1024;
    int t_num = (npw%t_size) ? (npw/t_size + 1) : (npw/t_size);

    dim3 npwgrid(((t_num%THREADS_PER_BLOCK) ? (t_num/THREADS_PER_BLOCK + 1) : (t_num/THREADS_PER_BLOCK)));


    hipLaunchKernelGGL(HIP_KERNEL_NAME(cal_force_npw<FPTYPE>), npwgrid, THREADS_PER_BLOCK,0,0,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(psiv),
        gv_x, gv_y, gv_z, rhocgigg_vec, force, pos_x, pos_y, pos_z,
        npw, omega, tpiba
    );

    return ;
}

template <typename FPTYPE>
FPTYPE cal_multi_dot_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const int& npw,
                                                                   const FPTYPE& fac,
                                                                   const FPTYPE* gk1,
                                                                   const FPTYPE* gk2,
                                                                   const FPTYPE* d_kfac,
                                                                   const std::complex<FPTYPE>* psi)
{
    FPTYPE* d_sum = nullptr;
    hipMalloc(&d_sum, sizeof(FPTYPE) * 1);
    hipMemset(d_sum, 0, sizeof(FPTYPE) * 1);
    int block = (npw + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    hipLaunchKernelGGL(HIP_KERNEL_NAME(cal_multi_dot<FPTYPE>), dim3(block), dim3(THREADS_PER_BLOCK), 0, 0,
        npw, fac, gk1, gk2, d_kfac, reinterpret_cast<const thrust::complex<FPTYPE>*>(psi), d_sum);
    FPTYPE sum;
    hipMemcpy(&sum, d_sum, sizeof(FPTYPE) * 1, hipMemcpyDeviceToHost);
    hipFree(d_sum);
    hipCheckOnDebug();
    return sum;
}

template struct cal_vq_op<double, base_device::DEVICE_GPU>;
template struct cal_vq_op<float, base_device::DEVICE_GPU>;

template struct cal_vq_deri_op<double, base_device::DEVICE_GPU>;
template struct cal_vq_deri_op<float, base_device::DEVICE_GPU>;

template struct cal_vkb_op<double, base_device::DEVICE_GPU>;
template struct cal_vkb_op<float, base_device::DEVICE_GPU>;

template struct cal_vkb_deri_op<double, base_device::DEVICE_GPU>;
template struct cal_vkb_deri_op<float, base_device::DEVICE_GPU>;

template struct cal_stress_drhoc_aux_op<double, base_device::DEVICE_GPU>;
template struct cal_stress_drhoc_aux_op<float, base_device::DEVICE_GPU>;

template <>
void pointer_array_malloc<base_device::DEVICE_GPU>::operator()(
        void **ptr,
        const int n
){
    hipErrcheck(hipMalloc(ptr, n * sizeof(void*)));
}

template struct pointer_array_malloc<base_device::DEVICE_GPU>;

template <>
void synchronize_ptrs<base_device::DEVICE_GPU>::operator()(
    void **ptr_out,
    const void **ptr_in,
    const int size)
{
    hipErrcheck(hipMemcpy(ptr_out, ptr_in, sizeof(void*) * size, hipMemcpyHostToDevice));
}

template struct synchronize_ptrs<base_device::DEVICE_GPU>;

template struct cal_stress_mgga_op<std::complex<float>,  base_device::DEVICE_GPU>;
template struct cal_stress_mgga_op<std::complex<double>, base_device::DEVICE_GPU>;

template struct cal_dbecp_noevc_nl_op<float, base_device::DEVICE_GPU>;
template struct cal_dbecp_noevc_nl_op<double, base_device::DEVICE_GPU>;

template struct cal_stress_nl_op<float, base_device::DEVICE_GPU>;
template struct cal_stress_nl_op<double, base_device::DEVICE_GPU>;

template struct cal_force_npw_op<double, base_device::DEVICE_GPU>;
template struct cal_force_npw_op<float, base_device::DEVICE_GPU>;

template struct cal_multi_dot_op<double, base_device::DEVICE_GPU>;
template struct cal_multi_dot_op<float, base_device::DEVICE_GPU>;
}  // namespace hamilt