#include "module_hamilt_pw/hamilt_pwdft/kernels/force_op.h"

#include <complex>

#include <thrust/complex.h>
#include <hip/hip_runtime.h>
#include <base/macros/macros.h>

#define THREADS_PER_BLOCK 256

namespace hamilt {

template <typename FPTYPE>
__global__ void cal_vkb1_nl(
        const int npwx,
        const int vkb_nc,
        const int nbasis,
        const int ipol,
        const thrust::complex<FPTYPE> NEG_IMAG_UNIT,
        const thrust::complex<FPTYPE> *vkb,
        const FPTYPE *gcar,
        thrust::complex<FPTYPE> *vkb1)
{
    thrust::complex<FPTYPE> *pvkb1 = vkb1 + blockIdx.x * npwx;
    const thrust::complex<FPTYPE> *pvkb = vkb + blockIdx.x * vkb_nc;
    for (int ig = threadIdx.x; ig < nbasis; ig += blockDim.x) {
        pvkb1[ig] = pvkb[ig] * NEG_IMAG_UNIT * gcar[ig * 3 + ipol];
    }
}

template <typename FPTYPE>
__global__ void cal_force_nl(
        const bool nondiagonal,
        const int ntype,
        const int spin,
        const int deeq_2,
        const int deeq_3,
        const int deeq_4,
        const int forcenl_nc,
        const int nbands,
        const int nkb,
        const int *atom_nh,
        const int *atom_na,
        const FPTYPE tpiba,
        const FPTYPE *d_wg,
        const bool occ,
        const FPTYPE* d_ekb,
        const FPTYPE* qq_nt,
        const FPTYPE *deeq,
        const thrust::complex<FPTYPE> *becp,
        const thrust::complex<FPTYPE> *dbecp,
        FPTYPE *force)
{
    const int ib = blockIdx.x / ntype;
    const int it = blockIdx.x % ntype;

    int iat = 0;
    int sum = 0;
    for (int ii = 0; ii < it; ii++) {
        iat += atom_na[ii];
        sum += atom_na[ii] * atom_nh[ii];
    }

    int nproj = atom_nh[it];
    FPTYPE fac;
    if(occ)
    {
        fac = d_wg[ib] * 2.0 * tpiba;
    }
    else
    {
        fac = d_wg[0] * 2.0 * tpiba;
    }
    FPTYPE ekb_now = 0.0;
    if (d_ekb != nullptr)
    {
        ekb_now = d_ekb[ib];
    }
    for (int ia = 0; ia < atom_na[it]; ia++) {
        for (int ip = threadIdx.x; ip < nproj; ip += blockDim.x) {
            FPTYPE ps_qq = 0;
            if(ekb_now != 0)
            {
                ps_qq = - ekb_now * qq_nt[it * deeq_3 * deeq_4 + ip * deeq_4 + ip];
            }
            // FPTYPE ps = GlobalC::ppcell.deeq[spin, iat, ip, ip];
            FPTYPE ps = deeq[((spin * deeq_2 + iat) * deeq_3 + ip) * deeq_4 + ip] + ps_qq;
            const int inkb = sum + ip;
            //out<<"\n ps = "<<ps;

            for (int ipol = 0; ipol < 3; ipol++) {
                const FPTYPE dbb = (conj(dbecp[ipol * nbands * nkb + ib * nkb + inkb]) *
                                    becp[ib * nkb + inkb]).real();
                // force[iat * forcenl_nc + ipol] -= ps * fac * dbb;
                atomicAdd(force + iat * forcenl_nc + ipol, -ps * fac * dbb);
                //cf[iat*3+ipol] += ps * fac * dbb;
            }

            if (nondiagonal) {
                //for (int ip2=0; ip2<nproj; ip2++)
                for (int ip2 = 0; ip2 < nproj; ip2++) {
                    if (ip != ip2) {
                        const int jnkb = sum + ip2;
                        FPTYPE ps_qq = 0;
                        if(ekb_now != 0)
                        {
                            ps_qq = - ekb_now * qq_nt[it * deeq_3 * deeq_4 + ip * deeq_4 + ip2];
                        }
                        ps = deeq[((spin * deeq_2 + iat) * deeq_3 + ip) * deeq_4 + ip2] + ps_qq;
                        for (int ipol = 0; ipol < 3; ipol++) {
                            const FPTYPE dbb = (conj(dbecp[ipol * nbands * nkb + ib * nkb + inkb]) *
                                                becp[ib * nkb + jnkb]).real();
                            atomicAdd(force + iat * forcenl_nc + ipol, -ps * fac * dbb);
                        }
                    }
                }
            }
        }
        iat += 1;
        sum += nproj;
    }
}

template <typename FPTYPE>
void cal_vkb1_nl_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* ctx,
                                                                 const int& nkb,
                                                                 const int& npwx,
                                                                 const int& vkb_nc,
                                                                 const int& nbasis,
                                                                 const int& ipol,
                                                                 const std::complex<FPTYPE>& NEG_IMAG_UNIT,
                                                                 const std::complex<FPTYPE>* vkb,
                                                                 const FPTYPE* gcar,
                                                                 std::complex<FPTYPE>* vkb1)
{
    hipLaunchKernelGGL(HIP_KERNEL_NAME(cal_vkb1_nl<FPTYPE>), dim3(nkb), dim3(THREADS_PER_BLOCK), 0, 0,
            npwx,
            vkb_nc,
            nbasis,
            ipol,
            static_cast<const thrust::complex<FPTYPE>>(NEG_IMAG_UNIT), // array of data
            reinterpret_cast<const thrust::complex<FPTYPE>*>(vkb),
            gcar,// array of data
            reinterpret_cast<thrust::complex<FPTYPE>*>(vkb1)); // array of data

    hipCheckOnDebug();
}

template <typename FPTYPE>
void cal_force_nl_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* ctx,
                                                                  const bool& nondiagonal,
                                                                  const int& nbands_occ,
                                                                  const int& ntype,
                                                                  const int& spin,
                                                                  const int& deeq_2,
                                                                  const int& deeq_3,
                                                                  const int& deeq_4,
                                                                  const int& forcenl_nc,
                                                                  const int& nbands,
                                                                  const int& nkb,
                                                                  const int* atom_nh,
                                                                  const int* atom_na,
                                                                  const FPTYPE& tpiba,
                                                                  const FPTYPE* d_wg,
                                                                  const bool& occ,
                                                                  const FPTYPE* d_ekb,
                                                                  const FPTYPE* qq_nt,
                                                                  const FPTYPE* deeq,
                                                                  const std::complex<FPTYPE>* becp,
                                                                  const std::complex<FPTYPE>* dbecp,
                                                                  FPTYPE* force)
{
    hipLaunchKernelGGL(HIP_KERNEL_NAME(cal_force_nl<FPTYPE>), dim3(nbands_occ * ntype), dim3(THREADS_PER_BLOCK), 0, 0,
            nondiagonal,
            ntype, spin,
            deeq_2, deeq_3, deeq_4,
            forcenl_nc, nbands, nkb,
            atom_nh, atom_na,
            tpiba,
            d_wg, occ, d_ekb, qq_nt, deeq,
            reinterpret_cast<const thrust::complex<FPTYPE>*>(becp),
            reinterpret_cast<const thrust::complex<FPTYPE>*>(dbecp),
            force);// array of data

    hipCheckOnDebug();
}

template <typename FPTYPE>
__global__ void cal_force_nl(
        const int ntype,
        const int deeq_2,
        const int deeq_3,
        const int deeq_4,
        const int forcenl_nc,
        const int nbands,
        const int nkb,
        const int *atom_nh,
        const int *atom_na,
        const FPTYPE tpiba,
        const FPTYPE *d_wg,
        const bool occ,
        const FPTYPE* d_ekb,
        const FPTYPE* qq_nt,
        const thrust::complex<FPTYPE> *deeq_nc,
        const thrust::complex<FPTYPE> *becp,
        const thrust::complex<FPTYPE> *dbecp,
        FPTYPE *force)
{
    const int ib = blockIdx.x / ntype;
    const int ib2 = ib * 2;
    const int it = blockIdx.x % ntype;

    int iat = 0;
    int sum = 0;
    for (int ii = 0; ii < it; ii++) {
        iat += atom_na[ii];
        sum += atom_na[ii] * atom_nh[ii];
    }

    int nproj = atom_nh[it];
    FPTYPE fac;
    if(occ)
    {
        fac = d_wg[ib] * 2.0 * tpiba;
    }
    else
    {
        fac = d_wg[0] * 2.0 * tpiba;
    }
    FPTYPE ekb_now = 0.0;
    if (d_ekb != nullptr)
    {
        ekb_now = d_ekb[ib];
    }
    for (int ia = 0; ia < atom_na[it]; ia++) {
        for (int ip = threadIdx.x; ip < nproj; ip += blockDim.x) {
            const int inkb = sum + ip;
            for (int ip2 = 0; ip2 < nproj; ip2++) 
            {
                // Effective values of the D-eS coefficients
                thrust::complex<FPTYPE> ps_qq = 0;
                if (ekb_now)
                {
                    ps_qq = thrust::complex<FPTYPE>(-ekb_now * qq_nt[it * deeq_3 * deeq_4 + ip * deeq_4 + ip2], 0.0);
                }
                const int jnkb = sum + ip2;
                const thrust::complex<FPTYPE> ps0 = deeq_nc[((0 * deeq_2 + iat) * deeq_3 + ip) * deeq_4 + ip2] + ps_qq;
                const thrust::complex<FPTYPE> ps1 = deeq_nc[((1 * deeq_2 + iat) * deeq_3 + ip) * deeq_4 + ip2];
                const thrust::complex<FPTYPE> ps2 = deeq_nc[((2 * deeq_2 + iat) * deeq_3 + ip) * deeq_4 + ip2];
                const thrust::complex<FPTYPE> ps3 = deeq_nc[((3 * deeq_2 + iat) * deeq_3 + ip) * deeq_4 + ip2] + ps_qq;

                for (int ipol = 0; ipol < 3; ipol++) {
                    const int index0 = ipol * nbands * 2 * nkb + ib2 * nkb + inkb;
                    const int index1 = ib2 * nkb + jnkb;
                    const thrust::complex<FPTYPE> dbb0 = conj(dbecp[index0]) * becp[index1];
                    const thrust::complex<FPTYPE> dbb1 = conj(dbecp[index0]) * becp[index1 + nkb];
                    const thrust::complex<FPTYPE> dbb2 = conj(dbecp[index0 + nkb]) * becp[index1];
                    const thrust::complex<FPTYPE> dbb3 = conj(dbecp[index0 + nkb]) * becp[index1 + nkb];
                    const FPTYPE tmp = - fac * (ps0 * dbb0 + ps1 * dbb1 + ps2 * dbb2 + ps3 * dbb3).real();
                    atomicAdd(force + iat * forcenl_nc + ipol, tmp);
                }
            }
        }
        iat += 1;
        sum += nproj;
    }
}

// interface for nspin=4 only
template <typename FPTYPE>
void cal_force_nl_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* ctx,
                    const int& nbands_occ,
                    const int& ntype,
                    const int& deeq_2,
                    const int& deeq_3,
                    const int& deeq_4,
                    const int& forcenl_nc,
                    const int& nbands,
                    const int& nkb,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE& tpiba,
                    const FPTYPE* d_wg,
                    const bool& occ,
                    const FPTYPE* d_ekb,
                    const FPTYPE* qq_nt,
                    const std::complex<FPTYPE>* deeq_nc,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* force)
{
    hipLaunchKernelGGL(HIP_KERNEL_NAME(cal_force_nl<FPTYPE>), dim3(nbands_occ * ntype), dim3(THREADS_PER_BLOCK), 0, 0,
            ntype,
            deeq_2, deeq_3, deeq_4,
            forcenl_nc, nbands, nkb,
            atom_nh, atom_na,
            tpiba,
            d_wg, occ, d_ekb, qq_nt, 
            reinterpret_cast<const thrust::complex<FPTYPE>*>(deeq_nc),
            reinterpret_cast<const thrust::complex<FPTYPE>*>(becp),
            reinterpret_cast<const thrust::complex<FPTYPE>*>(dbecp),
            force);// array of data

    hipCheckOnDebug();
}

template <typename FPTYPE>
__global__ void cal_force_onsite(int wg_nc,
                                  int ntype,
                                  int forcenl_nc,
                                  int nbands,
                                  int ik,
                                  int nkb,
                                  const int* atom_nh,
                                  const int* atom_na,
                                  int tpiba,
                                  const FPTYPE* d_wg,
                                  const thrust::complex<FPTYPE>* vu,
                                  const int* orbital_corr,
                                  const thrust::complex<FPTYPE>* becp,
                                  const thrust::complex<FPTYPE>* dbecp,
                                  FPTYPE* force)
{
    const int ib = blockIdx.x / ntype; // index of loop-nbands
    const int ib2 = ib * 2;
    const int it = blockIdx.x % ntype; // index of loop-ntype
    if (orbital_corr[it] == -1)
        return;
    const int orbital_l = orbital_corr[it];
    const int ip_begin = orbital_l * orbital_l;
    const int tlp1 = 2 * orbital_l + 1;
    const int tlp1_2 = tlp1 * tlp1;

    int iat = 0; // calculate the begin of atomic index
    int sum = 0; // calculate the begin of atomic-orbital index
    for (int ii = 0; ii < it; ii++)
    {
        iat += atom_na[ii];
        sum += atom_na[ii] * atom_nh[ii];
        vu += 4 * tlp1_2 * atom_na[ii]; // step for vu
    }

    const FPTYPE fac = d_wg[ik * wg_nc + ib] * 2.0 * tpiba;
    const int nprojs = atom_nh[it];
    for (int ia = 0; ia < atom_na[it]; ia++)
    {
        for (int mm = threadIdx.x; mm < tlp1_2; mm += blockDim.x)
        {
            const int m1 = mm / tlp1;
            const int m2 = mm % tlp1;
            const int ip1 = ip_begin + m1;
            const int ip2 = ip_begin + m2;
            const int inkb1 = sum + ip1 + ib2 * nkb;
            const int inkb2 = sum + ip2 + ib2 * nkb;
            thrust::complex<FPTYPE> ps[4] = {vu[mm], vu[mm + tlp1_2], vu[mm + 2 * tlp1_2], vu[mm + 3 * tlp1_2]};
            // out<<"\n ps = "<<ps;
            for (int ipol = 0; ipol < 3; ipol++)
            {
                const int inkb0 = ipol * nbands * 2 * nkb + inkb1;
                const thrust::complex<FPTYPE> dbb0 = conj(dbecp[inkb0]) * becp[inkb2];
                const thrust::complex<FPTYPE> dbb1 = conj(dbecp[inkb0]) * becp[inkb2 + nkb];
                const thrust::complex<FPTYPE> dbb2 = conj(dbecp[inkb0 + nkb]) * becp[inkb2];
                const thrust::complex<FPTYPE> dbb3 = conj(dbecp[inkb0 + nkb]) * becp[inkb2 + nkb];
                const FPTYPE tmp = -fac * (ps[0] * dbb0 + ps[1] * dbb1 + ps[2] * dbb2 + ps[3] * dbb3).real();
                atomicAdd(force + iat * forcenl_nc + ipol, tmp);
            }
        }
        ++iat;
        sum += nprojs;
        vu += 4 * tlp1_2;
    } // ia
}

template <typename FPTYPE>
__global__ void cal_force_onsite(int wg_nc,
                                 int ntype,
                                 int forcenl_nc,
                                 int nbands,
                                 int ik,
                                 int nkb,
                                 const int* atom_nh,
                                 const int* atom_na,
                                 int tpiba,
                                 const FPTYPE* d_wg,
                                 const FPTYPE* lambda,
                                 const thrust::complex<FPTYPE>* becp,
                                 const thrust::complex<FPTYPE>* dbecp,
                                 FPTYPE* force)
{
    const int ib = blockIdx.x / ntype; // index of loop-nbands
    const int ib2 = ib * 2;
    const int it = blockIdx.x % ntype; // index of loop-ntype

    int iat = 0; // calculate the begin of atomic index
    int sum = 0; // calculate the begin of atomic-orbital index
    for (int ii = 0; ii < it; ii++)
    {
        iat += atom_na[ii];
        sum += atom_na[ii] * atom_nh[ii];
    }

    const FPTYPE fac = d_wg[ik * wg_nc + ib] * 2.0 * tpiba;
    const int nprojs = atom_nh[it];
    for (int ia = 0; ia < atom_na[it]; ia++)
    {
        const thrust::complex<FPTYPE> coefficients0(lambda[iat * 3 + 2], 0.0);
        const thrust::complex<FPTYPE> coefficients1(lambda[iat * 3], lambda[iat * 3 + 1]);
        const thrust::complex<FPTYPE> coefficients2(lambda[iat * 3], -1 * lambda[iat * 3 + 1]);
        const thrust::complex<FPTYPE> coefficients3(-1 * lambda[iat * 3 + 2], 0.0);
        for (int ip = threadIdx.x; ip < nprojs; ip += blockDim.x)
        {
            const int inkb = sum + ip + ib2 * nkb;
            // out<<"\n ps = "<<ps;
            for (int ipol = 0; ipol < 3; ipol++)
            {
                const int inkb0 = ipol * nbands * 2 * nkb + inkb;
                const thrust::complex<FPTYPE> dbb0 = conj(dbecp[inkb0]) * becp[inkb];
                const thrust::complex<FPTYPE> dbb1 = conj(dbecp[inkb0]) * becp[inkb + nkb];
                const thrust::complex<FPTYPE> dbb2 = conj(dbecp[inkb0 + nkb]) * becp[inkb];
                const thrust::complex<FPTYPE> dbb3 = conj(dbecp[inkb0 + nkb]) * becp[inkb + nkb];
                const FPTYPE tmp
                    = -fac
                      * (coefficients0 * dbb0 + coefficients1 * dbb1 + coefficients2 * dbb2 + coefficients3 * dbb3)
                            .real();
                atomicAdd(force + iat * forcenl_nc + ipol, tmp);
            }
        }
        ++iat;
        sum += nprojs;
    } // ia
}

// kernel for DFTU force
template <typename FPTYPE>
void cal_force_nl_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* ctx,
                                                                  const int& nbands_occ,
                                                                  const int& wg_nc,
                                                                  const int& ntype,
                                                                  const int& forcenl_nc,
                                                                  const int& nbands,
                                                                  const int& ik,
                                                                  const int& nkb,
                                                                  const int* atom_nh,
                                                                  const int* atom_na,
                                                                  const FPTYPE& tpiba,
                                                                  const FPTYPE* d_wg,
                                                                  const std::complex<FPTYPE>* vu,
                                                                  const int* orbital_corr,
                                                                  const std::complex<FPTYPE>* becp,
                                                                  const std::complex<FPTYPE>* dbecp,
                                                                  FPTYPE* force)
{
    hipLaunchKernelGGL(HIP_KERNEL_NAME(cal_force_onsite<FPTYPE>),
                       dim3(nbands_occ * ntype),
                       dim3(THREADS_PER_BLOCK),
                       0,
                       0,
                       wg_nc,
                       ntype,
                       forcenl_nc,
                       nbands,
                       ik,
                       nkb,
                       atom_nh,
                       atom_na,
                       tpiba,
                       d_wg,
                       reinterpret_cast<const thrust::complex<FPTYPE>*>(vu),
                       orbital_corr,
                       reinterpret_cast<const thrust::complex<FPTYPE>*>(becp),
                       reinterpret_cast<const thrust::complex<FPTYPE>*>(dbecp),
                       force); // array of data

    hipCheckOnDebug();
}
// kernel for DeltaSpin force
template <typename FPTYPE>
void cal_force_nl_op<FPTYPE, base_device::DEVICE_GPU>::operator()(const base_device::DEVICE_GPU* ctx,
                                                                  const int& nbands_occ,
                                                                  const int& wg_nc,
                                                                  const int& ntype,
                                                                  const int& forcenl_nc,
                                                                  const int& nbands,
                                                                  const int& ik,
                                                                  const int& nkb,
                                                                  const int* atom_nh,
                                                                  const int* atom_na,
                                                                  const FPTYPE& tpiba,
                                                                  const FPTYPE* d_wg,
                                                                  const FPTYPE* lambda,
                                                                  const std::complex<FPTYPE>* becp,
                                                                  const std::complex<FPTYPE>* dbecp,
                                                                  FPTYPE* force)
{
    hipLaunchKernelGGL(HIP_KERNEL_NAME(cal_force_onsite<FPTYPE>),
                       dim3(nbands_occ * ntype),
                       dim3(THREADS_PER_BLOCK),
                       0,
                       0,
                       wg_nc,
                       ntype,
                       forcenl_nc,
                       nbands,
                       ik,
                       nkb,
                       atom_nh,
                       atom_na,
                       tpiba,
                       d_wg,
                       lambda,
                       reinterpret_cast<const thrust::complex<FPTYPE>*>(becp),
                       reinterpret_cast<const thrust::complex<FPTYPE>*>(dbecp),
                       force); // array of data

    hipCheckOnDebug();
}

template <typename FPTYPE>
__global__ void saveVkbValues_(
    const int *gcar_zero_ptrs, 
    const thrust::complex<FPTYPE> *vkb_ptr, 
    thrust::complex<FPTYPE> *vkb_save_ptr, 
    int nkb,  
    int npw,  
    int ipol,
    int npwx)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x; // global index
    int n_total_gcar_zeros = gcar_zero_ptrs[ipol * npwx];
    const int* gcar_zero_ptr = gcar_zero_ptrs + ipol * npwx + 1; // skip the first element
    int ikb = index / n_total_gcar_zeros;              // index of nkb
    int icount = index % n_total_gcar_zeros;           // index of n_total_gcar_zeros
    
    // check if the index is valid
    if(ikb < nkb )
    {
        int ig = gcar_zero_ptr[icount]; // get ig from gcar_zero_ptrs
        // use the flat index to get the saved position, pay attention to the relationship between ikb and npw,
        vkb_save_ptr[index] = vkb_ptr[ikb * npw + ig];    // save the value
    }
}

template <typename FPTYPE>
void saveVkbValues(
    const int *gcar_zero_ptrs, 
    const std::complex<FPTYPE> *vkb_ptr, 
    std::complex<FPTYPE> *vkb_save_ptr, 
    int nkb, 
    int gcar_zero_count,
    int npw,  
    int ipol,
    int npwx 
    )
{
    hipLaunchKernelGGL(HIP_KERNEL_NAME(saveVkbValues_<FPTYPE>), dim3(nkb*gcar_zero_count), dim3(THREADS_PER_BLOCK), 0, 0,
        gcar_zero_ptrs, 
        reinterpret_cast<const thrust::complex<FPTYPE>*>(vkb_ptr), 
        reinterpret_cast<thrust::complex<FPTYPE>*>(vkb_save_ptr), 
        nkb,  
        npw,  
        ipol,
        npwx);
    hipCheckOnDebug();
}

template <typename FPTYPE>
__global__ void revertVkbValues_(
    const int *gcar_zero_ptrs, 
    thrust::complex<FPTYPE> *vkb_ptr, 
    const thrust::complex<FPTYPE> *vkb_save_ptr, 
    int nkb, 
    int npw,  
    int ipol,
    int npwx,
    const thrust::complex<FPTYPE> coeff)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x; // global index
    int n_total_gcar_zeros = gcar_zero_ptrs[ipol * npwx];
    const int* gcar_zero_ptr = gcar_zero_ptrs + ipol * npwx + 1; // skip the first element
    int ikb = index / n_total_gcar_zeros;              // index of nkb
    int icount = index % n_total_gcar_zeros;           // index of n_total_gcar_zeros
    
    // check if the index is valid
    if(ikb < nkb)
    {
        int ig = gcar_zero_ptr[icount]; // get ig from gcar_zero_ptrs
        // use the flat index to get the saved position, pay attention to the relationship between ikb and npw,
        vkb_ptr[ikb * npw + ig] = vkb_save_ptr[index] * coeff;    // revert the values
    }
}

template <typename FPTYPE>
void revertVkbValues(
    const int *gcar_zero_ptrs, 
    std::complex<FPTYPE> *vkb_ptr, 
    const std::complex<FPTYPE> *vkb_save_ptr, 
    int nkb, 
    int gcar_zero_count,
    int npw,  
    int ipol,
    int npwx,
    const std::complex<FPTYPE> coeff)
{
    hipLaunchKernelGGL(HIP_KERNEL_NAME(revertVkbValues_<FPTYPE>), dim3(nkb*gcar_zero_count), dim3(THREADS_PER_BLOCK), 0, 0,
        gcar_zero_ptrs, 
        reinterpret_cast<thrust::complex<FPTYPE>*>(vkb_ptr), 
        reinterpret_cast<const thrust::complex<FPTYPE>*>(vkb_save_ptr), 
        nkb, 
        npw, 
        ipol,
        npwx,
        static_cast<const thrust::complex<FPTYPE>>(coeff));
    hipCheckOnDebug();
}

// for revertVkbValues functions instantiation
template void revertVkbValues<double>(const int *gcar_zero_ptrs, std::complex<double> *vkb_ptr, const std::complex<double> *vkb_save_ptr, int nkb, int gcar_zero_count, int npw, int ipol, int npwx, const std::complex<double> coeff);
// for saveVkbValues functions instantiation
template void saveVkbValues<double>(const int *gcar_zero_ptrs, const std::complex<double> *vkb_ptr, std::complex<double> *vkb_save_ptr, int nkb, int gcar_zero_count, int npw, int ipol, int npwx);

template struct cal_vkb1_nl_op<float, base_device::DEVICE_GPU>;
template struct cal_force_nl_op<float, base_device::DEVICE_GPU>;

template struct cal_vkb1_nl_op<double, base_device::DEVICE_GPU>;
template struct cal_force_nl_op<double, base_device::DEVICE_GPU>;

}  // namespace hamilt
