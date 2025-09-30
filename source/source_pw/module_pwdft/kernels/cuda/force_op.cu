#include "source_pw/module_pwdft/kernels/force_op.h"
// #include "source_psi/kernels/device.h"
#include "source_base/module_device/types.h"
#include "source_base/constants.h"

#include <complex>

#include <thrust/complex.h>
#include <cuda_runtime.h>
#include <base/macros/macros.h>
#include <source_base/module_device/device.h>

#define THREADS_PER_BLOCK 256
#define FULL_MASK 0xffffffff
#define WARP_SIZE 32

namespace hamilt {

template <typename FPTYPE>
__forceinline__
__device__
void warp_reduce(FPTYPE & val) {
    for (int offset = 16; offset > 0; offset >>= 1) {
        val += __shfl_down_sync(FULL_MASK, val, offset);
    }
}

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

    int iat = 0, sum = 0;
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
    cal_vkb1_nl<FPTYPE><<<nkb, THREADS_PER_BLOCK>>>(
            npwx,
            vkb_nc,
            nbasis,
            ipol,
            static_cast<const thrust::complex<FPTYPE>>(NEG_IMAG_UNIT), // array of data
            reinterpret_cast<const thrust::complex<FPTYPE>*>(vkb),
            gcar,// array of data
            reinterpret_cast<thrust::complex<FPTYPE>*>(vkb1)); // array of data

    cudaCheckOnDebug();
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
    cal_force_nl<FPTYPE><<<nbands_occ * ntype, THREADS_PER_BLOCK>>>(
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

    cudaCheckOnDebug();
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

    int iat = 0, sum = 0;
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
    cal_force_nl<FPTYPE><<<nbands_occ * ntype, THREADS_PER_BLOCK>>>(
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

    cudaCheckOnDebug();
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
    cal_force_onsite<FPTYPE>
        <<<nbands_occ * ntype, THREADS_PER_BLOCK>>>(wg_nc,
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

    cudaCheckOnDebug();
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
    cal_force_onsite<FPTYPE>
        <<<nbands_occ * ntype, THREADS_PER_BLOCK>>>(wg_nc,
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

    cudaCheckOnDebug();
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
    if(ikb < nkb)
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
    int npwx)
{
    saveVkbValues_<FPTYPE><<<nkb*gcar_zero_count , THREADS_PER_BLOCK>>>(
        gcar_zero_ptrs, 
        reinterpret_cast<const thrust::complex<FPTYPE>*>(vkb_ptr), 
        reinterpret_cast<thrust::complex<FPTYPE>*>(vkb_save_ptr), 
        nkb, 
        npw, 
        ipol,
        npwx);
    cudaCheckOnDebug();
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
    if(ikb < nkb && icount < n_total_gcar_zeros)
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
    revertVkbValues_<FPTYPE><<<nkb*gcar_zero_count , THREADS_PER_BLOCK>>>(
        gcar_zero_ptrs, 
        reinterpret_cast<thrust::complex<FPTYPE>*>(vkb_ptr), 
        reinterpret_cast<const thrust::complex<FPTYPE>*>(vkb_save_ptr), 
        nkb, 
        npw, 
        ipol,
        npwx,
        static_cast<const thrust::complex<FPTYPE>>(coeff));
    cudaCheckOnDebug();
}

template <typename FPTYPE>
__global__ void force_loc_kernel(
    const int nat,
    const int npw,
    const FPTYPE tpiba_omega,
    const int* iat2it,
    const int* ig2gg_d,
    const FPTYPE* gcar_d,
    const FPTYPE* tau_d,
    const thrust::complex<FPTYPE>* aux_d,
    const FPTYPE* vloc_d,
    const int vloc_nc,
    FPTYPE* forcelc_d)
{
    const int iat = blockIdx.x;
    const int tid = threadIdx.x;
    const int warp_id = tid / WARP_SIZE;
    const int lane_id = tid % WARP_SIZE;
    
    if (iat >= nat) return;
    const int it = iat2it[iat]; // get the type of atom

    // Initialize force components
    FPTYPE force_x = 0.0;
    FPTYPE force_y = 0.0;
    FPTYPE force_z = 0.0;
    
    const auto tau_x = tau_d[iat * 3 + 0];
    const auto tau_y = tau_d[iat * 3 + 1];
    const auto tau_z = tau_d[iat * 3 + 2];

    // Process all plane waves in chunks of blockDim.x
    for (int ig = tid; ig < npw; ig += blockDim.x) {
        const auto gcar_x = gcar_d[ig * 3 + 0];
        const auto gcar_y = gcar_d[ig * 3 + 1];
        const auto gcar_z = gcar_d[ig * 3 + 2];

        // Calculate phase factor
        const FPTYPE phase = ModuleBase::TWO_PI * (gcar_x * tau_x +
                                                   gcar_y * tau_y +
                                                   gcar_z * tau_z);
        FPTYPE sinp, cosp;
        sincos(phase, &sinp, &cosp);
        
        // Get vloc value
        const FPTYPE vloc_val = vloc_d[it * vloc_nc + ig2gg_d[ig]];
        
        // Calculate factor
        const auto aux_val = aux_d[ig];
        const FPTYPE factor = vloc_val * (cosp * aux_val.imag() + sinp * aux_val.real());
        
        // Multiply by gcar components
        force_x += gcar_x * factor;
        force_y += gcar_y * factor;
        force_z += gcar_z * factor;
    }
    
    // Warp-level reduction
    warp_reduce<FPTYPE>(force_x);
    warp_reduce<FPTYPE>(force_y);
    warp_reduce<FPTYPE>(force_z);
    
    // First thread in each warp writes to shared memory
    __shared__ FPTYPE warp_sums_x[THREADS_PER_BLOCK / WARP_SIZE]; // 256 threads / 32 = 8 warps
    __shared__ FPTYPE warp_sums_y[THREADS_PER_BLOCK / WARP_SIZE];
    __shared__ FPTYPE warp_sums_z[THREADS_PER_BLOCK / WARP_SIZE];
    
    if (lane_id == 0) {
        warp_sums_x[warp_id] = force_x;
        warp_sums_y[warp_id] = force_y;
        warp_sums_z[warp_id] = force_z;
    }
    
    __syncthreads();
    
    // Final reduction by first warp
    if (warp_id == 0) {
        FPTYPE final_x = (lane_id < blockDim.x/WARP_SIZE) ? warp_sums_x[lane_id] : 0.0;
        FPTYPE final_y = (lane_id < blockDim.x/WARP_SIZE) ? warp_sums_y[lane_id] : 0.0;
        FPTYPE final_z = (lane_id < blockDim.x/WARP_SIZE) ? warp_sums_z[lane_id] : 0.0;
        
        warp_reduce<FPTYPE>(final_x);
        warp_reduce<FPTYPE>(final_y);
        warp_reduce<FPTYPE>(final_z);
        
        if (lane_id == 0) {
            forcelc_d[iat * 3 + 0] = final_x * tpiba_omega;
            forcelc_d[iat * 3 + 1] = final_y * tpiba_omega;
            forcelc_d[iat * 3 + 2] = final_z * tpiba_omega;
        }
    }
}

template <typename FPTYPE>
__global__ void force_ew_kernel(
    const int nat,
    const int npw,
    const int ig_gge0,
    const int* iat2it,
    const FPTYPE* gcar_d,
    const FPTYPE* tau_d,
    const FPTYPE* it_fact_d,
    const thrust::complex<FPTYPE>* aux_d,
    FPTYPE* forceion_d)
{
    const int iat = blockIdx.x;
    const int tid = threadIdx.x;
    const int warp_id = tid / WARP_SIZE;
    const int lane_id = tid % WARP_SIZE;

    if( iat >= nat) return;
    const int it = iat2it[iat]; // get the type of atom
    const FPTYPE it_fact_val = it_fact_d[it]; // Get it_fact value

    // Initialize force components
    FPTYPE force_x = 0.0;
    FPTYPE force_y = 0.0;
    FPTYPE force_z = 0.0;

    const auto tau_x = tau_d[iat * 3 + 0];
    const auto tau_y = tau_d[iat * 3 + 1];
    const auto tau_z = tau_d[iat * 3 + 2];

    for (int ig = tid; ig < npw; ig += blockDim.x) {
        if(ig == ig_gge0)
        { continue; }
        const auto gcar_x = gcar_d[ig * 3 + 0];
        const auto gcar_y = gcar_d[ig * 3 + 1];
        const auto gcar_z = gcar_d[ig * 3 + 2];

        // Calculate phase factor
        const FPTYPE phase = ModuleBase::TWO_PI * (gcar_x * tau_x +
                                                   gcar_y * tau_y +
                                                   gcar_z * tau_z);
        FPTYPE sinp, cosp;
        sincos(phase, &sinp, &cosp);

        // Calculate force contribution
        const FPTYPE sumnb = -cosp * aux_d[ig].imag() + sinp * aux_d[ig].real();
        
        // Multiply by gcar components
        force_x += gcar_x * sumnb;
        force_y += gcar_y * sumnb;
        force_z += gcar_z * sumnb;
    }

    // Warp-level reduction
    warp_reduce<FPTYPE>(force_x);
    warp_reduce<FPTYPE>(force_y);
    warp_reduce<FPTYPE>(force_z);

    // First thread in each warp writes to shared memory
    __shared__ FPTYPE warp_sums_x[THREADS_PER_BLOCK / WARP_SIZE]; // 256 threads / 32 = 8 warps
    __shared__ FPTYPE warp_sums_y[THREADS_PER_BLOCK / WARP_SIZE];
    __shared__ FPTYPE warp_sums_z[THREADS_PER_BLOCK / WARP_SIZE];

    if (lane_id == 0) {
        warp_sums_x[warp_id] = force_x;
        warp_sums_y[warp_id] = force_y;
        warp_sums_z[warp_id] = force_z;
    }

    __syncthreads();

    // Final reduction by first warp
    if (warp_id == 0) {
        FPTYPE final_x = (lane_id < blockDim.x/WARP_SIZE) ? warp_sums_x[lane_id] : 0.0;
        FPTYPE final_y = (lane_id < blockDim.x/WARP_SIZE) ? warp_sums_y[lane_id] : 0.0;
        FPTYPE final_z = (lane_id < blockDim.x/WARP_SIZE) ? warp_sums_z[lane_id] : 0.0;

        warp_reduce<FPTYPE>(final_x);
        warp_reduce<FPTYPE>(final_y);
        warp_reduce<FPTYPE>(final_z);

        if (lane_id == 0) {
            forceion_d[iat * 3 + 0] = final_x * it_fact_val;
            forceion_d[iat * 3 + 1] = final_y * it_fact_val;
            forceion_d[iat * 3 + 2] = final_z * it_fact_val;
        }
    }
}

template <typename FPTYPE>
void cal_force_loc_op<FPTYPE, base_device::DEVICE_GPU>::operator()(
    const int nat,
    const int npw,
    const FPTYPE tpiba_omega,
    const int* iat2it,
    const int* ig2igg,
    const FPTYPE* gcar,
    const FPTYPE* tau,
    const std::complex<FPTYPE>* aux,
    const FPTYPE* vloc,
    const int vloc_nc,
    FPTYPE* forcelc)
{
    force_loc_kernel<FPTYPE>
        <<<nat, THREADS_PER_BLOCK>>>(nat,
                                     npw,
                                     tpiba_omega,
                                     iat2it,
                                     ig2igg,
                                     gcar,
                                     tau,
                                     reinterpret_cast<const thrust::complex<FPTYPE>*>(aux),
                                     vloc,
                                     vloc_nc,
                                     forcelc); // array of data

}

template <typename FPTYPE>
void cal_force_ew_op<FPTYPE, base_device::DEVICE_GPU>::operator()(
    const int nat,
    const int npw,
    const int ig_gge0,
    const int* iat2it,
    const FPTYPE* gcar,
    const FPTYPE* tau,
    const FPTYPE* it_fact,
    const std::complex<FPTYPE>* aux,
    FPTYPE* forceion)
{
    force_ew_kernel<FPTYPE>
        <<<nat, THREADS_PER_BLOCK>>>(nat,
                                     npw,
                                     ig_gge0,
                                     iat2it,
                                     gcar,
                                     tau,
                                     it_fact,
                                     reinterpret_cast<const thrust::complex<FPTYPE>*>(aux),
                                     forceion); // array of data
}


// for revertVkbValues functions instantiation
template void revertVkbValues<double>(const int *gcar_zero_ptrs, std::complex<double> *vkb_ptr, const std::complex<double> *vkb_save_ptr, int nkb, int gcar_zero_count, int npw, int ipol, int npwx, const std::complex<double> coeff);
// for saveVkbValues functions instantiation
template void saveVkbValues<double>(const int *gcar_zero_ptrs, const std::complex<double> *vkb_ptr, std::complex<double> *vkb_save_ptr, int nkb, int gcar_zero_count, int npw, int ipol, int npwx);

template struct cal_vkb1_nl_op<float, base_device::DEVICE_GPU>;
template struct cal_force_nl_op<float, base_device::DEVICE_GPU>;
template struct cal_force_loc_op<float, base_device::DEVICE_GPU>;
template struct cal_force_ew_op<float, base_device::DEVICE_GPU>;

template struct cal_vkb1_nl_op<double, base_device::DEVICE_GPU>;
template struct cal_force_nl_op<double, base_device::DEVICE_GPU>;
template struct cal_force_loc_op<double, base_device::DEVICE_GPU>;
template struct cal_force_ew_op<double, base_device::DEVICE_GPU>;
}  // namespace hamilt
