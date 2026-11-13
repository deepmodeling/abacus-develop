#include "source_pw/module_pwdft/kernels/vnl_op.h"
#include "source_base/module_device/memory_op.h"

#include <cstdint>
#include <cfloat>
#include <cmath>
#include "vnl_tools_cu.hpp"

#include <complex>

#include <thrust/complex.h>
#include <cuda_runtime.h>
#include <base/macros/macros.h>

#define THREADS_PER_BLOCK 256

namespace hamilt {

#if __CUDA || __UT_USE_CUDA
constexpr double range_guard_factor = 8.0; // cover small host/device rounding differences
__global__ void cal_gk_kernel(const int count,
                              const double* gcar,
                              const double* kvec,
                              const double shift_x,
                              const double shift_y,
                              const double shift_z,
                              double* gk)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < count)
    {
        const int component = i % 3;
        const double shift = component == 0 ? shift_x : (component == 1 ? shift_y : shift_z);
        gk[i] = gcar[i] + kvec[component] + shift;
    }
}

void cal_gk_gpu(int npw,
                const double* gcar,
                const double* kvec,
                double shift_x,
                double shift_y,
                double shift_z,
                double* gk)
{
    if (npw == 0)
    {
        return;
    }
    const int count = 3 * npw;
    cal_gk_kernel<<<(count + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK>>>(
        count, gcar, kvec, shift_x, shift_y, shift_z, gk);
}

__global__ void gk_range_kernel(const double* gk,
                                const int npw,
                                const double tpiba,
                                const double dq,
                                const double bound,
                                int* invalid)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < npw)
    {
        const double x = gk[3 * i];
        const double y = gk[3 * i + 1];
        const double z = gk[3 * i + 2];
        const double q = std::sqrt(x * x + y * y + z * z) * tpiba;
        const double index = q / dq;
        const double guarded = index * (1.0 + range_guard_factor * DBL_EPSILON);
        if (!std::isfinite(guarded) || guarded < 0.0 || std::floor(guarded) + 3.0 >= bound)
        {
            atomicOr(invalid, 1);
        }
    }
}

bool gk_range_valid_gpu(const double* gk, int npw, double tpiba, double dq, double bound)
{
    if (npw == 0)
    {
        return true;
    }
    int* invalid = nullptr;
    int host_invalid = 0;
    base_device::memory::resize_memory_op<int, base_device::DEVICE_GPU>()(invalid, 1);
    base_device::memory::set_memory_op<int, base_device::DEVICE_GPU>()(invalid, 0, 1);
    gk_range_kernel<<<(npw + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK>>>(
        gk, npw, tpiba, dq, bound, invalid);
    base_device::memory::synchronize_memory_op<int, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(
        &host_invalid, invalid, 1);
    base_device::memory::delete_memory_op<int, base_device::DEVICE_GPU>()(invalid);
    return host_invalid == 0;
}
#endif

template<typename FPTYPE>
__global__ void cal_vnl(
    const int ntype,
    const int npw,
    const int npwx,
    const int nhm,
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
            const std::int64_t g_offset = static_cast<std::int64_t>(ig) * 3;
            const FPTYPE gnorm = sqrt(gk[g_offset + 0] * gk[g_offset + 0] + gk[g_offset + 1] * gk[g_offset + 1] +
                                      gk[g_offset + 2] * gk[g_offset + 2]) * tpiba;

            vq = _polynomial_interpolation(
                    tab, it, nb, tab_2, tab_3, DQ, gnorm);

            // add spherical harmonic part
            for (int ih = 0; ih < nh; ih++) {
                if (nb == indv[static_cast<std::int64_t>(it) * nhm + ih]) {
                    const int lm = static_cast<int>(nhtolm[static_cast<std::int64_t>(it) * nhm + ih]);
                    vkb1[static_cast<std::int64_t>(ih) * npw + ig] = ylm[static_cast<std::int64_t>(lm) * npw + ig] * vq;
                }
            } // end ih
        } // end nbeta

        // vkb1 contains all betas including angular part for type nt
        // now add the structure factor and factor (-i)^l
        for (int ia = 0; ia < atom_na[it]; ia++) {
            for (int ih = 0; ih < nh; ih++) {
                thrust::complex<FPTYPE> pref = pow(NEG_IMAG_UNIT, nhtol[static_cast<std::int64_t>(it) * nhm + ih]);    //?
                thrust::complex<FPTYPE> *pvkb = vkb_in + static_cast<std::int64_t>(jkb) * npwx;
                pvkb[ig] = vkb1[static_cast<std::int64_t>(ih) * npw + ig] * sk[static_cast<std::int64_t>(iat) * npw + ig] * pref;
                ++jkb;
            } // end ih
            iat++;
        } // end ia
    } // enddo
}

template <typename FPTYPE>
void cal_vnl_op<FPTYPE, base_device::DEVICE_GPU>::operator() (
    const base_device::DEVICE_GPU *ctx,
    const int &ntype,
    const int &npw,
    const int &npwx,
    const int &nhm,
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
    int block = npw / THREADS_PER_BLOCK + (npw % THREADS_PER_BLOCK != 0);
    cal_vnl<FPTYPE><<<block, THREADS_PER_BLOCK>>>(
            ntype, npw, npwx, nhm, tab_2, tab_3,
            atom_na, atom_nb, atom_nh,
            DQ, tpiba,
            static_cast<thrust::complex<FPTYPE>>(NEG_IMAG_UNIT),
            gk, ylm, indv, nhtol, nhtolm, tab, vkb1,
            reinterpret_cast<const thrust::complex<FPTYPE>*>(sk),
            reinterpret_cast<thrust::complex<FPTYPE>*>(vkb_in));

    CHECK_CUDA_SYNC();
}

template struct cal_vnl_op<float, base_device::DEVICE_GPU>;
template struct cal_vnl_op<double, base_device::DEVICE_GPU>;

}  // namespace hamilt
