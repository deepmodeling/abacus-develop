#include "source_hamilt/module_xc/ncgga_sf_gpu.h"

#include "source_base/module_device/device_check.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/vector3.h"
#include "source_hamilt/module_xc/xc_math.h"

#include <base/macros/macros.h>
#include <cuda_runtime.h>

#include <algorithm>
#include <complex>

namespace ModuleXC
{
namespace
{
constexpr int threads = 256;
constexpr double vanishing = 1e-10;
constexpr double epsr = 1e-6;
constexpr double small = 1e-10;

__device__ double jacobian_value(const NcggaSpinMapPoint& point, int spin, int channel)
{
    if (channel == 0)
    {
        if (point.saturated)
        {
            return spin == 0 ? point.density_sign : 0.0;
        }
        return 0.5 * point.density_sign;
    }
    if (point.saturated)
    {
        return 0.0;
    }
    const double spin_sign = spin == 0 ? 0.5 : -0.5;
    return spin_sign * point.radial.gradient[channel - 1];
}

__device__ double dot3(const double* a, const double* b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

__device__ void local_pbe(double rho, double zeta, double& exc, double& vup, double& vdw)
{
    double ex = 0.0;
    double vxup = 0.0;
    double vxdw = 0.0;
    xc_math::slater_spin(rho, zeta, ex, vxup, vxdw);
    exc = ex;
    vup = vxup;
    vdw = vxdw;
    const double pi34 = 0.62035049089940;
    const double rs = pi34 / cbrt(rho);
    double ec = 0.0;
    double vcup = 0.0;
    double vcdw = 0.0;
    xc_math::pw_spin(rs, zeta, ec, vcup, vcdw);
    exc += ec;
    vup += vcup;
    vdw += vcdw;
}

__device__ void gcx_pbe(double rhoup, double rhodw, double grhoup2, double grhodw2,
                        double& sx, double& v1xup, double& v1xdw, double& v2xup, double& v2xdw)
{
    double sxup = 0.0;
    double sxdw = 0.0;
    if (rhoup + rhodw <= small)
    {
        sx = 0.0;
        v1xup = v1xdw = v2xup = v2xdw = 0.0;
        return;
    }
    if (rhoup > small && sqrt(fabs(grhoup2)) > small)
    {
        xc_math::pbex(2.0 * rhoup, 4.0 * grhoup2, 0, sxup, v1xup, v2xup);
    }
    if (rhodw > small && sqrt(fabs(grhodw2)) > small)
    {
        xc_math::pbex(2.0 * rhodw, 4.0 * grhodw2, 0, sxdw, v1xdw, v2xdw);
    }
    sx = 0.5 * (sxup + sxdw);
    v2xup = 2.0 * v2xup;
    v2xdw = 2.0 * v2xdw;
}

__device__ void gcc_pbe(double rho, double zeta, double grho, double& sc,
                        double& v1cup, double& v1cdw, double& v2c)
{
    sc = 0.0;
    v1cup = 0.0;
    v1cdw = 0.0;
    v2c = 0.0;
    if (fabs(zeta) - 1.0 > small || rho <= small || sqrt(fabs(grho)) <= small)
    {
        return;
    }
    const double clipped = fmin(fabs(zeta), 1.0 - epsr);
    zeta = zeta > 0.0 ? clipped : -clipped;
    xc_math::pbec_spin(rho, zeta, grho, 1, sc, v1cup, v1cdw, v2c);
}

__global__ void ncgga_local_kernel(int nrxx,
                                   const NcggaSpinMapPoint* spin_map,
                                   const double* rhotmp1,
                                   const double* rhotmp2,
                                   const double* gdr1,
                                   const double* gdr2,
                                   double e2,
                                   bool is_gga,
                                   double* v,
                                   double* h1,
                                   double* h2,
                                   double* partial)
{
    const int ir = blockIdx.x * blockDim.x + threadIdx.x;
    double energy = 0.0;
    if (ir < nrxx)
    {
        const NcggaSpinMapPoint point = spin_map[ir];
        for (int channel = 0; channel < 4; ++channel)
        {
            v[channel * nrxx + ir] = 0.0;
        }
        for (int component = 0; component < 3; ++component)
        {
            h1[3 * ir + component] = 0.0;
            h2[3 * ir + component] = 0.0;
        }
        const double arho = point.absolute_density;
        if (arho > vanishing)
        {
            double zeta = point.clipped_magnitude / arho;
            zeta = zeta > 1.0 ? 1.0 : (zeta < -1.0 ? -1.0 : zeta);
            double exc = 0.0;
            double vup = 0.0;
            double vdw = 0.0;
            local_pbe(arho, zeta, exc, vup, vdw);
            for (int channel = 0; channel < 4; ++channel)
            {
                v[channel * nrxx + ir] = e2 * (jacobian_value(point, 0, channel) * vup
                                               + jacobian_value(point, 1, channel) * vdw);
            }
            energy += e2 * exc * arho;
        }
        if (is_gga)
        {
            const double* g1 = gdr1 + 3 * ir;
            const double* g2 = gdr2 + 3 * ir;
            const double grho2a = dot3(g1, g1);
            const double grho2b = dot3(g2, g2);
            const double rh = rhotmp1[ir] + rhotmp2[ir];
            double sx = 0.0;
            double v1xup = 0.0;
            double v1xdw = 0.0;
            double v2xup = 0.0;
            double v2xdw = 0.0;
            gcx_pbe(rhotmp1[ir], rhotmp2[ir], grho2a, grho2b, sx, v1xup, v1xdw, v2xup, v2xdw);
            double sc = 0.0;
            double v1cup = 0.0;
            double v1cdw = 0.0;
            double v2c = 0.0;
            if (rh > epsr)
            {
                const double zeta_input = fabs((rhotmp1[ir] - rhotmp2[ir]) / rh);
                double zeta = zeta_input;
                const double sum[3] = {g1[0] + g2[0], g1[1] + g2[1], g1[2] + g2[2]};
                const double grh2 = dot3(sum, sum);
                gcc_pbe(rh, zeta, grh2, sc, v1cup, v1cdw, v2c);
                if (zeta_input > 1.0 - epsr)
                {
                    const double fixed = 0.5 * ((1.0 + zeta) * v1cup + (1.0 - zeta) * v1cdw);
                    v1cup = fixed;
                    v1cdw = fixed;
                }
            }
            const double vup = e2 * (v1xup + v1cup);
            const double vdw = e2 * (v1xdw + v1cdw);
            for (int component = 0; component < 3; ++component)
            {
                h1[3 * ir + component] = e2 * ((v2xup + v2c) * g1[component] + v2c * g2[component]);
                h2[3 * ir + component] = e2 * ((v2xdw + v2c) * g2[component] + v2c * g1[component]);
            }
            for (int channel = 0; channel < 4; ++channel)
            {
                v[channel * nrxx + ir] += jacobian_value(point, 0, channel) * vup
                                          + jacobian_value(point, 1, channel) * vdw;
            }
            energy += e2 * (sx + sc);
        }
    }
    __shared__ double shared[threads];
    shared[threadIdx.x] = energy;
    __syncthreads();
    for (int stride = threads / 2; stride > 0; stride /= 2)
    {
        if (threadIdx.x < stride)
        {
            shared[threadIdx.x] += shared[threadIdx.x + stride];
        }
        __syncthreads();
    }
    if (threadIdx.x == 0)
    {
        partial[blockIdx.x] = shared[0];
    }
}

bool allocate(void** pointer, size_t bytes)
{
    const auto status = cudaMalloc(pointer, bytes);
    if (status == cudaErrorMemoryAllocation)
    {
        cudaGetLastError();
        *pointer = nullptr;
        return false;
    }
    CHECK_CUDA(status);
    return true;
}
} // namespace

bool ncgga_sf_local_gpu(int nrxx,
                        const NcggaSpinMapPoint* spin_map,
                        const double* spin_density_up,
                        const double* spin_density_down,
                        const ModuleBase::Vector3<double>* gradient_up,
                        const ModuleBase::Vector3<double>* gradient_down,
                        double e2,
                        bool is_gga,
                        ModuleBase::matrix& potential,
                        std::vector<ModuleBase::Vector3<double>>& flux_up,
                        std::vector<ModuleBase::Vector3<double>>& flux_down,
                        double& etxc)
{
    if (nrxx <= 0)
    {
        return false;
    }
    const int blocks = (nrxx + threads - 1) / threads;
    const size_t points = static_cast<size_t>(nrxx);
    const size_t vectors = 3 * points;
    void* d_map = nullptr;
    void* d_rho1 = nullptr;
    void* d_rho2 = nullptr;
    void* d_g1 = nullptr;
    void* d_g2 = nullptr;
    void* d_v = nullptr;
    void* d_h1 = nullptr;
    void* d_h2 = nullptr;
    void* d_partial = nullptr;
    const size_t map_bytes = points * sizeof(NcggaSpinMapPoint);
    const size_t vec_bytes = vectors * sizeof(double);
    const size_t v_bytes = 4 * points * sizeof(double);
    const size_t partial_bytes = blocks * sizeof(double);
    if (!allocate(&d_map, map_bytes) || !allocate(&d_rho1, points * sizeof(double))
        || !allocate(&d_rho2, points * sizeof(double)) || !allocate(&d_g1, vec_bytes)
        || !allocate(&d_g2, vec_bytes) || !allocate(&d_v, v_bytes)
        || !allocate(&d_h1, vec_bytes) || !allocate(&d_h2, vec_bytes)
        || !allocate(&d_partial, partial_bytes))
    {
        cudaFree(d_map); cudaFree(d_rho1); cudaFree(d_rho2); cudaFree(d_g1);
        cudaFree(d_g2); cudaFree(d_v); cudaFree(d_h1); cudaFree(d_h2); cudaFree(d_partial);
        return false;
    }
    CHECK_CUDA(cudaMemcpy(d_map, spin_map, map_bytes, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_rho1, spin_density_up, points * sizeof(double), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_rho2, spin_density_down, points * sizeof(double), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_g1, gradient_up, vec_bytes, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_g2, gradient_down, vec_bytes, cudaMemcpyHostToDevice));
    ncgga_local_kernel<<<blocks, threads>>>(
        nrxx, static_cast<const NcggaSpinMapPoint*>(d_map),
        static_cast<const double*>(d_rho1), static_cast<const double*>(d_rho2),
        static_cast<const double*>(d_g1), static_cast<const double*>(d_g2), e2, is_gga,
        static_cast<double*>(d_v), static_cast<double*>(d_h1),
        static_cast<double*>(d_h2), static_cast<double*>(d_partial));
    CHECK_CUDA(cudaGetLastError());
    CHECK_CUDA(cudaDeviceSynchronize());
    CHECK_CUDA(cudaMemcpy(potential.c, d_v, v_bytes, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(flux_up.data(), d_h1, vec_bytes, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(flux_down.data(), d_h2, vec_bytes, cudaMemcpyDeviceToHost));
    std::vector<double> partial(blocks);
    CHECK_CUDA(cudaMemcpy(partial.data(), d_partial, partial_bytes, cudaMemcpyDeviceToHost));
    double sum = 0.0;
    for (double value : partial)
    {
        sum += value;
    }
    etxc += sum;
    cudaFree(d_map); cudaFree(d_rho1); cudaFree(d_rho2); cudaFree(d_g1);
    cudaFree(d_g2); cudaFree(d_v); cudaFree(d_h1); cudaFree(d_h2); cudaFree(d_partial);
    return true;
}

} // namespace ModuleXC
