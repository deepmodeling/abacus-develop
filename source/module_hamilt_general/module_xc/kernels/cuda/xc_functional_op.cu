#include <module_hamilt_general/module_xc/kernels/xc_functional_op.h>
#include "module_base/module_device/types.h"
#include "module_hamilt_general/module_xc/xc_math.h"
#include <thrust/complex.h>
#include <base/macros/macros.h>

#define THREADS_PER_BLOCK 256

namespace hamilt {

namespace {

__global__ void density_gradient(int npw, int component, const int* index,
                                const double* gcar, const thrust::complex<double>* rhog,
                                thrust::complex<double>* grid)
{
    const int ig = blockIdx.x * blockDim.x + threadIdx.x;
    if (ig < npw)
    {
        grid[index[ig]] = thrust::complex<double>(0.0, 1.0) * rhog[ig] * gcar[3 * ig + component];
    }
}

__global__ void extract_density_real(int size, int stride, int offset, double scale,
                                    const thrust::complex<double>* grid, double* result)
{
    const int ir = blockIdx.x * blockDim.x + threadIdx.x;
    if (ir < size)
    {
        result[stride * ir + offset] = grid[ir].real() * scale;
    }
}

__global__ void pack_density_component(int size, int component, const double* vectors,
                                      thrust::complex<double>* grid)
{
    const int ir = blockIdx.x * blockDim.x + threadIdx.x;
    if (ir < size)
    {
        grid[ir] = thrust::complex<double>(vectors[3 * ir + component], 0.0);
    }
}

__global__ void accumulate_density_divergence(int npw, int nxyz, int component,
                                            const int* index, const double* gcar,
                                            const thrust::complex<double>* grid,
                                            thrust::complex<double>* reciprocal)
{
    const int ig = blockIdx.x * blockDim.x + threadIdx.x;
    if (ig < npw)
    {
        const auto coefficient = grid[index[ig]] / double(nxyz);
        const auto term = thrust::complex<double>(0.0, 1.0) * coefficient * gcar[3 * ig + component];
        if (component == 0)
        {
            reciprocal[ig] = term;
        }
        else
        {
            reciprocal[ig] += term;
        }
    }
}

__device__ double magnetic_norm(int size, int ir, const double* rho)
{
    const double x = rho[size + ir];
    const double y = rho[2 * size + ir];
    const double z = rho[3 * size + ir];
    return sqrt(x * x + y * y + z * z);
}

__device__ double block_sum(double value, double* scratch)
{
    const int lane = threadIdx.x;
    scratch[lane] = value;
    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride /= 2)
    {
        if (lane < stride)
        {
            scratch[lane] += scratch[lane + stride];
        }
        __syncthreads();
    }
    const double result = scratch[0];
    __syncthreads();
    return result;
}

__global__ void pack_density(int size, const double* density, thrust::complex<double>* grid)
{
    const int ir = blockIdx.x * blockDim.x + threadIdx.x;
    if (ir < size)
    {
        grid[ir] = thrust::complex<double>(density[ir], 0.0);
    }
}

__global__ void gather_density(int npw, int nxyz, const int* index,
                               const thrust::complex<double>* grid, thrust::complex<double>* reciprocal)
{
    const int ig = blockIdx.x * blockDim.x + threadIdx.x;
    if (ig < npw)
    {
        reciprocal[ig] = grid[index[ig]] / double(nxyz);
    }
}

__global__ void add_reciprocal(int npw, const thrust::complex<double>* source,
                              thrust::complex<double>* destination)
{
    const int ig = blockIdx.x * blockDim.x + threadIdx.x;
    if (ig < npw)
    {
        destination[ig] += source[ig];
    }
}

__global__ void spin_gradient(int size, int spin, int component, double scale,
                             const thrust::complex<double>* grid, const double* rho, double* gradient)
{
    const int ir = blockIdx.x * blockDim.x + threadIdx.x;
    if (ir >= size) { return; }
    const double value = 0.5 * scale * grid[ir].real();
    if (spin == 0)
    {
        gradient[3 * ir + component] = value;
        gradient[3 * (size + ir) + component] = value;
    }
    else
    {
        const double norm = magnetic_norm(size, ir, rho);
        const double direction = norm > 1e-12 ? rho[spin * size + ir] / norm : 0.0;
        gradient[3 * ir + component] += value * direction;
        gradient[3 * (size + ir) + component] -= value * direction;
    }
}

__device__ void local_pbe_spin(int size, int ir, const double* rho, const double* core, double norm,
                              double* potential, double& energy, double& potential_energy)
{
    for (int spin = 0; spin < 4; ++spin) { potential[spin * size + ir] = 0.0; }
    const double density = rho[ir] + core[ir];
    const double absolute_density = fabs(density);
    if (absolute_density <= 1e-10) { return; }
    const double zeta = fmin(norm / absolute_density, 1.0);
    const double rs = 0.62035049089940 / pow(absolute_density, 1.0 / 3.0);
    double exchange, correlation, exchange_up, exchange_down, correlation_up, correlation_down;
    xc_math::slater_spin(absolute_density, zeta, exchange, exchange_up, exchange_down);
    xc_math::pw_spin(rs, zeta, correlation, correlation_up, correlation_down);
    const double up = exchange_up + correlation_up;
    const double down = exchange_down + correlation_down;
    energy = 2.0 * (exchange + correlation) * density;
    potential[ir] = 2.0 * (0.5 * (up + down));
    potential_energy = potential[ir] * rho[ir];
    const double spin_potential = 0.5 * (up - down);
    if (norm > 1e-10)
    {
        for (int spin = 1; spin < 4; ++spin)
        {
            potential[spin * size + ir] = 2.0 * spin_potential * rho[spin * size + ir] / norm;
            potential_energy += potential[spin * size + ir] * rho[spin * size + ir];
        }
    }
}

__global__ void pbe_spin(int size, const double* rho, const double* core, double* gradient_flux,
                        double* potential, double* sums, bool include_local)
{
    const int ir = blockIdx.x * blockDim.x + threadIdx.x;
    double energy = 0.0, potential_energy = 0.0;
    if (ir < size)
    {
        const double norm = magnetic_norm(size, ir, rho);
        if (include_local) { local_pbe_spin(size, ir, rho, core, norm, potential, energy, potential_energy); }
        const double up = 0.5 * (rho[ir] + norm) + 0.5 * core[ir];
        const double down = 0.5 * (rho[ir] - norm) + 0.5 * core[ir];
        const double total = up + down;
        double gu[3], gd[3];
        double gu2 = 0.0, gd2 = 0.0, total_gradient2 = 0.0;
        for (int component = 0; component < 3; ++component)
        {
            gu[component] = gradient_flux[3 * ir + component];
            gd[component] = gradient_flux[3 * (size + ir) + component];
            gu2 += gu[component] * gu[component];
            gd2 += gd[component] * gd[component];
            const double sum = gu[component] + gd[component];
            total_gradient2 += sum * sum;
        }
        double sxup = 0.0, sxdw = 0.0, v1xup = 0.0, v1xdw = 0.0, v2xup = 0.0, v2xdw = 0.0;
        if (total > 1e-10)
        {
            if (up > 1e-10 && sqrt(fabs(gu2)) > 1e-10)
            {
                xc_math::pbex(2.0 * up, 4.0 * gu2, 0, sxup, v1xup, v2xup);
            }
            if (down > 1e-10 && sqrt(fabs(gd2)) > 1e-10)
            {
                xc_math::pbex(2.0 * down, 4.0 * gd2, 0, sxdw, v1xdw, v2xdw);
            }
        }
        double sc = 0.0, v1cup = 0.0, v1cdw = 0.0, v2c = 0.0;
        if (total > 1e-6)
        {
            double zeta = fabs((up - down) / total);
            if (zeta - 1.0 <= 1e-10 && sqrt(fabs(total_gradient2)) > 1e-10)
            {
                zeta = fmin(zeta, 1.0 - 1e-6);
                xc_math::pbec_spin(total, zeta, total_gradient2, 1, sc, v1cup, v1cdw, v2c);
            }
        }
        const double vu = 2.0 * (v1xup + v1cup);
        const double vd = 2.0 * (v1xdw + v1cdw);
        potential[ir] += 0.5 * (vu + vd);
        for (int spin = 1; spin < 4; ++spin)
        {
            const double direction = norm > 1e-12 ? rho[spin * size + ir] / norm : 0.0;
            potential[spin * size + ir] += 0.5 * (vu - vd) * direction;
        }
        // The same storage holds gradients before evaluation and fluxes afterwards.
        for (int component = 0; component < 3; ++component)
        {
            gradient_flux[3 * ir + component] = 2.0 * ((2.0 * v2xup + v2c) * gu[component] + v2c * gd[component]);
            gradient_flux[3 * (size + ir) + component]
                = 2.0 * ((2.0 * v2xdw + v2c) * gd[component] + v2c * gu[component]);
        }
        energy += 2.0 * (0.5 * (sxup + sxdw) + sc);
        potential_energy += vu * (up - 0.5 * core[ir]) + vd * (down - 0.5 * core[ir]);
    }
    __shared__ double scratch[THREADS_PER_BLOCK];
    const double energy_sum = block_sum(energy, scratch);
    const double potential_sum = block_sum(potential_energy, scratch);
    if (threadIdx.x == 0)
    {
        sums[blockIdx.x] = energy_sum;
        sums[gridDim.x + blockIdx.x] = potential_sum;
    }
}

__global__ void pack_spin_flux(int size, int spin, int component, const double* rho,
                              const double* flux, thrust::complex<double>* grid)
{
    const int ir = blockIdx.x * blockDim.x + threadIdx.x;
    if (ir >= size) { return; }
    const double up = flux[3 * ir + component];
    const double down = flux[3 * (size + ir) + component];
    double value = 0.5 * (up + down);
    if (spin != 0)
    {
        const double norm = magnetic_norm(size, ir, rho);
        const double direction = norm > 1e-12 ? rho[spin * size + ir] / norm : 0.0;
        value = 0.5 * (up - down) * direction;
    }
    grid[ir] = thrust::complex<double>(value, 0.0);
}

__global__ void apply_spin_divergence(int size, int spin, double scale,
                                    const thrust::complex<double>* grid, const double* rho,
                                    double* potential, double* magnetic_divergence, double* sums)
{
    const int ir = blockIdx.x * blockDim.x + threadIdx.x;
    double potential_energy = 0.0;
    if (ir < size)
    {
        const double value = scale * grid[ir].real();
        if (spin == 0)
        {
            potential[ir] -= value;
            potential_energy = value * rho[ir];
            magnetic_divergence[ir] = 0.0;
        }
        else
        {
            const double norm = magnetic_norm(size, ir, rho);
            const double direction = norm > 1e-12 ? rho[spin * size + ir] / norm : 0.0;
            magnetic_divergence[ir] += value * direction;
        }
    }
    if (spin == 0)
    {
        __shared__ double scratch[THREADS_PER_BLOCK];
        const double sum = block_sum(potential_energy, scratch);
        if (threadIdx.x == 0) { sums[2 * gridDim.x + blockIdx.x] = sum; }
    }
}

__global__ void finish_spin_potential(int size, const double* rho, const double* magnetic_divergence,
                                    double* potential, double* sums)
{
    const int ir = blockIdx.x * blockDim.x + threadIdx.x;
    double potential_energy = 0.0;
    if (ir < size)
    {
        const double norm = magnetic_norm(size, ir, rho);
        const double value = magnetic_divergence[ir];
        for (int spin = 1; spin < 4; ++spin)
        {
            const double direction = norm > 1e-12 ? rho[spin * size + ir] / norm : 0.0;
            potential[spin * size + ir] -= value * direction;
        }
        potential_energy = value * norm;
    }
    __shared__ double scratch[THREADS_PER_BLOCK];
    const double sum = block_sum(potential_energy, scratch);
    if (threadIdx.x == 0) { sums[3 * gridDim.x + blockIdx.x] = sum; }
}

} // namespace

void xc_gradient_coefficients_gpu(int npw, int component, const int* grid_index,
                                 const double* gcar, const std::complex<double>* rhog,
                                 std::complex<double>* grid)
{
    density_gradient<<<(npw + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK>>>(
        npw, component, grid_index, gcar, reinterpret_cast<const thrust::complex<double>*>(rhog),
        reinterpret_cast<thrust::complex<double>*>(grid));
    cudaCheckOnDebug();
}

void xc_extract_real_gpu(int size, int stride, int offset, double scale,
                         const std::complex<double>* grid, double* result)
{
    extract_density_real<<<(size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK>>>(
        size, stride, offset, scale, reinterpret_cast<const thrust::complex<double>*>(grid), result);
    cudaCheckOnDebug();
}

void xc_pack_component_gpu(int size, int component, const double* vectors, std::complex<double>* grid)
{
    pack_density_component<<<(size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK>>>(
        size, component, vectors, reinterpret_cast<thrust::complex<double>*>(grid));
    cudaCheckOnDebug();
}

void xc_accumulate_divergence_gpu(int npw, int nxyz, int component,
                                const int* grid_index, const double* gcar,
                                const std::complex<double>* grid, std::complex<double>* reciprocal)
{
    accumulate_density_divergence<<<(npw + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK>>>(
        npw, nxyz, component, grid_index, gcar, reinterpret_cast<const thrust::complex<double>*>(grid),
        reinterpret_cast<thrust::complex<double>*>(reciprocal));
    cudaCheckOnDebug();
}

int xc_reduction_size_gpu(int size)
{
    return (size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
}

void xc_pack_density_gpu(int size, const double* density, std::complex<double>* grid)
{
    pack_density<<<xc_reduction_size_gpu(size), THREADS_PER_BLOCK>>>(
        size, density, reinterpret_cast<thrust::complex<double>*>(grid));
    cudaCheckOnDebug();
}

void xc_gather_density_gpu(int npw, int nxyz, const int* index, const std::complex<double>* grid, std::complex<double>* reciprocal)
{
    gather_density<<<xc_reduction_size_gpu(npw), THREADS_PER_BLOCK>>>(
        npw, nxyz, index, reinterpret_cast<const thrust::complex<double>*>(grid), reinterpret_cast<thrust::complex<double>*>(reciprocal));
    cudaCheckOnDebug();
}

void xc_add_reciprocal_gpu(int npw, const std::complex<double>* source, std::complex<double>* destination)
{
    add_reciprocal<<<xc_reduction_size_gpu(npw), THREADS_PER_BLOCK>>>(
        npw, reinterpret_cast<const thrust::complex<double>*>(source), reinterpret_cast<thrust::complex<double>*>(destination));
    cudaCheckOnDebug();
}

void xc_spin_gradient_gpu(int size, int spin, int component, double scale, const std::complex<double>* grid, const double* rho, double* gradient)
{
    spin_gradient<<<xc_reduction_size_gpu(size), THREADS_PER_BLOCK>>>(
        size, spin, component, scale, reinterpret_cast<const thrust::complex<double>*>(grid), rho, gradient);
    cudaCheckOnDebug();
}

void xc_pbe_spin_gpu(int size, const double* rho, const double* core, double* gradient_flux,
                     double* potential, double* sums, bool include_local)
{
    pbe_spin<<<xc_reduction_size_gpu(size), THREADS_PER_BLOCK>>>(
        size, rho, core, gradient_flux, potential, sums, include_local);
    cudaCheckOnDebug();
}

void xc_pack_spin_flux_gpu(int size, int spin, int component, const double* rho, const double* flux, std::complex<double>* grid)
{
    pack_spin_flux<<<xc_reduction_size_gpu(size), THREADS_PER_BLOCK>>>(
        size, spin, component, rho, flux, reinterpret_cast<thrust::complex<double>*>(grid));
    cudaCheckOnDebug();
}

void xc_apply_spin_divergence_gpu(int size, int spin, double scale, const std::complex<double>* grid, const double* rho, double* potential, double* magnetic_divergence, double* sums)
{
    apply_spin_divergence<<<xc_reduction_size_gpu(size), THREADS_PER_BLOCK>>>(
        size, spin, scale, reinterpret_cast<const thrust::complex<double>*>(grid), rho, potential, magnetic_divergence, sums);
    cudaCheckOnDebug();
}

void xc_finish_spin_potential_gpu(int size, const double* rho, const double* magnetic_divergence, double* potential, double* sums)
{
    finish_spin_potential<<<xc_reduction_size_gpu(size), THREADS_PER_BLOCK>>>(
        size, rho, magnetic_divergence, potential, sums);
    cudaCheckOnDebug();
}

template <typename T>
__global__ void xc_functional_grad_wfc(
    const int ik,
    const int pol,
    const int npw,
    const int npwx,
	const T tpiba,
    const T* gcar,
    const T* kvec_c,
    const thrust::complex<T>* rhog,
    thrust::complex<T>* porter)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= npw) { return; }
	// the formula is : rho(r)^prime = \int iG * rho(G)e^{iGr} dG
	// double kplusg = wfc_basis->getgpluskcar(ik,ig)[ipol] * tpiba;
    T kplusg = (gcar[(ik * npwx + idx) * 3 + pol] +
                   kvec_c[ik * 3 + pol]) * tpiba;
	// calculate the charge density gradient in reciprocal space.
	porter[idx] = thrust::complex<T>(0.0, kplusg) * rhog[idx];
}

template <typename T>
__global__ void xc_functional_grad_wfc(
    const int ipol,
    const int nrxx,
    const T* porter,
    T* grad)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= nrxx) { return; }
    grad[ipol * nrxx + idx] = porter[idx];
}

template <typename T, typename Device>
void xc_functional_grad_wfc_op<T, Device>::operator()(
    const int& ik,
    const int& pol,
    const int& npw,
    const int& npwx,
	const Real& tpiba,
    const Real * gcar,
    const Real * kvec_c,
    const T * rhog,
    T* porter)
{
    auto porter_ = reinterpret_cast<thrust::complex<Real>*>(porter);
    auto rhog_ = reinterpret_cast<const thrust::complex<Real>*>(rhog);
    const int block = (npw + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    xc_functional_grad_wfc<Real><<<block, THREADS_PER_BLOCK>>>(
        ik, pol, npw, npwx, tpiba, gcar, kvec_c, rhog_, porter_);

    cudaCheckOnDebug();
}

template <typename T, typename Device>
void xc_functional_grad_wfc_op<T, Device>::operator()(
    const int& ipol,
    const int& nrxx,
    const T * porter,
    T* grad)
{
    auto grad_ = reinterpret_cast<thrust::complex<Real>*>(grad);
    auto porter_ = reinterpret_cast<const thrust::complex<Real>*>(porter);
    const int block = (nrxx + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    xc_functional_grad_wfc<<<block, THREADS_PER_BLOCK>>>(
        ipol, nrxx, porter_, grad_);

    cudaCheckOnDebug();
}

template struct xc_functional_grad_wfc_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct xc_functional_grad_wfc_op<std::complex<double>, base_device::DEVICE_GPU>;

} // namespace hamilt
