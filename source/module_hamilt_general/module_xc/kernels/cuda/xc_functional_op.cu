#include <module_hamilt_general/module_xc/kernels/xc_functional_op.h>
#include "module_base/module_device/types.h"
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
