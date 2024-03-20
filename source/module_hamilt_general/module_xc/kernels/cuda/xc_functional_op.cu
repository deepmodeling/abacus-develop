#include <module_hamilt_general/module_xc/kernels/xc_functional_op.h>
#include <module_psi/kernels/device.h>
#include <thrust/complex.h>
#include <base/macros/macros.h>

#define THREADS_PER_BLOCK 256

namespace hamilt {

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
    
    cudaErrcheck(cudaGetLastError());
    cudaErrcheck(cudaDeviceSynchronize());
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
    
    cudaErrcheck(cudaGetLastError());
    cudaErrcheck(cudaDeviceSynchronize());
}

template struct xc_functional_grad_wfc_op<std::complex<float> , psi::DEVICE_GPU>;
template struct xc_functional_grad_wfc_op<std::complex<double>, psi::DEVICE_GPU>;

} // namespace hamilt   