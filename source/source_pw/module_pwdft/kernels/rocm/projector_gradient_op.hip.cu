#include "source_pw/module_pwdft/kernels/projector_gradient_op.h"
#include <hip/hip_runtime.h>
#include <thrust/complex.h>
#include <stdexcept>
namespace hamilt {
template <typename Real>
__global__ void gradient_kernel(int npw, int ld, int nkb, int nbeta, int nq, Real dq, Real tpiba,
    const int* metadata, const Real* q, const Real* tab, const Real* derivative,
    const thrust::complex<Real>* sk, thrust::complex<Real>* out)
{
    const int ig=blockIdx.x*blockDim.x+threadIdx.x;
    if(ig<npw) gradient_element(ig,static_cast<int>(blockIdx.y),npw,ld,nkb,nbeta,nq,dq,tpiba,metadata,q,tab,derivative,sk,out);
}
template <typename Real>
void projector_gradient_op<Real, base_device::DEVICE_GPU>::operator()(int npw, int ld, int nkb, int nbeta, int nq,
    Real dq, Real tpiba, const int* metadata, const Real* q, const Real* tab, const Real* derivative,
    const std::complex<Real>* sk, std::complex<Real>* out) const
{
    if(npw==0 || nkb==0) return;
    gradient_kernel<<<dim3((npw+255)/256,nkb),256>>>(npw,ld,nkb,nbeta,nq,dq,tpiba,metadata,q,tab,derivative,
        reinterpret_cast<const thrust::complex<Real>*>(sk),reinterpret_cast<thrust::complex<Real>*>(out));
    const hipError_t error=hipGetLastError();
    if(error!=hipSuccess) throw std::runtime_error(hipGetErrorString(error));
}
template class projector_gradient_op<float, base_device::DEVICE_GPU>;
template class projector_gradient_op<double, base_device::DEVICE_GPU>;
}
