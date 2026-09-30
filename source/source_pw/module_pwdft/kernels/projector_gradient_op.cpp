#include "source_pw/module_pwdft/kernels/projector_gradient_op.h"
namespace hamilt {
template <typename Real, typename Device>
void projector_gradient_op<Real, Device>::operator()(int npw, int ld, int nkb, int nbeta, int nq,
    Real dq, Real tpiba, const int* metadata, const Real* q, const Real* tab, const Real* derivative,
    const std::complex<Real>* sk, std::complex<Real>* out) const
{
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int p=0; p<nkb; ++p)
        for (int ig=0; ig<npw; ++ig)
            gradient_element(ig,p,npw,ld,nkb,nbeta,nq,dq,tpiba,metadata,q,tab,derivative,sk,out);
}
template class projector_gradient_op<float, base_device::DEVICE_CPU>;
template class projector_gradient_op<double, base_device::DEVICE_CPU>;
}
