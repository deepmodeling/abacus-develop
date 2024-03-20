#include <module_hamilt_general/module_xc/kernels/xc_functional_op.h>

namespace hamilt {

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
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for(int ig = 0; ig < npw; ig++) {
		// the formula is : rho(r)^prime = \int iG * rho(G)e^{iGr} dG
		// double kplusg = wfc_basis->getgpluskcar(ik,ig)[ipol] * tpiba;
        Real kplusg = (gcar[(ik * npwx + ig) * 3 + pol] +
                       kvec_c[ik * 3 + pol]) * tpiba;
                       
		// calculate the charge density gradient in reciprocal space.
		porter[ig] = T(0.0, kplusg) * rhog[ig];
	}
}

template <typename T, typename Device>
void xc_functional_grad_wfc_op<T, Device>::operator()(
    const int& ipol,
    const int& nrxx,
    const T * porter,
    T* grad)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for (int ir = 0; ir < nrxx; ++ir) {
        grad[ipol * nrxx + ir] = porter[ir];
	}
}

template struct xc_functional_grad_wfc_op<std::complex<float> , psi::DEVICE_CPU>;
template struct xc_functional_grad_wfc_op<std::complex<double>, psi::DEVICE_CPU>;

} // namespace hamilt   