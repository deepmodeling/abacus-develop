#ifndef source_hamilt_MODULE_XC_KERNELS_H_
#define source_hamilt_MODULE_XC_KERNELS_H_

#include <complex>
#include <source_psi/psi.h>
#include <source_base/macros.h>

namespace hamilt {

template <typename T, typename Device>
struct xc_functional_grad_wfc_op {
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const int& ik,
        const int& pol,
        const int& npw,
        const int& npwx,
        const Real& tpiba,
        const Real * gcar,
        const Real * kvec_c,
        const T * rhog,
        T* porter);
    
    void operator()(
        const int& ipol,
        const int& nrxx,
        const T * porter,
        T* grad);
};

} // namespace hamilt
#endif // source_hamilt_MODULE_XC_KERNELS_H_