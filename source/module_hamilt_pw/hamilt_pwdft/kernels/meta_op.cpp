#include "module_hamilt_pw/hamilt_pwdft/kernels/meta_op.h"

namespace hamilt{

template <typename FPTYPE>
struct meta_pw_op<FPTYPE, psi::DEVICE_CPU> {
    void operator() (
            const psi::DEVICE_CPU* /*ctx*/,
            const int& ik,
            const int& pol,
            const int& npw,
            const int& npwx,
            const FPTYPE& tpiba,
            const FPTYPE* gcar,
            const FPTYPE* kvec_c,
            const std::complex<FPTYPE>* in,
            std::complex<FPTYPE>* out,
            const bool add)
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int ig = 0; ig < npw; ig++) {
            FPTYPE fact = (gcar[(ik * npwx + ig) * 3 + pol] +
                           kvec_c[ik * 3 + pol]) * tpiba;
            if (add) {
                out[ig] -= in[ig] * std::complex<FPTYPE>(0.0, fact);
            }
            else {
                out[ig] = in[ig] * std::complex<FPTYPE>(0.0, fact);
            }
        }
    }
};

template struct meta_pw_op<float, psi::DEVICE_CPU>;
template struct meta_pw_op<double, psi::DEVICE_CPU>;

}  // namespace hamilt

