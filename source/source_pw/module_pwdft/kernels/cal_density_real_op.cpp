#include "source_pw/module_pwdft/kernels/cal_density_real_op.h"
#include "source_psi/psi.h"
namespace hamilt
{
template <typename T>
struct cal_density_real_op<T, base_device::DEVICE_CPU>
{
    void operator()(const T *in1, const T *in2, T *out, double omega, int nrxx)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ir = 0; ir < nrxx; ir++)
        {
            // assert(is_finite(psi_nk_real[ir]));
            // assert(is_finite(psi_mq_real[ir]));
            out[ir] = in1[ir] * std::conj(in2[ir]) / static_cast<T>(omega); // Phase e^(i(q-k)r)
        }
    }

};

template struct cal_density_real_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct cal_density_real_op<std::complex<double>, base_device::DEVICE_CPU>;
}