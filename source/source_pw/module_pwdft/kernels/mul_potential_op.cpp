#include "source_pw/module_pwdft/kernels/mul_potential_op.h"
#include "source_io/module_parameter/parameter.h"
#include "source_psi/psi.h"
#include "source_base/macros.h"
namespace hamilt {
template <typename FPTYPE>
struct mul_potential_op<std::complex<FPTYPE>, base_device::DEVICE_CPU>
{
    using T = std::complex<FPTYPE>;

    // Single operator
    void operator()(const FPTYPE *pot, T *density_recip, int npw, int nks, int ik, int iq)
    {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (int ig = 0; ig < npw; ig++)
        {
            // int ig_kq = ik * nks * npw + iq * npw + ig;
            density_recip[ig] *= pot[ig];

        }
    }

    // Batched operator - CPU fallback: process sequentially
    void operator_batch(
        const FPTYPE *pot,
        T *density_recip_batch,
        int npw,
        int batch_size)
    {
        // CPU: Process each batch element sequentially
        for (int ib = 0; ib < batch_size; ib++)
        {
            T* density_recip_ib = density_recip_batch + ib * npw;
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static)
            #endif
            for (int ig = 0; ig < npw; ig++)
            {
                density_recip_ib[ig] *= pot[ig];
            }
        }
    }
};

template struct mul_potential_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct mul_potential_op<std::complex<double>, base_device::DEVICE_CPU>;
} // hamilt