#include "force_sincos_op.h"
#include "module_base/libm/libm.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace hamilt
{

template <typename FPTYPE>
struct cal_force_loc_sincos_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* ctx,
                    const int& nat,
                    const int& npw,
                    const FPTYPE* gcar,
                    const FPTYPE* tau,
                    const int* iat2it,
                    const FPTYPE* vloc_factors,
                    const std::complex<FPTYPE>* aux,
                    const FPTYPE& scale_factor,
                    FPTYPE* force)
    {
        const FPTYPE TWO_PI = 2.0 * M_PI;

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int iat = 0; iat < nat; ++iat)
        {
            const int it = iat2it[iat];
            const FPTYPE tau_x = tau[iat * 3 + 0];
            const FPTYPE tau_y = tau[iat * 3 + 1];
            const FPTYPE tau_z = tau[iat * 3 + 2];

            FPTYPE local_force[3] = {0.0, 0.0, 0.0};

            for (int ig = 0; ig < npw; ig++)
            {
                const FPTYPE phase = TWO_PI * (gcar[ig * 3 + 0] * tau_x + 
                                               gcar[ig * 3 + 1] * tau_y + 
                                               gcar[ig * 3 + 2] * tau_z);
                FPTYPE sinp, cosp;
                ModuleBase::libm::sincos(phase, &sinp, &cosp);
                
                const FPTYPE factor = vloc_factors[ig] * 
                                     (cosp * aux[ig].imag() + sinp * aux[ig].real());
                
                local_force[0] += gcar[ig * 3 + 0] * factor;
                local_force[1] += gcar[ig * 3 + 1] * factor;
                local_force[2] += gcar[ig * 3 + 2] * factor;
            }

            force[iat * 3 + 0] += local_force[0] * scale_factor;
            force[iat * 3 + 1] += local_force[1] * scale_factor;
            force[iat * 3 + 2] += local_force[2] * scale_factor;
        }
    }
};

template <typename FPTYPE>
struct cal_force_ew_sincos_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* ctx,
                    const int& nat,
                    const int& npw,
                    const int& ig_gge0,
                    const FPTYPE* gcar,
                    const FPTYPE* tau,
                    const int* iat2it,
                    const FPTYPE* it_facts,
                    const std::complex<FPTYPE>* aux,
                    FPTYPE* force)
    {
        const FPTYPE TWO_PI = 2.0 * M_PI;

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int iat = 0; iat < nat; ++iat)
        {
            const int it = iat2it[iat];
            const FPTYPE tau_x = tau[iat * 3 + 0];
            const FPTYPE tau_y = tau[iat * 3 + 1];
            const FPTYPE tau_z = tau[iat * 3 + 2];
            const FPTYPE it_fact = it_facts[iat];

            FPTYPE local_force[3] = {0.0, 0.0, 0.0};

            for (int ig = 0; ig < npw; ig++)
            {
                // Skip G=0 term
                if (ig == ig_gge0) continue;

                const FPTYPE arg = TWO_PI * (gcar[ig * 3 + 0] * tau_x + 
                                             gcar[ig * 3 + 1] * tau_y + 
                                             gcar[ig * 3 + 2] * tau_z);
                FPTYPE sinp, cosp;
                ModuleBase::libm::sincos(arg, &sinp, &cosp);
                
                const FPTYPE sumnb = -cosp * aux[ig].imag() + sinp * aux[ig].real();
                
                local_force[0] += gcar[ig * 3 + 0] * sumnb;
                local_force[1] += gcar[ig * 3 + 1] * sumnb;
                local_force[2] += gcar[ig * 3 + 2] * sumnb;
            }

            force[iat * 3 + 0] += local_force[0] * it_fact;
            force[iat * 3 + 1] += local_force[1] * it_fact;
            force[iat * 3 + 2] += local_force[2] * it_fact;
        }
    }
};

// Template instantiations
template struct cal_force_loc_sincos_op<float, base_device::DEVICE_CPU>;
template struct cal_force_loc_sincos_op<double, base_device::DEVICE_CPU>;

template struct cal_force_ew_sincos_op<float, base_device::DEVICE_CPU>;
template struct cal_force_ew_sincos_op<double, base_device::DEVICE_CPU>;

} // namespace hamilt 