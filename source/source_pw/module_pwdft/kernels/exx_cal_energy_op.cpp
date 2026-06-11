#include "source_pw/module_pwdft/kernels/exx_cal_energy_op.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_psi/psi.h"

namespace hamilt
{

// #ifdef _OPENMP
// #pragma omp parallel for reduction(+:Eexx_ik_real)
// #endif
// for (int ig = 0; ig < rhopw_dev->npw; ig++)
// {
//     int nks = wfcpw->nks;
//     int npw = rhopw_dev->npw;
//     int nk = nks / nk_fac;
//     Real Fac = pot[ik * nks * npw + iq * npw + ig];

// Eexx_ik_real += Fac * (density_recip[ig] * std::conj(density_recip[ig])).real()
//                 * wg_iqb_real / nqs * wg_ikb_real / kv->wk[ik];
// }

template <typename FPTYPE>
struct exx_cal_energy_op<std::complex<FPTYPE>, base_device::DEVICE_CPU>
{
    using T = std::complex<FPTYPE>;
    FPTYPE operator()(const T *den, const FPTYPE *pot, FPTYPE scalar, int npw)
    {
        FPTYPE energy = 0.0;
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+:energy)
        #endif
        for (int ig = 0; ig < npw; ++ig)
        {
            // Calculate the energy contribution from each reciprocal lattice vector
            energy += (den[ig] * std::conj(den[ig])).real() * pot[ig];
        }
        // Scale the energy by the scalar factor
        return scalar * energy;
    }
};

template <typename FPTYPE>
struct exx_density_potential_mul_op<std::complex<FPTYPE>, base_device::DEVICE_CPU>
{
    using T = std::complex<FPTYPE>;
    FPTYPE operator()(const T *vector_in, FPTYPE *vector_buffer, const FPTYPE *pot, FPTYPE *vec_temp, FPTYPE *weights,
        int npw, int batch_idx)
    {
        const FPTYPE one{1}; // Scalar one for gemv
        const FPTYPE zero{0}; // Scalar one for gemv
        const int inc{1}; // Increment for gemv
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (int i = 0; i < npw * batch_idx; ++i)
        {
            // Calculate the energy contribution from each reciprocal lattice vector
            vector_buffer[i] = (vector_in[i] * std::conj(vector_in[i])).real();
        }
        ModuleBase::gemv_op<FPTYPE, base_device::DEVICE_CPU>()(
            'T',
            npw,          // m
            batch_idx,           // n
            &one,            // alpha
            vector_buffer,   // matrix (batch_idx, npw), memory layout as
            npw,         // lda
            pot,  // vector (npw, )
            inc,                   // incx
            &zero,            // beta
            vec_temp,          // batch_idx vector (batch_idx)
            inc);                  // incy
        // Obtain the energy
        double energy = 0.0;

        return ModuleBase::dot_real_op<FPTYPE, base_device::DEVICE_CPU>()(
            batch_idx,
            vec_temp,
            weights,
            false
        );
    }
};


template struct exx_cal_energy_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct exx_cal_energy_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct exx_density_potential_mul_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct exx_density_potential_mul_op<std::complex<double>, base_device::DEVICE_CPU>;

template <typename FPTYPE>
struct exx_cal_energy_accumulate_op<std::complex<FPTYPE>, base_device::DEVICE_CPU>
{
    using T = std::complex<FPTYPE>;
    void operator()(const T* den, const FPTYPE* pot, FPTYPE scalar, int npw, FPTYPE* result)
    {
        *result += exx_cal_energy_op<T, base_device::DEVICE_CPU>()(den, pot, scalar, npw);
    }
};

template struct exx_cal_energy_accumulate_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct exx_cal_energy_accumulate_op<std::complex<double>, base_device::DEVICE_CPU>;
} // namespace hamilt
