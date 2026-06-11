#include "source_pw/module_pwdft/kernels/vec_mul_vec_complex_op.h"
#include "source_io/module_parameter/parameter.h"
#include "source_psi/psi.h"
namespace hamilt{
template <typename FPTYPE>
struct vec_mul_vec_complex_op<std::complex<FPTYPE>, base_device::DEVICE_CPU>
{
    using T = std::complex<FPTYPE>;
    void operator()(const T *vec1, const T *vec2, T *out, int n)
    {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (int i = 0; i < n; i++)
        {
            out[i] = vec1[i] * vec2[i];
        }
    }

    // Batched operator - CPU fallback (sequential over batches with OpenMP per element)
    void operator_batch(
        const T *vec1_batch,
        const T *vec2_batch,
        T *out_batch,
        int n,
        int batch_size)
    {
        // Process each batch element sequentially with OpenMP per element
        for (int ib = 0; ib < batch_size; ib++)
        {
            const T* vec1_ib = vec1_batch + ib * n;
            const T* vec2_ib = vec2_batch + ib * n;
            T* out_ib = out_batch + ib * n;

            #ifdef _OPENMP
            #pragma omp parallel for schedule(static)
            #endif
            for (int i = 0; i < n; i++)
            {
                out_ib[i] = vec1_ib[i] * vec2_ib[i];
            }
        }
    }

};
template struct vec_mul_vec_complex_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct vec_mul_vec_complex_op<std::complex<double>, base_device::DEVICE_CPU>;
} // hamilt