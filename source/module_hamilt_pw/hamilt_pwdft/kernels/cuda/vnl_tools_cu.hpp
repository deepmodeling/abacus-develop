#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_PW_HAMILT_PWDFT_KERNELS_CUDA_VNL_TOOLS_CU_HPP
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_PW_HAMILT_PWDFT_KERNELS_CUDA_VNL_TOOLS_CU_HPP

#include <thrust/complex.h>
#ifdef __CUDA
#include <cuda_runtime.h>
#endif
#ifdef __ROCM
#include <hip/hip_runtime.h>
#endif

namespace hamilt
{

template <typename FPTYPE>
__device__ FPTYPE _polynomial_interpolation(const FPTYPE* table,
                                            const int& dim1,
                                            const int& dim2,
                                            const int& tab_2,
                                            const int& tab_3,
                                            const FPTYPE& table_interval,
                                            const FPTYPE& x)
{
    const FPTYPE position = x / table_interval;
    const int iq = static_cast<int>(position);

    const FPTYPE x0 = position - static_cast<FPTYPE>(iq);
    const FPTYPE x1 = 1.0 - x0;
    const FPTYPE x2 = 2.0 - x0;
    const FPTYPE x3 = 3.0 - x0;
    const FPTYPE y = table[(dim1 * tab_2 + dim2) * tab_3 + iq + 0] * x1 * x2 * x3 / 6.0
                     + table[(dim1 * tab_2 + dim2) * tab_3 + iq + 0 + 1] * x0 * x2 * x3 / 2.0
                     - table[(dim1 * tab_2 + dim2) * tab_3 + iq + 0 + 2] * x1 * x0 * x3 / 2.0
                     + table[(dim1 * tab_2 + dim2) * tab_3 + iq + 0 + 3] * x1 * x2 * x0 / 6.0;

    return y;
}

template <typename FPTYPE>
__device__ FPTYPE _polynomial_interpolation_nl(const FPTYPE* table,
                                               const int& dim1,
                                               const int& dim2,
                                               const int& tab_2,
                                               const int& tab_3,
                                               const FPTYPE& table_interval,
                                               const FPTYPE& x)
{
    const FPTYPE position = x / table_interval;
    const int iq = static_cast<int>(position);

    const FPTYPE x0 = position - static_cast<FPTYPE>(iq);
    const FPTYPE x1 = 1.0 - x0;
    const FPTYPE x2 = 2.0 - x0;
    const FPTYPE x3 = 3.0 - x0;
    const FPTYPE y = (table[(dim1 * tab_2 + dim2) * tab_3 + iq + 0] * (-x2 * x3 - x1 * x3 - x1 * x2) / 6.0
                      + table[(dim1 * tab_2 + dim2) * tab_3 + iq + 0 + 1] * (+x2 * x3 - x0 * x3 - x0 * x2) / 2.0
                      - table[(dim1 * tab_2 + dim2) * tab_3 + iq + 0 + 2] * (+x1 * x3 - x0 * x3 - x0 * x1) / 2.0
                      + table[(dim1 * tab_2 + dim2) * tab_3 + iq + 0 + 3] * (+x1 * x2 - x0 * x2 - x0 * x1) / 6.0)
                     / table_interval;

    return y;
}

} // namespace hamilt
#endif
