#include "module_base/kernels/math_kernel_op.h"

namespace ModuleBase
{

template <typename FPTYPE>
struct scal_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const int& N,
                    const std::complex<FPTYPE>* alpha,
                    std::complex<FPTYPE>* X,
                    const int& incx)
    {
        BlasConnector::scal(N, *alpha, X, incx);
    }
};

template <typename T>
struct vector_mul_real_op<T, base_device::DEVICE_CPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const int dim, T* result, const T* vector, const Real constant)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(Real))
#endif
        for (int i = 0; i < dim; i++)
        {
            result[i] = vector[i] * constant;
        }
    }
};

template <typename T>
struct vector_mul_vector_op<T, base_device::DEVICE_CPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const int& dim, T* result, const T* vector1, const Real* vector2, const bool& add)
    {
        if (add)
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(Real))
#endif
            for (int i = 0; i < dim; i++)
            {
                result[i] += vector1[i] * vector2[i];
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(Real))
#endif
            for (int i = 0; i < dim; i++)
            {
                result[i] = vector1[i] * vector2[i];
            }
        }
    }
};

template <typename T>
struct vector_div_constant_op<T, base_device::DEVICE_CPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const int& dim, T* result, const T* vector, const Real constant)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(Real))
#endif
        for (int i = 0; i < dim; i++)
        {
            result[i] = vector[i] / constant;
        }
    }
};

template <typename T>
struct vector_div_vector_op<T, base_device::DEVICE_CPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const int& dim, T* result, const T* vector1, const Real* vector2)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(Real))
#endif
        for (int i = 0; i < dim; i++)
        {
            result[i] = vector1[i] / vector2[i];
        }
    }
};

template <typename T>
struct axpy_op<T, base_device::DEVICE_CPU>
{
    void operator()(const int& dim,
                    const T* alpha,
                    const T* X,
                    const int& incX,
                    T* Y,
                    const int& incY)
    {
        BlasConnector::axpy(dim, *alpha, X, incX, Y, incY);
    }
};


template <typename T>
struct vector_add_vector_op<T, base_device::DEVICE_CPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const int& dim,
                    T* result,
                    const T* vector1,
                    const Real constant1,
                    const T* vector2,
                    const Real constant2)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 8192 / sizeof(T))
#endif
        for (int i = 0; i < dim; i++)
        {
            result[i] = vector1[i] * constant1 + vector2[i] * constant2;
        }
    }
};



template <typename FPTYPE>
struct dot_real_op<FPTYPE, base_device::DEVICE_CPU>
{
    FPTYPE operator()(const int& dim, const FPTYPE* psi_L, const FPTYPE* psi_R, const bool reduce)
    {
        FPTYPE result = BlasConnector::dot(dim, psi_L, 1, psi_R, 1);
        if (reduce)
        {
            Parallel_Reduce::reduce_pool(result);
        }
        return result;
    }
};

template <typename FPTYPE>
struct dot_real_op<std::complex<FPTYPE>, base_device::DEVICE_CPU>
{
    FPTYPE operator()(const int& dim,
                      const std::complex<FPTYPE>* psi_L,
                      const std::complex<FPTYPE>* psi_R,
                      const bool reduce)
    {
        // Note that  ddot_(2*dim,a,1,b,1) = REAL( zdotc_(dim,a,1,b,1) )
        const FPTYPE* pL = reinterpret_cast<const FPTYPE*>(psi_L);
        const FPTYPE* pR = reinterpret_cast<const FPTYPE*>(psi_R);
        FPTYPE result = BlasConnector::dot(2 * dim, pL, 1, pR, 1);
        if (reduce)
        {
            Parallel_Reduce::reduce_pool(result);
        }
        return result;
    }
};

template struct scal_op<float, base_device::DEVICE_CPU>;
template struct scal_op<double, base_device::DEVICE_CPU>;

template struct vector_mul_real_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct vector_mul_real_op<double, base_device::DEVICE_CPU>;
template struct vector_mul_real_op<std::complex<double>, base_device::DEVICE_CPU>;

template struct vector_mul_vector_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct vector_mul_vector_op<double, base_device::DEVICE_CPU>;
template struct vector_mul_vector_op<std::complex<double>, base_device::DEVICE_CPU>;

template struct vector_div_constant_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct vector_div_constant_op<double, base_device::DEVICE_CPU>;
template struct vector_div_constant_op<std::complex<double>, base_device::DEVICE_CPU>;

template struct vector_div_vector_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct vector_div_vector_op<double, base_device::DEVICE_CPU>;
template struct vector_div_vector_op<std::complex<double>, base_device::DEVICE_CPU>;

template struct axpy_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct axpy_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct axpy_op<double, base_device::DEVICE_CPU>;

template struct vector_add_vector_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct vector_add_vector_op<double, base_device::DEVICE_CPU>;
template struct vector_add_vector_op<std::complex<double>, base_device::DEVICE_CPU>;

template struct dot_real_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct dot_real_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct dot_real_op<double, base_device::DEVICE_CPU>;
} // namespace ModuleBase
