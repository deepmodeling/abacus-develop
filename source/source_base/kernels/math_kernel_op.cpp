#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_external/blas_connector.h"

#include <iomanip>
#include <iostream>

namespace ModuleBase
{

template <typename T>
struct gemv_op<T, base_device::DEVICE_CPU>
{
    void operator()(const char& trans,
                    const int& m,
                    const int& n,
                    const T* alpha,
                    const T* A,
                    const int& lda,
                    const T* X,
                    const int& incx,
                    const T* beta,
                    T* Y,
                    const int& incy)
    {
        BlasConnector::gemv(trans, m, n, *alpha, A, lda, X, incx, *beta, Y, incy);
    }
};

template <typename T>
struct gemm_op<T, base_device::DEVICE_CPU>
{
    void operator()(const char& transa,
                    const char& transb,
                    const int& m,
                    const int& n,
                    const int& k,
                    const T* alpha,
                    const T* a,
                    const int& lda,
                    const T* b,
                    const int& ldb,
                    const T* beta,
                    T* c,
                    const int& ldc)
    {
        BlasConnector::gemm(transb, transa, n, m, k, *alpha, b, ldb, a, lda, *beta, c, ldc);
    }
};

#ifdef __DSP
template <typename T>
struct gemm_op_mt<T, base_device::DEVICE_CPU>
{
    void operator()(const char& transa,
                    const char& transb,
                    const int& m,
                    const int& n,
                    const int& k,
                    const T* alpha,
                    const T* a,
                    const int& lda,
                    const T* b,
                    const int& ldb,
                    const T* beta,
                    T* c,
                    const int& ldc)
    {
        BlasConnector::gemm(transb, transa, n, m, k, *alpha, b, ldb, a, lda, *beta, c, ldc, base_device::AbacusDevice_t::DspDevice);
    }
};
#endif

template <typename T>
struct matrixTranspose_op<T, base_device::DEVICE_CPU>
{
    void operator()(const int& row,
                    const int& col,
                    const T* input_matrix,
                    T* output_matrix)
    {
        T* temp = nullptr;
        base_device::memory::resize_memory_op<T, base_device::DEVICE_CPU>()(temp, row * col, "MTransOp");
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 8192 / sizeof(T))
#endif
        for (int j = 0; j < col; j++)
        {
            for (int i = 0; i < row; i++)
            {
                temp[j * row + i] = input_matrix[i * col + j];
            }
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 8192 / sizeof(T))
#endif
        for (int i = 0; i < row * col; i++)
        {
            output_matrix[i] = temp[i];
        }
        base_device::memory::delete_memory_op<T, base_device::DEVICE_CPU>()(temp);
    }
};

template <typename T>
struct matrixCopy<T, base_device::DEVICE_CPU>
{
    void operator()(const int& n1, const int& n2, const T* A, const int& LDA, T* B, const int& LDB)
    {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 8192 / sizeof(T))
#endif
        for (int i = 0; i < n1; i++)
        {
            for (int j = 0; j < n2; j++)
            {
                B[i * LDB + j] = A[i * LDA + j];
            }
        }
    }
};

template <typename T>
struct matrix_mul_vector_op<T, base_device::DEVICE_CPU> {
    using Real = typename GetTypeReal<T>::type;
    void operator()(const int& m, const int &n,
                  T *a,
                  const int &lda,
                  const Real *b,
                  const Real alpha,
                  T *c,
                  const int &ldc){
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 8192 / sizeof(T))
#endif
        for (int j = 0; j < n; j++){
            for (int i = 0; i < m; i++){
                c[j * ldc + i] = a[j * lda + i] * b[j] * alpha;
            }
        }

    }
};

template struct gemv_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct gemv_op<float, base_device::DEVICE_CPU>;
template struct gemm_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct gemm_op<float, base_device::DEVICE_CPU>;
template struct matrixTranspose_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct matrixCopy<std::complex<float>, base_device::DEVICE_CPU>;
template struct matrix_mul_vector_op<std::complex<float>, base_device::DEVICE_CPU>;

template struct gemv_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct gemv_op<double, base_device::DEVICE_CPU>;
template struct gemm_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct gemm_op<double, base_device::DEVICE_CPU>;
template struct matrixTranspose_op<std::complex<double>, base_device::DEVICE_CPU>;
template struct matrixCopy<double, base_device::DEVICE_CPU>;
template struct matrixCopy<std::complex<double>, base_device::DEVICE_CPU>;
template struct matrix_mul_vector_op<double, base_device::DEVICE_CPU>;
template struct matrix_mul_vector_op<std::complex<double>, base_device::DEVICE_CPU>;

#ifdef __LCAO
template struct matrixTranspose_op<double, base_device::DEVICE_CPU>;
#endif
#ifdef __DSP
template struct gemm_op_mt<std::complex<float>, base_device::DEVICE_CPU>;
template struct gemm_op_mt<std::complex<double>, base_device::DEVICE_CPU>;
#endif
} // namespace hsolver