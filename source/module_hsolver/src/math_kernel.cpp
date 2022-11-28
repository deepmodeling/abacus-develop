#include "module_hsolver/include/math_kernel.h"

#include <iomanip>
#include <iostream>

namespace hsolver
{

// CPU specialization of actual computation.
template <typename FPTYPE> 
struct zdot_real_op<FPTYPE, psi::DEVICE_CPU> {
  FPTYPE operator() (
      const psi::DEVICE_CPU* d,
      const int& dim,
      const std::complex<FPTYPE>* psi_L,
      const std::complex<FPTYPE>* psi_R,
      const bool reduce) 
  {
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // qianrui modify 2021-3-14
    // Note that  ddot_(2*dim,a,1,b,1) = REAL( zdotc_(dim,a,1,b,1) )
    const FPTYPE* pL = reinterpret_cast<const FPTYPE*>(psi_L);
    const FPTYPE* pR = reinterpret_cast<const FPTYPE*>(psi_R);
    FPTYPE result = BlasConnector::dot(2 * dim, pL, 1, pR, 1);
    if (reduce) {
      Parallel_Reduce::reduce_double_pool(result);
    }
    return result;
  }
};

template <typename FPTYPE> struct vector_div_constant_op<FPTYPE, psi::DEVICE_CPU>
{
    void operator()(const psi::DEVICE_CPU* d,
                    const int dim,
                    std::complex<FPTYPE>* result,
                    const std::complex<FPTYPE>* vector,
                    const FPTYPE constant)
    {
        for (int i = 0; i < dim; i++)
        {
            result[i] = vector[i] / constant;
        }
    }
};

template <typename FPTYPE> struct vector_mul_vector_op<FPTYPE, psi::DEVICE_CPU>
{
    void operator()(const psi::DEVICE_CPU* d,
                    const int& dim,
                    std::complex<FPTYPE>* result,
                    const std::complex<FPTYPE>* vector1,
                    const FPTYPE* vector2)
    {
        for (int i = 0; i < dim; i++)
        {
            result[i] = vector1[i] * vector2[i];
        }
    }
};

template <typename FPTYPE> struct vector_div_vector_op<FPTYPE, psi::DEVICE_CPU>
{
    void operator()(const psi::DEVICE_CPU* d,
                    const int& dim,
                    std::complex<FPTYPE>* result,
                    const std::complex<FPTYPE>* vector1,
                    const FPTYPE* vector2)
    {
        for (int i = 0; i < dim; i++)
        {
            result[i] = vector1[i] / vector2[i];
        }
    }
};

template <typename FPTYPE> struct constantvector_addORsub_constantVector_op<FPTYPE, psi::DEVICE_CPU>
{
    void operator()(const psi::DEVICE_CPU* d,
                    const int& dim,
                    std::complex<FPTYPE>* result,
                    const std::complex<FPTYPE>* vector1,
                    const FPTYPE constant1,
                    const std::complex<FPTYPE>* vector2,
                    const FPTYPE constant2)
    {
        for (int i = 0; i < dim; i++)
        {
            result[i] = vector1[i] * constant1 + vector2[i] * constant2;
        }
    }
};

template <>
void scal_op<double, psi::DEVICE_CPU>::operator()(const psi::DEVICE_CPU* d,
                                                  const int& N,
                                                  const std::complex<double>* alpha,
                                                  std::complex<double>* X,
                                                  const int& incx)
{
    zscal_(&N, alpha, X, &incx);
}

template <>
void scal_op<float, psi::DEVICE_CPU>::operator()(const psi::DEVICE_CPU* d,
                                                 const int& N,
                                                 const std::complex<float>* alpha,
                                                 std::complex<float>* X,
                                                 const int& incx)

{
    cscal_(&N, alpha, X, &incx);
}

template <>
void gemv_op<float, psi::DEVICE_CPU>::operator()(const psi::DEVICE_CPU* d,
                                                 const char& trans,
                                                 const int& m,
                                                 const int& n,
                                                 const std::complex<float>* alpha,
                                                 const std::complex<float>* A,
                                                 const int& lda,
                                                 const std::complex<float>* X,
                                                 const int& incx,
                                                 const std::complex<float>* beta,
                                                 std::complex<float>* Y,
                                                 const int& incy)
{
    cgemv_(&trans, &m, &n, alpha, A, &lda, X, &incx, beta, Y, &incy);
}

template <>
void gemv_op<double, psi::DEVICE_CPU>::operator()(const psi::DEVICE_CPU* d,
                                                  const char& trans,
                                                  const int& m,
                                                  const int& n,
                                                  const std::complex<double>* alpha,
                                                  const std::complex<double>* A,
                                                  const int& lda,
                                                  const std::complex<double>* X,
                                                  const int& incx,
                                                  const std::complex<double>* beta,
                                                  std::complex<double>* Y,
                                                  const int& incy)
{
    zgemv_(&trans, &m, &n, alpha, A, &lda, X, &incx, beta, Y, &incy);
}

template <>
void axpy_op<double, psi::DEVICE_CPU>::operator()(const psi::DEVICE_CPU* d,
                                                  const int& dim,
                                                  const std::complex<double>* alpha,
                                                  const std::complex<double>* X,
                                                  const int& incX,
                                                  std::complex<double>* Y,
                                                  const int& incY)
{
    zaxpy_(&dim, alpha, X, &incX, Y, &incY);
}

template <>
void axpy_op<float, psi::DEVICE_CPU>::operator()(const psi::DEVICE_CPU* d,
                                                 const int& dim,
                                                 const std::complex<float>* alpha,
                                                 const std::complex<float>* X,
                                                 const int& incX,
                                                 std::complex<float>* Y,
                                                 const int& incY)
{
    caxpy_(&dim, alpha, X, &incX, Y, &incY);
}

template <>
void gemm_op<float, psi::DEVICE_CPU>::operator()(const psi::DEVICE_CPU* d,
                                                 const char& transa, 
                                                 const char& transb, 
                                                 const int& m, 
                                                 const int& n, 
                                                 const int& k,
                                                 const std::complex<float> *alpha, 
                                                 const std::complex<float> *a, 
                                                 const int& lda, 
                                                 const std::complex<float> *b, 
                                                 const int& ldb,
                                                 const std::complex<float> *beta, 
                                                 std::complex<float> *c, 
                                                 const int& ldc)
{
    cgemm_(&transa, &transb, &m, &n ,&k, alpha, a ,&lda, b, &ldb, beta, c, &ldc);
}

template <>
void gemm_op<double, psi::DEVICE_CPU>::operator()(const psi::DEVICE_CPU* d,
                                                 const char& transa, 
                                                 const char& transb, 
                                                 const int& m, 
                                                 const int& n, 
                                                 const int& k,
                                                 const std::complex<double> *alpha, 
                                                 const std::complex<double> *a, 
                                                 const int& lda, 
                                                 const std::complex<double> *b, 
                                                 const int& ldb,
                                                 const std::complex<double> *beta, 
                                                 std::complex<double> *c, 
                                                 const int& ldc)
{
    zgemm_(&transa, &transb, &m, &n ,&k, alpha, a ,&lda, b, &ldb, beta, c, &ldc);
}

template <typename FPTYPE> struct matrixTranspose_op<FPTYPE, psi::DEVICE_CPU>
{
    void operator()(const psi::DEVICE_CPU* d,
                    const int& row,
                    const int& col,
                    const std::complex<FPTYPE>* input_matrix,
                    std::complex<FPTYPE>* output_matrix)
    {
        std::complex<FPTYPE>* temp = nullptr;
        psi::memory::resize_memory_op<std::complex<FPTYPE>, psi::DEVICE_CPU>()(d, temp, row * col);
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                temp[j * row + i] = input_matrix[i * col + j];
            }
        }
        for (int i = 0; i < row * col; i++)
        {
            output_matrix[i] = temp[i];
        }
        psi::memory::delete_memory_op<std::complex<FPTYPE>, psi::DEVICE_CPU>()(d, temp);
    }
};

template <typename FPTYPE> struct matrixSetToAnother<FPTYPE, psi::DEVICE_CPU>
{
    void operator()(const psi::DEVICE_CPU* d,
                    const int& n,
                    const std::complex<FPTYPE>* A,
                    const int& LDA,
                    std::complex<FPTYPE>* B,
                    const int& LDB)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < LDA; j++)
            {
                B[i * LDB + j] = A[i * LDA + j];
            }
        }
    }
};



// Explicitly instantiate functors for the types of functor registered.
template struct zdot_real_op<double, psi::DEVICE_CPU>;
template struct vector_div_constant_op<double, psi::DEVICE_CPU>;
template struct vector_mul_vector_op<double, psi::DEVICE_CPU>;
template struct vector_div_vector_op<double, psi::DEVICE_CPU>;
template struct constantvector_addORsub_constantVector_op<double, psi::DEVICE_CPU>;

template struct matrixTranspose_op<double, psi::DEVICE_CPU>;
template struct matrixSetToAnother<double, psi::DEVICE_CPU>;


} // namespace hsolver