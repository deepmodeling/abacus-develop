#include "blas_op.h"

namespace container {
namespace op {

//// CPU specialization of actual computation.
//template <typename T>
//struct zdot_real_op<T, DEVICE_CPU> {
//    T operator() (
//            const DEVICE_CPU* d,
//            const int& dim,
//            const std::complex<T>* psi_L,
//            const std::complex<T>* psi_R,
//            const bool reduce)
//    {
//        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//        // qianrui modify 2021-3-14
//        // Note that  ddot_(2*dim,a,1,b,1) = REAL( zdotc_(dim,a,1,b,1) )
//        const T* pL = reinterpret_cast<const T*>(psi_L);
//        const T* pR = reinterpret_cast<const T*>(psi_R);
//        T result = BlasConnector::dot(2 * dim, pL, 1, pR, 1);
//        if (reduce) {
//            Parallel_Reduce::reduce_double_pool(result);
//        }
//        return result;
//    }
//};

template <typename T>
struct scal_op<T, DEVICE_CPU> {
    void operator()(
            const int &N,
            const std::complex<T> *alpha,
            std::complex<T> *X,
            const int &incx)
    {
        BlasConnector::scal(N, *alpha, X, incx);
    }
};

template <typename T>
struct gemv_op<T, DEVICE_CPU> {
    void operator()(
            const DEVICE_CPU *d,
            const char &trans,
            const int &m,
            const int &n,
            const std::complex<T> *alpha,
            const std::complex<T> *A,
            const int &lda,
            const std::complex<T> *X,
            const int &incx,
            const std::complex<T> *beta,
            std::complex<T> *Y,
            const int &incy)
    {
        BlasConnector::gemv(trans, m, n, *alpha, A, lda, X, incx, *beta, Y, incy);
    }
};

template <typename T>
struct axpy_op<T, DEVICE_CPU> {
    void operator()(
            const DEVICE_CPU * /*ctx*/,
            const int &dim,
            const std::complex<T> *alpha,
            const std::complex<T> *X,
            const int &incX,
            std::complex<T> *Y,
            const int &incY)
    {
        BlasConnector::axpy(dim, *alpha, X, incX, Y, incY);
    }
};

template <typename T>
struct gemm_op<T, DEVICE_CPU> {
    void operator()(
            const DEVICE_CPU * /*ctx*/,
            const char &transa,
            const char &transb,
            const int &m,
            const int &n,
            const int &k,
            const std::complex<T> *alpha,
            const std::complex<T> *a,
            const int &lda,
            const std::complex<T> *b,
            const int &ldb,
            const std::complex<T> *beta,
            std::complex<T> *c,
            const int &ldc)
    {
        BlasConnector::gemm(transb, transa, n, m, k, *alpha, b, ldb, a, lda, *beta, c, ldc);
    }
};

// Explicitly instantiate functors for the types of functor registered.
template struct scal_op<float, DEVICE_CPU>;
template struct axpy_op<float, DEVICE_CPU>;
template struct gemv_op<float, DEVICE_CPU>;
template struct gemm_op<float, DEVICE_CPU>;

template struct scal_op<double, DEVICE_CPU>;
template struct axpy_op<double, DEVICE_CPU>;
template struct gemv_op<double, DEVICE_CPU>;
template struct gemm_op<double, DEVICE_CPU>;

} // namespace op
} // namespace container