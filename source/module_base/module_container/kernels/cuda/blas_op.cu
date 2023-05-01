#include "../blas_op.h"

namespace container {
namespace op {

static cublasHandle_t cublas_handle = nullptr;

// static inline
// void xdot_wrapper(const int &n, const float * x, const int &incx, const float * y, const int &incy, float &result) {
//     cublasErrcheck(cublasSdot(cublas_handle, n, x, incx, y, incy, &result));
// }
//
// static inline
// void xdot_wrapper(const int &n, const double * x, const int &incx, const double * y, const int &incy, double &result) {
//     cublasErrcheck(cublasDdot(cublas_handle, n, x, incx, y, incy, &result));
// }

static inline
void xscal_wrapper(const int &n, const std::complex<float> * alpha, std::complex<float> * X, const int &incx) {
    cublasErrcheck(cublasCscal(cublas_handle, n, (float2*)alpha, (float2*)X, incx));
}

static inline
void xscal_wrapper(const int &n, const std::complex<double> * alpha, std::complex<double> * X, const int &incx) {
    cublasErrcheck(cublasZscal(cublas_handle, n, (double2*)alpha, (double2*)X, incx));
}

void createBlasHandle() {
    if (cublas_handle == nullptr) {
        cublasErrcheck(cublasCreate(&cublas_handle));
    }
}

void destoryBlasHandle() {
    if (cublas_handle != nullptr) {
        cublasErrcheck(cublasDestroy(cublas_handle));
        cublas_handle = nullptr;
    }
}

template <typename T>
struct scal_op<T, DEVICE_GPU> {
    void operator()(
            const int &N,
            const std::complex<T> *alpha,
            std::complex<T> *X,
            const int &incx)
    {
        xscal_wrapper(N, alpha, X, incx);
    }
};

template <>
struct axpy_op<float, DEVICE_GPU> {
    void operator()(
            const int &N,
            const std::complex<float> *alpha,
            const std::complex<float> *X,
            const int &incX,
            std::complex<float> *Y,
            const int &incY)
    {
        cublasErrcheck(cublasCaxpy(cublas_handle, N, (float2 *) alpha, (float2 *) X, incX, (float2 *) Y, incY));
    }
};

template <>
struct axpy_op<double, DEVICE_GPU> {
    void operator()(
            const int &N,
            const std::complex<double> *alpha,
            const std::complex<double> *X,
            const int &incX,
            std::complex<double> *Y,
            const int &incY)
    {
        cublasErrcheck(cublasZaxpy(cublas_handle, N, (double2 *) alpha, (double2 *) X, incX, (double2 *) Y, incY));
    }
};

template <>
struct gemv_op<float, DEVICE_GPU> {
    void operator()(
            const char &trans,
            const int &m,
            const int &n,
            const std::complex<float> *alpha,
            const std::complex<float> *A,
            const int &lda,
            const std::complex<float> *X,
            const int &incx,
            const std::complex<float> *beta,
            std::complex<float> *Y,
            const int &incy) {
        cublasOperation_t cutrans = {};
        if (trans == 'N') {
            cutrans = CUBLAS_OP_N;
        } else if (trans == 'T') {
            cutrans = CUBLAS_OP_T;
        } else if (trans == 'C') {
            cutrans = CUBLAS_OP_C;
        }
        cublasErrcheck(
                cublasCgemv(cublas_handle, cutrans, m, n, (float2 *) alpha, (float2 *) A, lda, (float2 *) X, incx,
                            (float2 *) beta, (float2 *) Y, incy));
    }
};

template <>
struct gemv_op<double, DEVICE_GPU> {
    void operator()(
            const char &trans,
            const int &m,
            const int &n,
            const std::complex<double> *alpha,
            const std::complex<double> *A,
            const int &lda,
            const std::complex<double> *X,
            const int &incx,
            const std::complex<double> *beta,
            std::complex<double> *Y,
            const int &incy)
    {
        cublasOperation_t cutrans = {};
        if (trans == 'N') {
            cutrans = CUBLAS_OP_N;
        } else if (trans == 'T') {
            cutrans = CUBLAS_OP_T;
        } else if (trans == 'C') {
            cutrans = CUBLAS_OP_C;
        }
        cublasErrcheck(
                cublasZgemv(cublas_handle, cutrans, m, n, (double2 *) alpha, (double2 *) A, lda, (double2 *) X, incx,
                            (double2 *) beta, (double2 *) Y, incy));
    }
};


template <>
struct gemm_op<float, DEVICE_GPU> {
    void operator()(
            const char &transa,
            const char &transb,
            const int &m,
            const int &n,
            const int &k,
            const std::complex<float> *alpha,
            const std::complex<float> *a,
            const int &lda,
            const std::complex<float> *b,
            const int &ldb,
            const std::complex<float> *beta,
            std::complex<float> *c,
            const int &ldc)
    {
        cublasOperation_t cutransA = {};
        cublasOperation_t cutransB = {};
        // cutransA
        if (transa == 'N') {
            cutransA = CUBLAS_OP_N;
        } else if (transa == 'T') {
            cutransA = CUBLAS_OP_T;
        } else if (transa == 'C') {
            cutransA = CUBLAS_OP_C;
        }
        // cutransB
        if (transb == 'N') {
            cutransB = CUBLAS_OP_N;
        } else if (transb == 'T') {
            cutransB = CUBLAS_OP_T;
        } else if (transb == 'C') {
            cutransB = CUBLAS_OP_C;
        }
        cublasErrcheck(cublasCgemm(cublas_handle, cutransA, cutransB, m, n, k, (float2 *) alpha, (float2 *) a, lda,
                                   (float2 *) b, ldb, (float2 *) beta, (float2 *) c, ldc));
    }
};

template <>
struct gemm_op<double, DEVICE_GPU> {
    void operator()(
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
        cublasOperation_t cutransA;
        cublasOperation_t cutransB;
        // cutransA
        if (transa == 'N'){
            cutransA = CUBLAS_OP_N;
        }
        else if (transa == 'T'){
            cutransA = CUBLAS_OP_T;
        }
        else if (transa == 'C'){
            cutransA = CUBLAS_OP_C;
        }
        // cutransB
        if (transb == 'N'){
            cutransB = CUBLAS_OP_N;
        }
        else if (transb == 'T'){
            cutransB = CUBLAS_OP_T;
        }
        else if (transb == 'C'){
            cutransB = CUBLAS_OP_C;
        }
        cublasErrcheck(cublasZgemm(cublas_handle, cutransA, cutransB, m, n ,k, (double2*)alpha, (double2*)a , lda, (double2*)b, ldb, (double2*)beta, (double2*)c, ldc));
    }
};

// Explicitly instantiate functors for the types of functor registered.
template struct axpy_op<float, DEVICE_GPU>;
template struct scal_op<float, DEVICE_GPU>;

template struct axpy_op<double, DEVICE_GPU>;
template struct scal_op<double, DEVICE_GPU>;

} // namespace op
} // namespace container