#include "module_hsolver/include/dngvd_op.h"
#include "src_pdiag/helper_cuda.h"

#include <cusolverDn.h>

#define cusolverErrcheck(res)                      \
    {                                              \
        cusolverAssert((res), __FILE__, __LINE__); \
    }

// cuSOLVER API errors
static const char* _cusolverGetErrorEnum(cusolverStatus_t error)
{
    switch (error)
    {
    case CUSOLVER_STATUS_SUCCESS:
        return "CUSOLVER_STATUS_SUCCESS";
    case CUSOLVER_STATUS_NOT_INITIALIZED:
        return "CUSOLVER_STATUS_NOT_INITIALIZED";
    case CUSOLVER_STATUS_ALLOC_FAILED:
        return "CUSOLVER_STATUS_ALLOC_FAILED";
    case CUSOLVER_STATUS_INVALID_VALUE:
        return "CUSOLVER_STATUS_INVALID_VALUE";
    case CUSOLVER_STATUS_ARCH_MISMATCH:
        return "CUSOLVER_STATUS_ARCH_MISMATCH";
    case CUSOLVER_STATUS_MAPPING_ERROR:
        return "CUSOLVER_STATUS_MAPPING_ERROR";
    case CUSOLVER_STATUS_EXECUTION_FAILED:
        return "CUSOLVER_STATUS_EXECUTION_FAILED";
    case CUSOLVER_STATUS_INTERNAL_ERROR:
        return "CUSOLVER_STATUS_INTERNAL_ERROR";
    case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
        return "CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
    case CUSOLVER_STATUS_NOT_SUPPORTED:
        return "CUSOLVER_STATUS_NOT_SUPPORTED ";
    case CUSOLVER_STATUS_ZERO_PIVOT:
        return "CUSOLVER_STATUS_ZERO_PIVOT";
    case CUSOLVER_STATUS_INVALID_LICENSE:
        return "CUSOLVER_STATUS_INVALID_LICENSE";
    }
    return "<unknown>";
}

inline void cusolverAssert(cusolverStatus_t code, const char* file, int line, bool abort = true)
{
    if (code != CUBLAS_STATUS_SUCCESS)
    {
        fprintf(stderr, "cuSOLVER Assert: %s %s %d\n", _cusolverGetErrorEnum(code), file, line);
        if (abort)
            exit(code);
    }
}

namespace hsolver
{

static cusolverDnHandle_t cusolver_H = nullptr;

void createCUSOLVERhandle()
{
    if (cusolver_H == nullptr)
    {
        cusolverErrcheck(cusolverDnCreate(&cusolver_H));
    }
}

void destoryCUSOLVERhandle()
{
    if (cusolver_H != nullptr)
    {
        cusolverErrcheck(cusolverDnDestroy(cusolver_H));
        cusolver_H = nullptr;
    }
}

static inline
void xhegvd_wrapper (
        const cublasFillMode_t& uplo,
        const int& n,
        std::complex<float> * A, const int& lda,
        std::complex<float> * B, const int& ldb,
        float * W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    int * devInfo = nullptr;
    int lwork = 0, info_gpu = 0;
    float2 * work = nullptr;
    checkCudaErrors(cudaMalloc((void**)&devInfo, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnChegvd_bufferSize(cusolver_H, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                                 reinterpret_cast<const float2 *>(A), lda,
                                                 reinterpret_cast<const float2 *>(B), ldb, W, &lwork));
    // allocate memery
    checkCudaErrors(cudaMalloc((void**)&work, sizeof(float2) * lwork));

    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnChegvd(cusolver_H, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                      reinterpret_cast<float2 *>(A), lda, reinterpret_cast<float2 *>(B), ldb, W, work, lwork, devInfo));

    checkCudaErrors(cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost));
    assert(0 == info_gpu);
    // free the buffer
    checkCudaErrors(cudaFree(work));
    checkCudaErrors(cudaFree(devInfo));
}

static inline
void xhegvd_wrapper (
        const cublasFillMode_t& uplo,
        const int& n,
        std::complex<double> * A, const int& lda,
        std::complex<double> * B, const int& ldb,
        double * W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    int * devInfo = nullptr;
    int lwork = 0, info_gpu = 0;
    double2 * work = nullptr;
    checkCudaErrors(cudaMalloc((void**)&devInfo, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnZhegvd_bufferSize(cusolver_H, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                                 reinterpret_cast<const double2 *>(A), lda,
                                                 reinterpret_cast<const double2 *>(B), ldb, W, &lwork));
    // allocate memery
    checkCudaErrors(cudaMalloc((void**)&work, sizeof(double2) * lwork));

    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnZhegvd(cusolver_H, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                      reinterpret_cast<double2 *>(A), lda, reinterpret_cast<double2 *>(B), ldb, W, work, lwork, devInfo));

    checkCudaErrors(cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost));
    assert(0 == info_gpu);
    // free the buffer
    checkCudaErrors(cudaFree(work));
    checkCudaErrors(cudaFree(devInfo));
}

static inline
void xheevd_wrapper (
        const cublasFillMode_t& uplo,
        const int& n,
        std::complex<float> * A, const int& lda,
        float * W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    int * devInfo = nullptr;
    int lwork = 0, info_gpu = 0;
    float2 * work = nullptr;
    cusolverDnHandle_t cusolverH = {};
    cusolverErrcheck(cusolverDnCreate(&cusolverH));
    checkCudaErrors(cudaMalloc((void**)&devInfo, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnCheevd_bufferSize(cusolverH, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                                 reinterpret_cast<const float2 *>(A), lda, W, &lwork));
    // allocate memery
    checkCudaErrors(cudaMalloc((void**)&work, sizeof(float2) * lwork));
    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnCheevd(cusolverH, CUSOLVER_EIG_MODE_VECTOR, uplo, n, reinterpret_cast<float2 *>(A), lda, W, work, lwork, devInfo));

    checkCudaErrors(cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost));
    assert(0 == info_gpu);
    checkCudaErrors(cudaFree(work));
    checkCudaErrors(cudaFree(devInfo));
    cusolverErrcheck(cusolverDnDestroy(cusolverH));
}

static inline
void xheevd_wrapper (
        const cublasFillMode_t& uplo,
        const int& n,
        std::complex<double> * A, const int& lda,
        double * W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    int * devInfo = nullptr;
    int lwork = 0, info_gpu = 0;
    double2 * work = nullptr;
    cusolverDnHandle_t cusolverH = {};
    cusolverErrcheck(cusolverDnCreate(&cusolverH));
    checkCudaErrors(cudaMalloc((void**)&devInfo, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnZheevd_bufferSize(cusolverH, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                                 reinterpret_cast<const double2 *>(A), lda, W, &lwork));
    // allocate memery
    checkCudaErrors(cudaMalloc((void**)&work, sizeof(double2) * lwork));
    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnZheevd(cusolverH, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                      reinterpret_cast<double2 *>(A), lda, W, work, lwork, devInfo));

    checkCudaErrors(cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost));
    assert(0 == info_gpu);
    checkCudaErrors(cudaFree(work));
    checkCudaErrors(cudaFree(devInfo));
    cusolverErrcheck(cusolverDnDestroy(cusolverH));
}

template <typename FPTYPE>
struct dngvx_op<FPTYPE, psi::DEVICE_GPU> {
    void operator()(
            const psi::DEVICE_GPU * d,
            const int nstart,
            const int ldh,
            const std::complex<FPTYPE> *A, // hcc
            const std::complex<FPTYPE> *B, // scc
            const int m, // nbands
            FPTYPE *W, // eigenvalue
            std::complex<FPTYPE> *V)
    {
        using transpose_op = matrixTranspose_op<FPTYPE, psi::DEVICE_GPU>;
        using matrixset_op = matrixSetToAnother<FPTYPE, psi::DEVICE_GPU>;
        // init A_eigenvectors, transpose_B and all_W
        std::complex<FPTYPE> * A_eigenvectors = nullptr, * transpose_B = nullptr;
        if (nstart == ldh) {
            checkCudaErrors(cudaMalloc((void **) &A_eigenvectors, sizeof(std::complex<FPTYPE>) * nstart * nstart));
            checkCudaErrors(cudaMalloc((void **) &transpose_B, sizeof(std::complex<FPTYPE>) * nstart * nstart));

            transpose_op()(d, nstart, nstart, A,A_eigenvectors);
            transpose_op()(d, nstart, nstart, B, transpose_B);
        } else if (nstart < ldh) {
            // nstart < ldh
            checkCudaErrors(cudaMalloc((void **) &A_eigenvectors, sizeof(std::complex<FPTYPE>) * nstart * nstart));
            checkCudaErrors(cudaMalloc((void **) &transpose_B, sizeof(std::complex<FPTYPE>) * nstart * nstart));

            matrixset_op()(d, nstart, A, ldh,  A_eigenvectors, nstart);
            matrixset_op()(d, nstart, B, ldh, transpose_B, nstart);

            transpose_op()(d,nstart,nstart,A_eigenvectors,A_eigenvectors);
            transpose_op()(d,nstart,nstart,transpose_B,transpose_B);
        } else if (nstart > ldh) {
            assert(nstart < ldh);
        }

        FPTYPE * all_W = nullptr;
        checkCudaErrors(cudaMalloc((void **) &all_W, sizeof(FPTYPE) * nstart));

        xhegvd_wrapper(CUBLAS_FILL_MODE_LOWER, nstart, A_eigenvectors, nstart,
                       transpose_B, nstart, all_W);

        // get eigenvalues and eigenvectors.  only m !
        checkCudaErrors(cudaMemcpy(W, all_W, sizeof(FPTYPE) * m, cudaMemcpyDeviceToDevice));

        if (ldh == nstart) {
            transpose_op()(d, nstart, nstart, V, V);
            checkCudaErrors(
                    cudaMemcpy(V, A_eigenvectors, sizeof(std::complex<FPTYPE>) * nstart * m, cudaMemcpyDeviceToDevice));
            transpose_op()(d, nstart, nstart, V, V);
        } else {
            transpose_op()(d, ldh, ldh, V, V);
            matrixset_op()(d, m, A_eigenvectors, nstart, V, ldh);
            transpose_op()(d, ldh, ldh, V, V);
        }
        // free resources and destroy
        checkCudaErrors(cudaFree(A_eigenvectors));
        checkCudaErrors(cudaFree(transpose_B));
        checkCudaErrors(cudaFree(all_W));
    }
};

template <>
void dngv_op<double, psi::DEVICE_GPU>::operator()(const psi::DEVICE_GPU* d,
                                                  const int nstart,
                                                  const int ldh,
                                                  const std::complex<double>* A,
                                                  const std::complex<double>* B,
                                                  double* W,
                                                  std::complex<double>* V)
{
    assert(nstart == ldh);
    // init A_eigenvectors & transpose_B
    double2 *A_eigenvectors, *transpose_B;
    checkCudaErrors(cudaMalloc((void**)&A_eigenvectors, sizeof(double2) * ldh * nstart));
    checkCudaErrors(cudaMalloc((void**)&transpose_B, sizeof(double2) * ldh * nstart));

    // transpose A, B  to A_eigenvectors, transpose_B
    matrixTranspose_op<double, psi::DEVICE_GPU>()(d, ldh, nstart, A, (std::complex<double>*)A_eigenvectors);
    matrixTranspose_op<double, psi::DEVICE_GPU>()(d, ldh, nstart, B, (std::complex<double>*)transpose_B);

    // init all_W
    double* all_W;
    checkCudaErrors(cudaMalloc((void**)&all_W, sizeof(double) * ldh));

    // prepare some values for cusolverDnZhegvd_bufferSize
    cusolverDnHandle_t cusolverH;
    cusolverErrcheck(cusolverDnCreate(&cusolverH));
    int* devInfo;
    checkCudaErrors(cudaMalloc((void**)&devInfo, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    int lwork = 0;
    cusolverErrcheck(cusolverDnZhegvd_bufferSize(
        cusolverH,
        CUSOLVER_EIG_TYPE_1, // itype = CUSOLVER_EIG_TYPE_1: A*x = (lambda)*B*x.
        CUSOLVER_EIG_MODE_VECTOR, // jobz = CUSOLVER_EIG_MODE_VECTOR : Compute eigenvalues and eigenvectors.
        CUBLAS_FILL_MODE_UPPER,
        ldh,
        A_eigenvectors,
        nstart,
        transpose_B,
        nstart,
        all_W,
        &lwork));

    // allocate memery
    cuDoubleComplex* d_work;
    checkCudaErrors(cudaMalloc((void**)&d_work, sizeof(cuDoubleComplex) * lwork));

    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnZhegvd(
        cusolverH,
        CUSOLVER_EIG_TYPE_1, // itype = CUSOLVER_EIG_TYPE_1: A*x = (lambda)*B*x.
        CUSOLVER_EIG_MODE_VECTOR, // jobz = CUSOLVER_EIG_MODE_VECTOR : Compute eigenvalues and eigenvectors.
        CUBLAS_FILL_MODE_UPPER,
        ldh,
        A_eigenvectors,
        nstart,
        transpose_B,
        nstart,
        all_W,
        d_work,
        lwork,
        devInfo));

    checkCudaErrors(cudaDeviceSynchronize());

    // get all eigenvalues and eigenvectors.
    checkCudaErrors(cudaMemcpy(W, all_W, sizeof(double) * ldh, cudaMemcpyDeviceToDevice));
    matrixTranspose_op<double, psi::DEVICE_GPU>()(d, ldh, nstart, (std::complex<double>*)A_eigenvectors, V);

    int info_gpu;
    checkCudaErrors(cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost));
    assert(0 == info_gpu);

    // free the buffer
    checkCudaErrors(cudaFree(d_work));
    // free resources and destroy
    checkCudaErrors(cudaFree(A_eigenvectors));
    checkCudaErrors(cudaFree(transpose_B));
    checkCudaErrors(cudaFree(all_W));
    checkCudaErrors(cudaFree(devInfo));
    cusolverErrcheck(cusolverDnDestroy(cusolverH));
}

template <typename FPTYPE>
struct dngvd_op<FPTYPE, psi::DEVICE_GPU> {
    void operator()(
            const psi::DEVICE_GPU *d,
            const int nstart,
            const int ldh,
            const std::complex<FPTYPE> *A, // hcc
            const std::complex<FPTYPE> *B, // scc
            FPTYPE *W, // eigenvalue
            std::complex<FPTYPE> *V)
    {
        assert(nstart == ldh);
        // A to V
        checkCudaErrors(cudaMemcpy(V, A, sizeof(std::complex<FPTYPE>) * ldh * nstart, cudaMemcpyDeviceToDevice));

        xhegvd_wrapper(CUBLAS_FILL_MODE_UPPER, nstart, V, ldh,
                       (std::complex<FPTYPE> *)B, ldh, W);
    }
};

template <typename FPTYPE>
struct dnevx_op<FPTYPE, psi::DEVICE_GPU> {
    void operator()(
            const psi::DEVICE_GPU *d,
            const int nstart,
            const int ldh,
            const std::complex<FPTYPE> *A, // hcc
            const int m,
            FPTYPE *W, // eigenvalue
            std::complex<FPTYPE> *V)
    {
        assert(nstart <= ldh);

        // A to V
        checkCudaErrors(cudaMemcpy(V, A, sizeof(double2) * nstart * ldh, cudaMemcpyDeviceToDevice));

        xheevd_wrapper(CUBLAS_FILL_MODE_LOWER, nstart, V, nstart, W);

        // get eigenvalues and eigenvectors.  only m !
        checkCudaErrors(cudaMemcpy(W, all_W, sizeof(FPTYPE) * m, cudaMemcpyDeviceToDevice));

        if (ldh == nstart) {
            transpose_op()(d, nstart, nstart, V, V);
            checkCudaErrors(
                    cudaMemcpy(V, A_eigenvectors, sizeof(std::complex<FPTYPE>) * nstart * m, cudaMemcpyDeviceToDevice));
            transpose_op()(d, nstart, nstart, V, V);
        } else {
            transpose_op()(d, ldh, ldh, V, V);
            matrixset_op()(d, m, A_eigenvectors, nstart, V, ldh);
            transpose_op()(d, ldh, ldh, V, V);
        }

        // free resources and destroy
        checkCudaErrors(cudaFree(A_eigenvectors));
        checkCudaErrors(cudaFree(all_W));
    }
};

template struct dngvx_op<float, psi::DEVICE_GPU>;
template struct dngvd_op<float, psi::DEVICE_GPU>;
template struct dnevx_op<float, psi::DEVICE_GPU>;
template struct dngvx_op<double, psi::DEVICE_GPU>;
template struct dngvd_op<double, psi::DEVICE_GPU>;
template struct dnevx_op<double, psi::DEVICE_GPU>;

} // namespace hsolver