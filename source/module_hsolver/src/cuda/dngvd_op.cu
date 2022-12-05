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

template <>
void dngvx_op<double, psi::DEVICE_GPU>::operator()(const psi::DEVICE_GPU* d,
                                                   const int nstart,
                                                   const int ldh,
                                                   const std::complex<double>* A,
                                                   const std::complex<double>* B,
                                                   const int m,
                                                   double* W,
                                                   std::complex<double>* V)
{
    // init A_eigenvectors, transpose_B and all_W
    double2 *A_eigenvectors, *transpose_B;
    if (nstart == ldh)
    {
        checkCudaErrors(cudaMalloc((void**)&A_eigenvectors, sizeof(double2) * nstart * nstart));
        checkCudaErrors(cudaMalloc((void**)&transpose_B, sizeof(double2) * nstart * nstart));

        matrixTranspose_op<double, psi::DEVICE_GPU>()(d, nstart, nstart, A, (std::complex<double>*)A_eigenvectors);
        matrixTranspose_op<double, psi::DEVICE_GPU>()(d, nstart, nstart, B, (std::complex<double>*)transpose_B);
    }
    else if (nstart < ldh)
    {
        // nstart < ldh
        checkCudaErrors(cudaMalloc((void**)&A_eigenvectors, sizeof(double2) * nstart * nstart));
        checkCudaErrors(cudaMalloc((void**)&transpose_B, sizeof(double2) * nstart * nstart));

        matrixSetToAnother<double, psi::DEVICE_GPU>()(d, nstart, A, ldh, (std::complex<double>*)A_eigenvectors, nstart);
        matrixSetToAnother<double, psi::DEVICE_GPU>()(d, nstart, B, ldh, (std::complex<double>*)transpose_B, nstart);

        matrixTranspose_op<double, psi::DEVICE_GPU>()(d,
                                                      nstart,
                                                      nstart,
                                                      (std::complex<double>*)A_eigenvectors,
                                                      (std::complex<double>*)A_eigenvectors);
        matrixTranspose_op<double, psi::DEVICE_GPU>()(d,
                                                      nstart,
                                                      nstart,
                                                      (std::complex<double>*)transpose_B,
                                                      (std::complex<double>*)transpose_B);
    }
    else if (nstart > ldh)
    {
        assert(nstart < ldh);
    }

    double* all_W;
    checkCudaErrors(cudaMalloc((void**)&all_W, sizeof(double) * nstart));

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
        CUBLAS_FILL_MODE_LOWER,
        nstart,
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
        CUBLAS_FILL_MODE_LOWER,
        nstart,
        A_eigenvectors,
        nstart,
        transpose_B,
        nstart,
        all_W,
        d_work,
        lwork,
        devInfo));

    checkCudaErrors(cudaDeviceSynchronize());

    // get eigenvalues and eigenvectors.  only m !
    checkCudaErrors(cudaMemcpy(W, all_W, sizeof(double) * m, cudaMemcpyDeviceToDevice));

    if (ldh == nstart)
    {
        matrixTranspose_op<double, psi::DEVICE_GPU>()(d, nstart, nstart, V, V);
        checkCudaErrors(
            cudaMemcpy(V, A_eigenvectors, sizeof(std::complex<double>) * nstart * m, cudaMemcpyDeviceToDevice));
        matrixTranspose_op<double, psi::DEVICE_GPU>()(d, nstart, nstart, V, V);
    }
    else
    {
        matrixTranspose_op<double, psi::DEVICE_GPU>()(d, ldh, ldh, V, V);
        matrixSetToAnother<double, psi::DEVICE_GPU>()(d, m, (std::complex<double>*)A_eigenvectors, nstart, V, ldh);
        matrixTranspose_op<double, psi::DEVICE_GPU>()(d, ldh, ldh, V, V);
    }

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

} // namespace hsolver