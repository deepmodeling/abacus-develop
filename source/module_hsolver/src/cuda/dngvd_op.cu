#include "module_hsolver/include/dngvd_op.h"


#include <cusolverDn.h>

namespace hsolver
{

template <>
void dngvx_op<double, psi::DEVICE_GPU>::operator()(const psi::DEVICE_GPU* d,
                                                   const int row,
                                                   const int col,
                                                   const std::complex<double>* A,
                                                   const std::complex<double>* B,
                                                   const int m,
                                                   double* W,
                                                   std::complex<double>* V)
{
    createBLAShandle();
    // init A_eigenvectors, transpose_B and all_W
    double2 *A_eigenvectors, *transpose_B;
    cudaMalloc((void**)&A_eigenvectors, sizeof(double2) * row * col);
    cudaMalloc((void**)&transpose_B, sizeof(double2) * row * col);

    matrixTranspose_op<double, psi::DEVICE_GPU>()(d, col, row, A, (std::complex<double>*)A_eigenvectors);
    matrixTranspose_op<double, psi::DEVICE_GPU>()(d, col, row, B, (std::complex<double>*)transpose_B);




    double* all_W;
    cudaMalloc((void**)&all_W, sizeof(double) * col);

    // prepare some values for cusolverDnZhegvd_bufferSize
    cusolverDnHandle_t cusolverH;
    cusolverDnCreate(&cusolverH);
    int* devInfo;
    cudaMalloc((void**)&devInfo, sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    int lwork = 0;
    cusolverDnZhegvd_bufferSize(
        cusolverH,
        CUSOLVER_EIG_TYPE_1, // itype = CUSOLVER_EIG_TYPE_1: A*x = (lambda)*B*x.
        CUSOLVER_EIG_MODE_VECTOR, // jobz = CUSOLVER_EIG_MODE_VECTOR : Compute eigenvalues and eigenvectors.
        CUBLAS_FILL_MODE_LOWER,
        col,
        A_eigenvectors,
        // (double2)A,
        row,
        transpose_B,
        // (double2)B,
        row,
        all_W,
        &lwork);

    // allocate memery
    cuDoubleComplex* d_work;
    cudaMalloc((void**)&d_work, sizeof(cuDoubleComplex) * lwork);

    // compute eigenvalues and eigenvectors.
    cusolverDnZhegvd(
        cusolverH,
        CUSOLVER_EIG_TYPE_1, // itype = CUSOLVER_EIG_TYPE_1: A*x = (lambda)*B*x.
        CUSOLVER_EIG_MODE_VECTOR, // jobz = CUSOLVER_EIG_MODE_VECTOR : Compute eigenvalues and eigenvectors.
        CUBLAS_FILL_MODE_LOWER,
        col,
        A_eigenvectors,
        // (double2)A,
        row,
        transpose_B,
        // (double2)B,
        row,
        all_W,
        d_work,
        lwork,
        devInfo);

    cudaDeviceSynchronize();

    // get eigenvalues and eigenvectors.  only m !
    cudaMemcpy(W, all_W, sizeof(double) * m, cudaMemcpyDeviceToDevice);

    cudaMemcpy(V, A_eigenvectors, sizeof(std::complex<double>)*col*m, cudaMemcpyDeviceToDevice);

    matrixTranspose_op<double, psi::DEVICE_GPU>()(d, col, row, V, V);

    int info_gpu;
    cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost);
    assert(0 == info_gpu);
    
    // free the buffer
    cudaFree(d_work);
    // free resources and destroy
    cudaFree(A_eigenvectors);
    cudaFree(all_W);
    cudaFree(devInfo);
    cusolverDnDestroy(cusolverH);
    destoryBLAShandle();
}

template <>
void dngv_op<double, psi::DEVICE_GPU>::operator()(const psi::DEVICE_GPU* d,
                                                  const int row,
                                                  const int col,
                                                  const std::complex<double>* A,
                                                  const std::complex<double>* B,
                                                  double* W,
                                                  std::complex<double>* V)
{
    createBLAShandle();
    // init A_eigenvectors & transpose_B
    double2 *A_eigenvectors, *transpose_B;
    cudaMalloc((void**)&A_eigenvectors, sizeof(double2) * row * col);
    cudaMalloc((void**)&transpose_B, sizeof(double2) * row * col);
    
    // transpose A to A_eigenvectors
    matrixTranspose_op<double, psi::DEVICE_GPU>()(d, row, col, A, (std::complex<double>*)A_eigenvectors);
    // transpose B to transpose_B
    matrixTranspose_op<double, psi::DEVICE_GPU>()(d, row, col, B, (std::complex<double>*)transpose_B);
    
    // init all_W
    double* all_W;
    cudaMalloc((void**)&all_W, sizeof(double) * row);

    // prepare some values for cusolverDnZhegvd_bufferSize
    cusolverDnHandle_t cusolverH;
    cusolverDnCreate(&cusolverH);
    int* devInfo;
    cudaMalloc((void**)&devInfo, sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    int lwork = 0;
    cusolverDnZhegvd_bufferSize(
        cusolverH,
        CUSOLVER_EIG_TYPE_1, // itype = CUSOLVER_EIG_TYPE_1: A*x = (lambda)*B*x.
        CUSOLVER_EIG_MODE_VECTOR, // jobz = CUSOLVER_EIG_MODE_VECTOR : Compute eigenvalues and eigenvectors.
        CUBLAS_FILL_MODE_UPPER,
        row,
        A_eigenvectors,
        col,
        transpose_B,
        col,
        all_W,
        &lwork);

    // allocate memery
    cuDoubleComplex* d_work;
    cudaMalloc((void**)&d_work, sizeof(cuDoubleComplex) * lwork);

    // compute eigenvalues and eigenvectors.
    cusolverDnZhegvd(
        cusolverH,
        CUSOLVER_EIG_TYPE_1, // itype = CUSOLVER_EIG_TYPE_1: A*x = (lambda)*B*x.
        CUSOLVER_EIG_MODE_VECTOR, // jobz = CUSOLVER_EIG_MODE_VECTOR : Compute eigenvalues and eigenvectors.
        CUBLAS_FILL_MODE_UPPER,
        row,
        A_eigenvectors,
        col,
        transpose_B,
        col,
        all_W,
        d_work,
        lwork,
        devInfo);

    cudaDeviceSynchronize();

    // get all eigenvalues and eigenvectors.
    cudaMemcpy(W, all_W, sizeof(double) * row, cudaMemcpyDeviceToDevice);
    matrixTranspose_op<double, psi::DEVICE_GPU>()(d, row, col, (std::complex<double>*)A_eigenvectors, V);

    int info_gpu;
    cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost);
    assert(0 == info_gpu);

    // free the buffer
    cudaFree(d_work);
    // free resources and destroy
    cudaFree(A_eigenvectors);
    cudaFree(all_W);
    cudaFree(devInfo);
    cusolverDnDestroy(cusolverH);
    destoryBLAShandle();
}

} // namespace hsolver