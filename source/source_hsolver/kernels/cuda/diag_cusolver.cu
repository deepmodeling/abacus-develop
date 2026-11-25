#include <assert.h>
#include "diag_cusolver.cuh"
#include "source_base/module_device/device_check.h"

#ifdef ABACUS_DM_DIRECT_PILOT
#include <chrono>
#include <cstdio>
#include <cstdlib>
#endif

Diag_Cusolver_gvd::Diag_Cusolver_gvd(){
// step 1: create cusolver/cublas handle
    cusolverH = NULL;
    CHECK_CUSOLVER( cusolverDnCreate(&cusolverH) );

    itype = CUSOLVER_EIG_TYPE_1; // A*x = (lambda)*B*x
    jobz = CUSOLVER_EIG_MODE_VECTOR; // compute eigenvalues and eigenvectors.
    // LCAO supplies the authoritative upper triangle of the Hermitian H/S matrices.
    uplo = CUBLAS_FILL_MODE_UPPER;

    d_A = NULL;
    d_B = NULL;
    d_work = NULL;

    d_A2 = NULL;
    d_B2 = NULL;
    d_work2 = NULL;

    d_W = NULL;
    devInfo = NULL;

    lwork = 0;
    info_gpu = 0;
    is_init = 0;
}

void Diag_Cusolver_gvd::finalize(){
// free resources and destroy
    if (d_A      ) {CHECK_CUDA( cudaFree(d_A) );  d_A  = NULL;}
    if (d_B      ) {CHECK_CUDA( cudaFree(d_B) );  d_B  = NULL;}
    if (d_A2     ) {CHECK_CUDA( cudaFree(d_A2) ); d_A2 = NULL;}
    if (d_B2     ) {CHECK_CUDA( cudaFree(d_B2) ); d_B2 = NULL;}
    if (d_W      ) {CHECK_CUDA( cudaFree(d_W) );  d_W  = NULL;}
    if (devInfo  ) {CHECK_CUDA( cudaFree(devInfo) );   devInfo = NULL;}
}

Diag_Cusolver_gvd::~Diag_Cusolver_gvd(){
    finalize();
    if (cusolverH) {CHECK_CUSOLVER( cusolverDnDestroy(cusolverH) );    cusolverH = NULL;}
    //CHECK_CUDA( cudaDeviceReset() );
}


void Diag_Cusolver_gvd::init_double(int N){
// step 2: Malloc A and B on device
    m = lda = N;
    CHECK_CUDA( cudaMalloc ((void**)&d_A, sizeof(double) * lda * m) );
    CHECK_CUDA( cudaMalloc ((void**)&d_B, sizeof(double) * lda * m) );
    CHECK_CUDA( cudaMalloc ((void**)&d_W, sizeof(double) * m) );
    CHECK_CUDA( cudaMalloc ((void**)&devInfo, sizeof(int)) );
}

void Diag_Cusolver_gvd::init_complex(int N){
// step 2: Malloc A and B on device
    m = lda = N;
    CHECK_CUDA( cudaMalloc ((void**)&d_A2, sizeof(cuDoubleComplex) * lda * m) );
    CHECK_CUDA( cudaMalloc ((void**)&d_B2, sizeof(cuDoubleComplex) * lda * m) );
    CHECK_CUDA( cudaMalloc ((void**)&d_W, sizeof(double) * m) );
    CHECK_CUDA( cudaMalloc ((void**)&devInfo, sizeof(int)) );
}
        
void Diag_Cusolver_gvd::Dngvd_double(int N, int M, double *A, double *B, double *W, double *V){
#ifdef ABACUS_DM_DIRECT_PILOT
    assert(N == M);
    Dngvd_double_retained(N, A, B, W, V, N);
}

void Diag_Cusolver_gvd::Dngvd_double_retained(int N, double* A, double* B,
                                             double* W, double* V, int host_columns)
{
    const int M = N;
    const std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now();
#endif

    // copy A, B to the GPU
        assert(N == M);
        if (M != m) {
            this->finalize();
            this->init_double(M);
        }
        CHECK_CUDA( cudaMemcpy(d_A, A, sizeof(double) * lda * m, cudaMemcpyHostToDevice) );
        CHECK_CUDA( cudaMemcpy(d_B, B, sizeof(double) * lda * m, cudaMemcpyHostToDevice) );

    // Query working space of sygvd
    // The helper functions below can calculate the sizes needed for pre-allocated buffer.
    // The S and D data types are real valued single and double precision, respectively.
    // The C and Z data types are complex valued single and double precision, respectively.
        CHECK_CUSOLVER(cusolverDnDsygvd_bufferSize(
            cusolverH,
            itype,
            jobz,
            uplo,
            m,
            d_A,
            lda,
            d_B,
            lda,
            d_W,
            &lwork
        ));
        CHECK_CUDA( cudaMalloc((void**)&d_work, sizeof(double)*lwork) );
#ifdef ABACUS_DM_DIRECT_PILOT
        std::size_t available = 0;
        std::size_t total = 0;
        CHECK_CUDA(cudaMemGetInfo(&available, &total));
        const std::size_t matrix_bytes = sizeof(double) * static_cast<std::size_t>(N) * N;
        const std::size_t workspace_bytes = sizeof(double) * static_cast<std::size_t>(lwork);
        std::printf("DM_DIRECT_SOLVER_MEMORY n=%d matrix_bytes=%zu workspace_bytes=%zu free_at_peak=%zu total=%zu\n",
                    N, matrix_bytes, workspace_bytes, available, total);
#endif

    // compute spectrum of (A,B)
        CHECK_CUSOLVER(cusolverDnDsygvd(
            cusolverH,
            itype,
            jobz,
            uplo,
            m,
            d_A,
            lda,
            d_B,
            lda,
            d_W,
            d_work,
            lwork,
            devInfo
        ));
        CHECK_CUDA( cudaDeviceSynchronize() );

    // copy (W, V) to the cpu root
        CHECK_CUDA( cudaMemcpy(W, d_W, sizeof(double)*m, cudaMemcpyDeviceToHost) );
#ifdef ABACUS_DM_DIRECT_PILOT
        const std::size_t vector_bytes = sizeof(double) * static_cast<std::size_t>(lda) * host_columns;
        CHECK_CUDA(cudaMemcpy(V, d_A, vector_bytes, cudaMemcpyDeviceToHost));
#else
        CHECK_CUDA( cudaMemcpy(V, d_A, sizeof(double)*lda*m, cudaMemcpyDeviceToHost) );
#endif
        CHECK_CUDA( cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost) );
        assert(0 == info_gpu);
    // free the buffer
        if (d_work ) CHECK_CUDA( cudaFree(d_work) );
#ifdef ABACUS_DM_DIRECT_PILOT
        d_work = nullptr;
        const std::chrono::steady_clock::time_point finished = std::chrono::steady_clock::now();
        const double seconds = std::chrono::duration<double>(finished - started).count();
        const std::size_t hs_bytes = 2 * matrix_bytes;
        const std::size_t eigenvalue_bytes = sizeof(double) * static_cast<std::size_t>(m);
        std::printf("DM_DIRECT_SOLVER n=%d host_columns=%d seconds=%.9e hs_h2d_calls=2 hs_h2d_bytes=%zu psi_d2h_calls=1 psi_d2h_bytes=%zu eigenvalue_d2h_calls=1 eigenvalue_d2h_bytes=%zu info_d2h_calls=1 info_d2h_bytes=4 info=%d\n",
                    N, host_columns, seconds, hs_bytes, vector_bytes, eigenvalue_bytes, info_gpu);
#endif

}

#ifdef ABACUS_DM_DIRECT_PILOT
namespace
{
__global__ void weight_retained_columns(const int rows, const int columns,
                                         const double* eigenvectors,
                                         const double* weights, double* weighted)
{
    const std::size_t count = static_cast<std::size_t>(rows) * columns;
    const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index < count)
    {
        const std::size_t band = index / rows;
        weighted[index] = eigenvectors[index] * weights[band];
    }
}
}

bool Diag_Cusolver_gvd::density_from_retained(int bands, const double* weights, double* density)
{
    const std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now();
    if (d_A == nullptr || d_B == nullptr || bands <= 0 || bands > m || info_gpu != 0)
    {
        return false;
    }
    const std::size_t matrix_bytes = sizeof(double) * static_cast<std::size_t>(m) * m;
    const std::size_t weight_bytes = sizeof(double) * static_cast<std::size_t>(bands);
    double* device_density = nullptr;
    double* device_weights = nullptr;
    cublasHandle_t handle = nullptr;
    bool success = cudaMalloc(&device_density, matrix_bytes) == cudaSuccess;
    if (success) success = cudaMalloc(&device_weights, weight_bytes) == cudaSuccess;
    if (success) success = cublasCreate(&handle) == CUBLAS_STATUS_SUCCESS;
    if (success) success = cublasSetMathMode(handle, CUBLAS_DEFAULT_MATH) == CUBLAS_STATUS_SUCCESS;
    if (success) success = cudaMemcpy(device_weights, weights, weight_bytes, cudaMemcpyHostToDevice) == cudaSuccess;
    if (success)
    {
        const std::size_t count = static_cast<std::size_t>(m) * bands;
        const int blocks = static_cast<int>((count + 255) / 256);
        // d_B's Cholesky factor is dead after diagonalization; reuse its storage.
        weight_retained_columns<<<blocks, 256>>>(m, bands, d_A, device_weights, d_B);
        success = cudaGetLastError() == cudaSuccess;
    }
    if (success)
    {
        const double alpha = 1.0;
        const double beta = 0.0;
        success = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, m, m, bands,
                              &alpha, d_B, lda, d_A, lda, &beta, device_density, lda) == CUBLAS_STATUS_SUCCESS;
    }
    if (success) success = cudaMemcpy(density, device_density, matrix_bytes, cudaMemcpyDeviceToHost) == cudaSuccess;
    if (success) success = cudaDeviceSynchronize() == cudaSuccess;
    if (device_density != nullptr) success = (cudaFree(device_density) == cudaSuccess) && success;
    if (device_weights != nullptr) success = (cudaFree(device_weights) == cudaSuccess) && success;
    if (handle != nullptr) success = (cublasDestroy(handle) == CUBLAS_STATUS_SUCCESS) && success;
    const std::chrono::steady_clock::time_point finished = std::chrono::steady_clock::now();
    const double seconds = std::chrono::duration<double>(finished - started).count();
    std::printf("DM_DIRECT_DEVICE n=%d bands=%d seconds=%.9e success=%d psi_h2d_calls=0 psi_h2d_bytes=0 weights_h2d_calls=1 weights_h2d_bytes=%zu dm_d2h_calls=1 dm_d2h_bytes=%zu\n",
                m, bands, seconds, success, weight_bytes, matrix_bytes);
    return success;
}
#endif


void Diag_Cusolver_gvd::Dngvd_complex(int N, int M, std::complex<double> *A, std::complex<double> *B, double *W, std::complex<double> *V){

    // copy A, B to the GPU
        assert(N == M);
        if (M != m) {
            this->finalize();
            this->init_complex(M);
        }
        CHECK_CUDA( cudaMemcpy(d_A2, A, sizeof(cuDoubleComplex) * lda * m, cudaMemcpyHostToDevice) );
        CHECK_CUDA( cudaMemcpy(d_B2, B, sizeof(cuDoubleComplex) * lda * m, cudaMemcpyHostToDevice) );

    // Query working space of Zhegvd
    // The helper functions below can calculate the sizes needed for pre-allocated buffer.
    // The S and D data types are real valued single and double precision, respectively.
    // The C and Z data types are complex valued single and double precision, respectively.
        CHECK_CUSOLVER(
            cusolverDnZhegvd_bufferSize(
                cusolverH,
                itype,
                jobz,
                uplo,
                m,
                d_A2,
                lda,
                d_B2,
                lda,
                d_W,
                &lwork)
        );
        CHECK_CUDA( cudaMalloc((void**)&d_work2, sizeof(cuDoubleComplex)*lwork) );

    // compute spectrum of (A,B)
        CHECK_CUSOLVER(
            cusolverDnZhegvd(
                cusolverH,
                itype,
                jobz,
                uplo,
                m,
                d_A2,
                lda,
                d_B2,
                lda,
                d_W,
                d_work2,
                lwork,
                devInfo)
        );
        CHECK_CUDA( cudaDeviceSynchronize() );

    // copy (W, V) to the cpu root
        CHECK_CUDA( cudaMemcpy(W, d_W, sizeof(double)*m, cudaMemcpyDeviceToHost) );
        CHECK_CUDA( cudaMemcpy(V, d_A2, sizeof(std::complex<double>)*lda*m, cudaMemcpyDeviceToHost) );
        CHECK_CUDA( cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost) );
        assert(0 == info_gpu);

    // free the buffer
        if (d_work2 ) CHECK_CUDA( cudaFree(d_work2) );
}
