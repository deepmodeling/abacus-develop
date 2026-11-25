#include "diag_cusolver.cuh"

#include "source_base/module_device/device_check.h"

#include <cstddef>
#include <cstdio>
#include <cstdlib>

#include <cuda_runtime.h>

#ifdef ABACUS_DM_DIRECT_PILOT
#include <chrono>
#endif

namespace
{
constexpr cusolverEigType_t problem_type = CUSOLVER_EIG_TYPE_1;
constexpr cusolverEigMode_t job = CUSOLVER_EIG_MODE_VECTOR;
constexpr cublasFillMode_t triangle = CUBLAS_FILL_MODE_UPPER;

void solver_error(
    const char *mode,
    const int n,
    const int nbands,
    const char *error,
    const int code)
{
    std::fprintf(stderr, "cuSOLVER %s failed (n=%d, nbands=%d): %s (%d)\n", mode, n, nbands, error, code);
    std::exit(EXIT_FAILURE);
}

void check_status(
    const cusolverStatus_t status,
    const char *mode,
    const int n,
    const int nbands)
{
    if (status != CUSOLVER_STATUS_SUCCESS)
    {
        solver_error(mode, n, nbands, _cusolverGetErrorString(status), status);
    }
}

void check_info(
    const int info,
    const char *mode,
    const int n,
    const int nbands)
{
    if (info < 0)
    {
        solver_error(mode, n, nbands, "invalid argument to eigensolver, devInfo", info);
    }
    if (info > n)
    {
        solver_error(mode, n, nbands, "overlap matrix is not positive definite, devInfo", info);
    }
    if (info > 0)
    {
        solver_error(mode, n, nbands, "eigensolver did not converge, devInfo", info);
    }
}

#if CUDART_VERSION >= 10010
bool subset_unsupported(const cusolverStatus_t status)
{
    return status == CUSOLVER_STATUS_NOT_SUPPORTED || status == CUSOLVER_STATUS_ARCH_MISMATCH;
}
#endif

template <typename T>
struct DeviceBuffers
{
    T *h = nullptr;
    T *s = nullptr;
    double *eigenvalues = nullptr;
    int *info = nullptr;
    T *work = nullptr;

    ~DeviceBuffers()
    {
        CHECK_CUDA(cudaFree(work));
        CHECK_CUDA(cudaFree(info));
        CHECK_CUDA(cudaFree(eigenvalues));
        CHECK_CUDA(cudaFree(s));
        CHECK_CUDA(cudaFree(h));
    }

    void upload(const int n, const void *host_h, const void *host_s)
    {
        const std::size_t matrix_bytes = sizeof(T) * static_cast<std::size_t>(n) * n;
        const std::size_t eigenvalue_bytes = sizeof(double) * n;
        CHECK_CUDA(cudaMalloc(&h, matrix_bytes));
        CHECK_CUDA(cudaMalloc(&s, matrix_bytes));
        CHECK_CUDA(cudaMalloc(&eigenvalues, eigenvalue_bytes));
        CHECK_CUDA(cudaMalloc(&info, sizeof(int)));
        CHECK_CUDA(cudaMemcpy(h, host_h, matrix_bytes, cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(s, host_s, matrix_bytes, cudaMemcpyHostToDevice));
    }

    int read_info() const
    {
        CHECK_CUDA(cudaDeviceSynchronize());
        int host_info = 0;
        CHECK_CUDA(cudaMemcpy(&host_info, info, sizeof(host_info), cudaMemcpyDeviceToHost));
        return host_info;
    }

    void download(
        const int n,
        const int nbands,
        double *w,
        void *v) const
    {
        const std::size_t eigenvalue_bytes = sizeof(double) * nbands;
        const std::size_t eigenvector_bytes = sizeof(T) * static_cast<std::size_t>(n) * nbands;
        CHECK_CUDA(cudaMemcpy(w, eigenvalues, eigenvalue_bytes, cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(v, h, eigenvector_bytes, cudaMemcpyDeviceToHost));
    }
};
#ifdef ABACUS_DM_DIRECT_PILOT
void report_solver_memory(const int n, const std::size_t workspace_bytes)
{
    std::size_t available = 0;
    std::size_t total = 0;
    CHECK_CUDA(cudaMemGetInfo(&available, &total));
    const std::size_t matrix_bytes = sizeof(double) * static_cast<std::size_t>(n) * n;
    std::printf("DM_DIRECT_SOLVER_MEMORY n=%d matrix_bytes=%zu workspace_bytes=%zu free_at_peak=%zu total=%zu\n",
                n, matrix_bytes, workspace_bytes, available, total);
}

void report_real_solver(const char* algorithm, const int n, const int bands,
                        const std::chrono::steady_clock::time_point& started)
{
    const std::chrono::steady_clock::time_point finished = std::chrono::steady_clock::now();
    const std::chrono::steady_clock::duration duration = finished - started;
    const double seconds = std::chrono::duration<double>(duration).count();
    const std::size_t hs_bytes = 2 * sizeof(double) * static_cast<std::size_t>(n) * n;
    const std::size_t vector_bytes = sizeof(double) * static_cast<std::size_t>(n) * bands;
    const std::size_t eigenvalue_bytes = sizeof(double) * static_cast<std::size_t>(bands);
    std::printf("DM_DIRECT_SOLVER algorithm=%s n=%d host_columns=%d seconds=%.9e hs_h2d_calls=2 hs_h2d_bytes=%zu "
                "psi_d2h_calls=1 psi_d2h_bytes=%zu eigenvalue_d2h_calls=1 eigenvalue_d2h_bytes=%zu "
                "info_d2h_calls=1 info_d2h_bytes=4 info=0\n",
                algorithm, n, bands, seconds, hs_bytes, vector_bytes, eigenvalue_bytes);
}
#endif
} // namespace

Diag_Cusolver_gvd::Diag_Cusolver_gvd() { CHECK_CUSOLVER(cusolverDnCreate(&handle_)); }

Diag_Cusolver_gvd::~Diag_Cusolver_gvd()
{
#ifdef ABACUS_DM_DIRECT_PILOT
    release_eigenvectors();
#endif
    CHECK_CUSOLVER(cusolverDnDestroy(handle_));
}

#ifdef ABACUS_DM_DIRECT_PILOT
void Diag_Cusolver_gvd::release_eigenvectors()
{
    CHECK_CUDA(cudaFree(eigenvectors_));
    CHECK_CUDA(cudaFree(weighted_vectors_));
    eigenvectors_ = nullptr;
    weighted_vectors_ = nullptr;
    rows_ = 0;
    bands_ = 0;
}
#endif

bool Diag_Cusolver_gvd::use_subset(const int n, const int nbands)
{
    return CUDART_VERSION >= 10010 && n >= 4096 && nbands > 0 && nbands <= n / 4;
}

#if CUDART_VERSION >= 10010
bool Diag_Cusolver_gvd::solve_subset(
    const int n,
    const int nbands,
    const double *h,
    const double *s,
    double *w,
    double *v)
{
#ifdef ABACUS_DM_DIRECT_PILOT
    release_eigenvectors();
    const std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now();
#endif
    if (nbands <= 0 || nbands > n)
    {
        solver_error("Dsygvdx", n, nbands, "invalid eigenpair count", nbands);
    }
    DeviceBuffers<double> device;
    device.upload(n, h, s);
    const int first = 1;
    const double unused_bound = 0.0;
    int found = 0;
    int lwork = 0;
    auto status = cusolverDnDsygvdx_bufferSize(
        handle_, problem_type, job, CUSOLVER_EIG_RANGE_I, triangle, n, device.h, n, device.s, n, unused_bound,
        unused_bound, first, nbands, &found, device.eigenvalues, &lwork);
    if (subset_unsupported(status) || status == CUSOLVER_STATUS_ALLOC_FAILED)
    {
        return false;
    }
    check_status(status, "Dsygvdx workspace query", n, nbands);
    const std::size_t work_bytes = sizeof(double) * static_cast<std::size_t>(lwork);
    const auto alloc_status = cudaMalloc(&device.work, work_bytes);
    if (alloc_status == cudaErrorMemoryAllocation)
    {
        cudaGetLastError();
        return false;
    }
    CHECK_CUDA(alloc_status);
#ifdef ABACUS_DM_DIRECT_PILOT
    report_solver_memory(n, work_bytes);
#endif

    status = cusolverDnDsygvdx(
        handle_, problem_type, job, CUSOLVER_EIG_RANGE_I, triangle, n, device.h, n, device.s, n, unused_bound,
        unused_bound, first, nbands, &found, device.eigenvalues, device.work, lwork, device.info);
    if (subset_unsupported(status) || status == CUSOLVER_STATUS_ALLOC_FAILED)
    {
        return false;
    }
    check_status(status, "Dsygvdx", n, nbands);

    const int info = device.read_info();
    if (info < 0 || info > n)
    {
        check_info(info, "Dsygvdx", n, nbands);
    }
    if (info != 0 || found != nbands)
    {
        return false;
    }
    device.download(n, nbands, w, v);
#ifdef ABACUS_DM_DIRECT_PILOT
    eigenvectors_ = device.h;
    weighted_vectors_ = device.s;
    rows_ = n;
    bands_ = nbands;
    device.h = nullptr;
    device.s = nullptr;
    report_real_solver("Dsygvdx", n, nbands, started);
#endif
    return true;
}
#else
bool Diag_Cusolver_gvd::solve_subset(
    int,
    int,
    const double *,
    const double *,
    double *,
    double *)
{
#ifdef ABACUS_DM_DIRECT_PILOT
    release_eigenvectors();
#endif
    return false;
}
#endif

void Diag_Cusolver_gvd::solve_full(
    const int n,
    const int nbands,
    const double *h,
    const double *s,
    double *w,
    double *v)
{
#ifdef ABACUS_DM_DIRECT_PILOT
    release_eigenvectors();
    const std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now();
#endif
    DeviceBuffers<double> device;
    device.upload(n, h, s);
    int lwork = 0;
    auto status = cusolverDnDsygvd_bufferSize(
        handle_, problem_type, job, triangle, n, device.h, n, device.s, n, device.eigenvalues, &lwork);
    check_status(status, "Dsygvd workspace query", n, nbands);
    const std::size_t work_bytes = sizeof(double) * static_cast<std::size_t>(lwork);
    CHECK_CUDA(cudaMalloc(&device.work, work_bytes));
#ifdef ABACUS_DM_DIRECT_PILOT
    report_solver_memory(n, work_bytes);
#endif
    status = cusolverDnDsygvd(
        handle_, problem_type, job, triangle, n, device.h, n, device.s, n, device.eigenvalues, device.work, lwork,
        device.info);
    check_status(status, "Dsygvd", n, nbands);
    const int info = device.read_info();
    check_info(info, "Dsygvd", n, nbands);
    device.download(n, nbands, w, v);
#ifdef ABACUS_DM_DIRECT_PILOT
    eigenvectors_ = device.h;
    weighted_vectors_ = device.s;
    rows_ = n;
    bands_ = nbands;
    device.h = nullptr;
    device.s = nullptr;
    report_real_solver("Dsygvd", n, nbands, started);
#endif
}

void Diag_Cusolver_gvd::solve_full(
    const int n,
    const int nbands,
    const std::complex<double> *h,
    const std::complex<double> *s,
    double *w,
    std::complex<double> *v)
{
#ifdef ABACUS_DM_DIRECT_PILOT
    release_eigenvectors();
#endif
    static_assert(sizeof(std::complex<double>) == sizeof(cuDoubleComplex), "Unexpected complex representation");
    DeviceBuffers<cuDoubleComplex> device;
    device.upload(n, h, s);
    int lwork = 0;
    auto status = cusolverDnZhegvd_bufferSize(
        handle_, problem_type, job, triangle, n, device.h, n, device.s, n, device.eigenvalues, &lwork);
    check_status(status, "Zhegvd workspace query", n, nbands);
    const std::size_t work_bytes = sizeof(cuDoubleComplex) * static_cast<std::size_t>(lwork);
    CHECK_CUDA(cudaMalloc(&device.work, work_bytes));
    status = cusolverDnZhegvd(
        handle_, problem_type, job, triangle, n, device.h, n, device.s, n, device.eigenvalues, device.work, lwork,
        device.info);
    check_status(status, "Zhegvd", n, nbands);
    const int info = device.read_info();
    check_info(info, "Zhegvd", n, nbands);
    device.download(n, nbands, w, v);
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
    if (eigenvectors_ == nullptr || weighted_vectors_ == nullptr || bands <= 0 || bands > bands_
        || weights == nullptr || density == nullptr)
    {
        return false;
    }
    const std::size_t matrix_bytes = sizeof(double) * static_cast<std::size_t>(rows_) * rows_;
    const std::size_t weight_bytes = sizeof(double) * static_cast<std::size_t>(bands);
    double* device_density = nullptr;
    double* device_weights = nullptr;
    cublasHandle_t handle = nullptr;
    int weight_calls = 0;
    int density_calls = 0;
    std::size_t transferred_weights = 0;
    std::size_t transferred_density = 0;
    bool success = cudaMalloc(&device_density, matrix_bytes) == cudaSuccess;
    if (success) success = cudaMalloc(&device_weights, weight_bytes) == cudaSuccess;
    if (success) success = cublasCreate(&handle) == CUBLAS_STATUS_SUCCESS;
    if (success) success = cublasSetMathMode(handle, CUBLAS_DEFAULT_MATH) == CUBLAS_STATUS_SUCCESS;
    if (success)
    {
        ++weight_calls;
        success = cudaMemcpy(device_weights, weights, weight_bytes, cudaMemcpyHostToDevice) == cudaSuccess;
        if (success) transferred_weights = weight_bytes;
    }
    if (success)
    {
        const std::size_t count = static_cast<std::size_t>(rows_) * bands;
        const std::size_t block_count = (count + 255) / 256;
        const int blocks = static_cast<int>(block_count);
        // weighted_vectors_'s Cholesky factor is dead after diagonalization; reuse its storage.
        weight_retained_columns<<<blocks, 256>>>(rows_, bands, eigenvectors_, device_weights, weighted_vectors_);
        success = cudaGetLastError() == cudaSuccess;
    }
    if (success)
    {
        const double alpha = 1.0;
        const double beta = 0.0;
        success = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, rows_, rows_, bands,
                              &alpha, weighted_vectors_, rows_, eigenvectors_, rows_, &beta, device_density, rows_) == CUBLAS_STATUS_SUCCESS;
    }
    if (success)
    {
        ++density_calls;
        success = cudaMemcpy(density, device_density, matrix_bytes, cudaMemcpyDeviceToHost) == cudaSuccess;
        if (success) transferred_density = matrix_bytes;
    }
    if (success) success = cudaDeviceSynchronize() == cudaSuccess;
    if (device_density != nullptr) success = (cudaFree(device_density) == cudaSuccess) && success;
    if (device_weights != nullptr) success = (cudaFree(device_weights) == cudaSuccess) && success;
    if (handle != nullptr) success = (cublasDestroy(handle) == CUBLAS_STATUS_SUCCESS) && success;
    const std::chrono::steady_clock::time_point finished = std::chrono::steady_clock::now();
    const std::chrono::steady_clock::duration duration = finished - started;
    const double seconds = std::chrono::duration<double>(duration).count();
    std::printf("DM_DIRECT_DEVICE n=%d bands=%d seconds=%.9e success=%d psi_h2d_calls=0 psi_h2d_bytes=0 "
                "weights_h2d_calls=%d weights_h2d_bytes=%zu dm_d2h_calls=%d dm_d2h_bytes=%zu\n",
                rows_, bands, seconds, success, weight_calls, transferred_weights, density_calls, transferred_density);
    return success;
}
#endif
