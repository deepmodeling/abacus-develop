#include "diag_cusolver.cuh"

#include "source_base/module_device/device_check.h"

#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <cstdlib>

#include <cuda_runtime.h>

#include <sstream>

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
} // namespace

Diag_Cusolver_gvd::Diag_Cusolver_gvd() { CHECK_CUSOLVER(cusolverDnCreate(&handle_)); }

Diag_Cusolver_gvd::~Diag_Cusolver_gvd()
{
    release_eigenvectors();
    CHECK_CUSOLVER(cusolverDnDestroy(handle_));
}

void Diag_Cusolver_gvd::release_eigenvectors()
{
    if (eigenvectors_ != nullptr) CHECK_CUDA(cudaFree(eigenvectors_));
    if (weighted_vectors_ != nullptr) CHECK_CUDA(cudaFree(weighted_vectors_));
    if (density_workspace_ != nullptr) CHECK_CUDA(cudaFree(density_workspace_));
    if (occupations_ != nullptr) CHECK_CUDA(cudaFree(occupations_));
    eigenvectors_ = nullptr;
    weighted_vectors_ = nullptr;
    density_workspace_ = nullptr;
    occupations_ = nullptr;
    rows_ = 0;
    bands_ = 0;
}

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
    double *v,
    const bool retain_vectors)
{
    release_eigenvectors();
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
    cusolverStatus_t status = cusolverDnDsygvdx_bufferSize(
        handle_, problem_type, job, CUSOLVER_EIG_RANGE_I, triangle, n, device.h, n, device.s, n, unused_bound,
        unused_bound, first, nbands, &found, device.eigenvalues, &lwork);
    if (subset_unsupported(status) || status == CUSOLVER_STATUS_ALLOC_FAILED)
    {
        return false;
    }
    check_status(status, "Dsygvdx workspace query", n, nbands);
    const std::size_t density_elements = retain_vectors ? static_cast<std::size_t>(n) * n : 0;
    const std::size_t solver_elements = static_cast<std::size_t>(lwork);
    const std::size_t work_elements = std::max(solver_elements, density_elements);
    const std::size_t work_bytes = sizeof(double) * work_elements;
    const cudaError_t alloc_status = cudaMalloc(&device.work, work_bytes);
    if (alloc_status == cudaErrorMemoryAllocation)
    {
        cudaGetLastError();
        return false;
    }
    CHECK_CUDA(alloc_status);

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
    if (retain_vectors)
    {
        eigenvectors_ = device.h;
        weighted_vectors_ = device.s;
        density_workspace_ = device.work;
        occupations_ = device.eigenvalues;
        rows_ = n;
        bands_ = nbands;
        device.h = nullptr;
        device.s = nullptr;
        device.work = nullptr;
        device.eigenvalues = nullptr;
    }
    return true;
}
#else
bool Diag_Cusolver_gvd::solve_subset(
    int,
    int,
    const double *,
    const double *,
    double *,
    double *,
    bool)
{
    release_eigenvectors();
    return false;
}
#endif

void Diag_Cusolver_gvd::solve_full(
    const int n,
    const int nbands,
    const double *h,
    const double *s,
    double *w,
    double *v,
    const bool retain_vectors)
{
    release_eigenvectors();
    DeviceBuffers<double> device;
    device.upload(n, h, s);
    int lwork = 0;
    cusolverStatus_t status = cusolverDnDsygvd_bufferSize(
        handle_, problem_type, job, triangle, n, device.h, n, device.s, n, device.eigenvalues, &lwork);
    check_status(status, "Dsygvd workspace query", n, nbands);
    const std::size_t density_elements = retain_vectors ? static_cast<std::size_t>(n) * n : 0;
    const std::size_t solver_elements = static_cast<std::size_t>(lwork);
    const std::size_t work_elements = std::max(solver_elements, density_elements);
    const std::size_t work_bytes = sizeof(double) * work_elements;
    CHECK_CUDA(cudaMalloc(&device.work, work_bytes));
    status = cusolverDnDsygvd(
        handle_, problem_type, job, triangle, n, device.h, n, device.s, n, device.eigenvalues, device.work, lwork,
        device.info);
    check_status(status, "Dsygvd", n, nbands);
    const int info = device.read_info();
    check_info(info, "Dsygvd", n, nbands);
    device.download(n, nbands, w, v);
    if (retain_vectors)
    {
        eigenvectors_ = device.h;
        weighted_vectors_ = device.s;
        density_workspace_ = device.work;
        occupations_ = device.eigenvalues;
        rows_ = n;
        bands_ = nbands;
        device.h = nullptr;
        device.s = nullptr;
        device.work = nullptr;
        device.eigenvalues = nullptr;
    }
}

void Diag_Cusolver_gvd::solve_full(
    const int n,
    const int nbands,
    const std::complex<double> *h,
    const std::complex<double> *s,
    double *w,
    std::complex<double> *v)
{
    release_eigenvectors();
    static_assert(sizeof(std::complex<double>) == sizeof(cuDoubleComplex), "Unexpected complex representation");
    DeviceBuffers<cuDoubleComplex> device;
    device.upload(n, h, s);
    int lwork = 0;
    cusolverStatus_t status = cusolverDnZhegvd_bufferSize(
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

namespace
{
void density_error(std::string& error, const char* stage, const char* detail, const int code,
                   const int rows, const int bands, const std::size_t bytes)
{
    // Cleanup must not replace the original cause of a failed operation.
    if (!error.empty()) return;
    std::ostringstream message;
    message << "Direct GPU density matrix failed during " << stage << ": " << detail
            << " (code=" << code << ", n=" << rows << ", bands=" << bands;
    if (bytes != 0) message << ", requested_bytes=" << bytes;
    message << "). For CPU execution, rerun with device cpu and a CPU ks_solver such as scalapack_gvx.";
    error = message.str();
}

bool density_cuda_status(const cudaError_t status, const char* stage, const int rows,
                         const int bands, const std::size_t bytes, std::string& error)
{
    if (status == cudaSuccess) return true;
    const char* detail = cudaGetErrorString(status);
    const int code = static_cast<int>(status);
    density_error(error, stage, detail, code, rows, bands, bytes);
    return false;
}

bool density_blas_status(const cublasStatus_t status, const char* stage, const int rows,
                         const int bands, std::string& error)
{
    if (status == CUBLAS_STATUS_SUCCESS) return true;
    const char* detail = _cublasGetErrorString(status);
    const int code = static_cast<int>(status);
    density_error(error, stage, detail, code, rows, bands, 0);
    return false;
}

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

bool Diag_Cusolver_gvd::density_from_retained(int bands, const double* weights, double* density, std::string& error)
{
    error.clear();
    if (eigenvectors_ == nullptr || weighted_vectors_ == nullptr || density_workspace_ == nullptr
        || occupations_ == nullptr || bands <= 0 || bands > bands_
        || weights == nullptr || density == nullptr)
    {
        density_error(error, "validate retained eigenvectors", "invalid eigenvector state or DM arguments",
                      -1, rows_, bands, 0);
        return false;
    }
    const std::size_t matrix_elements = static_cast<std::size_t>(rows_) * rows_;
    const std::size_t matrix_bytes = sizeof(double) * matrix_elements;
    const std::size_t weight_bytes = sizeof(double) * static_cast<std::size_t>(bands);
    // Reuse the completed solve's workspace and eigenvalue storage without new device arrays.
    double* device_density = density_workspace_;
    double* device_weights = occupations_;
    cudaError_t cuda_status = cudaSuccess;
    cublasHandle_t handle = nullptr;
    const cublasStatus_t create_status = cublasCreate(&handle);
    bool success = density_blas_status(create_status, "create cuBLAS handle", rows_, bands, error);
    if (success)
    {
        const cublasStatus_t status = cublasSetMathMode(handle, CUBLAS_DEFAULT_MATH);
        success = density_blas_status(status, "set cuBLAS math mode", rows_, bands, error);
    }
    if (success)
    {
        cuda_status = cudaMemcpy(device_weights, weights, weight_bytes, cudaMemcpyHostToDevice);
        success = density_cuda_status(cuda_status, "upload occupations", rows_, bands, weight_bytes, error);
    }
    if (success)
    {
        const std::size_t count = static_cast<std::size_t>(rows_) * bands;
        const std::size_t block_count = (count + 255) / 256;
        const int blocks = static_cast<int>(block_count);
        // weighted_vectors_'s Cholesky factor is dead after diagonalization; reuse its storage.
        weight_retained_columns<<<blocks, 256>>>(rows_, bands, eigenvectors_, device_weights, weighted_vectors_);
        cuda_status = cudaGetLastError();
        success = density_cuda_status(cuda_status, "weight retained eigenvectors", rows_, bands, 0, error);
    }
    if (success)
    {
        const double alpha = 1.0;
        const double beta = 0.0;
        const cublasStatus_t status = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, rows_, rows_, bands,
                                                 &alpha, weighted_vectors_, rows_, eigenvectors_, rows_,
                                                 &beta, device_density, rows_);
        success = density_blas_status(status, "form density matrix (cublasDgemm)", rows_, bands, error);
    }
    if (success)
    {
        cuda_status = cudaMemcpy(density, device_density, matrix_bytes, cudaMemcpyDeviceToHost);
        success = density_cuda_status(cuda_status, "download density matrix", rows_, bands, matrix_bytes, error);
    }
    if (success)
    {
        cuda_status = cudaDeviceSynchronize();
        success = density_cuda_status(cuda_status, "synchronize density matrix", rows_, bands, 0, error);
    }
    if (handle != nullptr)
    {
        const cublasStatus_t status = cublasDestroy(handle);
        const bool released = density_blas_status(status, "destroy cuBLAS handle", rows_, bands, error);
        success = released && success;
    }
    return success;
}
