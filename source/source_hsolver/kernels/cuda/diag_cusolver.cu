#include "diag_cusolver.cuh"

#include "source_base/module_device/device_check.h"

#include <cstddef>
#include <cstdio>
#include <cstdlib>

#include <cuda_runtime.h>

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

Diag_Cusolver_gvd::~Diag_Cusolver_gvd() { CHECK_CUSOLVER(cusolverDnDestroy(handle_)); }

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
    const std::size_t work_bytes = sizeof(double) * static_cast<std::size_t>(lwork);
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
    DeviceBuffers<double> device;
    device.upload(n, h, s);
    int lwork = 0;
    cusolverStatus_t status = cusolverDnDsygvd_bufferSize(
        handle_, problem_type, job, triangle, n, device.h, n, device.s, n, device.eigenvalues, &lwork);
    check_status(status, "Dsygvd workspace query", n, nbands);
    const std::size_t work_bytes = sizeof(double) * static_cast<std::size_t>(lwork);
    CHECK_CUDA(cudaMalloc(&device.work, work_bytes));
    status = cusolverDnDsygvd(
        handle_, problem_type, job, triangle, n, device.h, n, device.s, n, device.eigenvalues, device.work, lwork,
        device.info);
    check_status(status, "Dsygvd", n, nbands);
    const int info = device.read_info();
    check_info(info, "Dsygvd", n, nbands);
    device.download(n, nbands, w, v);
}

void Diag_Cusolver_gvd::solve_full(
    const int n,
    const int nbands,
    const std::complex<double> *h,
    const std::complex<double> *s,
    double *w,
    std::complex<double> *v)
{
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
