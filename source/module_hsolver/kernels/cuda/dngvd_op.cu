#include "module_hsolver/kernels/dngvd_op.h"
#include "module_hsolver/kernels/fixed_overlap.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "helper_cuda.h"

#include <base/macros/macros.h>

#include <cusolverDn.h>
#include <algorithm>
#include <stdexcept>
#include <string>

namespace hsolver
{

// The workspace shares the single-threaded solver handle lifetime. Retaining
// it avoids a device-wide cudaFree fence while other streams prefetch data.
static struct
{
    cusolverDnHandle_t handle = nullptr;
    double2* zhegvd_work = nullptr;
    int* zhegvd_info = nullptr;
    int zhegvd_capacity = 0;
} solver;

void createGpuSolverHandle()
{
    if (solver.handle == nullptr)
    {
        cusolverErrcheck(cusolverDnCreate(&solver.handle));
    }
}

void destroyGpuSolverHandle()
{
    if (solver.handle != nullptr)
    {
        cudaErrcheck(cudaFree(solver.zhegvd_work));
        cudaErrcheck(cudaFree(solver.zhegvd_info));
        cusolverErrcheck(cusolverDnDestroy(solver.handle));
        solver = {};
    }
}

namespace
{

__global__ void expand_hermitian_upper(const double2* input, double2* output, int size)
{
    const size_t index = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index >= size_t(size) * size) { return; }
    const int row = index % size, column = index / size;
    const auto value = row <= column ? input[index] : input[size_t(row) * size + column];
    output[index] = make_double2(value.x, row == column ? 0.0 : row < column ? value.y : -value.y);
}

} // namespace

FixedOverlapSubspace::FixedOverlapSubspace(int points, int bands, int projector_rows)
    : points_(points), bands_(bands), rows_(projector_rows)
{
    assert(points > 0 && bands > 0 && projector_rows > 0 && solver.handle != nullptr);
    const size_t matrix = size_t(bands) * bands;
    const size_t projections = size_t(bands) * projector_rows;
    const size_t scalar_storage = (size_t(bands) * sizeof(double) + sizeof(int) + sizeof(Complex) - 1)
                                  / sizeof(Complex);
    const size_t elements = size_t(points) * (matrix + projections) + 2 * matrix + scalar_storage;
    auto status = cudaMalloc(&storage_, elements * sizeof(Complex));
    if (status == cudaErrorMemoryAllocation)
    {
        cudaGetLastError();
        storage_ = nullptr;
        return;
    }
    cudaErrcheck(status);
    h_ = storage_;
    r_ = h_ + size_t(points) * matrix;
    cholesky_ = r_ + size_t(points) * projections;
    vectors_ = cholesky_ + matrix;
    eigenvalues_ = reinterpret_cast<double*>(vectors_ + matrix);
    info_ = reinterpret_cast<int*>(eigenvalues_ + bands);
    cusolverErrcheck(cusolverDnZpotrf_bufferSize(solver.handle, CUBLAS_FILL_MODE_UPPER,
        bands, reinterpret_cast<double2*>(cholesky_), bands, &factor_work_));
    cusolverErrcheck(cusolverDnZheevd_bufferSize(solver.handle, CUSOLVER_EIG_MODE_VECTOR,
        CUBLAS_FILL_MODE_UPPER, bands, reinterpret_cast<double2*>(vectors_), bands,
        eigenvalues_, &eigen_work_));
    const int work = std::max(factor_work_, eigen_work_);
    if (work <= solver.zhegvd_capacity) { return; }

    // The optional cache cannot displace the established solver workspace.
    // Reserve its larger workspace before releasing the existing allocation.
    double2* replacement = nullptr;
    status = cudaMalloc(&replacement, size_t(work) * sizeof(double2));
    if (status == cudaErrorMemoryAllocation)
    {
        cudaGetLastError();
        cudaErrcheck(cudaFree(storage_));
        storage_ = nullptr;
        return;
    }
    cudaErrcheck(status);
    cudaErrcheck(cudaFree(solver.zhegvd_work));
    solver.zhegvd_work = replacement;
    solver.zhegvd_capacity = work;
}

FixedOverlapSubspace::~FixedOverlapSubspace()
{
    cudaErrcheck(cudaFree(storage_));
}

const FixedOverlapSubspace::Complex* FixedOverlapSubspace::hamiltonian(int point) const
{
    assert(available() && point >= 0 && point < points_);
    return h_ + size_t(point) * bands_ * bands_;
}

const FixedOverlapSubspace::Complex* FixedOverlapSubspace::projectors(int point) const
{
    assert(available() && point >= 0 && point < points_);
    return r_ + size_t(point) * bands_ * rows_;
}

void FixedOverlapSubspace::require_success(const char* operation) const
{
    int info = 0;
    cudaErrcheck(cudaMemcpy(&info, info_, sizeof(int), cudaMemcpyDeviceToHost));
    if (info != 0)
    {
        throw std::runtime_error(std::string("Fixed-overlap ") + operation
            + " failed for dimension " + std::to_string(bands_) + ", cuSolver info=" + std::to_string(info));
    }
}

void FixedOverlapSubspace::prepare(int point, const Complex* h, const Complex* s, const Complex* projectors)
{
    assert(available() && point >= 0 && point < points_);
    const size_t matrix = size_t(bands_) * bands_;
    auto* transformed_h = h_ + size_t(point) * matrix;
    auto* transformed_r = r_ + size_t(point) * bands_ * rows_;
    cudaErrcheck(cudaMemcpy(cholesky_, s, matrix * sizeof(Complex), cudaMemcpyDeviceToDevice));
    cusolverErrcheck(cusolverDnZpotrf(solver.handle, CUBLAS_FILL_MODE_UPPER, bands_,
        reinterpret_cast<double2*>(cholesky_), bands_, solver.zhegvd_work, factor_work_, info_));
    require_success("Cholesky factorization");
    expand_hermitian_upper<<<(matrix + 255) / 256, 256>>>(
        reinterpret_cast<const double2*>(h), reinterpret_cast<double2*>(transformed_h), bands_);
    cudaErrcheck(cudaGetLastError());
    cudaErrcheck(cudaMemcpy(transformed_r, projectors, size_t(rows_) * bands_ * sizeof(Complex),
                           cudaMemcpyDeviceToDevice));
    whiten(transformed_h, transformed_r);
}

void FixedOverlapSubspace::diagonalize(const Complex* h, const Complex* projectors,
                                      Complex* rotated, double* eigenvalues)
{
    assert(available());
    cudaErrcheck(cudaMemcpy(vectors_, h, size_t(bands_) * bands_ * sizeof(Complex), cudaMemcpyDeviceToDevice));
    cusolverErrcheck(cusolverDnZheevd(solver.handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER,
        bands_, reinterpret_cast<double2*>(vectors_), bands_, eigenvalues_,
        solver.zhegvd_work, eigen_work_, info_));
    require_success("Hermitian diagonalization");
    cudaErrcheck(cudaMemcpy(eigenvalues, eigenvalues_, bands_ * sizeof(double), cudaMemcpyDeviceToHost));
    const Complex one(1, 0), zero(0, 0);
    gemm_op<Complex, base_device::DEVICE_GPU>()(nullptr, 'N', 'N', rows_, bands_, bands_,
        &one, projectors, rows_, vectors_, bands_, &zero, rotated, rows_);
}

static inline
void xhegvd_wrapper(
    const cublasFillMode_t& uplo,
    const int& n,
    double* A, const int& lda,
    double* B, const int& ldb,
    double* W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    int* devInfo = nullptr;
    int lwork = 0, info_gpu = 0;
    double* work = nullptr;
    cudaErrcheck(cudaMalloc((void**)&devInfo, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnDsygvd_bufferSize(solver.handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
        A, lda, B, ldb, W, &lwork));
    // allocate memery
    cudaErrcheck(cudaMalloc((void**)&work, sizeof(double) * lwork));

    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnDsygvd(solver.handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
        A, lda, B, ldb, W, work, lwork, devInfo));

    cudaErrcheck(cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost));
    assert(0 == info_gpu);
    // free the buffer
    cudaErrcheck(cudaFree(work));
    cudaErrcheck(cudaFree(devInfo));
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
    cudaErrcheck(cudaMalloc((void**)&devInfo, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnChegvd_bufferSize(solver.handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                                 reinterpret_cast<const float2 *>(A), lda,
                                                 reinterpret_cast<const float2 *>(B), ldb, W, &lwork));
    // allocate memery
    cudaErrcheck(cudaMalloc((void**)&work, sizeof(float2) * lwork));

    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnChegvd(solver.handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                      reinterpret_cast<float2 *>(A), lda, reinterpret_cast<float2 *>(B), ldb, W, work, lwork, devInfo));

    cudaErrcheck(cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost));
    assert(0 == info_gpu);
    // free the buffer
    cudaErrcheck(cudaFree(work));
    cudaErrcheck(cudaFree(devInfo));
}

static inline
void xhegvd_wrapper (
        const cublasFillMode_t& uplo,
        const int& n,
        std::complex<double> * A, const int& lda,
        std::complex<double> * B, const int& ldb,
        double * W)
{
    int lwork = 0, info_gpu = 0;
    if (solver.zhegvd_info == nullptr)
    {
        cudaErrcheck(cudaMalloc(reinterpret_cast<void**>(&solver.zhegvd_info), sizeof(int)));
    }

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnZhegvd_bufferSize(solver.handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                                 reinterpret_cast<const double2 *>(A), lda,
                                                 reinterpret_cast<const double2 *>(B), ldb, W, &lwork));
    if (lwork > solver.zhegvd_capacity)
    {
        cudaErrcheck(cudaFree(solver.zhegvd_work));
        cudaErrcheck(cudaMalloc(reinterpret_cast<void**>(&solver.zhegvd_work), sizeof(double2) * lwork));
        solver.zhegvd_capacity = lwork;
    }

    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnZhegvd(solver.handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                      reinterpret_cast<double2 *>(A), lda, reinterpret_cast<double2 *>(B), ldb, W, solver.zhegvd_work, lwork, solver.zhegvd_info));

    cudaErrcheck(cudaMemcpy(&info_gpu, solver.zhegvd_info, sizeof(int), cudaMemcpyDeviceToHost));
    assert(0 == info_gpu);
}

static inline
void xheevd_wrapper(
    const cublasFillMode_t& uplo,
    const int& n,
    double* A, const int& lda,
    double* W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    int* devInfo = nullptr;
    int lwork = 0, info_gpu = 0;
    double* work = nullptr;
    cudaErrcheck(cudaMalloc((void**)&devInfo, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnDsyevd_bufferSize(solver.handle, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
        A, lda, W, &lwork));
    // allocate memery
    cudaErrcheck(cudaMalloc((void**)&work, sizeof(double) * lwork));
    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnDsyevd(solver.handle, CUSOLVER_EIG_MODE_VECTOR, uplo, n, A, lda, W, work, lwork, devInfo));

    cudaErrcheck(cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost));
    assert(0 == info_gpu);
    cudaErrcheck(cudaFree(work));
    cudaErrcheck(cudaFree(devInfo));
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
    cudaErrcheck(cudaMalloc((void**)&devInfo, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnCheevd_bufferSize(solver.handle, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                                 reinterpret_cast<const float2 *>(A), lda, W, &lwork));
    // allocate memery
    cudaErrcheck(cudaMalloc((void**)&work, sizeof(float2) * lwork));
    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnCheevd(solver.handle, CUSOLVER_EIG_MODE_VECTOR, uplo, n, reinterpret_cast<float2 *>(A), lda, W, work, lwork, devInfo));

    cudaErrcheck(cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost));
    assert(0 == info_gpu);
    cudaErrcheck(cudaFree(work));
    cudaErrcheck(cudaFree(devInfo));
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
    cudaErrcheck(cudaMalloc((void**)&devInfo, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnZheevd_bufferSize(solver.handle, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                                 reinterpret_cast<const double2 *>(A), lda, W, &lwork));
    // allocate memery
    cudaErrcheck(cudaMalloc((void**)&work, sizeof(double2) * lwork));
    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnZheevd(solver.handle, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                      reinterpret_cast<double2 *>(A), lda, W, work, lwork, devInfo));

    cudaErrcheck(cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost));
    assert(0 == info_gpu);
    cudaErrcheck(cudaFree(work));
    cudaErrcheck(cudaFree(devInfo));
}

template <typename T>
struct dngvd_op<T, base_device::DEVICE_GPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const base_device::DEVICE_GPU* d,
                    const int nstart,
                    const int ldh,
                    const T* A, // hcc
                    const T* B, // scc
                    Real* W,    // eigenvalue
                    T* V)
    {
        assert(nstart == ldh);
        // A to V
        cudaErrcheck(cudaMemcpy(V, A, sizeof(T) * ldh * nstart, cudaMemcpyDeviceToDevice));
        xhegvd_wrapper(CUBLAS_FILL_MODE_UPPER, nstart, V, ldh,
            (T*)B, ldh, W);
    }
};

template <typename T>
struct dnevx_op<T, base_device::DEVICE_GPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const base_device::DEVICE_GPU* d,
                    const int nstart,
                    const int ldh,
                    const T* A, // hcc
                    const int m,
                    Real* W, // eigenvalue
                    T* V)
    {
        assert(nstart <= ldh);
        // A to V
        cudaErrcheck(cudaMemcpy(V, A, sizeof(T) * nstart * ldh, cudaMemcpyDeviceToDevice));
        xheevd_wrapper(CUBLAS_FILL_MODE_LOWER, nstart, V, ldh, W);
    }
};

template <typename T>
struct dngvx_op<T, base_device::DEVICE_GPU>
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const base_device::DEVICE_GPU* d,
                    const int nbase,
                    const int ldh,
                    T* hcc,
                    T* scc,
                    const int m,
                    Real* eigenvalue,
                    T* vcc)
    {

    }
};

template struct dngvd_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct dnevx_op<std::complex<float>, base_device::DEVICE_GPU>;
template struct dngvx_op<std::complex<float>, base_device::DEVICE_GPU>;

template struct dngvd_op<std::complex<double>, base_device::DEVICE_GPU>;
template struct dnevx_op<std::complex<double>, base_device::DEVICE_GPU>;
template struct dngvx_op<std::complex<double>, base_device::DEVICE_GPU>;

#ifdef __LCAO
template struct dngvd_op<double, base_device::DEVICE_GPU>;
template struct dnevx_op<double, base_device::DEVICE_GPU>;
template struct dngvx_op<double, base_device::DEVICE_GPU>;
#endif

} // namespace hsolver
