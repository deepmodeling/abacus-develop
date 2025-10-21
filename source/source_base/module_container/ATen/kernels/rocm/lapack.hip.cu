#include <vector>
#include <ATen/kernels/lapack.h>
#include <base/third_party/lapack.h>

#include <hip/hip_runtime.h>
#include <thrust/complex.h>
#include <hipsolver/hipsolver.h>

namespace container {
namespace kernels {


static hipsolverHandle_t hipsolver_handle = nullptr;

void createGpuSolverHandle() {
    if (hipsolver_handle == nullptr) {
        hipsolverErrcheck(hipsolverCreate(&hipsolver_handle));
    }
}

void destroyGpuSolverHandle() {
    if (hipsolver_handle != nullptr) {
        hipsolverErrcheck(hipsolverDestroy(hipsolver_handle));
        hipsolver_handle = nullptr;
    }
}

template <typename T>
__global__ void set_matrix_kernel(
    const char uplo,
    T* A,
    const int dim)
{
    int bid = blockIdx.x;
    int tid = threadIdx.x;

    for (int ii = tid; ii < bid + 1; ii += THREADS_PER_BLOCK) {
        if (uplo == 'L') {
            A[ii * dim + bid + 1] = static_cast<T>(0);
        }
        else {
            A[(bid + 1) * dim + ii] = static_cast<T>(0);
        }
    }
}

template <typename T>
struct set_matrix<T, DEVICE_GPU> {
    using Type = typename GetTypeThrust<T>::type;
    void operator() (
        const char& uplo,
        T* A,
        const int& dim)
    {
        set_matrix_kernel<Type><<<dim - 1, THREADS_PER_BLOCK>>>(
            uplo, reinterpret_cast<Type*>(A), dim);
    }
};

template <typename T>
struct lapack_trtri<T, DEVICE_GPU> {
    void operator()(
        const char& uplo,
        const char& diag,
        const int& dim,
        T* Mat,
        const int& lda)
    {
        // TODO: trtri is not implemented in this method yet
        // Cause the trtri in cuSolver is not stable for ABACUS!
        // hipSolverConnector::trtri(hipsolver_handle, uplo, diag, dim, Mat, lda);
        // hipSolverConnector::potri(hipsolver_handle, uplo, diag, dim, Mat, lda);
        std::vector<T> H_Mat(dim * dim, static_cast<T>(0.0));
        hipMemcpy(H_Mat.data(), Mat, sizeof(T) * H_Mat.size(), hipMemcpyDeviceToHost);
        lapack_trtri<T, DEVICE_CPU>()(uplo, diag, dim, H_Mat.data(), lda);
        hipMemcpy(Mat, H_Mat.data(), sizeof(T) * H_Mat.size(), hipMemcpyHostToDevice);
    }
};

template <typename T>
struct lapack_potrf<T, DEVICE_GPU> {
    void operator()(
        const char& uplo,
        const int& dim,
        T* Mat,
        const int& lda)
    {
        // hipSolverConnector::potrf(hipsolver_handle, uplo, dim, Mat, dim);
        std::vector<T> H_Mat(dim * dim, static_cast<T>(0.0));
        hipMemcpy(H_Mat.data(), Mat, sizeof(T) * H_Mat.size(), hipMemcpyDeviceToHost);
        lapack_potrf<T, DEVICE_CPU>()(uplo, dim, H_Mat.data(), lda);
        hipMemcpy(Mat, H_Mat.data(), sizeof(T) * H_Mat.size(), hipMemcpyHostToDevice);
    }
};

template <typename T>
struct lapack_heevd<T, DEVICE_GPU> {
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const char& jobz,
        const char& uplo,
        T* Mat,
        const int& dim,
        Real* eigen_val)
    {
        // hipSolverConnector::heevd(hipsolver_handle, jobz, uplo, dim, Mat, dim, eigen_val);
        std::vector<T> H_Mat(dim * dim, static_cast<T>(0.0));
        std::vector<Real> H_eigen_val(dim, static_cast<Real>(0.0));
        hipMemcpy(H_Mat.data(), Mat, sizeof(T) * H_Mat.size(), hipMemcpyDeviceToHost);
        hipMemcpy(H_eigen_val.data(), eigen_val, sizeof(Real) * H_eigen_val.size(), hipMemcpyDeviceToHost);
        lapack_heevd<T, DEVICE_CPU>()(jobz, uplo, H_Mat.data(), dim, H_eigen_val.data());
        hipMemcpy(Mat, H_Mat.data(), sizeof(T) * H_Mat.size(), hipMemcpyHostToDevice);
        hipMemcpy(eigen_val, H_eigen_val.data(), sizeof(Real) * H_eigen_val.size(), hipMemcpyHostToDevice);
    }
};

template <typename T>
struct lapack_hegvd<T, DEVICE_GPU> {
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const int dim,
        const int lda,
        T* Mat_A,
        T* Mat_B,
        Real* eigen_val,
        T *eigen_vec)
    {
        const int itype = 1;
        const char jobz = 'V';
        const char uplo = 'U';
        hipErrcheck(hipMemcpy(eigen_vec, Mat_A, sizeof(T) * dim * lda, hipMemcpyDeviceToDevice));
        hipSolverConnector::hegvd(hipsolver_handle, itype, jobz, uplo, dim, Mat_A, lda, Mat_B, lda, eigen_val);
    }
};

template struct set_matrix<float,  DEVICE_GPU>;
template struct set_matrix<double, DEVICE_GPU>;
template struct set_matrix<std::complex<float>,  DEVICE_GPU>;
template struct set_matrix<std::complex<double>, DEVICE_GPU>;

template struct lapack_trtri<float,  DEVICE_GPU>;
template struct lapack_trtri<double, DEVICE_GPU>;
template struct lapack_trtri<std::complex<float>,  DEVICE_GPU>;
template struct lapack_trtri<std::complex<double>, DEVICE_GPU>;

template struct lapack_potrf<float,  DEVICE_GPU>;
template struct lapack_potrf<double, DEVICE_GPU>;
template struct lapack_potrf<std::complex<float>,  DEVICE_GPU>;
template struct lapack_potrf<std::complex<double>, DEVICE_GPU>;

template struct lapack_heevd<float,  DEVICE_GPU>;
template struct lapack_heevd<double, DEVICE_GPU>;
template struct lapack_heevd<std::complex<float>,  DEVICE_GPU>;
template struct lapack_heevd<std::complex<double>, DEVICE_GPU>;

template struct lapack_hegvd<float,  DEVICE_GPU>;
template struct lapack_hegvd<double, DEVICE_GPU>;
template struct lapack_hegvd<std::complex<float>,  DEVICE_GPU>;
template struct lapack_hegvd<std::complex<double>, DEVICE_GPU>;

} // namespace kernels
} // namespace container
