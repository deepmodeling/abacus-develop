#include <ATen/kernels/lapack.h>
#include <base/third_party/lapack.h>

#include <cublas_v2.h>
#include <cusolverDn.h>
#include <cuda_runtime.h>
#include <thrust/complex.h>

#include <cassert>


namespace container {
namespace kernels {


static cusolverDnHandle_t cusolver_handle = nullptr;

void createGpuSolverHandle() {
    if (cusolver_handle == nullptr) {
        cusolverErrcheck(cusolverDnCreate(&cusolver_handle));
    }
}

void destroyGpuSolverHandle() {
    if (cusolver_handle != nullptr) {
        cusolverErrcheck(cusolverDnDestroy(cusolver_handle));
        cusolver_handle = nullptr;
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

        cudaCheckOnDebug();
    }
};



// --- 1. Matrix Decomposition ---
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
        cuSolverConnector::trtri(cusolver_handle, uplo, diag, dim, Mat, lda);
        // cuSolverConnector::potri(cusolver_handle, uplo, diag, dim, Mat, lda);
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
        cuSolverConnector::potrf(cusolver_handle, uplo, dim, Mat, dim);
    }
};

template <typename T>
struct lapack_getrf<T, DEVICE_GPU> {
    void operator()(
        const int& m,
        const int& n,
        T* Mat,
        const int& lda,
        int* ipiv)
    {
        cuSolverConnector::getrf(cusolver_handle, m, n, Mat, lda, ipiv);
    }
};

template <typename T>
struct lapack_getri<T, DEVICE_GPU> {
    void operator()(
        const int& n,
        T* Mat,
        const int& lda,
        const int* ipiv,
        T* work,
        const int& lwork)
    {
        throw std::runtime_error("cuSOLVER does not provide LU-based matrix inversion interface (getri). To compute the inverse on GPU, use getrs instead.");
    }
};


template <typename T>
struct lapack_geqrf_inplace<T, DEVICE_GPU> {
    void operator()(
        const int m,
        const int n,
        T *d_A,
        const int lda)
    {
        const int k = std::min(m, n);

        // Allocate tau on device
        T *d_tau;
        cudaErrcheck(cudaMalloc(&d_tau, sizeof(T) * k));

        cuSolverConnector::geqrf(cusolver_handle, m, n, d_A, lda, d_tau);

        cuSolverConnector::orgqr(cusolver_handle, m, n, k, d_A, lda, d_tau);

        cudaErrcheck(cudaFree(d_tau));

        // // geqrf: workspace query

        // // In practice, we use helper function to get lwork
        // // Or use magma for better interface
        // // Let's assume we have a way to get lwork
        // // For now, do a dummy call to get it
        // size_t workspaceInBytes = 0;
        // cusolverErrcheck(cusolverDnXgeqrf_bufferSize(
        //     cusolverH, m, n,
        //     getCudaDataType<T>::type, d_A, lda,
        //     getCudaDataType<T>::type, // for tau
        //     CUDA_R_32F, // numerical precision
        //     CUSOLVER_WORKSPACE_QUERY_USE_MAX, &workspaceInBytes));

        // lwork = static_cast<int>(workspaceInBytes / sizeof(T));

        // // Allocate workspace
        // T *d_work;
        // cudaErrcheck(cudaMalloc(&d_work, sizeof(T) * lwork));

        // // 3. Perform geqrf
        // cusolverErrcheck(cusolverDnXgeqrf(
        //     cusolverH, m, n,
        //     getCudaDataType<T>::type, d_A, lda,
        //     d_tau,
        //     getCudaDataType<T>::type,
        //     d_work, lwork * sizeof(T),
        //     d_info));

        // int info;
        // cudaErrcheck(cudaMemcpy(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
        // if (info != 0) {
        //     throw std::runtime_error("cuSOLVER geqrf failed with info = " + std::to_string(info));
        // }

        // // 4. Generate Q using orgqr
        // // Query workspace for orgqr
        // cusolverErrcheck(cusolverDnXorgqr_bufferSize(
        //     cusolverH, m, n, k,
        //     getCudaDataType<T>::type, d_A, lda,
        //     getCudaDataType<T>::type, d_tau,
        //     CUDA_R_32F,
        //     CUSOLVER_WORKSPACE_QUERY_USE_MAX, &workspaceInBytes));

        // lwork = static_cast<int>(workspaceInBytes / sizeof(T));
        // cudaErrcheck(cudaRealloc(&d_work, sizeof(T) * lwork)); // or realloc

        // // orgqr: generate Q
        // cusolverErrcheck(cusolverDnXorgqr(
        //     cusolverH, m, n, k,
        //     getCudaDataType<T>::type, d_A, lda,
        //     getCudaDataType<T>::type, d_tau,
        //     d_work, lwork * sizeof(T),
        //     d_info));

        // cudaErrcheck(cudaMemcpy(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
        // if (info != 0) {
        //     throw std::runtime_error("cuSOLVER orgqr failed with info = " + std::to_string(info));
        // }

        // // Clean up
        // cudaErrcheck(cudaFree(d_tau));
        // cudaErrcheck(cudaFree(d_work));
        // cudaErrcheck(cudaFree(d_info));
    }
};

// --- 2. Linear System Solvers ---
template <typename T>
struct lapack_getrs<T, DEVICE_GPU> {
    void operator()(
        const char& trans,
        const int& n,
        const int& nrhs,
        T* A,
        const int& lda,
        const int* ipiv,
        T* B,
        const int& ldb)
    {
        cuSolverConnector::getrs(cusolver_handle, trans, n, nrhs, A, lda, ipiv, B, ldb);
    }
};


// --- 3. Standard & Generalized Eigenvalue ---
template <typename T>
struct lapack_heevd<T, DEVICE_GPU> {
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const int dim,
        T* Mat,
        const int lda,
        Real* eigen_val)
    {
        char jobz = 'V';        // Compute eigenvalues and eigenvectors
        char uplo = 'U';
        cuSolverConnector::heevd(cusolver_handle, jobz, uplo, dim, Mat, lda, eigen_val);
    }
};

template <typename T>
struct lapack_heevx<T, DEVICE_GPU> {
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const int n,
        const int lda,
        const T *d_Mat,
        const int neig,
        Real *d_eigen_val,
        T *d_eigen_vec)
    {
        assert(n <= lda);
        // copy d_Mat to d_eigen_vec, and results will be overwritten into d_eigen_vec
        // by cuSolver
        cudaErrcheck(cudaMemcpy(d_eigen_vec, d_Mat, sizeof(T) * n * lda, cudaMemcpyDeviceToDevice));

        int meig = 0;

        cuSolverConnector::heevdx(
            cusolver_handle,
            n,
            lda,
            d_eigen_vec,
            'V',        // jobz: compute vectors
            'L',        // uplo: lower triangle
            'I',        // range: by index
            1, neig,    // il, iu
            Real(0), Real(0), // vl, vu (unused)
            d_eigen_val,
            &meig
        );

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
        const char uplo = 'L';
        cudaErrcheck(cudaMemcpy(eigen_vec, Mat_A, sizeof(T) * dim * lda, cudaMemcpyDeviceToDevice));

        // prevent B from being overwritten by Cholesky
        T *d_B_backup = nullptr;
        cudaErrcheck(cudaMalloc(&d_B_backup, sizeof(T) * dim * lda));
        cudaErrcheck(cudaMemcpy(d_B_backup, Mat_B, sizeof(T) * dim * lda, cudaMemcpyDeviceToDevice));

        cuSolverConnector::hegvd(cusolver_handle, itype, jobz, uplo, dim,
                eigen_vec, lda,
                d_B_backup, lda,
                eigen_val);
        cudaErrcheck(cudaFree(d_B_backup));
    }
};

template <typename T>
struct lapack_hegvx<T, DEVICE_GPU> {
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const int n,
        const int lda,
        T *A,
        T *B,
        const int m,
        Real *eigen_val,
        T *eigen_vec)
    {
        const int itype = 1;
        const char jobz = 'V';
        const char range = 'I';
        const char uplo = 'U';
        int meig = 0;

        // this hegvdx will protect the input A, B from being overwritten
        // and write the eigenvectors into eigen_vec.
        cuSolverConnector::hegvdx(cusolver_handle,
            itype, jobz, range, uplo,
            n, lda, A, B,
            Real(0), Real(0),
            1, m, &meig,
            eigen_val, eigen_vec);
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


template struct lapack_getrs<float,  DEVICE_GPU>;
template struct lapack_getrs<double, DEVICE_GPU>;
template struct lapack_getrs<std::complex<float>,  DEVICE_GPU>;
template struct lapack_getrs<std::complex<double>, DEVICE_GPU>;


template struct lapack_heevd<float,  DEVICE_GPU>;
template struct lapack_heevd<double, DEVICE_GPU>;
template struct lapack_heevd<std::complex<float>,  DEVICE_GPU>;
template struct lapack_heevd<std::complex<double>, DEVICE_GPU>;

template struct lapack_heevx<float, DEVICE_GPU>;
template struct lapack_heevx<double, DEVICE_GPU>;
template struct lapack_heevx<std::complex<float>, DEVICE_GPU>;
template struct lapack_heevx<std::complex<double>, DEVICE_GPU>;

template struct lapack_hegvd<float,  DEVICE_GPU>;
template struct lapack_hegvd<double, DEVICE_GPU>;
template struct lapack_hegvd<std::complex<float>,  DEVICE_GPU>;
template struct lapack_hegvd<std::complex<double>, DEVICE_GPU>;

template struct lapack_hegvx<float,  DEVICE_GPU>;
template struct lapack_hegvx<double, DEVICE_GPU>;
template struct lapack_hegvx<std::complex<float>,  DEVICE_GPU>;
template struct lapack_hegvx<std::complex<double>, DEVICE_GPU>;

template struct lapack_getrf<float,  DEVICE_GPU>;
template struct lapack_getrf<double, DEVICE_GPU>;
template struct lapack_getrf<std::complex<float>,  DEVICE_GPU>;
template struct lapack_getrf<std::complex<double>, DEVICE_GPU>;

template struct lapack_getri<float,  DEVICE_GPU>;
template struct lapack_getri<double, DEVICE_GPU>;
template struct lapack_getri<std::complex<float>,  DEVICE_GPU>;
template struct lapack_getri<std::complex<double>, DEVICE_GPU>;

template struct lapack_geqrf_inplace<float,  DEVICE_GPU>;
template struct lapack_geqrf_inplace<double, DEVICE_GPU>;
template struct lapack_geqrf_inplace<std::complex<float>,  DEVICE_GPU>;
template struct lapack_geqrf_inplace<std::complex<double>, DEVICE_GPU>;

} // namespace kernels
} // namespace container
