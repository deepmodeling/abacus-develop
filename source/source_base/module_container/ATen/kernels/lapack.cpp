#include <ATen/kernels/lapack.h>

#include <base/third_party/lapack.h>

// #include <cstring> // std::memcpy
#include <algorithm> // std::copy

namespace container {
namespace kernels {

template <typename T>
struct set_matrix<T, DEVICE_CPU> {
    void operator() (
        const char& uplo,
        T* A,
        const int& dim)
    {
        if (uplo == 'L') {
            for (int ii = 0; ii < dim; ii++) {
                for (int jj = ii + 1; jj < dim; jj++) {
                    A[ii * dim + jj] = 0;
                }
            }
        }
        else if (uplo == 'U') {
            for (int ii = 0; ii < dim; ii++) {
                for (int jj = 0; jj < ii; jj++) {
                    A[ii * dim + jj] = 0;
                }
            }
        }
    }
};

template <typename T>
struct lapack_trtri<T, DEVICE_CPU> {
    void operator()(
        const char& uplo,
        const char& diag,
        const int& dim,
        T* Mat,
        const int& lda)
    {
        int info = 0;
        lapackConnector::trtri(uplo, diag, dim, Mat, lda, info);
        if (info != 0) {
            throw std::runtime_error("potrf failed with info = " + std::to_string(info));
        }
    }
};

template <typename T>
struct lapack_potrf<T, DEVICE_CPU> {
    void operator()(
        const char& uplo,
        const int& dim,
        T* Mat,
        const int& lda)
    {
        int info = 0;
        lapackConnector::potrf(uplo, dim, Mat, dim, info);
        if (info != 0) {
            throw std::runtime_error("potrf failed with info = " + std::to_string(info));
        }
    }
};

template <typename T>
struct lapack_heevd<T, DEVICE_CPU> {
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const char& jobz,
        const char& uplo,
        T* Mat,
        const int& dim,
        Real* eigen_val)
    {
        int info = 0;
        int lwork = std::max(2 * dim + dim * dim, 1 + 6 * dim + 2 * dim * dim);
        Tensor work(DataTypeToEnum<T>::value, DeviceType::CpuDevice, {lwork});
        work.zero();

        int lrwork = 1 + 5 * dim + 2 * dim * dim;
        Tensor rwork(DataTypeToEnum<Real>::value, DeviceType::CpuDevice, {lrwork});
        rwork.zero();

        int liwork = 3 + 5 * dim;
        Tensor iwork(DataTypeToEnum<int>::value, DeviceType::CpuDevice, {liwork});
        iwork.zero();

        lapackConnector::heevd(jobz, uplo, dim, Mat, dim, eigen_val, work.data<T>(), lwork, rwork.data<Real>(), lrwork, iwork.data<int>(), liwork, info);
        if (info != 0) {
            throw std::runtime_error("heevd failed with info = " + std::to_string(info));
        }
    }
};

template <typename T>
struct lapack_hegvd<T, DEVICE_CPU> {
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const int dim,
        const int lda,
        T *Mat_A,
        T *Mat_B,
        Real *eigen_val,
        T *eigen_vec)
    {
        // first copy Mat_A to eigen_vec
        // then pass as argument "A" in lapack hegvd
        // and this block of memory will be overwritten by eigenvectors
        // for (int i = 0; i < dim * lda; ++i){
        //     eigen_vec[i] = Mat_A[i];
        // }
        // std::memcpy(eigen_vec, Mat_A, sizeof(T) * dim * lda);
        // eigen_vec = Mat_A
        std::copy(Mat_A, Mat_A + dim*lda, eigen_vec);

        const int itype = 1;
        const char jobz = 'V';
        const char uplo = 'U';
        int info = 0;
        int lwork = std::max(2 * dim + dim * dim, 1 + 6 * dim + 2 * dim * dim);
        Tensor work(DataTypeToEnum<T>::value, DeviceType::CpuDevice, {lwork});
        work.zero();

        int lrwork = 1 + 5 * dim + 2 * dim * dim;
        Tensor rwork(DataTypeToEnum<Real>::value, DeviceType::CpuDevice, {lrwork});
        rwork.zero();

        int liwork = 3 + 5 * dim;
        Tensor iwork(DataType::DT_INT, DeviceType::CpuDevice, {liwork});
        iwork.zero();

        // After this, eigen_vec will contain the matrix Z of eigenvectors
        lapackConnector::hegvd(itype, jobz, uplo, dim, eigen_vec, lda, Mat_B, lda, eigen_val, work.data<T>(), lwork, rwork.data<Real>(), lrwork, iwork.data<int>(), liwork, info);
        if (info != 0) {
            throw std::runtime_error("hegvd failed with info = " + std::to_string(info));
        }
    }
};


// template <typename T>
// struct lapack_hegvx<T, DEVICE_CPU> {
//     using Real = typename GetTypeReal<T>::type;
//     void operator()(
//         const int n,
//         const int lda,
//         T *A,
//         T *B,
//         const int m,
//         Real *eigen_val,
//         T *eigen_vec)
//     {
//         int info = 0;

//         int mm = m;

//         int lwork = -1;

//         T *work = new T[1];
//         Real *rwork = new Real[7 * n];
//         int *iwork = new int[5 * n];
//         int *ifail = new int[n];

//         // set lwork = -1 to query optimal work size
//         lapackConnector::hegvx(1, // ITYPE = 1:  A*x = (lambda)*B*x
//             'V', 'I', 'U',
//             n, A, lda, B, lda,
//             0.0, 0.0,
//             1, m, 0.0, mm,
//             eigen_val, eigen_vec, lda,
//             work,
//             lwork,      // lwork = 1, query optimal size.
//             rwork, iwork, ifail,
//             info);

//         // !>  If LWORK = -1, then a workspace query is assumed; the routine
//         // !>  only calculates the optimal size of the WORK array, returns
//         // !>  this value as the first entry of the WORK array.
//         lwork = int(get_real(work[0]));
//         delete[] work;
//         work = new T[lwork];

//         lapackConnector::hegvx(
//             1, // ITYPE = 1:  A*x = (lambda)*B*x
//             'V',    // JOBZ = 'V':  Compute eigenvalues and eigenvectors.
//             'I',    // RANGE = 'I': the IL-th through IU-th eigenvalues will be found.
//             'U',    // UPLO = 'U':  Upper triangles of A and B are stored.
//             n,      // order of the matrices A and B.
//             A,      // A is COMPLEX*16 array  dimension (LDA, N)
//             lda,    // leading dimension of the array A.
//             B,      // B is COMPLEX*16 array, dimension (LDB, N)
//             lda,    // assume that leading dimension of B is the same as A.
//             0.0,    // VL, Not referenced if RANGE = 'A' or 'I'.
//             0.0,    // VU, Not referenced if RANGE = 'A' or 'I'.
//             1,      // IL: If RANGE='I', the index of the smallest eigenvalue to be returned. 1 <= IL <= IU <= N,
//             m,      // IU: If RANGE='I', the index of the largest eigenvalue to be returned. 1 <= IL <= IU <= N,
//             0.0,    // ABSTOL
//             mm,     // M: The total number of eigenvalues found.  0 <= M <= N. if RANGE = 'I', M = IU-IL+1.
//             eigen_val,  // W store eigenvalues
//             eigen_vec,  // Z store eigenvector
//             lda,        // LDZ: The leading dimension of the array Z.
//             work,
//             lwork,
//             rwork,
//             iwork,
//             ifail,
//             info);

//         delete[] work;
//         delete[] rwork;
//         delete[] iwork;
//         delete[] ifail;
//     }
// };

template <typename T>
struct lapack_getrf<T, DEVICE_CPU> {
    void operator()(
        const int& m,
        const int& n,
        T* Mat,
        const int& lda,
        int* ipiv)
    {
        int info = 0;
        lapackConnector::getrf(m, n, Mat, lda, ipiv, info);
        if (info != 0) {
            throw std::runtime_error("getrf failed with info = " + std::to_string(info));
        }
    }
};

template <typename T>
struct lapack_getri<T, DEVICE_CPU> {
    void operator()(
        const int& n,
        T* Mat,
        const int& lda,
        const int* ipiv,
        T* work,
        const int& lwork)
    {
        int info = 0;
        lapackConnector::getri(n, Mat, lda, ipiv, work, lwork, info);
        if (info != 0) {
            throw std::runtime_error("getri failed with info = " + std::to_string(info));
        }
    }
};

template <typename T>
struct lapack_getrs<T, DEVICE_CPU> {
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
        int info = 0;
        lapackConnector::getrs(trans, n, nrhs, A, lda, ipiv, B, ldb, info);
        if (info != 0) {
            throw std::runtime_error("getrs failed with info = " + std::to_string(info));
        }
    }
};

template struct set_matrix<float,  DEVICE_CPU>;
template struct set_matrix<double, DEVICE_CPU>;
template struct set_matrix<std::complex<float>,  DEVICE_CPU>;
template struct set_matrix<std::complex<double>, DEVICE_CPU>;

template struct lapack_potrf<float,  DEVICE_CPU>;
template struct lapack_potrf<double, DEVICE_CPU>;
template struct lapack_potrf<std::complex<float>,  DEVICE_CPU>;
template struct lapack_potrf<std::complex<double>, DEVICE_CPU>;

template struct lapack_trtri<float,  DEVICE_CPU>;
template struct lapack_trtri<double, DEVICE_CPU>;
template struct lapack_trtri<std::complex<float>,  DEVICE_CPU>;
template struct lapack_trtri<std::complex<double>, DEVICE_CPU>;

template struct lapack_heevd<float,  DEVICE_CPU>;
template struct lapack_heevd<double, DEVICE_CPU>;
template struct lapack_heevd<std::complex<float>,  DEVICE_CPU>;
template struct lapack_heevd<std::complex<double>, DEVICE_CPU>;

template struct lapack_hegvd<float,  DEVICE_CPU>;
template struct lapack_hegvd<double, DEVICE_CPU>;
template struct lapack_hegvd<std::complex<float>,  DEVICE_CPU>;
template struct lapack_hegvd<std::complex<double>, DEVICE_CPU>;

template struct lapack_getrf<float,  DEVICE_CPU>;
template struct lapack_getrf<double, DEVICE_CPU>;
template struct lapack_getrf<std::complex<float>,  DEVICE_CPU>;
template struct lapack_getrf<std::complex<double>, DEVICE_CPU>;

template struct lapack_getri<float, DEVICE_CPU>;
template struct lapack_getri<double, DEVICE_CPU>;
template struct lapack_getri<std::complex<float>, DEVICE_CPU>;
template struct lapack_getri<std::complex<double>, DEVICE_CPU>;

template struct lapack_getrs<float, DEVICE_CPU>;
template struct lapack_getrs<double, DEVICE_CPU>;
template struct lapack_getrs<std::complex<float>, DEVICE_CPU>;
template struct lapack_getrs<std::complex<double>, DEVICE_CPU>;

} // namespace kernels
} // namespace container
