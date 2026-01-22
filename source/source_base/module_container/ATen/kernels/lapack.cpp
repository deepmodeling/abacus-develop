#include "source_base/module_device/types.h"
#include <ATen/kernels/lapack.h>

#include <base/third_party/lapack.h>

#include <algorithm> // std::copy
#include <complex>
#include <stdexcept>
#include <string>

namespace container {
namespace kernels {

inline double get_real(const std::complex<double> &x) { return x.real(); }
inline float get_real(const std::complex<float> &x) { return x.real(); }
inline double get_real(const double &x) { return x; }
inline float get_real(const float &x) { return x; }

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

// --- 1. Matrix Decomposition ---
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
struct lapack_geqrf_inplace<T, DEVICE_CPU> {
    void operator()(
        const int m,
        const int n,
        T *A,
        const int lda)
    {
        // Tensor or vector?
        // 1. tau for storing the Householder reflectors
        // tau should be dimension min(m, n)
        int k = std::min(m, n);
        Tensor tau(DataTypeToEnum<T>::value, DeviceType::CpuDevice, {k});
        tau.zero();

        int info = 0;

        // 2. query for workspace size
        int lwork = -1;
        T work_query;
        lapackConnector::geqrf(m, n, A, lda, tau.data<T>(), &work_query, lwork, info);
        if (info != 0) {
            throw std::runtime_error("geqrf workspace query failed with info = " + std::to_string(info));
        }
        // allocate workspace
        lwork = static_cast<int>(get_real(work_query));
        Tensor work(DataTypeToEnum<T>::value, DeviceType::CpuDevice, {lwork});
        work.zero();

        // 3. perform QR decomposition
        // and A is overwritten with upper R.
        // Lower A + tau => Q
        lapackConnector::geqrf(m, n, A, lda, tau.data<T>(), work.data<T>(), lwork, info);
        if (info != 0) {
            throw std::runtime_error("geqrf failed with info = " + std::to_string(info));
        }

        // 4. use orgqr to compute Q
        // workspace query
        lwork = -1;
        lapackConnector::orgqr(m, n, k, A, lda, tau.data<T>(), &work_query, lwork, info);
        if (info != 0) {
            throw std::runtime_error("orgqr workspace query failed with info = " + std::to_string(info));
        }
        // allocate workspace
        lwork = static_cast<int>(get_real(work_query));
        work.resize({lwork});

        // compute Q
        lapackConnector::orgqr(m, n, k, A, lda, tau.data<T>(), work.data<T>(), lwork, info);
        if (info != 0) {
            throw std::runtime_error("orgqr failed with info = " + std::to_string(info));
        }

        // now, A should be overwritten with Q, columns orthogonal

    }
};

// --- 2. Linear System Solvers ---
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


// --- 3. Standard & Generalized Eigenvalue ---
template <typename T>
struct lapack_heevd<T, DEVICE_CPU> {
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const int dim,
        T* Mat,
        const int lda,
        Real* eigen_val)
    {
        char jobz = 'V';        // Compute eigenvalues and eigenvectors
        char uplo = 'U';
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

        lapackConnector::heevd(jobz, uplo, dim, Mat, lda, eigen_val, work.data<T>(), lwork, rwork.data<Real>(), lrwork, iwork.data<int>(), liwork, info);
        if (info != 0) {
            throw std::runtime_error("heevd failed with info = " + std::to_string(info));
        }
    }
};

template <typename T>
struct lapack_heevx<T, DEVICE_CPU> {
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const int n,
        const int lda,
        const T *Mat,
        const int neig,
        Real *eigen_val,
        T *eigen_vec)
    {
        // copy Mat to aux, solve heevx(aux, eigen_val, eigen_vec)
        // input Mat is not referenced in actual heevx LAPACK routines, and aux is destroyed.
        Tensor aux(DataTypeToEnum<T>::value, DeviceType::CpuDevice, {n * lda});
        // Copy Mat to aux since heevx will destroy it
        // aux = Mat
        std::copy(Mat, Mat + n * lda, aux.data<T>());

        char jobz = 'V';        // Compute eigenvalues and eigenvectors
        char range = 'I';       // Find eigenvalues in index range [il, iu]
        char uplo = 'L';        // Use Lower triangle
        int info = 0;
        int found = 0;          // Number of eigenvalues found
        // found should be iu - il + 1, i.e. found = neig
        const int il = 1;
        const int iu = neig;
        Real abstol = 0.0;

        // Workspace query first
        int lwork = -1;
        T work_query;
        Real rwork_query;
        int iwork_query;
        int ifail_query;

        // Dummy call to get optimal workspace size
        // when lwork = -1
        lapackConnector::heevx(
            jobz, range, uplo, n,
            aux.data<T>(), lda,
            0.0, 0.0, il, iu,   // vl, vu not used when range='I'
            abstol,
            found,
            eigen_val,
            eigen_vec, lda,
            &work_query, lwork,
            &rwork_query,
            &iwork_query,
            &ifail_query,
            info);

        if (info != 0) {
            throw std::runtime_error("heevx workspace query failed with info = " + std::to_string(info));
        }

        lwork = static_cast<int>(get_real(work_query));

        // Allocate buffers using Tensor (RAII)
        Tensor work(DataTypeToEnum<T>::value, DeviceType::CpuDevice, {lwork});
        work.zero();

        Tensor rwork(DataTypeToEnum<Real>::value, DeviceType::CpuDevice, {7 * n});
        rwork.zero();

        Tensor iwork(DataType::DT_INT, DeviceType::CpuDevice, {5 * n});
        iwork.zero();

        Tensor ifail(DataType::DT_INT, DeviceType::CpuDevice, {n});
        ifail.zero();

        // Actual call to heevx
        lapackConnector::heevx(
            jobz, range, uplo, n,
            aux.data<T>(), lda,
            0.0, 0.0, il, iu,
            abstol,
            found,
            eigen_val,
            eigen_vec, lda,
            work.data<T>(), lwork,
            rwork.data<Real>(),
            iwork.data<int>(),
            ifail.data<int>(),
            info);

        if (info != 0) {
            throw std::runtime_error("heevx failed with info = " + std::to_string(info));
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
        // eigen_vec = Mat_A
        std::copy(Mat_A, Mat_A + dim*lda, eigen_vec);

        Tensor aux_B(DataTypeToEnum<T>::value, DeviceType::CpuDevice, {dim * lda});
        std::copy(Mat_B, Mat_B + dim * lda, aux_B.data<T>());

        const int itype = 1;
        const char jobz = 'V';
        const char uplo = 'L';
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
        lapackConnector::hegvd(itype, jobz, uplo, dim, eigen_vec, lda, aux_B.data<T>(), lda, eigen_val, work.data<T>(), lwork, rwork.data<Real>(), lrwork, iwork.data<int>(), liwork, info);
        if (info != 0) {
            throw std::runtime_error("hegvd failed with info = " + std::to_string(info));
        }
    }
};


template <typename T>
struct lapack_hegvx<T, DEVICE_CPU> {
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const int n,
        const int lda,
        T *Mat_A,
        T *Mat_B,
        const int m,
        Real *eigen_val,
        T *eigen_vec)
    {
        // first copy Mat_A and Mat_B to auxiliary memory
        // to avoid the origin block being overwritten by hegvx
        Tensor aux_A(DataTypeToEnum<T>::value, DeviceType::CpuDevice, {n * lda});
        std::copy(Mat_A, Mat_A + n * lda, aux_A.data<T>());
        Tensor aux_B(DataTypeToEnum<T>::value, DeviceType::CpuDevice, {n * lda});
        std::copy(Mat_B, Mat_B + n * lda, aux_B.data<T>());

        const int itype = 1;    // ITYPE = 1:  A*x = (lambda)*B*x
        const char jobz = 'V';// JOBZ = 'V':  Compute eigenvalues and eigenvectors.
        const char range = 'I'; // RANGE = 'I': the IL-th through IU-th eigenvalues will be found.
        const char uplo = 'L'; // UPLO = 'L':  Lower triangles of A and B are stored.

        const int il = 1;
        const int iu = m;
        int found = m; // Found, should be iu - il + 1
        int info = 0;

        int lwork = -1;

        T work_query;
        Real rwork_query;

        // set lwork = -1 to query optimal work size
        lapackConnector::hegvx(
                    itype, jobz, range, uplo,
                    n,
                    aux_A.data<T>(), lda,        // A (in/out)
                    aux_B.data<T>(), lda,        // B (in/out)
                    0.0, 0.0,                    // VL, VU (not used)
                    il, iu,                      // IL, IU
                    Real(0.0),                   // ABSTOL
                    found,                       // M (output)
                    eigen_val,                   // W (output)
                    eigen_vec, lda,              // Z (output)
                    &work_query,                 // WORK (query)
                    lwork,
                    &rwork_query,                // RWORK (query)
                    static_cast<int*>(nullptr),  // IWORK (query)
                    static_cast<int*>(nullptr),  // IFAIL (query)
                    info);

        // !>  If LWORK = -1, then a workspace query is assumed; the routine
        // !>  only calculates the optimal size of the WORK array, returns
        // !>  this value as the first entry of the WORK array.
        lwork = static_cast<int>(get_real(work_query));
        lwork = std::max(lwork, 1);

        // work space
        Tensor work(DataTypeToEnum<T>::value, DeviceType::CpuDevice, {lwork});
        work.zero();

        const int lrwork = 7 * n;
        Tensor rwork(DataTypeToEnum<Real>::value, DeviceType::CpuDevice, {lrwork});
        rwork.zero();

        const int liwork = 5 * n;
        Tensor iwork(DataType::DT_INT, DeviceType::CpuDevice, {liwork});
        iwork.zero();

        std::vector<int> ifail(n);

        lapackConnector::hegvx(
                    itype, jobz, range, uplo,
                    n,
                    aux_A.data<T>(), lda,        // A
                    aux_B.data<T>(), lda,        // B
                    0.0, 0.0,                    // VL, VU
                    il, iu,                      // IL, IU
                    Real(0.0),                   // ABSTOL
                    found,                       // M (output)
                    eigen_val,                   // W
                    eigen_vec, lda,              // Z (output)
                    work.data<T>(),              // WORK
                    lwork,
                    rwork.data<Real>(),          // RWORK
                    iwork.data<int>(),           // IWORK
                    ifail.data(),                // IFAIL
                    info);

        if (info < 0) {
            throw std::runtime_error("hegvx failed: illegal argument #" + std::to_string(-info));
        }
        if (info > 0) {
            throw std::runtime_error("hegvx failed to converge. Number of converged eigenvalues: " + std::to_string(info));
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

template struct lapack_geqrf_inplace<float,  DEVICE_CPU>;
template struct lapack_geqrf_inplace<double, DEVICE_CPU>;
template struct lapack_geqrf_inplace<std::complex<float>,  DEVICE_CPU>;
template struct lapack_geqrf_inplace<std::complex<double>, DEVICE_CPU>;

template struct lapack_heevd<float,  DEVICE_CPU>;
template struct lapack_heevd<double, DEVICE_CPU>;
template struct lapack_heevd<std::complex<float>,  DEVICE_CPU>;
template struct lapack_heevd<std::complex<double>, DEVICE_CPU>;

template struct lapack_heevx<float, DEVICE_CPU>;
template struct lapack_heevx<double, DEVICE_CPU>;
template struct lapack_heevx<std::complex<float>, DEVICE_CPU>;
template struct lapack_heevx<std::complex<double>, DEVICE_CPU>;

template struct lapack_hegvd<float,  DEVICE_CPU>;
template struct lapack_hegvd<double, DEVICE_CPU>;
template struct lapack_hegvd<std::complex<float>,  DEVICE_CPU>;
template struct lapack_hegvd<std::complex<double>, DEVICE_CPU>;

template struct lapack_hegvx<float,  DEVICE_CPU>;
template struct lapack_hegvx<double, DEVICE_CPU>;
template struct lapack_hegvx<std::complex<float>,  DEVICE_CPU>;
template struct lapack_hegvx<std::complex<double>, DEVICE_CPU>;

} // namespace kernels
} // namespace container
