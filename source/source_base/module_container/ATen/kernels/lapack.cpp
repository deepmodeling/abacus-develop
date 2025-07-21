#include <ATen/kernels/lapack.h>

#include <base/third_party/lapack.h>

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
        lapackConnector::trtri(uplo, diag, dim, Mat, lda);
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
        lapackConnector::potrf(uplo, dim, Mat, dim);
    }
};

template <typename T>
struct lapack_dnevd<T, DEVICE_CPU> {
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const char& jobz,
        const char& uplo,
        T* Mat,
        const int& dim,
        Real* eigen_val)
    {
        lapackConnector::dnevd(jobz, uplo, dim, Mat, dim, eigen_val);
    }
};

template <typename T>
struct lapack_dngvd<T, DEVICE_CPU> {
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const int& itype,
        const char& jobz,
        const char& uplo,
        T* Mat_A,
        T* Mat_B,
        const int& dim,
        Real* eigen_val)
    {
        lapackConnector::dngvd(itype, jobz, uplo, dim, Mat_A, dim, Mat_B, dim, eigen_val);
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
        lapackConnector::getrf(m, n, Mat, lda, ipiv);
    }
};

template <typename T>
struct lapack_getri<T, DEVICE_CPU> {
    void operator()(
        const int& n,
        T* Mat,
        const int& lda,
        int* ipiv,
        T* work,
        const int& lwork)
    {
        lapackConnector::getri(n, Mat, lda, ipiv);
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
        int* ipiv,
        T* B,
        const int& ldb)
    {
        lapackConnector::getrs(trans, n, nrhs, A, lda, ipiv, B, ldb);
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

template struct lapack_dnevd<float,  DEVICE_CPU>;
template struct lapack_dnevd<double, DEVICE_CPU>;
template struct lapack_dnevd<std::complex<float>,  DEVICE_CPU>;
template struct lapack_dnevd<std::complex<double>, DEVICE_CPU>;

template struct lapack_dngvd<float,  DEVICE_CPU>;
template struct lapack_dngvd<double, DEVICE_CPU>;
template struct lapack_dngvd<std::complex<float>,  DEVICE_CPU>;
template struct lapack_dngvd<std::complex<double>, DEVICE_CPU>;

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