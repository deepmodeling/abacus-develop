#include <ATen/kernels/lapack_op.h>
#include <base/third_party/lapack.h>

namespace container {
namespace op {

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
struct lapack_dnevd<T, DEVICE_CPU> {
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

        lapackConnector::dnevd(jobz, uplo, dim, Mat, dim, eigen_val,  work.data<T>(), lwork, rwork.data<Real>(), lrwork, iwork.data<int>(), liwork, info);
        if (info != 0) {
            throw std::runtime_error("dnevd failed with info = " + std::to_string(info));
        }
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

        lapackConnector::dngvd(itype, jobz, uplo, dim, Mat_A, dim, Mat_B, dim, eigen_val, work.data<T>(), lwork, rwork.data<Real>(), lrwork, iwork.data<int>(), liwork, info);
        if (info != 0) {
            throw std::runtime_error("dngvd failed with info = " + std::to_string(info));
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

template struct lapack_dnevd<float,  DEVICE_CPU>;
template struct lapack_dnevd<double, DEVICE_CPU>;
template struct lapack_dnevd<std::complex<float>,  DEVICE_CPU>;
template struct lapack_dnevd<std::complex<double>, DEVICE_CPU>;

template struct lapack_dngvd<float,  DEVICE_CPU>;
template struct lapack_dngvd<double, DEVICE_CPU>;
template struct lapack_dngvd<std::complex<float>,  DEVICE_CPU>;
template struct lapack_dngvd<std::complex<double>, DEVICE_CPU>;

} // namespace op
} // namespace container