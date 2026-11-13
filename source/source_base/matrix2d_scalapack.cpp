#include "source_base/matrix2d_scalapack.h"

#ifdef __MPI
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/matrix2d_view.h"

#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

namespace matrix2d_impl
{
namespace scalapack
{
void transpose(Matrix2dView<double> result,
               Matrix2dView<const double> source)
{
    const int m = source.nrglb();
    const int n = source.ncglb();
    const double alpha = 1.0;
    const double* source_data = source.data();
    const int ia = 1;
    const int ja = 1;
    const int* source_descriptor = source.desc();
    const double beta = 0.0;
    double* result_data = result.data();
    const int ic = 1;
    const int jc = 1;
    const int* result_descriptor = result.desc();
    ScalapackConnector::tranu(
        /* m= */ n,
        /* n= */ m,
        /* alpha= */ alpha,
        /* a= */ source_data,
        /* ia= */ ia,
        /* ja= */ ja,
        /* desca= */ source_descriptor,
        /* beta= */ beta,
        /* c= */ result_data,
        /* ic= */ ic,
        /* jc= */ jc,
        /* descc= */ result_descriptor);
}

void conjugate_transpose(Matrix2dView<double> result,
                         Matrix2dView<const double> source)
{
    transpose(result, source);
}

void transpose(Matrix2dView<std::complex<double>> result,
               Matrix2dView<const std::complex<double>> source)
{
    const int m = source.nrglb();
    const int n = source.ncglb();
    const std::complex<double> alpha(1.0, 0.0);
    const std::complex<double>* source_data = source.data();
    const int ia = 1;
    const int ja = 1;
    const int* source_descriptor = source.desc();
    const std::complex<double> beta(0.0, 0.0);
    std::complex<double>* result_data = result.data();
    const int ic = 1;
    const int jc = 1;
    const int* result_descriptor = result.desc();
    ScalapackConnector::tranu(
        /* m= */ n,
        /* n= */ m,
        /* alpha= */ alpha,
        /* a= */ source_data,
        /* ia= */ ia,
        /* ja= */ ja,
        /* desca= */ source_descriptor,
        /* beta= */ beta,
        /* c= */ result_data,
        /* ic= */ ic,
        /* jc= */ jc,
        /* descc= */ result_descriptor);
}

void conjugate_transpose(
    Matrix2dView<std::complex<double>> result,
    Matrix2dView<const std::complex<double>> source)
{
    const int m = source.nrglb();
    const int n = source.ncglb();
    const std::complex<double> alpha(1.0, 0.0);
    const std::complex<double>* source_data = source.data();
    const int ia = 1;
    const int ja = 1;
    const int* source_descriptor = source.desc();
    const std::complex<double> beta(0.0, 0.0);
    std::complex<double>* result_data = result.data();
    const int ic = 1;
    const int jc = 1;
    const int* result_descriptor = result.desc();
    ScalapackConnector::tranc(
        /* m= */ n,
        /* n= */ m,
        /* alpha= */ alpha,
        /* a= */ source_data,
        /* ia= */ ia,
        /* ja= */ ja,
        /* desca= */ source_descriptor,
        /* beta= */ beta,
        /* c= */ result_data,
        /* ic= */ ic,
        /* jc= */ jc,
        /* descc= */ result_descriptor);
}

namespace
{
void check_scalapack_error(const char* operation, const int info)
{
    if (info < 0)
    {
        throw std::invalid_argument(std::string(operation)
                                    + " received an invalid argument");
    }
    if (info > 0)
    {
        throw std::runtime_error(std::string(operation)
                                 + " failed because the matrix is singular");
    }
}

template <typename T>
int queried_workspace_size(const T& value)
{
    return static_cast<int>(value);
}

int queried_workspace_size(const std::complex<double>& value)
{
    return static_cast<int>(value.real());
}
} // namespace

void inv(Matrix2dView<double> matrix)
{
    const int n = matrix.nrglb();
    const int ia = 1;
    const int ja = 1;
    const int* descriptor = matrix.desc();
    const int pivot_size = matrix.nrloc() + descriptor[4];
    std::vector<int> ipiv(static_cast<std::size_t>(pivot_size));
    int info = 0;
    ScalapackConnector::getrf(n,
                              n,
                              matrix.data(),
                              ia,
                              ja,
                              descriptor,
                              ipiv.data(),
                              &info);
    check_scalapack_error("ScaLAPACK getrf", info);

    int lwork = -1;
    int liwork = -1;
    std::vector<double> work(1);
    std::vector<int> iwork(1);
    ScalapackConnector::getri(n,
                              matrix.data(),
                              ia,
                              ja,
                              descriptor,
                              ipiv.data(),
                              work.data(),
                              &lwork,
                              iwork.data(),
                              &liwork,
                              &info);
    check_scalapack_error("ScaLAPACK getri workspace query", info);
    lwork = queried_workspace_size(work[0]);
    liwork = iwork[0];
    work.resize(static_cast<std::size_t>(lwork));
    iwork.resize(static_cast<std::size_t>(liwork));
    ScalapackConnector::getri(n,
                              matrix.data(),
                              ia,
                              ja,
                              descriptor,
                              ipiv.data(),
                              work.data(),
                              &lwork,
                              iwork.data(),
                              &liwork,
                              &info);
    check_scalapack_error("ScaLAPACK getri", info);
}

void inv(Matrix2dView<std::complex<double>> matrix)
{
    const int n = matrix.nrglb();
    const int ia = 1;
    const int ja = 1;
    const int* descriptor = matrix.desc();
    const int pivot_size = matrix.nrloc() + descriptor[4];
    std::vector<int> ipiv(static_cast<std::size_t>(pivot_size));
    int info = 0;
    ScalapackConnector::getrf(n,
                              n,
                              matrix.data(),
                              ia,
                              ja,
                              descriptor,
                              ipiv.data(),
                              &info);
    check_scalapack_error("ScaLAPACK getrf", info);

    int lwork = -1;
    int liwork = -1;
    std::vector<std::complex<double>> work(1);
    std::vector<int> iwork(1);
    ScalapackConnector::getri(n,
                              matrix.data(),
                              ia,
                              ja,
                              descriptor,
                              ipiv.data(),
                              work.data(),
                              &lwork,
                              iwork.data(),
                              &liwork,
                              &info);
    check_scalapack_error("ScaLAPACK getri workspace query", info);
    lwork = queried_workspace_size(work[0]);
    liwork = iwork[0];
    work.resize(static_cast<std::size_t>(lwork));
    iwork.resize(static_cast<std::size_t>(liwork));
    ScalapackConnector::getri(n,
                              matrix.data(),
                              ia,
                              ja,
                              descriptor,
                              ipiv.data(),
                              work.data(),
                              &lwork,
                              iwork.data(),
                              &liwork,
                              &info);
    check_scalapack_error("ScaLAPACK getri", info);
}

void solve(Matrix2dView<double> matrix,
           Matrix2dView<double> rhs)
{
    const int n = matrix.nrglb();
    const int nrhs = rhs.ncglb();
    const int ia = 1;
    const int ja = 1;
    const int ib = 1;
    const int jb = 1;
    const int pivot_size = matrix.nrloc() + matrix.desc()[4];
    std::vector<int> ipiv(static_cast<std::size_t>(pivot_size));
    int info = 0;
    ScalapackConnector::gesv(n,
                             nrhs,
                             matrix.data(),
                             ia,
                             ja,
                             matrix.desc(),
                             ipiv.data(),
                             rhs.data(),
                             ib,
                             jb,
                             rhs.desc(),
                             &info);
    check_scalapack_error("ScaLAPACK gesv", info);
}

void solve(Matrix2dView<std::complex<double>> matrix,
           Matrix2dView<std::complex<double>> rhs)
{
    const int n = matrix.nrglb();
    const int nrhs = rhs.ncglb();
    const int ia = 1;
    const int ja = 1;
    const int ib = 1;
    const int jb = 1;
    const int pivot_size = matrix.nrloc() + matrix.desc()[4];
    std::vector<int> ipiv(static_cast<std::size_t>(pivot_size));
    int info = 0;
    ScalapackConnector::gesv(n,
                             nrhs,
                             matrix.data(),
                             ia,
                             ja,
                             matrix.desc(),
                             ipiv.data(),
                             rhs.data(),
                             ib,
                             jb,
                             rhs.desc(),
                             &info);
    check_scalapack_error("ScaLAPACK gesv", info);
}

void matmul(Matrix2dView<double> result,
            Matrix2dView<const double> left,
            Matrix2dView<const double> right)
{
    const char transa = 'N';
    const char transb = 'N';
    const int m = left.nrglb();
    const int n = right.ncglb();
    const int k = left.ncglb();
    const double alpha = 1.0;
    const double* left_data = left.data();
    const int ia = 1;
    const int ja = 1;
    const int* left_descriptor = left.desc();
    const double* right_data = right.data();
    const int ib = 1;
    const int jb = 1;
    const int* right_descriptor = right.desc();
    const double beta = 0.0;
    double* result_data = result.data();
    const int ic = 1;
    const int jc = 1;
    const int* result_descriptor = result.desc();
    ScalapackConnector::gemm(
        /* transa= */ transa,
        /* transb= */ transb,
        /* M= */ m,
        /* N= */ n,
        /* K= */ k,
        /* alpha= */ alpha,
        /* A= */ left_data,
        /* IA= */ ia,
        /* JA= */ ja,
        /* DESCA= */ left_descriptor,
        /* B= */ right_data,
        /* IB= */ ib,
        /* JB= */ jb,
        /* DESCB= */ right_descriptor,
        /* beta= */ beta,
        /* C= */ result_data,
        /* IC= */ ic,
        /* JC= */ jc,
        /* DESCC= */ result_descriptor);
}

void matmul(Matrix2dView<std::complex<double>> result,
            Matrix2dView<const std::complex<double>> left,
            Matrix2dView<const std::complex<double>> right)
{
    const char transa = 'N';
    const char transb = 'N';
    const int m = left.nrglb();
    const int n = right.ncglb();
    const int k = left.ncglb();
    const std::complex<double> alpha(1.0, 0.0);
    const std::complex<double>* left_data = left.data();
    const int ia = 1;
    const int ja = 1;
    const int* left_descriptor = left.desc();
    const std::complex<double>* right_data = right.data();
    const int ib = 1;
    const int jb = 1;
    const int* right_descriptor = right.desc();
    const std::complex<double> beta(0.0, 0.0);
    std::complex<double>* result_data = result.data();
    const int ic = 1;
    const int jc = 1;
    const int* result_descriptor = result.desc();
    ScalapackConnector::gemm(
        /* transa= */ transa,
        /* transb= */ transb,
        /* M= */ m,
        /* N= */ n,
        /* K= */ k,
        /* alpha= */ alpha,
        /* A= */ left_data,
        /* IA= */ ia,
        /* JA= */ ja,
        /* DESCA= */ left_descriptor,
        /* B= */ right_data,
        /* IB= */ ib,
        /* JB= */ jb,
        /* DESCB= */ right_descriptor,
        /* beta= */ beta,
        /* C= */ result_data,
        /* IC= */ ic,
        /* JC= */ jc,
        /* DESCC= */ result_descriptor);
}

} // namespace scalapack
} // namespace matrix2d_impl
#endif
