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
void assert_compatible_layout(const int* left_desc,
                              const int* right_desc)
{
    if (left_desc == nullptr
        || right_desc == nullptr
        || left_desc[1] != right_desc[1]
        || left_desc[4] != right_desc[4]
        || left_desc[5] != right_desc[5]
        || left_desc[6] != right_desc[6]
        || left_desc[7] != right_desc[7])
    {
        throw std::invalid_argument(
            "matrix operations require compatible ScaLAPACK descriptors");
    }
}

void transpose(Matrix2dView<double> result,
               Matrix2dView<const double> source)
{
    const int m = source.nrglb();
    const int n = source.ncglb();
    const double alpha = 1.0;
    const double* source_data = source.data();
    const int ia = 1;
    const int ja = 1;
    const int* desc_src = source.desc();
    const double beta = 0.0;
    double* result_data = result.data();
    const int ic = 1;
    const int jc = 1;
    const int* desc_res = result.desc();
    ScalapackConnector::tranu(
        /* m= */ n,
        /* n= */ m,
        /* alpha= */ alpha,
        /* a= */ source_data,
        /* ia= */ ia,
        /* ja= */ ja,
        /* desca= */ desc_src,
        /* beta= */ beta,
        /* c= */ result_data,
        /* ic= */ ic,
        /* jc= */ jc,
        /* descc= */ desc_res);
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
    const int* desc_src = source.desc();
    const std::complex<double> beta(0.0, 0.0);
    std::complex<double>* result_data = result.data();
    const int ic = 1;
    const int jc = 1;
    const int* desc_res = result.desc();
    ScalapackConnector::tranu(
        /* m= */ n,
        /* n= */ m,
        /* alpha= */ alpha,
        /* a= */ source_data,
        /* ia= */ ia,
        /* ja= */ ja,
        /* desca= */ desc_src,
        /* beta= */ beta,
        /* c= */ result_data,
        /* ic= */ ic,
        /* jc= */ jc,
        /* descc= */ desc_res);
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
    const int* desc_src = source.desc();
    const std::complex<double> beta(0.0, 0.0);
    std::complex<double>* result_data = result.data();
    const int ic = 1;
    const int jc = 1;
    const int* desc_res = result.desc();
    ScalapackConnector::tranc(
        /* m= */ n,
        /* n= */ m,
        /* alpha= */ alpha,
        /* a= */ source_data,
        /* ia= */ ia,
        /* ja= */ ja,
        /* desca= */ desc_src,
        /* beta= */ beta,
        /* c= */ result_data,
        /* ic= */ ic,
        /* jc= */ jc,
        /* descc= */ desc_res);
}

namespace
{
void throw_scalapack_error(const char* operation, const int info)
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
    const int* desc = matrix.desc();
    const int pivot_size = matrix.nrloc() + desc[4];
    std::vector<int> ipiv(static_cast<std::size_t>(pivot_size));
    int info = 0;
    ScalapackConnector::getrf(n,
                              n,
                              matrix.data(),
                              ia,
                              ja,
                              desc,
                              ipiv.data(),
                              &info);
    throw_scalapack_error("ScaLAPACK getrf", info);

    int lws = -1;
    int liws = -1;
    std::vector<double> ws(1);
    std::vector<int> iws(1);
    ScalapackConnector::getri(n,
                              matrix.data(),
                              ia,
                              ja,
                              desc,
                              ipiv.data(),
                              ws.data(),
                              &lws,
                              iws.data(),
                              &liws,
                              &info);
    throw_scalapack_error("ScaLAPACK getri workspace query", info);
    lws = queried_workspace_size(ws[0]);
    liws = iws[0];
    ws.resize(static_cast<std::size_t>(lws));
    iws.resize(static_cast<std::size_t>(liws));
    ScalapackConnector::getri(n,
                              matrix.data(),
                              ia,
                              ja,
                              desc,
                              ipiv.data(),
                              ws.data(),
                              &lws,
                              iws.data(),
                              &liws,
                              &info);
    throw_scalapack_error("ScaLAPACK getri", info);
}

void inv(Matrix2dView<std::complex<double>> matrix)
{
    const int n = matrix.nrglb();
    const int ia = 1;
    const int ja = 1;
    const int* desc = matrix.desc();
    const int pivot_size = matrix.nrloc() + desc[4];
    std::vector<int> ipiv(static_cast<std::size_t>(pivot_size));
    int info = 0;
    ScalapackConnector::getrf(n,
                              n,
                              matrix.data(),
                              ia,
                              ja,
                              desc,
                              ipiv.data(),
                              &info);
    throw_scalapack_error("ScaLAPACK getrf", info);

    int lws = -1;
    int liws = -1;
    std::vector<std::complex<double>> ws(1);
    std::vector<int> iws(1);
    ScalapackConnector::getri(n,
                              matrix.data(),
                              ia,
                              ja,
                              desc,
                              ipiv.data(),
                              ws.data(),
                              &lws,
                              iws.data(),
                              &liws,
                              &info);
    throw_scalapack_error("ScaLAPACK getri workspace query", info);
    lws = queried_workspace_size(ws[0]);
    liws = iws[0];
    ws.resize(static_cast<std::size_t>(lws));
    iws.resize(static_cast<std::size_t>(liws));
    ScalapackConnector::getri(n,
                              matrix.data(),
                              ia,
                              ja,
                              desc,
                              ipiv.data(),
                              ws.data(),
                              &lws,
                              iws.data(),
                              &liws,
                              &info);
    throw_scalapack_error("ScaLAPACK getri", info);
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
    throw_scalapack_error("ScaLAPACK gesv", info);
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
    throw_scalapack_error("ScaLAPACK gesv", info);
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
    const int* desc_lhs = left.desc(); // lhs = left operand
    const double* right_data = right.data();
    const int ib = 1;
    const int jb = 1;
    const int* desc_rhs = right.desc(); // rhs = right operand
    const double beta = 0.0;
    double* result_data = result.data();
    const int ic = 1;
    const int jc = 1;
    const int* desc_res = result.desc();
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
        /* DESCA= */ desc_lhs,
        /* B= */ right_data,
        /* IB= */ ib,
        /* JB= */ jb,
        /* DESCB= */ desc_rhs,
        /* beta= */ beta,
        /* C= */ result_data,
        /* IC= */ ic,
        /* JC= */ jc,
        /* DESCC= */ desc_res);
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
    const int* desc_lhs = left.desc();
    const std::complex<double>* right_data = right.data();
    const int ib = 1;
    const int jb = 1;
    const int* desc_rhs = right.desc();
    const std::complex<double> beta(0.0, 0.0);
    std::complex<double>* result_data = result.data();
    const int ic = 1;
    const int jc = 1;
    const int* desc_res = result.desc();
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
        /* DESCA= */ desc_lhs,
        /* B= */ right_data,
        /* IB= */ ib,
        /* JB= */ jb,
        /* DESCB= */ desc_rhs,
        /* beta= */ beta,
        /* C= */ result_data,
        /* IC= */ ic,
        /* JC= */ jc,
        /* DESCC= */ desc_res);
}

} // namespace scalapack
} // namespace matrix2d_impl
#endif
