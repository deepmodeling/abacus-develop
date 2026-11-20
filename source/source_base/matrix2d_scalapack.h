#ifndef ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_SCALAPACK_H
#define ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_SCALAPACK_H

#include <complex>

template <typename T>
class Matrix2dView;

#ifdef __MPI
namespace matrix2d_impl
{
namespace scalapack
{
void transpose(Matrix2dView<double> result,
               Matrix2dView<const double> source);
void transpose(Matrix2dView<std::complex<double>> result,
               Matrix2dView<const std::complex<double>> source);
void conjugate_transpose(Matrix2dView<double> result,
                         Matrix2dView<const double> source);
void conjugate_transpose(
    Matrix2dView<std::complex<double>> result,
    Matrix2dView<const std::complex<double>> source);
void inv(Matrix2dView<double> matrix);
void inv(Matrix2dView<std::complex<double>> matrix);
void solve(Matrix2dView<double> matrix,
           Matrix2dView<double> rhs);
void solve(Matrix2dView<std::complex<double>> matrix,
           Matrix2dView<std::complex<double>> rhs);
void matmul(Matrix2dView<double> result,
            Matrix2dView<const double> left,
            Matrix2dView<const double> right);
void matmul(Matrix2dView<std::complex<double>> result,
            Matrix2dView<const std::complex<double>> left,
            Matrix2dView<const std::complex<double>> right);

} // namespace scalapack
} // namespace matrix2d_impl
#endif

#endif
