#ifndef ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_LINALG_H
#define ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_LINALG_H

#include "source_base/matrix2d_view.h"

#include <complex>

namespace matrix2d_linalg
{
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

void linear_combination(Matrix2dView<double> result,
                        const double alpha,
                        Matrix2dView<const double> left,
                        const double beta,
                        Matrix2dView<const double> right);
void linear_combination(
    Matrix2dView<std::complex<double>> result,
    const std::complex<double> alpha,
    Matrix2dView<const std::complex<double>> left,
    const std::complex<double> beta,
    Matrix2dView<const std::complex<double>> right);

void scale(Matrix2dView<double> matrix, const double factor);
void scale(Matrix2dView<std::complex<double>> matrix,
           const std::complex<double> factor);

void set_identity(Matrix2dView<double> matrix);
void set_identity(Matrix2dView<std::complex<double>> matrix);
} // namespace matrix2d_linalg

#endif
