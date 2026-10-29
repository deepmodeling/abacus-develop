#include "source_hsolver/diago_cg.h"
#include "source_base/parallel_global.h"
#include "source_hsolver/hs_operator.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>

namespace
{
using Complex = std::complex<double>;
template <typename Scalar>
class DiagonalCGOperator : public hsolver::HSOperator<std::complex<Scalar>>
{
  public:
    using Value = std::complex<Scalar>;
    explicit DiagonalCGOperator(const std::vector<Scalar>& diagonal) : diagonal_(diagonal) {}
    void update_k(const int) override {}
    void hpsi(const Value* x, Value* hx, const int ld, const int columns) const override
    {
        for (int column = 0; column < columns; ++column)
        {
            for (int row = 0; row < static_cast<int>(diagonal_.size()); ++row)
            {
                hx[column * ld + row] = diagonal_[row] * x[column * ld + row];
            }
        }
    }
    void spsi(const Value* x, Value* sx, const int ld, const int columns) const override
    {
        const Value* end = x + ld * columns;
        std::copy(x, end, sx);
    }
  private:
    const std::vector<Scalar>& diagonal_;
};

class ShiftedTwoByTwoCGOperator : public hsolver::HSOperator<Complex>
{
  public:
    ShiftedTwoByTwoCGOperator(const double shift, const double coupling)
        : shift_(shift), coupling_(coupling) {}
    void update_k(const int) override {}
    void hpsi(const Complex* x, Complex* hx, const int ld, const int columns) const override
    {
        for (int column = 0; column < columns; ++column)
        {
            const int first = column * ld;
            const int second = first + 1;
            hx[first] = shift_ * x[first] + coupling_ * x[second];
            hx[second] = coupling_ * x[first] + shift_ * x[second];
        }
    }
    void spsi(const Complex* x, Complex* sx, const int ld, const int columns) const override
    {
        const Complex* end = x + ld * columns;
        std::copy(x, end, sx);
    }
  private:
    const double shift_;
    const double coupling_;
};

TEST(DiagoCGAccuracy, LargeCommonEnergyChoosesLowerRoot)
{
    const int dimension = 2;
    const int bands = 1;
    const double shift = std::ldexp(1.0, 30);
    const double coupling = std::ldexp(1.0, -30);
    ShiftedTwoByTwoCGOperator op(shift, coupling);
    std::vector<Complex> psi{Complex(1.0, 0.0), Complex(0.0, 0.0)};
    std::vector<double> eigen(bands, 0.0);
    std::vector<double> threshold(bands, 1.0e-20);
    std::vector<double> prec(dimension, 1.0);
#ifdef __MPI
    const hsolver::diag_comm_info comm(MPI_COMM_SELF, 0, 1);
#else
    const hsolver::diag_comm_info comm(0, 1);
#endif
    hsolver::DiagoCG<Complex> solver("pw", "scf", false, comm, 1.0e-20, 100);
    solver.diag(op, dimension, bands, dimension, psi.data(), eigen.data(), threshold, prec.data());
    // Both eigenvalues round to shift; the exact lower mode still has opposite components.
    const Complex upper_component = psi[0] + psi[1];
    const double error = std::abs(upper_component);
    ASSERT_TRUE(std::isfinite(error));
    const double norm = std::norm(psi[0]) + std::norm(psi[1]);
    ASSERT_TRUE(std::isfinite(norm));
    EXPECT_NEAR(norm, 1.0, 1.0e-12);
    EXPECT_DOUBLE_EQ(eigen[0], shift);

    EXPECT_LT(error, 1.0e-12);
}

class GeneralizedCGOperator : public DiagonalCGOperator<double>
{
  public:
    GeneralizedCGOperator(const std::vector<double>& h, const std::vector<double>& metric)
        : DiagonalCGOperator<double>(h), metric_(metric) {}
    void spsi(const Complex* x, Complex* sx, const int ld, const int columns) const override
    {
        for (int column = 0; column < columns; ++column)
        {
            for (int row = 0; row < static_cast<int>(metric_.size()); ++row)
            {
                sx[column * ld + row] = metric_[row] * x[column * ld + row];
            }
        }
    }
  private:
    const std::vector<double>& metric_;
};

TEST(DiagoCGAccuracy, GeneralizedPreconditionedNegativeRoot)
{
    const int dimension = 2;
    const int bands = 1;
    std::vector<double> h{-6.0, -3.0};
    std::vector<double> metric{2.0, 3.0};
    GeneralizedCGOperator op(h, metric);
    const double first = 1.0 / std::sqrt(2.0);
    const double second = 1.0e-9 / std::sqrt(6.0);
    std::vector<Complex> psi{Complex(first, 0.0), Complex(second, second)};
    std::vector<double> eigen(bands, 0.0);
    std::vector<double> threshold(bands, 1.0e-20);
    std::vector<double> prec{2.0, 10.0};
#ifdef __MPI
    const hsolver::diag_comm_info comm(MPI_COMM_SELF, 0, 1);
#else
    const hsolver::diag_comm_info comm(0, 1);
#endif
    hsolver::DiagoCG<Complex> solver("pw", "scf", false, comm, 1.0e-20, 100);
    solver.diag(op, dimension, bands, dimension, psi.data(), eigen.data(), threshold, prec.data());
    double residual_squared = 0.0;
    double norm = 0.0;
    for (int row = 0; row < dimension; ++row)
    {
        const Complex residual = (h[row] - eigen[0] * metric[row]) * psi[row];
        residual_squared += std::norm(residual);
        norm += metric[row] * std::norm(psi[row]);
    }
    const double residual = std::sqrt(residual_squared);
    ASSERT_TRUE(std::isfinite(residual));
    ASSERT_TRUE(std::isfinite(norm));
    EXPECT_NEAR(eigen[0], -3.0, 1.0e-14);
    EXPECT_NEAR(norm, 1.0, 1.0e-14);
    EXPECT_LT(residual, 1.0e-12);
}

TEST(DiagoCGAccuracy, FloatNearlyConvergedDirection)
{
    using FloatComplex = std::complex<float>;
    const int dimension = 2;
    const int bands = 1;
    std::vector<float> diagonal{1.0f, 2.0f};
    DiagonalCGOperator<float> op(diagonal);
    std::vector<FloatComplex> psi{FloatComplex(1.0f, 0.0f), FloatComplex(1.0e-4f, 0.0f)};
    std::vector<float> eigen(bands, 0.0f);
    std::vector<double> threshold(bands, 1.0e-12);
    std::vector<float> prec(dimension, 1.0f);
#ifdef __MPI
    const hsolver::diag_comm_info comm(MPI_COMM_SELF, 0, 1);
#else
    const hsolver::diag_comm_info comm(0, 1);
#endif
    hsolver::DiagoCG<FloatComplex> solver("pw", "scf", false, comm, 1.0e-12f, 100);
    solver.diag(op, dimension, bands, dimension, psi.data(), eigen.data(), threshold, prec.data());
    double residual_squared = 0.0;
    for (int row = 0; row < dimension; ++row)
    {
        const FloatComplex residual = (diagonal[row] - eigen[0]) * psi[row];
        residual_squared += std::norm(residual);
    }
    const double residual = std::sqrt(residual_squared);
    ASSERT_TRUE(std::isfinite(residual));
    const double norm = std::norm(psi[0]) + std::norm(psi[1]);
    ASSERT_TRUE(std::isfinite(norm));
    EXPECT_NEAR(norm, 1.0, 1.0e-6);
    EXPECT_NEAR(eigen[0], 1.0, 1.0e-6);

    EXPECT_LT(residual, 1.0e-6);
}

TEST(DiagoCGAccuracy, NearlyConvergedDirection)
{
    const int dimension = 2;
    const int bands = 1;
    std::vector<double> diagonal{1.0, 2.0};
    DiagonalCGOperator<double> op(diagonal);
    std::vector<Complex> psi{Complex(1.0, 0.0), Complex(1.0e-9, 0.0)};
    std::vector<double> eigen(bands, 0.0);
    std::vector<double> threshold(bands, 1.0e-20);
    std::vector<double> prec(dimension, 1.0);
#ifdef __MPI
    const hsolver::diag_comm_info comm(MPI_COMM_SELF, 0, 1);
#else
    const hsolver::diag_comm_info comm(0, 1);
#endif
    hsolver::DiagoCG<Complex> solver("pw", "scf", false, comm, 1.0e-20, 100);
    solver.diag(op, dimension, bands, dimension, psi.data(), eigen.data(), threshold, prec.data());
    double residual_squared = 0.0;
    for (int row = 0; row < dimension; ++row)
    {
        const Complex residual = (diagonal[row] - eigen[0]) * psi[row];
        residual_squared += std::norm(residual);
    }
    const double residual = std::sqrt(residual_squared);
    std::cout << std::setprecision(17) << "CG_PROBE eigen=" << eigen[0] << " residual=" << residual << std::endl;
    ASSERT_TRUE(std::isfinite(residual));
    const double norm = std::norm(psi[0]) + std::norm(psi[1]);
    ASSERT_TRUE(std::isfinite(norm));
    EXPECT_NEAR(norm, 1.0, 1.0e-12);
    EXPECT_NEAR(eigen[0], 1.0, 1.0e-14);

    EXPECT_LT(residual, 1.0e-12);
}
TEST(DiagoCGAccuracy, TightMultiBandOrthogonality)
{
    const int dimension = 48;
    const int bands = 6;
    std::vector<double> diagonal(dimension);
    for (int row = 0; row < dimension; ++row)
    {
        diagonal[row] = 1.0 + 0.1 * row;
    }
    DiagonalCGOperator<double> op(diagonal);
    std::vector<Complex> psi(dimension * bands);
    for (int column = 0; column < bands; ++column)
    {
        for (int row = 0; row < dimension; ++row)
        {
            const double real_angle = (row + 1) * (column + 2) * 0.75;
            const double imag_angle = (row + 3) * (column + 1) * 0.37;
            const double a = std::sin(real_angle);
            const double b = std::cos(imag_angle);
            psi[row + dimension * column] = Complex(a, b);
        }
    }
    std::vector<double> eigen(bands, 0.0);
    std::vector<double> threshold(bands, 1.0e-20);
    std::vector<double> prec(dimension, 1.0);
#ifdef __MPI
    const hsolver::diag_comm_info comm(MPI_COMM_SELF, 0, 1);
#else
    const hsolver::diag_comm_info comm(0, 1);
#endif
    hsolver::DiagoCG<Complex> solver("pw", "scf", true, comm, 1.0e-20, 500);
    solver.diag(op, dimension, bands, dimension, psi.data(), eigen.data(), threshold, prec.data());
    double max_residual = 0.0;
    double max_gram = 0.0;
    for (int column = 0; column < bands; ++column)
    {
        double residual_squared = 0.0;
        for (int row = 0; row < dimension; ++row)
        {
            const Complex residual = (diagonal[row] - eigen[column]) * psi[row + column * dimension];
            residual_squared += std::norm(residual);
        }
        const double residual = std::sqrt(residual_squared);
        ASSERT_TRUE(std::isfinite(residual));
        max_residual = std::max(max_residual, residual);
        const double eigen_error = std::abs(eigen[column] - diagonal[column]);
        EXPECT_LT(eigen_error, 1.0e-8);
        for (int other = 0; other < bands; ++other)
        {
            Complex gram(0.0, 0.0);
            for (int row = 0; row < dimension; ++row)
            {
                gram += std::conj(psi[row + column * dimension]) * psi[row + other * dimension];
            }
            if (column == other) { gram -= Complex(1.0, 0.0); }
            const double error = std::abs(gram);
            ASSERT_TRUE(std::isfinite(error));
            max_gram = std::max(max_gram, error);
        }
    }
    std::cout << std::setprecision(17) << "CG_MULTI residual=" << max_residual << " gram=" << max_gram << std::endl;
    EXPECT_LT(max_residual, 1.0e-8);
    EXPECT_LT(max_gram, 1.0e-10);
}

}

int main(int argc, char** argv)
{
    int nproc = 1;
    int threads = 1;
    int rank = 0;
    Parallel_Global::read_pal_param(argc, argv, nproc, threads, rank);
#ifdef __MPI
    POOL_WORLD = MPI_COMM_WORLD;
    KP_WORLD = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
#endif
    if (nproc != 1)
    {
        std::cerr << "CG accuracy tests require one process." << std::endl;
#ifdef __MPI
        Parallel_Global::finalize_mpi();
#endif
        return 2;
    }
    testing::InitGoogleTest(&argc, argv);
    const int status = RUN_ALL_TESTS();
    std::cout << "CG_ACCURACY_COMPLETE status=" << status << std::endl;
#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif
    return status;
}
