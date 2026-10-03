#include "source_base/parallel_device.h"
#include "source_base/parallel_global.h"
#include "source_hsolver/hsolver_linear.h"
#include "source_hsolver/test/linear_test_utils.h"

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <type_traits>

namespace
{
template <typename T>
class IdentityOperator final : public hsolver::LinearOperator<T>
{
  private:
    const int dim_;

  public:
    explicit IdentityOperator(const int dim) : dim_(dim)
    {
    }
    bool is_identity() const override
    {
        return true;
    }
    void apply(const T* x, T* y, const int ld, const int nvec) const override
    {
        for (int band = 0; band < nvec; ++band)
        {
            for (int i = 0; i < dim_; ++i)
            {
                y[band * ld + i] = x[band * ld + i];
            }
        }
    }
};

template <typename T>
class DenseOperator final : public hsolver::LinearOperator<T>
{
  public:
    DenseOperator(const std::vector<T>& matrix, const int n, const int start, const int dim, const hsolver::diag_comm_info& comm)
        : matrix_(matrix), n_(n), start_(start), dim_(dim), comm_(comm)
    {
    }

    void apply(const T* x, T* y, const int ld, const int nvec) const override
    {
        std::vector<T> full(n_ * nvec, T(0));
        for (int band = 0; band < nvec; ++band)
        {
            for (int i = 0; i < dim_; ++i)
            {
                full[band * n_ + start_ + i] = x[band * ld + i];
            }
        }
#ifdef __MPI
        const int count = full.size();
        Parallel_Common::reduce_data(full.data(), count, comm_.comm);
#endif
        const std::vector<T> product = linear_test::multiply(matrix_, full, n_, nvec);
        for (int band = 0; band < nvec; ++band)
        {
            for (int i = 0; i < dim_; ++i)
            {
                y[band * ld + i] = product[band * n_ + start_ + i];
            }
        }
    }

  private:
    const std::vector<T>& matrix_;
    const int n_;
    const int start_;
    const int dim_;
    const hsolver::diag_comm_info comm_;
};

template <typename T>
class RightPreconditioner final : public hsolver::LinearOperator<T>
{
  public:
    RightPreconditioner(const std::vector<T>& matrix, const int n, const int start, const int dim)
        : matrix_(matrix), n_(n), start_(start), dim_(dim)
    {
    }
    void apply(const T* x, T* y, const int ld, const int nvec) const override
    {
        for (int band = 0; band < nvec; ++band)
        {
            for (int i = 0; i < dim_; ++i)
            {
                const int row = start_ + i;
                y[band * ld + i] = x[band * ld + i] / matrix_[row * n_ + row];
            }
        }
    }

  private:
    const std::vector<T>& matrix_;
    const int n_;
    const int start_;
    const int dim_;
};

template <typename Real>
class LinearSolveTest : public testing::Test
{
};
using Precisions = testing::Types<double, float>;
TYPED_TEST_SUITE(LinearSolveTest, Precisions);

TYPED_TEST(LinearSolveTest, DenseReferenceAndReusedBatches)
{
    using T = std::complex<TypeParam>;
    const hsolver::diag_comm_info comm = linear_test::world_comm();
    const double tolerance = std::is_same<TypeParam, double>::value ? 1e-12 : 2e-6;
    const double solution_tol = std::is_same<TypeParam, double>::value ? 1e-10 : 1e-4;
    for (const hsolver::LinearMethod method: {hsolver::LinearMethod::bicgstab, hsolver::LinearMethod::cgs})
    {
        for (const bool preconditioned: {false, true})
        {
            hsolver::LinearSolveOptions options;
            options.method = method;
            options.tolerance = tolerance;
            hsolver::HSolverLinear<T> solver(options, comm);
            // Shrink, grow, and use zero local rows without replacing the solver.
            for (const int n: {11, 5, 17, 1})
            {
                const int start = n * comm.rank / comm.nproc;
                const int end = n * (comm.rank + 1) / comm.nproc;
                const int dim = end - start;
                const int ld = dim + 3;
                const int nvec = n == 5 ? 3 : 5;
                for (const bool cn_matrix: {true, false})
                {
                    SCOPED_TRACE(testing::Message() << "n=" << n << " method=" << static_cast<int>(method)
                                                    << " preconditioned=" << preconditioned << " cn=" << cn_matrix);
                    std::vector<T> a = linear_test::hermitian<T>(n);
                    for (int j = 0; j < n; ++j)
                    {
                        for (int i = 0; i < n; ++i)
                        {
                            if (cn_matrix)
                            {
                                a[j * n + i] *= T(0, 0.35);
                                if (i == j)
                                {
                                    a[j * n + i] += T(1);
                                }
                            }
                            else if (i > j)
                            {
                                a[j * n + i] *= T(0.3, 0.6);
                            }
                        }
                    }
                    std::vector<T> exact(n * nvec, T(0));
                    for (int band = 0; band < nvec; ++band)
                    {
                        for (int i = 0; i < n; ++i)
                        {
                            exact[band * n + i] = T(0.2 + 0.03 * i, 0.1 * (band + 1) / (i + 1));
                        }
                    }
                    // Middle columns retire before the first: exercise swaps and restoration.
                    std::fill(exact.begin() + n, exact.begin() + 2 * n, T(0));
                    std::fill(exact.begin() + 2 * n, exact.begin() + 3 * n, T(0));
                    exact[2 * n] = T(0.7, -0.2);
                    const std::vector<T> b = linear_test::multiply(a, exact, n, nvec);
                    const std::vector<linear_test::Complex> reference = linear_test::lapack_solve(a, b, n, nvec);
                    const T sentinel(123, -45);
                    std::vector<T> rhs(ld * nvec, sentinel);
                    std::vector<T> x(ld * nvec, sentinel);
                    for (int band = 0; band < nvec; ++band)
                    {
                        for (int i = 0; i < dim; ++i)
                        {
                            rhs[band * ld + i] = b[band * n + start + i];
                            x[band * ld + i] = band == 3 ? exact[band * n + start + i] : T(0);
                        }
                    }
                    const std::vector<T> rhs_before = rhs;
                    DenseOperator<T> op(a, n, start, dim, comm);
                    RightPreconditioner<T> precond(a, n, start, dim);
                    hsolver::LinearSolveResult result;
                    if (preconditioned)
                    {
                        result = solver.solve(op, precond, ld, nvec, dim, x.data(), rhs.data());
                    }
                    else
                    {
                        const IdentityOperator<T> identity(dim);
                        result = solver.solve(op, identity, ld, nvec, dim, x.data(), rhs.data());
                    }
                    EXPECT_EQ(result.status, hsolver::LinearSolveStatus::converged);
                    EXPECT_EQ(result.failed_band, -1);
                    EXPECT_EQ(rhs, rhs_before);
                    std::vector<linear_test::Complex> full(n * nvec, 0.0);
                    for (int band = 0; band < nvec; ++band)
                    {
                        for (int i = 0; i < dim; ++i)
                        {
                            full[band * n + start + i] = x[band * ld + i];
                            EXPECT_LT(std::abs(full[band * n + start + i] - reference[band * n + start + i]), solution_tol);
                        }
                        for (int i = dim; i < ld; ++i)
                        {
                            EXPECT_EQ(x[band * ld + i], sentinel);
                        }
                    }
#ifdef __MPI
                    const int count = full.size();
                    Parallel_Common::reduce_data(full.data(), count, comm.comm);
#endif
                    const std::vector<linear_test::Complex> a_double(a.begin(), a.end());
                    const std::vector<linear_test::Complex> ax = linear_test::multiply(a_double, full, n, nvec);
                    double max_residual = 0.0;
                    for (int band = 0; band < nvec; ++band)
                    {
                        double residual2 = 0.0;
                        double bnorm2 = 0.0;
                        for (int i = 0; i < n; ++i)
                        {
                            const linear_test::Complex value = b[band * n + i];
                            residual2 += std::norm(ax[band * n + i] - value);
                            bnorm2 += std::norm(value);
                        }
                        double scale = std::sqrt(bnorm2);
                        if (method == hsolver::LinearMethod::bicgstab)
                        {
                            scale = std::max(1.0, scale);
                        }
                        else if (scale == 0.0)
                        {
                            scale = 1.0;
                        }
                        const double residual = std::sqrt(residual2);
                        // Allow rounding in the independent double-precision residual of a float solve.
                        EXPECT_LE(residual, 1.2 * tolerance * scale);
                        max_residual = std::max(max_residual, residual);
                    }
                    EXPECT_NEAR(result.max_residual, max_residual, tolerance);
                }
            }
        }
    }
}

TEST(LinearSolveFailure, ReportsOriginalColumnAndTrueResidual)
{
    using T = linear_test::Complex;
    const hsolver::diag_comm_info comm = linear_test::world_comm();
    const int n = 7;
    const int start = n * comm.rank / comm.nproc;
    const int dim = n * (comm.rank + 1) / comm.nproc - start;
    const int ld = dim + 1;
    for (const hsolver::LinearMethod method: {hsolver::LinearMethod::bicgstab, hsolver::LinearMethod::cgs})
    {
        for (const bool singular: {false, true})
        {
            std::vector<T> a = linear_test::hermitian<T>(n);
            if (singular)
            {
                std::fill(a.begin(), a.end(), T(0));
            }
            DenseOperator<T> op(a, n, start, dim, comm);
            hsolver::LinearSolveOptions options;
            options.method = method;
            options.tolerance = 1e-13;
            options.max_iterations = singular ? 10 : 1;
            hsolver::HSolverLinear<T> solver(options, comm);
            std::vector<T> x(3 * ld, T(0));
            std::vector<T> b(3 * ld, T(0));
            // Retire column 0; column 2 is swapped forward but failure must still name column 1.
            for (int i = 0; i < dim; ++i)
            {
                b[ld + i] = T(1, 0.1 * (start + i));
            }
            const IdentityOperator<T> identity(dim);
            const hsolver::LinearSolveResult result = solver.solve(op, identity, ld, 3, dim, x.data(), b.data());
            const hsolver::LinearSolveStatus expected
                = singular ? hsolver::LinearSolveStatus::breakdown : hsolver::LinearSolveStatus::max_iterations;
            EXPECT_EQ(result.status, expected);
            EXPECT_EQ(result.failed_band, 1);
            std::vector<T> ax(3 * ld, T(0));
            op.apply(x.data(), ax.data(), ld, 3);
            double residual2 = 0.0;
            for (int i = 0; i < dim; ++i)
            {
                residual2 += std::norm(ax[ld + i] - b[ld + i]);
            }
#ifdef __MPI
            Parallel_Common::reduce_data(&residual2, 1, comm.comm);
#endif
            EXPECT_NEAR(result.max_residual, std::sqrt(residual2), 1e-12);
            EXPECT_GT(result.max_residual, options.tolerance);
        }
    }
}
} // namespace

int main(int argc, char** argv)
{
    int nproc = 1;
    int nthread = 1;
    int rank = 0;
    Parallel_Global::read_pal_param(argc, argv, nproc, nthread, rank);
#ifdef __MPI
    // This numerical test creates no application pools; the shared finalizer owns these handles.
    POOL_WORLD = MPI_COMM_NULL;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif
    return result;
}
