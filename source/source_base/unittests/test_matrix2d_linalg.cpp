#include "source_base/matrix2d_linalg.h"
#include "source_base/matrix2d.h"

#ifdef __MPI
#include "source_base/parallel_common.h"
#endif

#include "gtest/gtest.h"

#include <cmath>
#include <complex>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace
{
void expect_near(const double actual, const double expected)
{
    EXPECT_NEAR(actual, expected, 1.0e-12);
}

void expect_near(const std::complex<double>& actual,
                 const std::complex<double>& expected)
{
    EXPECT_NEAR(actual.real(), expected.real(), 1.0e-12);
    EXPECT_NEAR(actual.imag(), expected.imag(), 1.0e-12);
}

template <typename T>
void expect_host_matrices_near(const Matrix2d<T>& actual,
                               const Matrix2d<T>& expected)
{
    ASSERT_EQ(actual.nrglb(), expected.nrglb());
    ASSERT_EQ(actual.ncglb(), expected.ncglb());
    for (int i = 0; i < expected.nrglb(); ++i)
    {
        for (int j = 0; j < expected.ncglb(); ++j)
        {
            expect_near(actual.atloc(i, j), expected.atloc(i, j));
        }
    }
}

Matrix2d<double> make_real_system()
{
    Matrix2d<double> matrix(2, 2);
    matrix.atloc(0, 0) = 4.0;
    matrix.atloc(0, 1) = 1.0;
    matrix.atloc(1, 0) = 2.0;
    matrix.atloc(1, 1) = 3.0;
    return matrix;
}

Matrix2d<double> make_real_inverse_reference()
{
    Matrix2d<double> inverse(2, 2);
    inverse.atloc(0, 0) = 0.3;
    inverse.atloc(0, 1) = -0.1;
    inverse.atloc(1, 0) = -0.2;
    inverse.atloc(1, 1) = 0.4;
    return inverse;
}

Matrix2d<std::complex<double>> make_complex_system()
{
    Matrix2d<std::complex<double>> matrix(2, 2);
    matrix.atloc(0, 0) = std::complex<double>(1.0, 1.0);
    matrix.atloc(0, 1) = std::complex<double>(2.0, 0.0);
    matrix.atloc(1, 0) = std::complex<double>(3.0, 0.0);
    matrix.atloc(1, 1) = std::complex<double>(4.0, -1.0);
    return matrix;
}

Matrix2d<std::complex<double>> make_complex_solution()
{
    Matrix2d<std::complex<double>> solution(2, 2);
    solution.atloc(0, 0) = std::complex<double>(1.0, 0.0);
    solution.atloc(0, 1) = std::complex<double>(0.0, 2.0);
    solution.atloc(1, 0) = std::complex<double>(3.0, 0.0);
    solution.atloc(1, 1) = std::complex<double>(4.0, 0.0);
    return solution;
}

#ifdef __MPI
std::vector<std::pair<int, int>> available_process_grids()
{
    const int nproc = Parallel_Common::comm_size(MPI_COMM_WORLD);
    std::vector<std::pair<int, int>> grids;
    if (nproc == 1)
    {
        grids.push_back(std::make_pair(1, 1));
    }
    if (nproc == 2)
    {
        grids.push_back(std::make_pair(1, 2));
        grids.push_back(std::make_pair(2, 1));
    }
    if (nproc == 4)
    {
        grids.push_back(std::make_pair(1, 4));
        grids.push_back(std::make_pair(2, 2));
        grids.push_back(std::make_pair(4, 1));
    }
    return grids;
}

Matrix2d<double> make_distributed_real_system(const int size)
{
    Matrix2d<double> matrix(size, size);
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            if (i == j)
            {
                matrix.atloc(i, j) = 6.0 + i;
            }
            else if (std::abs(i - j) == 1)
            {
                matrix.atloc(i, j) = 0.5;
            }
            else
            {
                matrix.atloc(i, j) = 0.0;
            }
        }
    }
    return matrix;
}

Matrix2d<double> make_distributed_real_solution(const int nr,
                                                const int nc)
{
    Matrix2d<double> solution(nr, nc);
    for (int i = 0; i < nr; ++i)
    {
        for (int j = 0; j < nc; ++j)
        {
            solution.atloc(i, j)
                = static_cast<double>(1 + i + 2 * j);
        }
    }
    return solution;
}

Matrix2d<std::complex<double>> make_distributed_complex_system(const int size)
{
    Matrix2d<std::complex<double>> matrix(size, size);
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            if (i == j)
            {
                matrix.atloc(i, j)
                    = std::complex<double>(6.0 + i, 0.5);
            }
            else if (j == i + 1)
            {
                matrix.atloc(i, j)
                    = std::complex<double>(0.5, 0.25);
            }
            else if (j + 1 == i)
            {
                matrix.atloc(i, j)
                    = std::complex<double>(0.5, -0.25);
            }
            else
            {
                matrix.atloc(i, j) = std::complex<double>(0.0, 0.0);
            }
        }
    }
    return matrix;
}

Matrix2d<std::complex<double>> make_distributed_complex_solution(
    const int nr,
    const int nc)
{
    Matrix2d<std::complex<double>> solution(nr, nc);
    for (int i = 0; i < nr; ++i)
    {
        for (int j = 0; j < nc; ++j)
        {
            solution.atloc(i, j)
                = std::complex<double>(1 + i + 2 * j, i - j);
        }
    }
    return solution;
}

template <typename T>
void check_distributed_solve(const Matrix2d<T>& host_matrix,
                             const Matrix2d<T>& host_solution,
                             const int nrproc,
                             const int ncproc,
                             const int nrblk)
{
    const Matrix2d<T> host_rhs = host_matrix.matmul(host_solution);
    Matrix2d<T> matrix(host_matrix.nrglb(),
                       host_matrix.ncglb(),
                       nrproc,
                       ncproc,
                       nrblk,
                       MPI_COMM_WORLD);
    const std::shared_ptr<const Parallel_Blkcyclic> solution_distrib
        = make_borrowed_distribution(host_solution.nrglb(),
                                     host_solution.ncglb(),
                                     nrblk,
                                     matrix.ictxt());
    const std::shared_ptr<const Parallel_Blkcyclic> rhs_distrib
        = make_borrowed_distribution(host_rhs.nrglb(),
                                     host_rhs.ncglb(),
                                     nrblk,
                                     matrix.ictxt());
    Matrix2d<T> solution(solution_distrib);
    Matrix2d<T> rhs(rhs_distrib);
    matrix.scatter(host_matrix);
    solution.scatter(host_solution);
    const Matrix2d<T>& matrix_input = matrix;
    const Matrix2d<T>& solution_input = solution;
    matrix2d_linalg::matmul(rhs.view(),
                            matrix_input.view(),
                            solution_input.view());

    matrix2d_linalg::solve(matrix.view(), rhs.view());

    expect_host_matrices_near(rhs.gather(), host_solution);
}

template <typename T>
void check_distributed_inverse(const Matrix2d<T>& host_matrix,
                               const int nrproc,
                               const int ncproc,
                               const int nrblk)
{
    Matrix2d<T> expected = host_matrix;
    matrix2d_linalg::inv(expected.view());

    Matrix2d<T> matrix(host_matrix.nrglb(),
                       host_matrix.ncglb(),
                       nrproc,
                       ncproc,
                       nrblk,
                       MPI_COMM_WORLD);
    matrix.scatter(host_matrix);
    matrix2d_linalg::inv(matrix.view());

    expect_host_matrices_near(matrix.gather(), expected);
}
#endif
} // namespace

TEST(Matrix2dLinalgTest, HostRealInverseAndSolve)
{
    Matrix2d<double> inverse = make_real_system();
    matrix2d_linalg::inv(inverse.view());
    expect_host_matrices_near(inverse, make_real_inverse_reference());

    Matrix2d<double> rhs(2, 2);
    rhs.atloc(0, 0) = 1.0;
    rhs.atloc(0, 1) = 2.0;
    rhs.atloc(1, 0) = 3.0;
    rhs.atloc(1, 1) = 4.0;
    const Matrix2d<double> expected = make_real_inverse_reference().matmul(rhs);
    Matrix2d<double> matrix = make_real_system();
    matrix2d_linalg::solve(matrix.view(), rhs.view());
    expect_host_matrices_near(rhs, expected);
}

TEST(Matrix2dLinalgTest, ExplicitInverseSatisfiesInverseIdentity)
{
    const Matrix2d<double> real_matrix = make_real_system();
    Matrix2d<double> real_inverse = real_matrix;
    Matrix2dView<double> real_inverse_view = real_inverse.view();
    matrix2d_linalg::inv(real_inverse_view);
    Matrix2d<double> real_identity(2, 2);
    Matrix2dView<double> real_identity_view = real_identity.view();
    matrix2d_linalg::set_identity(real_identity_view);
    const Matrix2d<double> real_product = real_matrix.matmul(real_inverse);
    expect_host_matrices_near(real_product, real_identity);

    const Matrix2d<std::complex<double>> complex_matrix
        = make_complex_system();
    Matrix2d<std::complex<double>> complex_inverse = complex_matrix;
    Matrix2dView<std::complex<double>> complex_inverse_view
        = complex_inverse.view();
    matrix2d_linalg::inv(complex_inverse_view);
    Matrix2d<std::complex<double>> complex_identity(2, 2);
    Matrix2dView<std::complex<double>> complex_identity_view
        = complex_identity.view();
    matrix2d_linalg::set_identity(complex_identity_view);
    const Matrix2d<std::complex<double>> complex_product
        = complex_matrix.matmul(complex_inverse);
    expect_host_matrices_near(complex_product, complex_identity);
}

TEST(Matrix2dLinalgTest, HostComplexSolve)
{
    Matrix2d<std::complex<double>> matrix = make_complex_system();
    const Matrix2d<std::complex<double>> expected = make_complex_solution();
    Matrix2d<std::complex<double>> rhs = matrix.matmul(expected);

    matrix2d_linalg::solve(matrix.view(), rhs.view());

    expect_host_matrices_near(rhs, expected);

    Matrix2d<std::complex<double>> inverse = make_complex_system();
    matrix2d_linalg::inv(inverse.view());
    Matrix2d<std::complex<double>> inverse_reference(2, 2);
    inverse_reference.atloc(0, 0) = std::complex<double>(-0.7, -1.1);
    inverse_reference.atloc(0, 1) = std::complex<double>(0.2, 0.6);
    inverse_reference.atloc(1, 0) = std::complex<double>(0.3, 0.9);
    inverse_reference.atloc(1, 1) = std::complex<double>(0.2, -0.4);
    expect_host_matrices_near(inverse, inverse_reference);
}

TEST(Matrix2dLinalgTest, SerialDistributionUsesLocalBackend)
{
    const std::shared_ptr<const Parallel_Blkcyclic> distrib
        = make_frozen_distribution(2, 2);

    Matrix2d<double> matrix(distrib);
    matrix.atloc(0, 0) = 4.0;
    matrix.atloc(0, 1) = 1.0;
    matrix.atloc(1, 0) = 2.0;
    matrix.atloc(1, 1) = 3.0;

    Matrix2d<double> rhs(distrib);
    rhs.atloc(0, 0) = 1.0;
    rhs.atloc(0, 1) = 2.0;
    rhs.atloc(1, 0) = 3.0;
    rhs.atloc(1, 1) = 4.0;
    const Matrix2d<double> expected = make_real_inverse_reference().matmul(rhs);

    matrix2d_linalg::solve(matrix.view(), rhs.view());
    expect_host_matrices_near(rhs, expected);

    Matrix2d<double> inverse(distrib);
    inverse.atloc(0, 0) = 4.0;
    inverse.atloc(0, 1) = 1.0;
    inverse.atloc(1, 0) = 2.0;
    inverse.atloc(1, 1) = 3.0;
    matrix2d_linalg::inv(inverse.view());
    expect_host_matrices_near(inverse, make_real_inverse_reference());
}

TEST(Matrix2dLinalgTest, ExternalBufferViewDoesNotCopy)
{
    std::vector<double> values = {4.0, 1.0, 2.0, 3.0};
    Matrix2dView<double> matrix(values.data(), 2, 2);
    Matrix2d<double> rhs(2, 1);
    rhs.atloc(0, 0) = 1.0;
    rhs.atloc(1, 0) = 2.0;
    matrix2d_linalg::solve(matrix, rhs.view());

    expect_near(rhs.atloc(0, 0), 0.1);
    expect_near(rhs.atloc(1, 0), 0.6);
    EXPECT_EQ(matrix.data(), values.data());
}

TEST(Matrix2dLinalgTest, HostViewMatmulAndAxpby)
{
    const Matrix2d<double> left = make_real_system();
    const Matrix2d<double> right = make_real_inverse_reference();
    Matrix2d<double> product(2, 2);
    matrix2d_linalg::matmul(product.view(), left.view(), right.view());
    expect_host_matrices_near(product, left.matmul(right));

    Matrix2d<double> combined(2, 2);
    matrix2d_linalg::axpby(
        combined.view(), 2.0, left.view(), -0.5, right.view());
    Matrix2d<double> combined_reference(2, 2);
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            combined_reference.atloc(i, j)
                = 2.0 * left.atloc(i, j)
                    - 0.5 * right.atloc(i, j);
        }
    }
    expect_host_matrices_near(combined, combined_reference);

    Matrix2d<double> in_place = left;
    const Matrix2d<double> in_place_original = in_place;
    const Matrix2d<double>& in_place_input = in_place;
    matrix2d_linalg::axpby(in_place.view(),
                           2.0,
                           in_place_input.view(),
                           -0.5,
                           right.view());
    Matrix2d<double> in_place_reference(2, 2);
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            in_place_reference.atloc(i, j)
                = 2.0 * in_place_original.atloc(i, j)
                    - 0.5 * right.atloc(i, j);
        }
    }
    expect_host_matrices_near(in_place, in_place_reference);

    matrix2d_linalg::scale(combined.view(), 3.0);
    matrix2d_linalg::set_identity(product.view());
    const Matrix2d<double> generated_identity
        = matrix2d_linalg::eye<double>(2);
    const Matrix2d<std::complex<double>> generated_complex_identity
        = matrix2d_linalg::eye<std::complex<double>>(2);
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            const double identity_value = i == j ? 1.0 : 0.0;
            expect_near(combined.atloc(i, j),
                        3.0 * combined_reference.atloc(i, j));
            expect_near(product.atloc(i, j), identity_value);
            expect_near(generated_identity.atloc(i, j), identity_value);
            expect_near(generated_complex_identity.atloc(i, j).real(),
                        identity_value);
            expect_near(generated_complex_identity.atloc(i, j).imag(), 0.0);
        }
    }
}

TEST(Matrix2dLinalgTest, RejectsSingularAndMismatchedMatrices)
{
    Matrix2d<double> singular(2, 2);
    singular.atloc(0, 0) = 1.0;
    singular.atloc(0, 1) = 2.0;
    singular.atloc(1, 0) = 2.0;
    singular.atloc(1, 1) = 4.0;
    EXPECT_THROW(matrix2d_linalg::inv(singular.view()), std::runtime_error);

    Matrix2d<double> matrix = make_real_system();
    Matrix2d<double> wrong_rhs(1, 1);
    EXPECT_THROW(matrix2d_linalg::solve(matrix.view(), wrong_rhs.view()),
                 std::invalid_argument);
}

TEST(Matrix2dLinalgTest, RejectsSerialDistributionAgainstHostLayout)
{
    const std::shared_ptr<const Parallel_Blkcyclic> distrib
        = make_frozen_distribution(2, 2);
    Matrix2d<double> distributed(distrib);
    Matrix2d<double> host(2, 2);

    ASSERT_TRUE(distributed.has_distribution());
    ASSERT_FALSE(distributed.is_active());
    EXPECT_THROW(
        matrix2d_linalg::axpby(host.view(),
                               1.0,
                               distributed.view(),
                               0.0,
                               host.view()),
        std::invalid_argument);
}

#ifdef __MPI
TEST(Matrix2dLinalgTest, DistributedRealAndComplexSolve)
{
    const int size = 5;
    const int nrhs = 3;
    const Matrix2d<double> real_matrix = make_distributed_real_system(size);
    const Matrix2d<double> real_solution
        = make_distributed_real_solution(size, nrhs);
    const Matrix2d<std::complex<double>> complex_matrix
        = make_distributed_complex_system(size);
    const Matrix2d<std::complex<double>> complex_solution
        = make_distributed_complex_solution(size, nrhs);
    const std::vector<std::pair<int, int>> grids = available_process_grids();
    const int nrblks[] = {2, 3};
    for (std::size_t igrid = 0; igrid < grids.size(); ++igrid)
    {
        for (std::size_t iblk = 0; iblk < 2; ++iblk)
        {
            check_distributed_solve(real_matrix,
                                    real_solution,
                                    grids[igrid].first,
                                    grids[igrid].second,
                                    nrblks[iblk]);
            check_distributed_solve(complex_matrix,
                                    complex_solution,
                                    grids[igrid].first,
                                    grids[igrid].second,
                                    nrblks[iblk]);
        }
    }
}

TEST(Matrix2dLinalgTest, DistributedRealAndComplexInverse)
{
    const int size = 5;
    const Matrix2d<double> real_matrix = make_distributed_real_system(size);
    const Matrix2d<std::complex<double>> complex_matrix
        = make_distributed_complex_system(size);
    const std::vector<std::pair<int, int>> grids = available_process_grids();
    const int nrblks[] = {2, 3};
    for (std::size_t igrid = 0; igrid < grids.size(); ++igrid)
    {
        for (std::size_t iblk = 0; iblk < 2; ++iblk)
        {
            check_distributed_inverse(real_matrix,
                                      grids[igrid].first,
                                      grids[igrid].second,
                                      nrblks[iblk]);
            check_distributed_inverse(complex_matrix,
                                      grids[igrid].first,
                                      grids[igrid].second,
                                      nrblks[iblk]);
        }
    }
}

TEST(Matrix2dLinalgTest, DistributedIdentityAndAxpby)
{
    const int size = 5;
    const int nrblks[] = {2, 3};
    const std::vector<std::pair<int, int>> grids = available_process_grids();
    for (std::size_t igrid = 0; igrid < grids.size(); ++igrid)
    {
        for (std::size_t iblk = 0; iblk < 2; ++iblk)
        {
            Matrix2d<double> identity(size,
                                      size,
                                      grids[igrid].first,
                                      grids[igrid].second,
                                      nrblks[iblk],
                                      MPI_COMM_WORLD);
            matrix2d_linalg::set_identity(identity.view());
            Matrix2d<double> identity_expected(size, size);
            matrix2d_linalg::set_identity(identity_expected.view());
            expect_host_matrices_near(identity.gather(), identity_expected);

            const Matrix2d<double> host_left = make_distributed_real_system(size);
            Matrix2d<double> left(identity.distribution());
            Matrix2d<double> combined(identity.distribution());
            left.scatter(host_left);
            matrix2d_linalg::axpby(combined.view(),
                                   2.0,
                                   left.view(),
                                   0.5,
                                   identity.view());
            Matrix2d<double> combined_expected(size, size);
            for (int i = 0; i < size; ++i)
            {
                for (int j = 0; j < size; ++j)
                {
                    combined_expected.atloc(i, j)
                        = 2.0 * host_left.atloc(i, j)
                            + 0.5 * identity_expected.atloc(i, j);
                }
            }
            expect_host_matrices_near(combined.gather(), combined_expected);
        }
    }
}
#endif

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
