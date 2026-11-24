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
#include <type_traits>
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

template <typename Left, typename Right>
class IsTwoArgumentSolveCallable
{
  private:
    template <typename L, typename R>
    static auto test(int) -> decltype(
        matrix2d_linalg::solve(std::declval<L>(), std::declval<R>()),
        std::true_type());

    template <typename, typename>
    static std::false_type test(...);

  public:
    static const bool value
        = decltype(test<Left, Right>(0))::value;
};

template <typename Left, typename Right>
class IsSolveInplaceCallable
{
  private:
    template <typename L, typename R>
    static auto test(int) -> decltype(
        matrix2d_linalg::solve_inplace(std::declval<L>(), std::declval<R>()),
        std::true_type());

    template <typename, typename>
    static std::false_type test(...);

  public:
    static const bool value
        = decltype(test<Left, Right>(0))::value;
};

static_assert(
    std::is_same<
        decltype(matrix2d_linalg::solve(
            std::declval<const Matrix2d<double>&>(),
            std::declval<const Matrix2d<double>&>())),
        Matrix2d<double>>::value,
    "two-argument solve must return Matrix2d by value");
static_assert(
    std::is_same<
        decltype(matrix2d_linalg::solve(
            std::declval<const Matrix2d<double>&>(),
            std::declval<const Matrix2d<double>&>(),
            std::declval<Matrix2d<double>&>())),
        void>::value,
    "three-argument solve must be the output-parameter overload");
static_assert(
    std::is_same<
        decltype(matrix2d_linalg::matmul(
            std::declval<const Matrix2d<double>&>(),
            std::declval<const Matrix2d<double>&>())),
        Matrix2d<double>>::value,
    "two-argument matmul must return Matrix2d by value");
static_assert(
    std::is_same<
        decltype(matrix2d_linalg::inv(
            std::declval<const Matrix2d<double>&>())),
        Matrix2d<double>>::value,
    "one-argument inv must return Matrix2d by value");
static_assert(!IsTwoArgumentSolveCallable<Matrix2dView<double>,
                                          Matrix2dView<double>>::value,
              "views must use the explicitly named solve_inplace primitive");
static_assert(IsSolveInplaceCallable<Matrix2dView<double>,
                                     Matrix2dView<double>>::value,
              "mutable views must bind to solve_inplace");
static_assert(!IsSolveInplaceCallable<Matrix2dView<const double>,
                                      Matrix2dView<double>>::value,
              "solve_inplace must not accept a const coefficient matrix");
static_assert(!IsSolveInplaceCallable<Matrix2dView<double>,
                                      Matrix2dView<const double>>::value,
              "solve_inplace must not accept a const right-hand side");

Matrix2d<double> make_real_system()
{
    Matrix2d<double> matrix(2, 2);
    matrix.atloc(0, 0) = 4.0;
    matrix.atloc(0, 1) = 1.0;
    matrix.atloc(1, 0) = 2.0;
    matrix.atloc(1, 1) = 3.0;
    return matrix;
}

Matrix2d<double> make_real_rhs()
{
    Matrix2d<double> rhs(2, 2);
    rhs.atloc(0, 0) = 1.0;
    rhs.atloc(0, 1) = 2.0;
    rhs.atloc(1, 0) = 3.0;
    rhs.atloc(1, 1) = 4.0;
    return rhs;
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
void expect_distributed_context(const Matrix2d<T>& actual,
                                const Matrix2d<T>& source)
{
    ASSERT_TRUE(actual.has_distribution());
    ASSERT_TRUE(source.has_distribution());
    EXPECT_TRUE(actual.is_active());
    EXPECT_EQ(actual.nrproc(), source.nrproc());
    EXPECT_EQ(actual.ncproc(), source.ncproc());
    EXPECT_EQ(actual.ictxt(), source.ictxt());
    EXPECT_EQ(actual.distribution()->nb(), source.distribution()->nb());
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
    Matrix2d<T> solution(solution_distrib);
    matrix.scatter(host_matrix);
    solution.scatter(host_solution);
    const Matrix2d<T>& matrix_input = matrix;
    const Matrix2d<T>& solution_input = solution;
    const Matrix2d<T> rhs
        = matrix2d_linalg::matmul(matrix_input, solution_input);
    expect_distributed_context(rhs, matrix_input);
    expect_host_matrices_near(rhs.gather(), host_rhs);
    expect_host_matrices_near(matrix.gather(), host_matrix);
    expect_host_matrices_near(solution.gather(), host_solution);

    const Matrix2d<T> solution_value
        = matrix2d_linalg::solve(matrix_input, rhs);
    EXPECT_EQ(solution_value.distribution(), rhs.distribution());
    expect_distributed_context(solution_value, rhs);
    expect_host_matrices_near(solution_value.gather(), host_solution);
    expect_host_matrices_near(matrix.gather(), host_matrix);
    expect_host_matrices_near(rhs.gather(), host_rhs);

    Matrix2d<T> solution_output(1, 1);
    matrix2d_linalg::solve(matrix_input, rhs, solution_output);
    EXPECT_EQ(solution_output.distribution(), rhs.distribution());
    expect_distributed_context(solution_output, rhs);
    expect_host_matrices_near(solution_output.gather(), host_solution);
    expect_host_matrices_near(matrix.gather(), host_matrix);
    expect_host_matrices_near(rhs.gather(), host_rhs);

    Matrix2d<T> product_output(1, 1);
    matrix2d_linalg::matmul(matrix_input,
                            solution_input,
                            product_output);
    expect_distributed_context(product_output, matrix_input);
    expect_host_matrices_near(product_output.gather(), host_rhs);
    expect_host_matrices_near(matrix.gather(), host_matrix);
    expect_host_matrices_near(solution.gather(), host_solution);
}

template <typename T>
void check_distributed_inverse(const Matrix2d<T>& host_matrix,
                               const int nrproc,
                               const int ncproc,
                               const int nrblk)
{
    const Matrix2d<T> host_matrix_original = host_matrix;
    const Matrix2d<T> expected = matrix2d_linalg::inv(host_matrix);
    expect_host_matrices_near(host_matrix, host_matrix_original);

    Matrix2d<T> matrix(host_matrix.nrglb(),
                       host_matrix.ncglb(),
                       nrproc,
                       ncproc,
                       nrblk,
                       MPI_COMM_WORLD);
    matrix.scatter(host_matrix);
    const Matrix2d<T> inverse_value = matrix2d_linalg::inv(matrix);
    EXPECT_EQ(inverse_value.distribution(), matrix.distribution());
    expect_distributed_context(inverse_value, matrix);
    expect_host_matrices_near(inverse_value.gather(), expected);
    expect_host_matrices_near(matrix.gather(), host_matrix_original);

    Matrix2d<T> inverse_output(1, 1);
    matrix2d_linalg::inv(matrix, inverse_output);
    EXPECT_EQ(inverse_output.distribution(), matrix.distribution());
    expect_distributed_context(inverse_output, matrix);
    expect_host_matrices_near(inverse_output.gather(), expected);
    expect_host_matrices_near(matrix.gather(), host_matrix_original);
}
#endif
} // namespace

TEST(Matrix2dLinalgTest, HostRealInverseAndSolve)
{
    const Matrix2d<double> matrix = make_real_system();
    const Matrix2d<double> matrix_original = matrix;
    const Matrix2d<double> inverse_reference = make_real_inverse_reference();
    const Matrix2d<double> inverse = matrix2d_linalg::inv(matrix);
    expect_host_matrices_near(inverse, inverse_reference);

    Matrix2d<double> inverse_output(1, 1);
    matrix2d_linalg::inv(matrix, inverse_output);
    expect_host_matrices_near(inverse_output, inverse_reference);
    expect_host_matrices_near(matrix, matrix_original);

    const Matrix2d<double> rhs = make_real_rhs();
    const Matrix2d<double> rhs_original = rhs;
    const Matrix2d<double> expected
        = matrix2d_linalg::matmul(inverse_reference, rhs);

    const Matrix2d<double> solution = matrix2d_linalg::solve(matrix, rhs);
    expect_host_matrices_near(solution, expected);
    expect_host_matrices_near(rhs, rhs_original);
    expect_host_matrices_near(matrix, matrix_original);

    Matrix2d<double> solution_output(1, 1);
    matrix2d_linalg::solve(matrix, rhs, solution_output);
    expect_host_matrices_near(solution_output, expected);
    expect_host_matrices_near(rhs, rhs_original);
    expect_host_matrices_near(matrix, matrix_original);
}

TEST(Matrix2dLinalgTest, OutputOverloadsMayAliasInputs)
{
    const Matrix2d<double> matrix = make_real_system();
    const Matrix2d<double> rhs = make_real_rhs();
    const Matrix2d<double> expected = matrix2d_linalg::solve(matrix, rhs);

    Matrix2d<double> aliased_rhs = rhs;
    matrix2d_linalg::solve(matrix, aliased_rhs, aliased_rhs);
    expect_host_matrices_near(aliased_rhs, expected);

    Matrix2d<double> aliased_matrix = matrix;
    matrix2d_linalg::solve(aliased_matrix, rhs, aliased_matrix);
    expect_host_matrices_near(aliased_matrix, expected);

    Matrix2d<double> aliased_inverse = matrix;
    matrix2d_linalg::inv(aliased_inverse, aliased_inverse);
    const Matrix2d<double> inverse_reference = make_real_inverse_reference();
    expect_host_matrices_near(aliased_inverse, inverse_reference);

    Matrix2d<double> product = make_real_system();
    matrix2d_linalg::matmul(matrix, product, product);
    const Matrix2d<double> expected_product
        = matrix2d_linalg::matmul(matrix, matrix);
    expect_host_matrices_near(product, expected_product);
}

TEST(Matrix2dLinalgTest, HostOutputAndValueForms)
{
    const Matrix2d<double> left = make_real_system();
    const Matrix2d<double> right = make_real_inverse_reference();
    const Matrix2d<double> product = matrix2d_linalg::matmul(left, right);
    Matrix2d<double> product_output(1, 1);
    matrix2d_linalg::matmul(left, right, product_output);
    const Matrix2d<double> product_reference = left.matmul(right);
    expect_host_matrices_near(product, product_reference);
    expect_host_matrices_near(product_output, product);

    const Matrix2d<double> combination
        = matrix2d_linalg::axpby(2.0, left, -0.5, right);
    Matrix2d<double> combination_output(1, 1);
    matrix2d_linalg::axpby(2.0, left, -0.5, right, combination_output);
    expect_host_matrices_near(combination_output, combination);

    const Matrix2d<double> scaled = matrix2d_linalg::scale(left, 3.0);
    Matrix2d<double> scaled_output(1, 1);
    matrix2d_linalg::scale(left, 3.0, scaled_output);
    expect_host_matrices_near(scaled_output, scaled);
    const Matrix2d<double> left_original = make_real_system();
    expect_host_matrices_near(left, left_original);

    const Matrix2d<double> transposed = matrix2d_linalg::transpose(left);
    const Matrix2d<double> dagger = matrix2d_linalg::dagger(left);
    const Matrix2d<double> transposed_reference = left.transpose();
    const Matrix2d<double> dagger_reference = left.dagger();
    expect_host_matrices_near(transposed, transposed_reference);
    expect_host_matrices_near(dagger, dagger_reference);
}

TEST(Matrix2dLinalgTest, ExplicitInverseSatisfiesInverseIdentity)
{
    const Matrix2d<double> real_matrix = make_real_system();
    Matrix2d<double> real_inverse = real_matrix;
    Matrix2dView<double> real_inverse_view = real_inverse.view();
    matrix2d_linalg::inv_inplace(real_inverse_view);
    Matrix2d<double> real_identity(2, 2);
    Matrix2dView<double> real_identity_view = real_identity.view();
    matrix2d_linalg::set_identity_inplace(real_identity_view);
    const Matrix2d<double> real_product
        = matrix2d_linalg::matmul(real_matrix, real_inverse);
    expect_host_matrices_near(real_product, real_identity);

    const Matrix2d<std::complex<double>> complex_matrix
        = make_complex_system();
    const Matrix2d<std::complex<double>> complex_inverse
        = matrix2d_linalg::inv(complex_matrix);
    const Matrix2d<std::complex<double>> complex_identity
        = matrix2d_linalg::eye<std::complex<double>>(2);
    const Matrix2d<std::complex<double>> complex_product
        = matrix2d_linalg::matmul(complex_matrix, complex_inverse);
    expect_host_matrices_near(complex_product, complex_identity);
}

TEST(Matrix2dLinalgTest, HostComplexSolve)
{
    const Matrix2d<std::complex<double>> matrix = make_complex_system();
    const Matrix2d<std::complex<double>> matrix_original = matrix;
    const Matrix2d<std::complex<double>> expected = make_complex_solution();
    const Matrix2d<std::complex<double>> rhs
        = matrix2d_linalg::matmul(matrix, expected);
    const Matrix2d<std::complex<double>> rhs_original = rhs;

    const Matrix2d<std::complex<double>> solution
        = matrix2d_linalg::solve(matrix, rhs);
    expect_host_matrices_near(solution, expected);
    expect_host_matrices_near(rhs, rhs_original);
    expect_host_matrices_near(matrix, matrix_original);

    Matrix2d<std::complex<double>> solution_output(1, 1);
    matrix2d_linalg::solve(matrix, rhs, solution_output);
    expect_host_matrices_near(solution_output, expected);

    const Matrix2d<std::complex<double>> inverse
        = matrix2d_linalg::inv(matrix);
    Matrix2d<std::complex<double>> inverse_reference(2, 2);
    inverse_reference.atloc(0, 0) = std::complex<double>(-0.7, -1.1);
    inverse_reference.atloc(0, 1) = std::complex<double>(0.2, 0.6);
    inverse_reference.atloc(1, 0) = std::complex<double>(0.3, 0.9);
    inverse_reference.atloc(1, 1) = std::complex<double>(0.2, -0.4);
    expect_host_matrices_near(inverse, inverse_reference);
    expect_host_matrices_near(matrix, matrix_original);
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
    const Matrix2d<double> matrix_original = matrix;

    Matrix2d<double> rhs(distrib);
    rhs.atloc(0, 0) = 1.0;
    rhs.atloc(0, 1) = 2.0;
    rhs.atloc(1, 0) = 3.0;
    rhs.atloc(1, 1) = 4.0;
    const Matrix2d<double> rhs_original = rhs;
    const Matrix2d<double> inverse_reference = make_real_inverse_reference();
    const Matrix2d<double> expected = inverse_reference.matmul(rhs);
    const Matrix2d<double> solution = matrix2d_linalg::solve(matrix, rhs);
    expect_host_matrices_near(solution, expected);
    EXPECT_EQ(solution.distribution(), rhs.distribution());
    expect_host_matrices_near(matrix, matrix_original);
    expect_host_matrices_near(rhs, rhs_original);

    Matrix2d<double> solution_output(1, 1);
    matrix2d_linalg::solve(matrix, rhs, solution_output);
    expect_host_matrices_near(solution_output, expected);
    EXPECT_EQ(solution_output.distribution(), rhs.distribution());
    expect_host_matrices_near(matrix, matrix_original);
    expect_host_matrices_near(rhs, rhs_original);

    const Matrix2d<double> inverse = matrix2d_linalg::inv(matrix);
    expect_host_matrices_near(inverse, inverse_reference);
    EXPECT_EQ(inverse.distribution(), matrix.distribution());
    expect_host_matrices_near(matrix, matrix_original);

    const Matrix2d<double> combination
        = matrix2d_linalg::axpby(1.0, matrix, 0.0, matrix);
    EXPECT_EQ(combination.distribution(), matrix.distribution());
    expect_near(combination.dot(matrix), matrix.dot(matrix));

    Matrix2d<double> accumulated = combination;
    accumulated += matrix;
    Matrix2d<double> accumulated_reference(2, 2);
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            accumulated_reference.atloc(i, j)
                = 2.0 * matrix_original.atloc(i, j);
        }
    }
    expect_host_matrices_near(accumulated.gather(), accumulated_reference);
}

TEST(Matrix2dLinalgTest, ExternalBufferViewDoesNotCopy)
{
    std::vector<double> values = {4.0, 1.0, 2.0, 3.0};
    Matrix2dView<double> matrix(values.data(), 2, 2);
    Matrix2d<double> rhs(2, 1);
    rhs.atloc(0, 0) = 1.0;
    rhs.atloc(1, 0) = 2.0;
    matrix2d_linalg::solve_inplace(matrix, rhs.view());

    expect_near(rhs.atloc(0, 0), 0.1);
    expect_near(rhs.atloc(1, 0), 0.6);
    EXPECT_EQ(matrix.data(), values.data());
}

TEST(Matrix2dLinalgTest, HostViewMatmulAndAxpby)
{
    const Matrix2d<double> left = make_real_system();
    const Matrix2d<double> right = make_real_inverse_reference();
    Matrix2d<double> product(2, 2);
    matrix2d_linalg::matmul_inplace(product.view(),
                                    left.view(),
                                    right.view());
    expect_host_matrices_near(product, left.matmul(right));

    Matrix2d<double> combined(2, 2);
    matrix2d_linalg::axpby_inplace(
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
    matrix2d_linalg::axpby_inplace(in_place.view(),
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

    matrix2d_linalg::scale_inplace(combined.view(), 3.0);
    matrix2d_linalg::set_identity_inplace(product.view());
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
    EXPECT_THROW(matrix2d_linalg::inv(singular), std::runtime_error);
    EXPECT_THROW(matrix2d_linalg::inv_inplace(singular.view()),
                 std::runtime_error);

    const Matrix2d<double> matrix = make_real_system();
    const Matrix2d<double> wrong_rhs(1, 1);
    EXPECT_THROW(matrix2d_linalg::solve(matrix, wrong_rhs),
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
        matrix2d_linalg::axpby_inplace(host.view(),
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
            matrix2d_linalg::set_identity_inplace(identity.view());
            Matrix2d<double> identity_expected(size, size);
            matrix2d_linalg::set_identity_inplace(identity_expected.view());
            expect_host_matrices_near(identity.gather(), identity_expected);

            const Matrix2d<double> host_left = make_distributed_real_system(size);
            Matrix2d<double> left(identity.distribution());
            left.scatter(host_left);
            const Matrix2d<double> combined
                = matrix2d_linalg::axpby(2.0,
                                         left,
                                         0.5,
                                         identity);
            EXPECT_EQ(combined.distribution(), left.distribution());
            expect_distributed_context(combined, left);

            Matrix2d<double> combined_output(1, 1);
            matrix2d_linalg::axpby(2.0,
                                   left,
                                   0.5,
                                   identity,
                                   combined_output);
            EXPECT_EQ(combined_output.distribution(), left.distribution());
            expect_distributed_context(combined_output, left);

            const Matrix2d<double> scaled
                = matrix2d_linalg::scale(left, 3.0);
            EXPECT_EQ(scaled.distribution(), left.distribution());
            expect_distributed_context(scaled, left);

            Matrix2d<double> combined_expected(size, size);
            Matrix2d<double> scaled_expected(size, size);
            for (int i = 0; i < size; ++i)
            {
                for (int j = 0; j < size; ++j)
                {
                    combined_expected.atloc(i, j)
                        = 2.0 * host_left.atloc(i, j)
                            + 0.5 * identity_expected.atloc(i, j);
                    scaled_expected.atloc(i, j)
                        = 3.0 * host_left.atloc(i, j);
                }
            }
            expect_host_matrices_near(combined.gather(), combined_expected);
            expect_host_matrices_near(combined_output.gather(),
                                      combined_expected);
            expect_host_matrices_near(scaled.gather(), scaled_expected);
            expect_host_matrices_near(left.gather(), host_left);
            expect_host_matrices_near(identity.gather(), identity_expected);
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
