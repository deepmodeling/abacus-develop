#ifdef NDEBUG
#undef NDEBUG
#endif

#include "source_base/matrix2d.h"

#ifdef __MPI
#include "source_base/parallel_common.h"
#endif

#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

static_assert(std::is_same<Matrix2d<double>, NDArray<double, 2>>::value,
              "Matrix2d<T> must be the same type as NDArray<T, 2>");
static_assert(std::is_same<NDArray<double>, NDArray<double, 0>>::value,
              "NDArray<T> must remain the rank-zero specialization");
static_assert(!std::is_same<NDArray<double>, NDArray<double, 2>>::value,
              "Rank-zero and rank-two NDArray types must be distinct");
static_assert(!std::is_base_of<NDArray<double>, Matrix2d<double>>::value,
              "Matrix2d must not inherit from the rank-zero NDArray");
static_assert(std::is_copy_constructible<Matrix2d<double>>::value,
              "Matrix2d must be copy constructible");
static_assert(std::is_copy_assignable<Matrix2d<double>>::value,
              "Matrix2d must be copy assignable");
static_assert(std::is_nothrow_move_constructible<Matrix2d<double>>::value,
              "Matrix2d move construction must be noexcept");
static_assert(std::is_nothrow_move_assignable<Matrix2d<double>>::value,
              "Matrix2d move assignment must be noexcept");

namespace
{
template <typename T>
class HasAt
{
    template <typename U>
    static auto test(int) -> decltype(std::declval<U&>().at(0), std::true_type());

    template <typename>
    static std::false_type test(...);

  public:
    static constexpr bool value = decltype(test<T>(0))::value;
};

template <typename T>
class HasBegin
{
    template <typename U>
    static auto test(int) -> decltype(std::declval<U&>().begin(), std::true_type());

    template <typename>
    static std::false_type test(...);

  public:
    static constexpr bool value = decltype(test<T>(0))::value;
};

template <typename T>
class HasReshape
{
    template <typename U>
    static auto test(int) -> decltype(std::declval<U&>().reshape(1, 1), std::true_type());

    template <typename>
    static std::false_type test(...);

  public:
    static constexpr bool value = decltype(test<T>(0))::value;
};

template <typename T>
class HasShape
{
    template <typename U>
    static auto test(int) -> decltype(std::declval<U&>().shape(), std::true_type());

    template <typename>
    static std::false_type test(...);

  public:
    static constexpr bool value = decltype(test<T>(0))::value;
};

template <typename T>
class HasSize
{
    template <typename U>
    static auto test(int) -> decltype(std::declval<U&>().size(), std::true_type());

    template <typename>
    static std::false_type test(...);

  public:
    static constexpr bool value = decltype(test<T>(0))::value;
};

static_assert(!HasAt<Matrix2d<double>>::value,
              "Matrix2d must not expose the rank-zero at() interface");
static_assert(!HasBegin<Matrix2d<double>>::value,
              "Matrix2d must not expose rank-zero iterators");
static_assert(!HasReshape<Matrix2d<double>>::value,
              "Matrix2d must not expose rank-zero reshape()");
static_assert(!HasShape<Matrix2d<double>>::value,
              "Matrix2d must not expose a generic shape() container");
static_assert(HasSize<Matrix2d<double>>::value,
              "Matrix2d must expose the local element count as size()");

void check_host_layout()
{
    Matrix2d<double> matrix(2, 3);

    EXPECT_FALSE(matrix.is_distributed());
    EXPECT_EQ(matrix.distribution(), nullptr);
    EXPECT_EQ(matrix.nrglb(), 2);
    EXPECT_EQ(matrix.ncglb(), 3);
    EXPECT_EQ(matrix.nrloc(), 2);
    EXPECT_EQ(matrix.ncloc(), 3);
    EXPECT_EQ(matrix.size(), 6);
    EXPECT_EQ(matrix.nrproc(), 1);
    EXPECT_EQ(matrix.ncproc(), 1);
    EXPECT_EQ(matrix.irproc(), 0);
    EXPECT_EQ(matrix.icproc(), 0);
    EXPECT_EQ(matrix.lld(), 3);
    EXPECT_EQ(matrix.desc(), nullptr);
    EXPECT_EQ(matrix.ictxt(), -1);
    EXPECT_EQ(matrix.data(), &matrix.at_local(0, 0));

    double value = 1.0;
    for (int row = 0; row < 2; ++row)
    {
        for (int col = 0; col < 3; ++col)
        {
            const std::size_t expected_index = static_cast<std::size_t>(row) * 3U
                + static_cast<std::size_t>(col);
            EXPECT_EQ(&matrix.at_local(row, col), matrix.data() + expected_index);
            EXPECT_EQ(&matrix.at_global(row, col), &matrix.at_local(row, col));
            EXPECT_TRUE(matrix.owns_global_index(row, col));
            matrix.at_local(row, col) = value;
            value += 1.0;
        }
    }

    EXPECT_EQ(matrix.data()[0], 1.0);
    EXPECT_EQ(matrix.data()[1], 2.0);
    EXPECT_EQ(matrix.data()[2], 3.0);
    EXPECT_EQ(matrix.data()[3], 4.0);
    EXPECT_EQ(matrix.data()[4], 5.0);
    EXPECT_EQ(matrix.data()[5], 6.0);

    const Matrix2d<double>& const_matrix = matrix;
    EXPECT_EQ(const_matrix.data(), matrix.data());
    EXPECT_EQ(&const_matrix.at_local(1, 2), matrix.data() + 5);
    EXPECT_EQ(&const_matrix.at_global(1, 2), matrix.data() + 5);
}

void check_distribution_layout(Matrix2d<double>& matrix,
                               const int global_rows,
                               const int global_cols)
{
    const std::shared_ptr<const Parallel_Blkcyclic>& distribution
        = matrix.distribution();
    ASSERT_NE(distribution, nullptr);
    EXPECT_EQ(matrix.nrglb(), global_rows);
    EXPECT_EQ(matrix.ncglb(), global_cols);
    EXPECT_EQ(matrix.nrloc(), distribution->nrloc());
    EXPECT_EQ(matrix.ncloc(), distribution->ncloc());
    EXPECT_EQ(matrix.size(),
              distribution->size());
    EXPECT_EQ(matrix.nrproc(), distribution->nrproc());
    EXPECT_EQ(matrix.ncproc(), distribution->ncproc());
    EXPECT_EQ(matrix.irproc(), distribution->irproc());
    EXPECT_EQ(matrix.icproc(), distribution->icproc());
    EXPECT_EQ(matrix.lld(), distribution->lld());

    const int local_rows = matrix.nrloc();
    const int local_cols = matrix.ncloc();
    const std::size_t leading_dimension
        = static_cast<std::size_t>(matrix.lld());
    for (int local_row = 0; local_row < local_rows; ++local_row)
    {
        for (int local_col = 0; local_col < local_cols; ++local_col)
        {
            const std::size_t expected_index
                = static_cast<std::size_t>(local_col) * leading_dimension
                + static_cast<std::size_t>(local_row);
            EXPECT_EQ(&matrix.at_local(local_row, local_col), matrix.data() + expected_index);
        }
    }

    for (int global_row = 0; global_row < global_rows; ++global_row)
    {
        for (int global_col = 0; global_col < global_cols; ++global_col)
        {
            const bool expected_owner = distribution->in_this_processor(global_row, global_col);
            EXPECT_EQ(matrix.owns_global_index(global_row, global_col), expected_owner);
            if (!expected_owner)
            {
                continue;
            }

            const int local_row = distribution->irloc(global_row);
            const int local_col = distribution->icloc(global_col);
            ASSERT_NE(local_row, -1);
            ASSERT_NE(local_col, -1);
            EXPECT_EQ(&matrix.at_global(global_row, global_col),
                      &matrix.at_local(local_row, local_col));
            const double value = static_cast<double>(global_row * global_cols + global_col);
            matrix.at_global(global_row, global_col) = value;
            EXPECT_EQ(matrix.at_local(local_row, local_col), value);
        }
    }

    const Matrix2d<double>& const_matrix = matrix;
    EXPECT_EQ(const_matrix.data(), matrix.data());
}

#ifdef __MPI
void check_mpi_descriptor(Matrix2d<double>& matrix,
                          const int global_rows,
                          const int global_cols,
                          const int block_size)
{
    constexpr int DESCRIPTOR_MATRIX_TYPE_INDEX = 0;
    constexpr int DESCRIPTOR_BLACS_CONTEXT_INDEX = 1;
    constexpr int DESCRIPTOR_GLOBAL_ROWS_INDEX = 2;
    constexpr int DESCRIPTOR_GLOBAL_COLS_INDEX = 3;
    constexpr int DESCRIPTOR_BLOCK_ROWS_INDEX = 4;
    constexpr int DESCRIPTOR_BLOCK_COLS_INDEX = 5;
    constexpr int DESCRIPTOR_LEADING_DIM_INDEX = 8;

    const int* const descriptor = matrix.desc();
    ASSERT_NE(descriptor, nullptr);
    EXPECT_EQ(descriptor[DESCRIPTOR_MATRIX_TYPE_INDEX], 1);
    EXPECT_EQ(descriptor[DESCRIPTOR_BLACS_CONTEXT_INDEX], matrix.ictxt());
    EXPECT_EQ(descriptor[DESCRIPTOR_GLOBAL_ROWS_INDEX], global_rows);
    EXPECT_EQ(descriptor[DESCRIPTOR_GLOBAL_COLS_INDEX], global_cols);
    EXPECT_EQ(descriptor[DESCRIPTOR_BLOCK_ROWS_INDEX], block_size);
    EXPECT_EQ(descriptor[DESCRIPTOR_BLOCK_COLS_INDEX], block_size);
    EXPECT_EQ(descriptor[DESCRIPTOR_LEADING_DIM_INDEX], matrix.lld());
    EXPECT_GE(matrix.ictxt(), 0);
    EXPECT_EQ(matrix.nrproc() * matrix.ncproc(),
              Parallel_Common::comm_size(MPI_COMM_WORLD));
}
#endif
} // namespace

namespace
{
constexpr double DOUBLE_RELATIVE_TOLERANCE = 1.0e-12;
constexpr double COMPLEX_RELATIVE_TOLERANCE = 1.0e-6;

void report_reduction_tolerances()
{
    std::cout << "Matrix2d reduction tolerance: double relative <= "
              << DOUBLE_RELATIVE_TOLERANCE << ", complex<double> relative <= "
              << COMPLEX_RELATIVE_TOLERANCE << std::endl;
}

template <typename T>
double reduction_tolerance()
{
    return DOUBLE_RELATIVE_TOLERANCE;
}

template <>
double reduction_tolerance<std::complex<double>>()
{
    return COMPLEX_RELATIVE_TOLERANCE;
}

template <typename T>
T test_value(const int row, const int col, const double offset)
{
    return static_cast<T>(offset + 0.5 + 1.25 * row - 0.75 * col);
}

template <>
std::complex<double> test_value<std::complex<double>>(const int row,
                                                       const int col,
                                                       const double offset)
{
    return std::complex<double>(offset + 0.5 + 1.25 * row - 0.75 * col,
                                0.25 - 0.5 * row + 0.875 * col);
}

template <typename T>
T constant_value(const double real_part)
{
    return static_cast<T>(real_part);
}

template <>
std::complex<double> constant_value<std::complex<double>>(const double real_part)
{
    return std::complex<double>(real_part, -0.25 * real_part);
}

template <typename T>
double real_part(const T& value)
{
    return static_cast<double>(value);
}

template <>
double real_part<std::complex<double>>(const std::complex<double>& value)
{
    return value.real();
}

template <typename T>
double imaginary_part(const T&)
{
    return 0.0;
}

template <>
double imaginary_part<std::complex<double>>(const std::complex<double>& value)
{
    return value.imag();
}

template <typename T>
T conjugate_value(const T& value)
{
    return value;
}

template <>
std::complex<double> conjugate_value<std::complex<double>>(
    const std::complex<double>& value)
{
    return std::conj(value);
}

template <typename T>
void expect_scalar_near(const T& actual, const T& expected)
{
    const double tolerance = reduction_tolerance<T>();
    const double scale = std::max(1.0, std::abs(expected));
    EXPECT_NEAR(actual, expected, tolerance * scale);
}

template <>
void expect_scalar_near<std::complex<double>>(
    const std::complex<double>& actual,
    const std::complex<double>& expected)
{
    const double tolerance = reduction_tolerance<std::complex<double>>();
    const double scale = std::max(1.0, std::abs(expected));
    EXPECT_NEAR(actual.real(), expected.real(), tolerance * scale);
    EXPECT_NEAR(actual.imag(), expected.imag(), tolerance * scale);
}

template <typename T>
void fill_global_values(Matrix2d<T>& matrix, const double offset = 0.0)
{
    for (int row = 0; row < matrix.nrglb(); ++row)
    {
        for (int col = 0; col < matrix.ncglb(); ++col)
        {
            if (matrix.owns_global_index(row, col))
            {
                matrix.at_global(row, col) = test_value<T>(row, col, offset);
            }
        }
    }
}

template <typename T>
Matrix2d<T> make_host_reference(const int rows,
                                const int cols,
                                const double offset = 0.0)
{
    Matrix2d<T> reference(rows, cols);
    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < cols; ++col)
        {
            reference.at_local(row, col) = test_value<T>(row, col, offset);
        }
    }
    return reference;
}

template <typename T>
void verify_global_matrix(const Matrix2d<T>& actual, const Matrix2d<T>& expected)
{
    ASSERT_FALSE(actual.is_distributed());
    ASSERT_EQ(actual.nrglb(), expected.nrglb());
    ASSERT_EQ(actual.ncglb(), expected.ncglb());
    for (int row = 0; row < expected.nrglb(); ++row)
    {
        for (int col = 0; col < expected.ncglb(); ++col)
        {
            const std::size_t c_order_index
                = static_cast<std::size_t>(row)
                    * static_cast<std::size_t>(expected.ncglb())
                + static_cast<std::size_t>(col);
            expect_scalar_near(actual.data()[c_order_index],
                               expected.at_local(row, col));
        }
    }
}

template <typename T>
void check_elementwise_operations(Matrix2d<T>& matrix, Matrix2d<T>& other)
{
    const int rows = matrix.nrglb();
    const int cols = matrix.ncglb();

    Matrix2d<T> filled = matrix;
    const T fill_value = constant_value<T>(3.5);
    filled.fill(fill_value);
    Matrix2d<T> expected(rows, cols);
    expected.fill(fill_value);
    verify_global_matrix(filled.gather(), expected);

    filled = matrix;
    filled.scale(constant_value<T>(2.5));
    expected = make_host_reference<T>(rows, cols);
    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < cols; ++col)
        {
            expected.at_local(row, col) *= constant_value<T>(2.5);
        }
    }
    verify_global_matrix(filled.gather(), expected);

    filled = matrix;
    filled += other;
    expected = make_host_reference<T>(rows, cols);
    Matrix2d<T> expected_other = make_host_reference<T>(rows, cols, 10.0);
    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < cols; ++col)
        {
            expected.at_local(row, col) += expected_other.at_local(row, col);
        }
    }
    verify_global_matrix(filled.gather(), expected);

    filled = matrix;
    filled -= other;
    expected = make_host_reference<T>(rows, cols);
    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < cols; ++col)
        {
            expected.at_local(row, col) -= expected_other.at_local(row, col);
        }
    }
    verify_global_matrix(filled.gather(), expected);

    filled = matrix;
    filled *= other;
    expected = make_host_reference<T>(rows, cols);
    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < cols; ++col)
        {
            expected.at_local(row, col) *= expected_other.at_local(row, col);
        }
    }
    verify_global_matrix(filled.gather(), expected);

    filled = matrix;
    filled.conjugate();
    expected = make_host_reference<T>(rows, cols);
    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < cols; ++col)
        {
            expected.at_local(row, col)
                = conjugate_value(expected.at_local(row, col));
        }
    }
    verify_global_matrix(filled.gather(), expected);

    Matrix2d<double> expected_real(rows, cols);
    Matrix2d<double> expected_imag(rows, cols);
    const Matrix2d<T> matrix_reference = make_host_reference<T>(rows, cols);
    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < cols; ++col)
        {
            const T& value = matrix_reference.at_local(row, col);
            expected_real.at_local(row, col) = real_part(value);
            expected_imag.at_local(row, col) = imaginary_part(value);
        }
    }

    Matrix2d<double> real_output = matrix.is_distributed()
        ? Matrix2d<double>(matrix.distribution())
        : Matrix2d<double>(rows, cols);
    Matrix2d<double> imaginary_output = matrix.is_distributed()
        ? Matrix2d<double>(matrix.distribution())
        : Matrix2d<double>(rows, cols);
    matrix.real(real_output);
    matrix.imag(imaginary_output);
    verify_global_matrix(real_output.gather(), expected_real);
    verify_global_matrix(imaginary_output.gather(), expected_imag);
    verify_global_matrix(matrix.real().gather(), expected_real);
    verify_global_matrix(matrix.imag().gather(), expected_imag);
}

template <typename T>
void check_reductions(Matrix2d<T>& matrix, Matrix2d<T>& other)
{
    const Matrix2d<T> matrix_reference = make_host_reference<T>(
        matrix.nrglb(),
        matrix.ncglb());
    const Matrix2d<T> other_reference = make_host_reference<T>(
        matrix.nrglb(),
        matrix.ncglb(),
        10.0);
    T reference_sum = T();
    T reference_dot = T();
    double reference_sum_squares = 0.0;
    double reference_max_abs = 0.0;
    for (int row = 0; row < matrix.nrglb(); ++row)
    {
        for (int col = 0; col < matrix.ncglb(); ++col)
        {
            const T& value = matrix_reference.at_local(row, col);
            const T& other_value = other_reference.at_local(row, col);
            reference_sum += value;
            reference_dot += conjugate_value(value) * other_value;
            reference_sum_squares += real_part(value) * real_part(value)
                + imaginary_part(value) * imaginary_part(value);
            reference_max_abs = std::max(reference_max_abs,
                                         static_cast<double>(std::abs(value)));
        }
    }

    expect_scalar_near(matrix.sum(), reference_sum);
    expect_scalar_near(matrix.dot(other), reference_dot);
    const T expected_self_dot = static_cast<T>(matrix.sum_squares());
    expect_scalar_near(matrix.dot(matrix), expected_self_dot);
    expect_scalar_near(matrix.sum_squares(), reference_sum_squares);
    expect_scalar_near(matrix.norm2(), std::sqrt(reference_sum_squares));
    expect_scalar_near(matrix.max_abs(), reference_max_abs);

    const std::vector<T> diagonal = matrix.diagonal();
    const int diagonal_size = std::min(matrix.nrglb(), matrix.ncglb());
    ASSERT_EQ(diagonal.size(), static_cast<std::size_t>(diagonal_size));
    for (int i = 0; i < diagonal_size; ++i)
    {
        expect_scalar_near(diagonal[static_cast<std::size_t>(i)],
                           matrix_reference.at_local(i, i));
    }
}

template <typename T>
void check_gather_scatter(Matrix2d<T>& matrix)
{
    const int rows = matrix.nrglb();
    const int cols = matrix.ncglb();
    const Matrix2d<T> expected = make_host_reference<T>(rows, cols);
    const Matrix2d<T> gathered = matrix.gather();
    verify_global_matrix(gathered, expected);
    EXPECT_EQ(std::memcmp(gathered.data(),
                          expected.data(),
                          static_cast<std::size_t>(rows * cols) * sizeof(T)),
              0);
    const Matrix2d<T> gathered_again = matrix.gather();
    EXPECT_EQ(std::memcmp(gathered.data(),
                          gathered_again.data(),
                          static_cast<std::size_t>(rows * cols) * sizeof(T)),
              0);

    Matrix2d<T> restored(matrix.distribution());
    restored.scatter(gathered);
    verify_global_matrix(restored.gather(), expected);
    verify_global_matrix(matrix.gather(), expected);

#ifdef __MPI
    const int rank = Parallel_Common::comm_rank(MPI_COMM_WORLD);
    const int root_rank = Parallel_Common::comm_size(MPI_COMM_WORLD) - 1;
    const Matrix2d<T> root_result = matrix.gather_on_root(root_rank);
    if (rank == root_rank)
    {
        verify_global_matrix(root_result, expected);
    }
    else
    {
        EXPECT_EQ(root_result.nrglb(), 0);
        EXPECT_EQ(root_result.ncglb(), 0);
        EXPECT_EQ(root_result.size(), 0);
        EXPECT_EQ(root_result.lld(), 0);
        EXPECT_FALSE(root_result.is_distributed());
    }
#endif

    Matrix2d<T> wrong_shape(rows + 1, cols);
    EXPECT_THROW(restored.scatter(wrong_shape), std::invalid_argument);
}

template <typename T>
std::size_t set_padding_sentinel(Matrix2d<T>& matrix, const T& sentinel)
{
    const int local_rows = matrix.nrloc();
    const int local_cols = matrix.ncloc();
    const std::size_t leading_dimension
        = static_cast<std::size_t>(std::max(local_rows, 1));
    const std::size_t column_major_elements
        = leading_dimension * static_cast<std::size_t>(local_cols);
    const std::size_t row_major_elements
        = static_cast<std::size_t>(local_rows)
        * std::max(static_cast<std::size_t>(local_cols), leading_dimension);
    const std::size_t storage_elements = std::max(column_major_elements, row_major_elements);
    const std::size_t logical_elements
        = static_cast<std::size_t>(local_rows) * static_cast<std::size_t>(local_cols);
    for (std::size_t index = logical_elements; index < storage_elements; ++index)
    {
        matrix.data()[index] = sentinel;
    }
    return storage_elements - logical_elements;
}

template <typename T>
void expect_padding_sentinel(const Matrix2d<T>& matrix,
                             const std::size_t padding_elements,
                             const T& sentinel)
{
    const int local_rows = matrix.nrloc();
    const int local_cols = matrix.ncloc();
    const std::size_t logical_elements
        = static_cast<std::size_t>(local_rows) * static_cast<std::size_t>(local_cols);
    for (std::size_t index = 0; index < padding_elements; ++index)
    {
        expect_scalar_near(matrix.data()[logical_elements + index], sentinel);
    }
}

std::vector<std::pair<int, int>> available_process_grids()
{
    std::vector<std::pair<int, int>> grids;
#ifdef __MPI
    const int number_of_ranks = Parallel_Common::comm_size(MPI_COMM_WORLD);
#else
    const int number_of_ranks = 1;
#endif
    if (number_of_ranks == 1)
    {
        grids.push_back(std::make_pair(1, 1));
    }
    if (number_of_ranks == 2)
    {
        grids.push_back(std::make_pair(1, 2));
        grids.push_back(std::make_pair(2, 1));
    }
    if (number_of_ranks == 4)
    {
        grids.push_back(std::make_pair(1, 4));
        grids.push_back(std::make_pair(2, 2));
        grids.push_back(std::make_pair(4, 1));
    }
    return grids;
}

template <typename T>
Matrix2d<T> make_transpose_reference(const Matrix2d<T>& source)
{
    const int source_rows = source.nrglb();
    const int source_cols = source.ncglb();
    Matrix2d<T> reference(source_cols, source_rows);
    for (int row = 0; row < source_rows; ++row)
    {
        for (int col = 0; col < source_cols; ++col)
        {
            reference.at_local(col, row)
                = source.at_local(row, col);
        }
    }
    return reference;
}

template <typename T>
Matrix2d<T> make_dagger_reference(const Matrix2d<T>& source)
{
    const int source_rows = source.nrglb();
    const int source_cols = source.ncglb();
    Matrix2d<T> reference(source_cols, source_rows);
    for (int row = 0; row < source_rows; ++row)
    {
        for (int col = 0; col < source_cols; ++col)
        {
            reference.at_local(col, row)
                = conjugate_value(source.at_local(row, col));
        }
    }
    return reference;
}

template <typename T>
Matrix2d<T> make_matmul_reference(const Matrix2d<T>& left,
                                  const Matrix2d<T>& right)
{
    const int rows = left.nrglb();
    const int shared = left.ncglb();
    const int cols = right.ncglb();
    Matrix2d<T> reference(rows, cols);
    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < cols; ++col)
        {
            T value = T();
            for (int inner = 0; inner < shared; ++inner)
            {
                value += left.at_local(row, inner) * right.at_local(inner, col);
            }
            reference.at_local(row, col) = value;
        }
    }
    return reference;
}

template <typename T>
Matrix2d<T> make_host_transpose_reference(const int rows,
                                          const int cols)
{
    return make_transpose_reference(make_host_reference<T>(rows, cols));
}

template <typename T>
Matrix2d<T> make_host_matmul_reference(const int left_rows,
                                       const int shared_cols,
                                       const int right_cols)
{
    return make_matmul_reference(make_host_reference<T>(left_rows, shared_cols),
                                 make_host_reference<T>(shared_cols, right_cols, 10.0));
}

template <typename T>
void check_host_transpose(const int rows, const int cols)
{
    const Matrix2d<T> source = make_host_reference<T>(rows, cols);
    const Matrix2d<T> transposed = source.transpose();
    const Matrix2d<T> dagger = source.dagger();
    const Matrix2d<T> expected_transpose = make_transpose_reference(source);
    const Matrix2d<T> expected_dagger = make_dagger_reference(source);

    verify_global_matrix(transposed, expected_transpose);
    verify_global_matrix(dagger, expected_dagger);
    EXPECT_FALSE(transposed.is_distributed());
    EXPECT_FALSE(dagger.is_distributed());
    EXPECT_EQ(transposed.nrglb(), cols);
    EXPECT_EQ(transposed.ncglb(), rows);
    EXPECT_EQ(dagger.nrglb(), cols);
    EXPECT_EQ(dagger.ncglb(), rows);
}

template <typename T>
void check_host_matmul(const int left_rows,
                       const int shared_cols,
                       const int right_cols)
{
    const Matrix2d<T> left = make_host_reference<T>(left_rows, shared_cols);
    const Matrix2d<T> right = make_host_reference<T>(shared_cols, right_cols, 10.0);
    const Matrix2d<T> product = left.matmul(right);
    const Matrix2d<T> expected = make_matmul_reference(left, right);

    verify_global_matrix(product, expected);
    EXPECT_FALSE(product.is_distributed());
}

#ifdef __MPI
template <typename T>
void check_mpi_descriptor_for_type(const Matrix2d<T>& matrix,
                                   const int global_rows,
                                   const int global_cols,
                                   const int block_size)
{
    const int* const descriptor = matrix.desc();
    ASSERT_NE(descriptor, nullptr);
    EXPECT_EQ(descriptor[0], 1);
    EXPECT_EQ(descriptor[1], matrix.ictxt());
    EXPECT_EQ(descriptor[2], global_rows);
    EXPECT_EQ(descriptor[3], global_cols);
    EXPECT_EQ(descriptor[4], block_size);
    EXPECT_EQ(descriptor[5], block_size);
    EXPECT_EQ(descriptor[6], 0);
    EXPECT_EQ(descriptor[7], 0);
    EXPECT_GT(descriptor[8], 0);
}

template <typename T>
void check_distributed_transpose(const int rows,
                                 const int cols,
                                 const int grid_rows,
                                 const int grid_cols,
                                 const int block_size)
{
    Matrix2d<T> source(rows,
                       cols,
                       grid_rows,
                       grid_cols,
                       block_size,
                       MPI_COMM_WORLD);
    fill_global_values(source);
    const int* const source_descriptor_before = source.desc();

    const Matrix2d<T> transposed = source.transpose();
    const Matrix2d<T> dagger = source.dagger();
    const int* const source_descriptor_after = source.desc();
    const int* const result_descriptor = transposed.desc();
    const int* const dagger_descriptor = dagger.desc();

    ASSERT_TRUE(transposed.is_distributed());
    ASSERT_TRUE(dagger.is_distributed());
    EXPECT_EQ(source_descriptor_after, source_descriptor_before);
    EXPECT_NE(result_descriptor, nullptr);
    EXPECT_NE(dagger_descriptor, nullptr);
    EXPECT_EQ(source.ictxt(), transposed.ictxt());
    EXPECT_EQ(source.ictxt(), dagger.ictxt());
    check_mpi_descriptor_for_type(transposed, cols, rows, block_size);
    check_mpi_descriptor_for_type(dagger, cols, rows, block_size);
    verify_global_matrix(transposed.gather(),
                         make_host_transpose_reference<T>(rows, cols));
    verify_global_matrix(dagger.gather(),
                         make_dagger_reference(source.gather()));
    EXPECT_EQ(source.desc(), source_descriptor_after);
    EXPECT_EQ(transposed.desc(), result_descriptor);
    EXPECT_EQ(dagger.desc(), dagger_descriptor);
}

template <typename T>
Matrix2d<T> make_transpose_without_source(const int rows,
                                          const int cols,
                                          const int grid_rows,
                                          const int grid_cols,
                                          const int block_size)
{
    Matrix2d<T> source(rows,
                       cols,
                       grid_rows,
                       grid_cols,
                       block_size,
                       MPI_COMM_WORLD);
    fill_global_values(source);
    return source.transpose();
}

template <typename T>
void check_distributed_matmul(const int left_rows,
                              const int shared_cols,
                              const int right_cols,
                              const int grid_rows,
                              const int grid_cols,
                              const int block_size)
{
    Matrix2d<T> left(left_rows,
                     shared_cols,
                     grid_rows,
                     grid_cols,
                     block_size,
                     MPI_COMM_WORLD);
    const int left_blacs_context = left.ictxt();
    const std::shared_ptr<const Parallel_Blkcyclic> right_distribution
        = make_borrowed_distribution(shared_cols,
                                     right_cols,
                                     block_size,
                                     left_blacs_context);
    Matrix2d<T> right(right_distribution);
    fill_global_values(left);
    fill_global_values(right, 10.0);

    const int* const left_descriptor_before = left.desc();
    const int* const right_descriptor_before = right.desc();
    const Matrix2d<T> product = left.matmul(right);
    const int* const left_descriptor_after = left.desc();
    const int* const right_descriptor_after = right.desc();
    const int* const result_descriptor = product.desc();

    ASSERT_TRUE(product.is_distributed());
    EXPECT_EQ(left_descriptor_after, left_descriptor_before);
    EXPECT_EQ(right_descriptor_after, right_descriptor_before);
    EXPECT_NE(result_descriptor, nullptr);
    EXPECT_EQ(left.ictxt(), right.ictxt());
    EXPECT_EQ(left.ictxt(), product.ictxt());
    check_mpi_descriptor_for_type(product, left_rows, right_cols, block_size);
    verify_global_matrix(product.gather(),
                         make_host_matmul_reference<T>(left_rows,
                                                       shared_cols,
                                                       right_cols));
    EXPECT_EQ(left.desc(), left_descriptor_after);
    EXPECT_EQ(right.desc(), right_descriptor_after);
    EXPECT_EQ(product.desc(), result_descriptor);
}
#endif
} // namespace

TEST(Matrix2dTest, HostTransposeAndDaggerMatchReferences)
{
    check_host_transpose<double>(5, 3);
    check_host_transpose<double>(3, 5);
    check_host_transpose<std::complex<double>>(5, 3);
    check_host_transpose<std::complex<double>>(3, 5);
}

TEST(Matrix2dTest, HostMatmulMatchesTripleLoopReference)
{
    check_host_matmul<double>(5, 4, 6);
    check_host_matmul<std::complex<double>>(5, 4, 6);
}

TEST(Matrix2dTest, SerialDegenerateDistributionTransposeAndMatmul)
{
    const std::shared_ptr<const Parallel_Blkcyclic> transpose_distribution
        = make_frozen_distribution(5, 3);
    Matrix2d<double> transpose_source(transpose_distribution);
    fill_global_values(transpose_source);
    verify_global_matrix(transpose_source.transpose().gather(),
                         make_host_transpose_reference<double>(5, 3));
    verify_global_matrix(transpose_source.dagger().gather(),
                         make_host_transpose_reference<double>(5, 3));

    const std::shared_ptr<const Parallel_Blkcyclic> left_distribution
        = make_frozen_distribution(5, 4);
    const std::shared_ptr<const Parallel_Blkcyclic> right_distribution
        = make_frozen_distribution(4, 6);
    Matrix2d<double> left(left_distribution);
    Matrix2d<double> right(right_distribution);
    fill_global_values(left);
    fill_global_values(right, 10.0);
    verify_global_matrix(left.matmul(right).gather(),
                         make_host_matmul_reference<double>(5, 4, 6));

    const std::shared_ptr<const Parallel_Blkcyclic> complex_transpose_distribution
        = make_frozen_distribution(5, 3);
    Matrix2d<std::complex<double>> complex_transpose_source(complex_transpose_distribution);
    fill_global_values(complex_transpose_source);
    verify_global_matrix(
        complex_transpose_source.transpose().gather(),
        make_host_transpose_reference<std::complex<double>>(5, 3));
    verify_global_matrix(
        complex_transpose_source.dagger().gather(),
        make_dagger_reference(complex_transpose_source.gather()));

    const std::shared_ptr<const Parallel_Blkcyclic> complex_left_distribution
        = make_frozen_distribution(5, 4);
    const std::shared_ptr<const Parallel_Blkcyclic> complex_right_distribution
        = make_frozen_distribution(4, 6);
    Matrix2d<std::complex<double>> complex_left(complex_left_distribution);
    Matrix2d<std::complex<double>> complex_right(complex_right_distribution);
    fill_global_values(complex_left);
    fill_global_values(complex_right, 10.0);
    verify_global_matrix(
        complex_left.matmul(complex_right).gather(),
        make_host_matmul_reference<std::complex<double>>(5, 4, 6));
}

#ifdef __MPI
TEST(Matrix2dTest, TransposeAndDaggerDistributedMatchHost)
{
    const std::vector<std::pair<int, int>> grids = available_process_grids();
    for (std::size_t grid_index = 0; grid_index < grids.size(); ++grid_index)
    {
        check_distributed_transpose<double>(5, 3,
                                            grids[grid_index].first,
                                            grids[grid_index].second,
                                            3);
        check_distributed_transpose<double>(3, 5,
                                            grids[grid_index].first,
                                            grids[grid_index].second,
                                            3);
        check_distributed_transpose<std::complex<double>>(5, 3,
                                                          grids[grid_index].first,
                                                          grids[grid_index].second,
                                                          3);
        check_distributed_transpose<std::complex<double>>(3, 5,
                                                          grids[grid_index].first,
                                                          grids[grid_index].second,
                                                          3);
    }
}

TEST(Matrix2dTest, DistributedResultRetainsBlacsContext)
{
    const std::vector<std::pair<int, int>> grids = available_process_grids();
    if (grids.empty())
    {
        return;
    }

    const int grid_rows = grids.front().first;
    const int grid_cols = grids.front().second;
    const Matrix2d<double> transposed
        = make_transpose_without_source<double>(5, 3, grid_rows, grid_cols, 3);

    ASSERT_TRUE(transposed.is_distributed());
    ASSERT_NE(transposed.distribution(), nullptr);
    EXPECT_NE(transposed.distribution()->comm(), MPI_COMM_NULL);
    verify_global_matrix(transposed.gather(),
                         make_host_transpose_reference<double>(5, 3));
}

TEST(Matrix2dTest, MatmulDistributedMatchesHost)
{
    const std::vector<std::pair<int, int>> grids = available_process_grids();
    for (std::size_t grid_index = 0; grid_index < grids.size(); ++grid_index)
    {
        check_distributed_matmul<double>(5, 4, 6,
                                         grids[grid_index].first,
                                         grids[grid_index].second,
                                         3);
        check_distributed_matmul<std::complex<double>>(5, 4, 6,
                                                       grids[grid_index].first,
                                                       grids[grid_index].second,
                                                       3);
    }
}

TEST(Matrix2dTest, UnsupportedDistributedTypesThrow)
{
    const std::vector<std::pair<int, int>> grids = available_process_grids();
    const std::pair<int, int>& grid = grids[0];
    {
        Matrix2d<int> source(5, 3, grid.first, grid.second, 3, MPI_COMM_WORLD);
        EXPECT_THROW(source.transpose(), std::invalid_argument);
    }
    {
        Matrix2d<float> source(5, 3, grid.first, grid.second, 3, MPI_COMM_WORLD);
        EXPECT_THROW(source.transpose(), std::invalid_argument);
    }
    {
        Matrix2d<std::complex<float>> source(5, 3, grid.first, grid.second, 3, MPI_COMM_WORLD);
        EXPECT_THROW(source.transpose(), std::invalid_argument);
    }
    {
        Matrix2d<int> left(5, 4, grid.first, grid.second, 3, MPI_COMM_WORLD);
        Matrix2d<int> right(4, 6, grid.first, grid.second, 3, MPI_COMM_WORLD);
        EXPECT_THROW(left.matmul(right), std::invalid_argument);
    }
    {
        Matrix2d<float> left(5, 4, grid.first, grid.second, 3, MPI_COMM_WORLD);
        Matrix2d<float> right(4, 6, grid.first, grid.second, 3, MPI_COMM_WORLD);
        EXPECT_THROW(left.matmul(right), std::invalid_argument);
    }
    {
        Matrix2d<std::complex<float>> left(5, 4, grid.first, grid.second, 3, MPI_COMM_WORLD);
        Matrix2d<std::complex<float>> right(4, 6, grid.first, grid.second, 3, MPI_COMM_WORLD);
        EXPECT_THROW(left.matmul(right), std::invalid_argument);
    }
}

TEST(Matrix2dTest, IncompatibleDistributionThrows)
{
    const std::vector<std::pair<int, int>> grids = available_process_grids();
    const std::pair<int, int>& grid = grids[0];
    const int block = 3;
    {
        Matrix2d<double> first(8, 6, grid.first, grid.second, block, MPI_COMM_WORLD);
        Matrix2d<double> different_block(6, 8, grid.first, grid.second, block + 1, MPI_COMM_WORLD);
        fill_global_values(first);
        fill_global_values(different_block, 10.0);
        EXPECT_THROW(first.matmul(different_block), std::invalid_argument);
    }
    if (grids.size() > 1)
    {
        const std::pair<int, int>& other_grid = grids[1];
        Matrix2d<double> first(8, 6, grid.first, grid.second, block, MPI_COMM_WORLD);
        Matrix2d<double> different_grid(6, 8, other_grid.first, other_grid.second, block, MPI_COMM_WORLD);
        fill_global_values(first);
        fill_global_values(different_grid, 10.0);
        EXPECT_THROW(first.matmul(different_grid), std::invalid_argument);
    }
    {
        Matrix2d<double> first(8, 6, grid.first, grid.second, block, MPI_COMM_WORLD);
        Matrix2d<double> different_context(6, 8, grid.first, grid.second, block, MPI_COMM_WORLD);
        fill_global_values(first);
        fill_global_values(different_context, 10.0);
        EXPECT_THROW(first.matmul(different_context), std::invalid_argument);
    }
}

TEST(Matrix2dTest, DescriptorAddressStableAcrossOperators)
{
    const std::vector<std::pair<int, int>> grids = available_process_grids();
    const std::pair<int, int>& grid = grids[0];
    const int block = 3;
    Matrix2d<double> left(5, 4, grid.first, grid.second, block, MPI_COMM_WORLD);
    const int left_blacs_context = left.ictxt();
    const std::shared_ptr<const Parallel_Blkcyclic> right_distribution
        = make_borrowed_distribution(4, 6, block, left_blacs_context);
    Matrix2d<double> right(right_distribution);
    fill_global_values(left);
    fill_global_values(right, 10.0);

    const int* const left_descriptor = left.desc();
    const int* const right_descriptor = right.desc();

    const Matrix2d<double> transposed = left.transpose();
    const Matrix2d<double> product = left.matmul(right);

    EXPECT_EQ(left.desc(), left_descriptor);
    EXPECT_EQ(right.desc(), right_descriptor);
    EXPECT_EQ(left.ictxt(), transposed.ictxt());
    EXPECT_EQ(left.ictxt(), product.ictxt());

    verify_global_matrix(product.gather(),
                         make_host_matmul_reference<double>(5, 4, 6));
}
#endif

TEST(Matrix2dTest, HostConstructionAndRowMajorLayout)
{
    check_host_layout();
}

TEST(Matrix2dTest, HostCopyDeepCopiesStorage)
{
    Matrix2d<double> source(2, 2);
    source.at_local(0, 0) = 1.0;
    source.at_local(0, 1) = 2.0;
    source.at_local(1, 0) = 3.0;
    source.at_local(1, 1) = 4.0;

    Matrix2d<double> copy(source);
    EXPECT_EQ(copy.distribution(), nullptr);
    EXPECT_NE(copy.data(), source.data());
    EXPECT_EQ(copy.at_local(1, 1), 4.0);
    copy.at_local(1, 1) = 9.0;
    EXPECT_EQ(source.at_local(1, 1), 4.0);

    Matrix2d<double> assigned(1, 1);
    assigned = source;
    EXPECT_NE(assigned.data(), source.data());
    EXPECT_EQ(assigned.at_local(0, 1), 2.0);
}

TEST(Matrix2dTest, SerialDistributionHandleUsesColumnMajorLayout)
{
    const std::shared_ptr<const Parallel_Blkcyclic> distribution
        = make_frozen_distribution(3, 2);
    Matrix2d<double> matrix(distribution);

    EXPECT_TRUE(matrix.is_distributed());
    EXPECT_EQ(matrix.distribution().get(), distribution.get());
    EXPECT_EQ(matrix.ictxt(), -1);
    check_distribution_layout(matrix, 3, 2);
}

TEST(Matrix2dTest, CopiesShareDistributionAndDeepCopyStorage)
{
    const std::shared_ptr<const Parallel_Blkcyclic> distribution
        = make_frozen_distribution(3, 2);
    Matrix2d<double> source(distribution);
    source.at_local(0, 0) = 7.0;

    Matrix2d<double> copy(source);
    EXPECT_EQ(copy.distribution().get(), distribution.get());
    EXPECT_EQ(source.distribution().get(), copy.distribution().get());
    EXPECT_NE(copy.data(), source.data());
    EXPECT_EQ(copy.at_local(0, 0), 7.0);
    copy.at_local(0, 0) = 8.0;
    EXPECT_EQ(source.at_local(0, 0), 7.0);
}

TEST(Matrix2dTest, MatricesFromSameHandleShareContext)
{
    const std::shared_ptr<const Parallel_Blkcyclic> distribution
        = make_frozen_distribution(5, 4);
    Matrix2d<double> first(distribution);
    Matrix2d<double> second(distribution);

    EXPECT_EQ(first.distribution().get(), second.distribution().get());
    EXPECT_EQ(first.ictxt(), second.ictxt());
    EXPECT_NE(first.data(), second.data());
}

TEST(Matrix2dTest, RejectsNullDistribution)
{
    const std::shared_ptr<const Parallel_Blkcyclic> distribution;
    EXPECT_THROW(
        {
            Matrix2d<double> matrix(distribution);
            static_cast<void>(matrix);
        },
        std::invalid_argument);
}

#ifdef __MPI
TEST(Matrix2dTest, BlockSizeConstructorAtLocalAndAtGlobal)
{
    const int global_rows = 8;
    const int global_cols = 6;
    const int block_size = 2;
    Matrix2d<double> matrix(global_rows, global_cols, block_size, MPI_COMM_WORLD);

    check_mpi_descriptor(matrix, global_rows, global_cols, block_size);
    check_distribution_layout(matrix, global_rows, global_cols);
}

TEST(Matrix2dTest, UnitBlockConstructorAtLocalAndAtGlobal)
{
    const int global_rows = 5;
    const int global_cols = 3;
    Matrix2d<double> matrix(global_rows, global_cols, MPI_COMM_WORLD);

    check_mpi_descriptor(matrix, global_rows, global_cols, 1);
    check_distribution_layout(matrix, global_rows, global_cols);
}

TEST(Matrix2dTest, FrozenHandleConstructorSharesBlacsContext)
{
    const int global_rows = 8;
    const int global_cols = 6;
    const int block_size = 2;
    const std::shared_ptr<const Parallel_Blkcyclic> distribution
        = make_frozen_distribution(global_rows, global_cols, block_size, MPI_COMM_WORLD);
    Matrix2d<double> first(distribution);
    Matrix2d<double> second(distribution);

    const int distribution_context = distribution->ictxt();
    EXPECT_GE(distribution_context, 0);
    EXPECT_EQ(first.distribution().get(), distribution.get());
    EXPECT_EQ(second.distribution().get(), distribution.get());
    EXPECT_EQ(first.ictxt(), distribution_context);
    EXPECT_EQ(second.ictxt(), distribution_context);
    check_mpi_descriptor(first, global_rows, global_cols, block_size);
}
#endif

TEST(Matrix2dTest, DistributedReshapeIsRejectedAtCompileTime)
{
    EXPECT_FALSE(HasReshape<Matrix2d<double>>::value);
}

TEST(Matrix2dTest, HostElementwiseReductionsAndGatherScatter)
{
    Matrix2d<double> matrix(3, 4);
    Matrix2d<double> other(3, 4);
    fill_global_values(matrix);
    fill_global_values(other, 10.0);
    check_elementwise_operations(matrix, other);
    check_reductions(matrix, other);

    const Matrix2d<double> gathered = matrix.gather();
    EXPECT_NE(gathered.data(), matrix.data());
    verify_global_matrix(gathered, make_host_reference<double>(3, 4));

    Matrix2d<double> restored(3, 4);
    restored.scatter(gathered);
    verify_global_matrix(restored, make_host_reference<double>(3, 4));

    const Matrix2d<double> root_result = matrix.gather_on_root(0);
    verify_global_matrix(root_result, make_host_reference<double>(3, 4));
}

TEST(Matrix2dTest, HostComplexElementwiseAndReductions)
{
    Matrix2d<std::complex<double>> matrix(3, 4);
    Matrix2d<std::complex<double>> other(3, 4);
    fill_global_values(matrix);
    fill_global_values(other, 10.0);
    report_reduction_tolerances();
    check_elementwise_operations(matrix, other);
    check_reductions(matrix, other);
}

TEST(Matrix2dTest, SerialDistributionGatherReturnsHostRowMajorCopy)
{
    const std::shared_ptr<const Parallel_Blkcyclic> distribution
        = make_frozen_distribution(3, 2);
    Matrix2d<double> matrix(distribution);
    fill_global_values(matrix);

    const Matrix2d<double> gathered = matrix.gather();
    EXPECT_FALSE(gathered.is_distributed());
    verify_global_matrix(gathered, make_host_reference<double>(3, 2));

    const Matrix2d<double> root_gathered = matrix.gather_on_root(0);
    EXPECT_FALSE(root_gathered.is_distributed());
    verify_global_matrix(root_gathered, make_host_reference<double>(3, 2));
}

#ifdef __MPI
TEST(Matrix2dTest, DistributedDoubleAcrossBlocksAndProcessGrids)
{
    report_reduction_tolerances();
    const std::vector<std::pair<int, int>> grids = available_process_grids();
    const int blocks[] = {2, 3, 4};
    for (std::size_t grid_index = 0; grid_index < grids.size(); ++grid_index)
    {
        for (std::size_t block_index = 0; block_index < 3; ++block_index)
        {
            Matrix2d<double> matrix(9,
                                    7,
                                    grids[grid_index].first,
                                    grids[grid_index].second,
                                    blocks[block_index],
                                    MPI_COMM_WORLD);
            Matrix2d<double> other(matrix.distribution());
            fill_global_values(matrix);
            fill_global_values(other, 10.0);
            check_elementwise_operations(matrix, other);
            check_reductions(matrix, other);
            check_gather_scatter(matrix);
        }
    }
}

TEST(Matrix2dTest, DistributedComplexAcrossBlocksAndProcessGrids)
{
    report_reduction_tolerances();
    const std::vector<std::pair<int, int>> grids = available_process_grids();
    const int blocks[] = {2, 3, 4};
    for (std::size_t grid_index = 0; grid_index < grids.size(); ++grid_index)
    {
        for (std::size_t block_index = 0; block_index < 3; ++block_index)
        {
            Matrix2d<std::complex<double>> matrix(9,
                                                  7,
                                                  grids[grid_index].first,
                                                  grids[grid_index].second,
                                                  blocks[block_index],
                                                  MPI_COMM_WORLD);
            Matrix2d<std::complex<double>> other(matrix.distribution());
            fill_global_values(matrix);
            fill_global_values(other, 10.0);
            check_elementwise_operations(matrix, other);
            check_reductions(matrix, other);
            check_gather_scatter(matrix);
        }
    }
}

TEST(Matrix2dTest, PaddingDoesNotAffectOperatorsOrGather)
{
    const int number_of_ranks = Parallel_Common::comm_size(MPI_COMM_WORLD);
    Matrix2d<double> matrix(9, 4, 1, number_of_ranks, 2, MPI_COMM_WORLD);
    Matrix2d<double> other(matrix.distribution());
    fill_global_values(matrix);
    fill_global_values(other, 10.0);

    const double sentinel = -98765.25;
    const std::size_t padding_elements = set_padding_sentinel(matrix, sentinel);
    ASSERT_GT(padding_elements, 0U);
    set_padding_sentinel(other, 12345.5);

    matrix += other;
    Matrix2d<double> expected = make_host_reference<double>(9, 4);
    const Matrix2d<double> other_reference = make_host_reference<double>(9, 4, 10.0);
    for (int row = 0; row < 9; ++row)
    {
        for (int col = 0; col < 4; ++col)
        {
            expected.at_local(row, col) += other_reference.at_local(row, col);
        }
    }
    verify_global_matrix(matrix.gather(), expected);
    expect_padding_sentinel(matrix, padding_elements, sentinel);

    matrix.fill(2.0);
    expect_padding_sentinel(matrix, padding_elements, sentinel);
    matrix.scale(3.0);
    expect_padding_sentinel(matrix, padding_elements, sentinel);

    const Matrix2d<double> scatter_reference = make_host_reference<double>(9, 4);
    matrix.scatter(scatter_reference);
    verify_global_matrix(matrix.gather(), scatter_reference);
    expect_padding_sentinel(matrix, padding_elements, sentinel);
}

TEST(Matrix2dTest, RejectsIncompatibleDistributions)
{
    const int number_of_ranks = Parallel_Common::comm_size(MPI_COMM_WORLD);
    const std::vector<std::pair<int, int>> grids = available_process_grids();
    Matrix2d<double> first(8,
                           6,
                           grids.front().first,
                           grids.front().second,
                           2,
                           MPI_COMM_WORLD);
    Matrix2d<double> same_layout_other_handle(8,
                                              6,
                                              grids.front().first,
                                              grids.front().second,
                                              2,
                                              MPI_COMM_WORLD);
    EXPECT_THROW(first += same_layout_other_handle, std::invalid_argument);

    Matrix2d<double> different_block(8,
                                     6,
                                     grids.front().first,
                                     grids.front().second,
                                     3,
                                     MPI_COMM_WORLD);
    EXPECT_THROW(first += different_block, std::invalid_argument);

    Matrix2d<double> host(8, 6);
    EXPECT_THROW(first += host, std::invalid_argument);
    EXPECT_THROW(first.dot(host), std::invalid_argument);

    if (number_of_ranks > 1)
    {
        const std::pair<int, int>& other_grid = grids[1];
        Matrix2d<double> different_grid(8,
                                        6,
                                        other_grid.first,
                                        other_grid.second,
                                        2,
                                        MPI_COMM_WORLD);
        EXPECT_THROW(first += different_grid, std::invalid_argument);
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
