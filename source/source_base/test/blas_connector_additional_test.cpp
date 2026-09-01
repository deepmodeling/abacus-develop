#include "source_base/module_external/blas_connector.h"

#include "gtest/gtest.h"
#include <cmath>
#include <complex>

namespace
{

template <typename T>
void expect_near(const T& actual, const T& expected)
{
    EXPECT_NEAR(std::abs(actual - expected), 0.0, 1.0e-6);
}

template <typename T>
void exercise_level_one(const T& alpha)
{
    const int one = 1;
    const T x[1] = {static_cast<T>(2)};
    T y[1] = {static_cast<T>(3)};
    BlasConnector::axpy(one, alpha, x, one, y, one);
    expect_near(y[0], alpha * x[0] + static_cast<T>(3));

    T scaled[1] = {static_cast<T>(2)};
    BlasConnector::scal(one, alpha, scaled, one);
    expect_near(scaled[0], alpha * static_cast<T>(2));

    T copied[1] = {};
    BlasConnector::copy(one, x, one, copied, one);
    expect_near(copied[0], x[0]);
}

template <typename T>
void exercise_matrix_operations(const T& alpha)
{
    const T beta = static_cast<T>(1);
    const T a[1] = {static_cast<T>(2)};
    const T b[1] = {static_cast<T>(3)};
    T c[1] = {static_cast<T>(4)};

    BlasConnector::gemm('N', 'N', 1, 1, 1, alpha, a, 1, b, 1, beta, c, 1);
    expect_near(c[0], alpha * a[0] * b[0] + beta * static_cast<T>(4));

    c[0] = static_cast<T>(4);
    BlasConnector::gemm_cm('N', 'N', 1, 1, 1, alpha, a, 1, b, 1, beta, c, 1);
    expect_near(c[0], alpha * a[0] * b[0] + beta * static_cast<T>(4));

    c[0] = static_cast<T>(4);
    BlasConnector::symm_cm('L', 'U', 1, 1, alpha, a, 1, b, 1, beta, c, 1);
    expect_near(c[0], alpha * a[0] * b[0] + beta * static_cast<T>(4));

    T hermitian[1] = {static_cast<T>(2)};
    T right[1] = {static_cast<T>(3)};
    c[0] = static_cast<T>(4);
    BlasConnector::hemm_cm('L', 'U', 1, 1, alpha, hermitian, 1, right, 1, beta, c, 1);
    expect_near(c[0], alpha * hermitian[0] * right[0] + beta * static_cast<T>(4));

    c[0] = static_cast<T>(4);
    BlasConnector::gemv('N', 1, 1, alpha, a, 1, b, 1, beta, c, 1);
    expect_near(c[0], alpha * a[0] * b[0] + beta * static_cast<T>(4));
}

template <typename T, typename Operand>
void exercise_elementwise_operations(const T& left, const Operand& right)
{
    const int one = 1;
    T result[1] = {};
    const T left_values[1] = {left};
    const Operand right_values[1] = {right};

    BlasConnector::vector_mul_vector<T, Operand>(one,
                                                 result,
                                                 left_values,
                                                 right_values,
                                                 base_device::AbacusDevice_t::CpuDevice);
    expect_near(result[0], left * right);

    BlasConnector::vector_div_vector<T, Operand>(one,
                                                 result,
                                                 left_values,
                                                 right_values,
                                                 base_device::AbacusDevice_t::CpuDevice);
    expect_near(result[0], left / right);

    BlasConnector::vector_add_vector<T, Operand>(one,
                                                 result,
                                                 left_values,
                                                 static_cast<Operand>(2),
                                                 left_values,
                                                 static_cast<Operand>(3),
                                                 base_device::AbacusDevice_t::CpuDevice);
    expect_near(result[0], left * static_cast<Operand>(5));
}

} // namespace

TEST(BlasConnectorAdditionalTest, CoversLevelOneOverloads)
{
    exercise_level_one<float>(2.0F);
    exercise_level_one<double>(2.0);
    exercise_level_one<std::complex<float>>(std::complex<float>(2.0F, 0.0F));
    exercise_level_one<std::complex<double>>(std::complex<double>(2.0, 0.0));

    const float float_values[2] = {3.0F, 4.0F};
    const double double_values[2] = {3.0, 4.0};
    const std::complex<float> complex_float_values[2] = {{3.0F, 1.0F}, {4.0F, -1.0F}};
    const std::complex<double> complex_double_values[2] = {{3.0, 1.0}, {4.0, -1.0}};
    EXPECT_FLOAT_EQ(BlasConnector::dot(2, float_values, 1, float_values, 1), 25.0F);
    EXPECT_FLOAT_EQ(BlasConnector::dotu(2, float_values, 1, float_values, 1), 25.0F);
    EXPECT_DOUBLE_EQ(BlasConnector::dotu(2, double_values, 1, double_values, 1), 25.0);
    expect_near(BlasConnector::dotu(2, complex_float_values, 1, complex_float_values, 1),
                std::complex<float>(23.0F, -2.0F));
    expect_near(BlasConnector::dotu(2, complex_double_values, 1, complex_double_values, 1),
                std::complex<double>(23.0, -2.0));
    EXPECT_FLOAT_EQ(BlasConnector::dotc(2, float_values, 1, float_values, 1), 25.0F);
    EXPECT_DOUBLE_EQ(BlasConnector::dotc(2, double_values, 1, double_values, 1), 25.0);
    expect_near(BlasConnector::dotc(2, complex_float_values, 1, complex_float_values, 1),
                std::complex<float>(27.0F, 0.0F));
    expect_near(BlasConnector::dotc(2, complex_double_values, 1, complex_double_values, 1),
                std::complex<double>(27.0, 0.0));
    EXPECT_FLOAT_EQ(BlasConnector::nrm2(2, float_values, 1), 5.0F);
    EXPECT_DOUBLE_EQ(BlasConnector::nrm2(2, complex_double_values, 1), std::sqrt(27.0));
}

TEST(BlasConnectorAdditionalTest, CoversMatrixOverloads)
{
    exercise_matrix_operations<float>(2.0F);
    exercise_matrix_operations<double>(2.0);
    exercise_matrix_operations<std::complex<float>>(std::complex<float>(2.0F, 0.0F));
    exercise_matrix_operations<std::complex<double>>(std::complex<double>(2.0, 0.0));
}

TEST(BlasConnectorAdditionalTest, CoversElementwiseInstantiations)
{
    exercise_elementwise_operations<float, float>(2.0F, 4.0F);
    exercise_elementwise_operations<double, double>(2.0, 4.0);
    exercise_elementwise_operations<std::complex<float>, float>(std::complex<float>(2.0F, 1.0F), 4.0F);
    exercise_elementwise_operations<std::complex<double>, double>(std::complex<double>(2.0, 1.0), 4.0);
    exercise_elementwise_operations<std::complex<float>, std::complex<float>>(std::complex<float>(2.0F, 1.0F),
                                                                              std::complex<float>(4.0F, -1.0F));
}
