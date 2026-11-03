#include "../invmat_ri.h"
#include "../lri_cv_tools.h"
#include "source_base/parallel_global.h"
#include <gtest/gtest.h>
#include <complex>

namespace
{
template<typename T>
void check_derivative(const double cutoff, const T phase)
{
    const double angle = 0.31;
    const double step = 1e-5;
    auto metric = [phase](const double x) {
        RI::Tensor<T> q({3, 3});
        q(0, 0) = std::cos(x);
        q(0, 2) = -std::sin(x) * phase;
        q(1, 1) = 1;
        q(2, 0) = std::sin(x) * RI::Global_Func::get_conj(phase);
        q(2, 2) = std::cos(x);
        RI::Tensor<T> diagonal({3, 3});
        diagonal(0, 0) = 2;
        diagonal(1, 1) = 2;
        diagonal(2, 2) = 0.01;
        return q * diagonal * q.conjugate().transpose();
    };
    RI::Tensor<T> generator({3, 3});
    generator(0, 2) = -phase;
    generator(2, 0) = RI::Global_Func::get_conj(phase);
    const auto v = metric(angle);
    const auto dv = generator * v - v * generator;
    const std::array<RI::Tensor<T>, 3> directions = {dv, RI::Tensor<T>({3, 3}), RI::Tensor<T>({3, 3})};
    Inverse_Matrix<T> inverse;
    inverse.input(v);
    const auto derivatives = inverse.cal_inverse_with_derivatives(directions, cutoff);
    Inverse_Matrix<T> plus;
    const auto displaced_plus = metric(angle + step);
    plus.input(displaced_plus);
    plus.cal_inverse(Inverse_Matrix<T>::Method::syev, cutoff);
    Inverse_Matrix<T> minus;
    const auto displaced_minus = metric(angle - step);
    minus.input(displaced_minus);
    minus.cal_inverse(Inverse_Matrix<T>::Method::syev, cutoff);
    const auto fd = (plus.output() - minus.output()) * T(0.5 / step);
    const auto response = inverse.output() * v;
    for (std::size_t i = 0; i < 9; ++i)
    {
        EXPECT_NEAR(std::abs(derivatives[0].ptr()[i] - fd.ptr()[i]), 0, 2e-6);
        EXPECT_NEAR(std::abs(derivatives[1].ptr()[i]), 0, 1e-12);
        EXPECT_NEAR(std::abs(derivatives[2].ptr()[i]), 0, 1e-12);
    }
    const double expected_trace = cutoff > 0 ? 2 : 3;
    const T trace = response(0, 0) + response(1, 1) + response(2, 2);
    EXPECT_NEAR(std::abs(trace - T(expected_trace)), 0, 1e-10);
}
}

TEST(SpectralInverse, RealRetainedDiscardedRotation) { check_derivative<double>(0.1, 1.0); }
TEST(SpectralInverse, ComplexRetainedDiscardedRotation)
{
    const auto phase = std::polar(1.0, 0.7);
    check_derivative<std::complex<double>>(0.1, phase);
}
TEST(SpectralInverse, RealFullInverse) { check_derivative<double>(0.0, 1.0); }
TEST(SpectralInverse, ComplexFullInverse)
{
    const auto phase = std::polar(1.0, 0.7);
    check_derivative<std::complex<double>>(0.0, phase);
}

namespace
{
template<typename T>
LRI_CV_Tools::RIBlockMatrix<T> two_center_metric(const std::array<double, 3>& x, const T phase)
{
    LRI_CV_Tools::RIBlockMatrix<T> metric = {
        {RI::Tensor<T>({2, 2}), RI::Tensor<T>({2, 1})},
        {RI::Tensor<T>({1, 2}), RI::Tensor<T>({1, 1})}};
    metric[0][0](0, 0) = 3.0;
    metric[0][0](0, 1) = 0.2;
    metric[0][0](1, 0) = 0.2;
    metric[0][0](1, 1) = 1.8;
    metric[1][1](0, 0) = 0.8;
    metric[0][1](0, 0) = 0.35 + 0.11 * x[0] - 0.07 * x[1] + 0.04 * x[2];
    metric[0][1](1, 0) = phase * (0.18 - 0.06 * x[0] + 0.05 * x[1] + 0.09 * x[2]);
    metric[1][0] = metric[0][1].conjugate().transpose();
    return metric;
}

template<typename T>
std::vector<RI::Tensor<T>> three_center_derivative(const std::size_t axis, const T phase)
{
    const std::array<double, 3> slope = {0.1, -0.08, 0.05};
    std::vector<RI::Tensor<T>> derivative = {RI::Tensor<T>({2, 1, 2}), RI::Tensor<T>({1, 1, 2})};
    for (std::size_t block = 0; block < 2; ++block)
    {
        for (std::size_t i = 0; i < derivative[block].shape[0]; ++i)
        {
            const std::size_t row = block == 0 ? i : 2;
            for (std::size_t j = 0; j < 2; ++j)
            {
                derivative[block](i, 0, j) = phase * (slope[axis] * (row + 1) * (j + 0.5));
            }
        }
    }
    return derivative;
}

template<typename T>
std::vector<RI::Tensor<T>> three_center(const std::array<double, 3>& x, const T phase)
{
    std::vector<RI::Tensor<T>> a = {RI::Tensor<T>({2, 1, 2}), RI::Tensor<T>({1, 1, 2})};
    for (std::size_t block = 0; block < 2; ++block)
    {
        for (std::size_t i = 0; i < a[block].get_shape_all(); ++i)
        {
            a[block].ptr()[i] = T(0.3 + 0.1 * i + 0.2 * block);
        }
    }
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
        const auto derivative = three_center_derivative(axis, phase);
        for (std::size_t block = 0; block < 2; ++block)
        {
            a[block] += derivative[block] * T(x[axis]);
        }
    }
    return a;
}

template<typename T>
void check_coefficients(const double cutoff, const T phase)
{
    const double step = 1e-5;
    const std::array<double, 3> x = {0.13, -0.11, 0.07};
    const auto metric = two_center_metric(x, phase);
    const auto a = three_center(x, phase);
    const std::array<double, 3> first_slope = {0.11, -0.07, 0.04};
    const std::array<double, 3> second_slope = {-0.06, 0.05, 0.09};
    std::array<RI::Tensor<T>, 3> derivative01;
    std::array<RI::Tensor<T>, 3> derivative10;
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
        derivative01[axis] = RI::Tensor<T>({2, 1});
        derivative01[axis](0, 0) = first_slope[axis];
        derivative01[axis](1, 0) = phase * second_slope[axis];
        derivative10[axis] = derivative01[axis].conjugate().transpose();
    }
    const auto inverse = LRI_CV_Tools::cal_I_dI(metric, derivative01, derivative10, cutoff);
    const auto reference_inverse = LRI_CV_Tools::cal_I(metric, Inverse_Matrix<T>::Method::syev, cutoff);
    for (std::size_t row = 0; row < 2; ++row)
    {
        for (std::size_t col = 0; col < 2; ++col)
        {
            EXPECT_EQ(inverse.first[row][col].shape[0], row == 0 ? 2 : 1);
            EXPECT_EQ(inverse.first[row][col].shape[1], col == 0 ? 2 : 1);
            for (std::size_t i = 0; i < inverse.first[row][col].get_shape_all(); ++i)
            {
                EXPECT_NEAR(std::abs(inverse.first[row][col].ptr()[i]
                                     - reference_inverse[row][col].ptr()[i]), 0, 1e-12);
            }
        }
    }
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
        auto plus_x = x;
        auto minus_x = x;
        plus_x[axis] += step;
        minus_x[axis] -= step;
        const auto plus_metric = two_center_metric(plus_x, phase);
        const auto minus_metric = two_center_metric(minus_x, phase);
        const auto plus_a = three_center(plus_x, phase);
        const auto minus_a = three_center(minus_x, phase);
        const auto plus_inverse = LRI_CV_Tools::cal_I(plus_metric, Inverse_Matrix<T>::Method::syev, cutoff);
        const auto minus_inverse = LRI_CV_Tools::cal_I(minus_metric, Inverse_Matrix<T>::Method::syev, cutoff);
        const auto plus_c = LRI_CV_Tools::mul2(plus_inverse, plus_a);
        const auto minus_c = LRI_CV_Tools::mul2(minus_inverse, minus_a);
        const auto da = three_center_derivative(axis, phase);
        const auto first_term = LRI_CV_Tools::mul2(inverse.first, da);
        const auto second_term = LRI_CV_Tools::mul2(inverse.second[axis], a);
        for (std::size_t block = 0; block < 2; ++block)
        {
            const auto analytic = first_term[block] + second_term[block];
            const auto fd = (plus_c[block] - minus_c[block]) * T(0.5 / step);
            for (std::size_t i = 0; i < analytic.get_shape_all(); ++i)
            {
                EXPECT_NEAR(std::abs(analytic.ptr()[i] - fd.ptr()[i]), 0, 2e-8);
            }
        }
        for (std::size_t row = 0; row < 2; ++row)
        {
            for (std::size_t col = 0; col < 2; ++col)
            {
                const auto fd = (plus_inverse[row][col] - minus_inverse[row][col]) * T(0.5 / step);
                for (std::size_t i = 0; i < fd.get_shape_all(); ++i)
                {
                    EXPECT_NEAR(std::abs(inverse.second[axis][row][col].ptr()[i] - fd.ptr()[i]), 0, 2e-8);
                }
            }
        }
        // Fixed diagonal blocks and changing off-diagonal blocks vary the eigenvalues.
        const double plus_norm = std::norm(plus_metric[0][1](0, 0)) + std::norm(plus_metric[0][1](1, 0));
        const double minus_norm = std::norm(minus_metric[0][1](0, 0)) + std::norm(minus_metric[0][1](1, 0));
        EXPECT_GT(std::abs(plus_norm - minus_norm), 1e-9);
    }
}
}

TEST(SpectralInverse, RealUnequalBlocksWithMovingThreeCenter)
{
    check_coefficients<double>(0.3, 1.0);
}
TEST(SpectralInverse, ComplexUnequalBlocksWithMovingThreeCenter)
{
    const auto phase = std::polar(1.0, 0.7);
    check_coefficients<std::complex<double>>(0.3, phase);
}
TEST(SpectralInverse, RealUnequalBlocksWithoutTruncation)
{
    check_coefficients<double>(0.0, 1.0);
}
TEST(SpectralInverse, ComplexUnequalBlocksWithoutTruncation)
{
    const auto phase = std::polar(1.0, 0.7);
    check_coefficients<std::complex<double>>(0.0, phase);
}

int main(int argc, char** argv)
{
    int processes = 1;
    int threads = 1;
    int rank = 0;
    Parallel_Global::read_pal_param(argc, argv, processes, threads, rank);
    POOL_WORLD = MPI_COMM_NULL;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    Parallel_Global::finalize_mpi();
    return result;
}
