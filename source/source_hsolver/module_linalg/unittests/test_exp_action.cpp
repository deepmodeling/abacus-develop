#include "source_base/module_container/base/third_party/lapack.h"
#include "source_hsolver/module_linalg/exp_action.h"
#include "source_hsolver/module_linalg/unittests/linear_test_utils.h"

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <type_traits>

namespace
{
using Wide = std::complex<double>;
using Device = base_device::DEVICE_CPU;

// A full complex Hermitian matrix, with spin components separated by padding.
template <typename T>
class DenseExp final : public hsolver::HSOperator<T>
{
  public:
    std::vector<T> matrix;
    int npol;
    mutable int calls = 0;
    mutable int cols = 0;

    explicit DenseExp(int spins) : npol(spins)
    {
        matrix = linear_test::hermitian<T>(6);
    }
    int row(int i, int ld) const
    {
        return i / (6 / npol) * (ld / npol) + i % (6 / npol);
    }
    void update_k(int) override
    {
    }
    void hpsi(const T* x, T* y, int ld, int bands) const override
    {
        ++calls;
        cols += bands;
        std::fill_n(y, ld * bands, T(0));
        for (int b = 0; b < bands; ++b)
        {
            for (int i = 0; i < 6; ++i)
            {
                for (int j = 0; j < 6; ++j)
                {
                    y[b * ld + row(i, ld)] += matrix[j * 6 + i] * x[b * ld + row(j, ld)];
                }
            }
        }
    }
    void spsi(const T*, T*, int, int) const override
    {
        FAIL() << "PW exponential must not request an overlap operator";
    }
};

template <typename T>
std::vector<T> pack_input(const DenseExp<T>& op, const hsolver::OrbLayout& layout)
{
    std::vector<T> input(layout.ld * layout.bands, T(30, -20));
    for (int b = 0; b < layout.bands; ++b)
    {
        for (int i = 0; i < 6; ++i)
        {
            input[b * layout.ld + op.row(i, layout.ld)] = T(std::sin(i + 2 * b + 1), std::cos(2 * i + b));
        }
    }
    return input;
}

// Diagonalize the entire matrix, independently of the Krylov approximation.
template <typename T>
std::vector<Wide> exact_exp(const std::vector<T>& matrix, const std::vector<Wide>& input, double alpha)
{
    const int n = 6;
    const int lwork = 64;
    const char job = 'V';
    const char upper = 'U';
    std::vector<Wide> vectors(matrix.begin(), matrix.end());
    std::vector<Wide> work(lwork);
    std::vector<double> rwork(3 * n);
    std::vector<double> values(n);
    int status = 0;
    ct::lapackConnector::heev(job, upper, n, vectors.data(), n, values.data(), work.data(), lwork, rwork.data(), status);
    EXPECT_EQ(status, 0);
    std::vector<Wide> output(n, Wide(0));
    for (int k = 0; k < n; ++k)
    {
        Wide weight(0);
        for (int i = 0; i < n; ++i)
        {
            weight += std::conj(vectors[k * n + i]) * input[i];
        }
        weight *= std::exp(Wide(0, -alpha * values[k]));
        for (int i = 0; i < n; ++i)
        {
            output[i] += vectors[k * n + i] * weight;
        }
    }
    return output;
}

template <typename T>
void check_output(const DenseExp<T>& op,
                  const hsolver::OrbLayout& layout,
                  const std::vector<T>& input,
                  const std::vector<T>& output,
                  double alpha,
                  int order,
                  double tol)
{
    for (int b = 0; b < layout.bands; ++b)
    {
        std::vector<Wide> src(6);
        for (int i = 0; i < 6; ++i)
        {
            src[i] = input[b * layout.ld + op.row(i, layout.ld)];
        }
        std::vector<Wide> expected;
        if (order == 0)
        {
            expected = exact_exp(op.matrix, src, alpha);
        }
        else
        {
            // Direct dense polynomial in double precision, without batched kernels.
            expected = src;
            std::vector<Wide> term = src;
            const std::vector<Wide> matrix(op.matrix.begin(), op.matrix.end());
            for (int k = 1; k <= order; ++k)
            {
                term = linear_test::multiply(matrix, term, 6, 1);
                for (int i = 0; i < 6; ++i)
                {
                    term[i] *= Wide(0, -alpha / k);
                    expected[i] += term[i];
                }
            }
        }
        double norm_in = 0;
        double norm_out = 0;
        for (int i = 0; i < 6; ++i)
        {
            const Wide actual = output[b * layout.ld + op.row(i, layout.ld)];
            EXPECT_NEAR(std::abs(actual - expected[i]), 0.0, tol);
            norm_in += std::norm(src[i]);
            norm_out += std::norm(actual);
        }
        if (order == 0)
        {
            EXPECT_NEAR(norm_out, norm_in, 20 * tol);
        }
        for (int spin = 0; spin < layout.npol; ++spin)
        {
            for (int i = layout.npw; i < layout.ld / layout.npol; ++i)
            {
                EXPECT_EQ(output[b * layout.ld + spin * (layout.ld / layout.npol) + i], T(0));
            }
        }
    }
}

template <typename Real>
class ExpTest : public testing::Test
{
};
using Precisions = testing::Types<double, float>;
TYPED_TEST_SUITE(ExpTest, Precisions);

TYPED_TEST(ExpTest, LanczosMatchesFullExponential)
{
    using T = std::complex<TypeParam>;
    const double tol = std::is_same<TypeParam, float>::value ? 2e-5 : 2e-10;
    for (int npol: {1, 2})
    {
        SCOPED_TRACE(npol);
        DenseExp<T> op(npol);
        const hsolver::OrbLayout layout{10, 6 / npol, npol, 10};
        const std::vector<T> input = pack_input(op, layout);
        std::vector<T> output(input.size());
        hsolver::ExpControls exp_ctrl;
        std::unique_ptr<hsolver::ExpAction<T, Device>> exp = hsolver::make_exp_action<T, Device>(exp_ctrl, linear_test::world_comm());
        const hsolver::ExpResult result = exp->apply(op, layout, 0.4, input.data(), output.data());
        ASSERT_EQ(result.state, hsolver::ExpState::success);
        EXPECT_TRUE(result.has_err_est);
        EXPECT_LE(result.rel_err_est, exp->rel_tol());
        EXPECT_EQ(result.hpsi_calls, op.calls);
        EXPECT_EQ(result.hpsi_cols, op.cols);
        check_output(op, layout, input, output, 0.4, 0, tol);
    }
}

TYPED_TEST(ExpTest, TaylorPairMatchesBothPolynomials)
{
    using T = std::complex<TypeParam>;
    const double tol = std::is_same<TypeParam, float>::value ? 2e-6 : 2e-14;
    for (int npol: {1, 2})
    {
        SCOPED_TRACE(npol);
        DenseExp<T> op(npol);
        const hsolver::OrbLayout layout{10, 6 / npol, npol, 3};
        const std::vector<T> input = pack_input(op, layout);
        std::vector<T> half(input.size());
        std::vector<T> full(input.size());
        hsolver::ExpControls exp_ctrl;
        exp_ctrl.expansion = hsolver::ExpMethod::taylor;
        exp_ctrl.taylor_order = 6;
        std::unique_ptr<hsolver::ExpAction<T, Device>> exp = hsolver::make_exp_action<T, Device>(exp_ctrl, linear_test::world_comm());
        const hsolver::ExpResult result = exp->apply_pair(op, layout, 0.15, input.data(), half.data(), full.data());
        ASSERT_EQ(result.state, hsolver::ExpState::success);
        EXPECT_FALSE(result.has_err_est);
        EXPECT_EQ(op.calls, exp_ctrl.taylor_order);
        EXPECT_EQ(result.hpsi_calls, op.calls);
        EXPECT_EQ(result.hpsi_cols, exp_ctrl.taylor_order * layout.bands);
        check_output(op, layout, input, half, 0.15, exp_ctrl.taylor_order, tol);
        check_output(op, layout, input, full, 0.3, exp_ctrl.taylor_order, tol);
    }
}

TEST(ExpAction, LanczosSubdividesOrReportsExhaustion)
{
    DenseExp<Wide> op(1);
    const hsolver::OrbLayout layout{8, 6, 1, 2};
    const std::vector<Wide> input = pack_input(op, layout);
    std::vector<Wide> output(input.size());
    hsolver::ExpControls exp_ctrl;
    exp_ctrl.krylov_max = 4;
    exp_ctrl.rel_tol = 1e-8;
    std::unique_ptr<hsolver::ExpAction<Wide, Device>> exp = hsolver::make_exp_action<Wide, Device>(exp_ctrl, linear_test::world_comm());
    const hsolver::ExpResult result = exp->apply(op, layout, 1.0, input.data(), output.data());
    ASSERT_EQ(result.state, hsolver::ExpState::success);
    EXPECT_GT(result.bisect_count, 0);
    EXPECT_LE(result.rel_err_est, exp_ctrl.rel_tol);
    EXPECT_EQ(result.hpsi_calls, op.calls);
    EXPECT_EQ(result.hpsi_cols, op.cols);
    check_output(op, layout, input, output, 1.0, 0, 5e-8);

    exp_ctrl.krylov_max = 2;
    exp_ctrl.rel_tol = 1e-14;
    exp = hsolver::make_exp_action<Wide, Device>(exp_ctrl, linear_test::world_comm());
    const hsolver::ExpResult failed = exp->apply(op, layout, 1.0, input.data(), output.data());
    EXPECT_EQ(failed.state, hsolver::ExpState::failed);
    EXPECT_GT(failed.bisect_count, 0);
    EXPECT_GE(failed.failed_band, 0);
    EXPECT_LT(failed.failed_band, layout.bands);
}

TEST(ExpAction, NonfiniteInputReportsOriginalBand)
{
    DenseExp<Wide> op(2);
    const hsolver::OrbLayout layout{10, 3, 2, 10};
    std::vector<Wide> input = pack_input(op, layout);
    input[9 * layout.ld] = Wide(std::numeric_limits<double>::quiet_NaN(), 0);
    std::vector<Wide> output(input.size());
    for (hsolver::ExpMethod method: {hsolver::ExpMethod::lanczos, hsolver::ExpMethod::taylor})
    {
        hsolver::ExpControls exp_ctrl;
        exp_ctrl.expansion = method;
        std::unique_ptr<hsolver::ExpAction<Wide, Device>> exp = hsolver::make_exp_action<Wide, Device>(exp_ctrl, linear_test::world_comm());
        const hsolver::ExpResult result = exp->apply(op, layout, 0.1, input.data(), output.data());
        EXPECT_EQ(result.state, hsolver::ExpState::failed);
        EXPECT_EQ(result.failed_band, 9);
    }
}
} // namespace
