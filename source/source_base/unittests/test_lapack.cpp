#include "source_base/macros.h"
#include "source_base/module_container/base/third_party/lapack.h"

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

template <typename T>
class LapackGvxOutput : public ::testing::Test
{
};

using GvxTypes = ::testing::Types<float, double, std::complex<float>, std::complex<double>>;
TYPED_TEST_SUITE(LapackGvxOutput, GvxTypes);

namespace
{
template <typename T>
void check_gvx(const char range, const bool positive)
{
    using Real = typename GetTypeReal<T>::type;
    const int n = 4;
    const int ld = 6;
    const int first = 2;
    const int last = 3;
    const Real lower = 1.5;
    const Real upper = 3.5;
    const Real tolerance = 0;
    const T sentinel(731);
    std::vector<T> a(ld * n, sentinel);
    std::vector<T> b(ld * n, sentinel);
    std::vector<T> vectors(ld * n, sentinel);
    std::vector<Real> values(n, Real(-999));
    std::vector<Real> rwork(7 * n);
    std::vector<int> iwork(5 * n);
    std::vector<int> ifail(n);
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            a[i + j * ld] = i == j ? T(2 * (j + 1)) : T(0);
            b[i + j * ld] = i == j ? T(2) : T(0);
        }
    }
    if (!positive)
    {
        b[0] = T(-2);
    }
    int found = -123;
    int info = -456;
    T query(0);
    const int query_size = -1;
    container::lapackConnector::hegvx(1, 'V', range, 'U', n, a.data(), ld, b.data(), ld,
                                    lower, upper, first, last, tolerance, &found,
                                    values.data(), vectors.data(), ld, &query, query_size,
                                    rwork.data(), iwork.data(), ifail.data(), info);
    ASSERT_EQ(info, 0);
    ASSERT_TRUE(std::isfinite(std::real(query)));
    const int lwork = std::max(8 * n, static_cast<int>(std::real(query)));
    std::vector<T> work(lwork);
    found = -123; // The explicit-count overload must replace this sentinel.
    info = -456;
    container::lapackConnector::hegvx(1, 'V', range, 'U', n, a.data(), ld, b.data(), ld,
                                    lower, upper, first, last, tolerance, &found,
                                    values.data(), vectors.data(), ld, work.data(), lwork,
                                    rwork.data(), iwork.data(), ifail.data(), info);
    if (!positive)
    {
        EXPECT_EQ(info, n + 1); // The leading principal minor of B is not positive.
        return;
    }
    ASSERT_EQ(info, 0);
    const int expected = range == 'A' ? n : 2;
    ASSERT_EQ(found, expected);
    const double error = 100 * std::numeric_limits<Real>::epsilon();
    for (int j = 0; j < found; ++j)
    {
        const double eigenvalue = range == 'A' ? j + 1 : j + 2;
        EXPECT_NEAR(values[j], eigenvalue, error);
        for (int i = n; i < ld; ++i)
        {
            EXPECT_EQ(vectors[i + j * ld], sentinel);
        }
    }
}
} // namespace

TYPED_TEST(LapackGvxOutput, ReportsCountValuesAndOutputBounds)
{
    for (const char range : {'I', 'V', 'A'})
    {
        SCOPED_TRACE(range);
        check_gvx<TypeParam>(range, true);
    }
}

TYPED_TEST(LapackGvxOutput, ReturnsNonPositiveMetricInfo)
{
    check_gvx<TypeParam>('I', false);
}

TYPED_TEST(LapackGvxOutput, LegacyByValueInterfaceKeepsCallerCount)
{
    using T = TypeParam;
    using Real = typename GetTypeReal<T>::type;
    T a = T(6);
    T b = T(2);
    T vector = T(0);
    Real value = 0;
    int found = -123;
    int info = -1;
    T work[64] = {};
    Real rwork[7] = {};
    int iwork[5] = {};
    int ifail[1] = {};
    const Real zero = 0;
    container::lapackConnector::hegvx(1, 'V', 'I', 'U', 1, &a, 1, &b, 1,
                                    zero, zero, 1, 1, zero, found, &value, &vector, 1,
                                    work, 64, rwork, iwork, ifail, info);
    EXPECT_EQ(info, 0);
    EXPECT_EQ(found, -123);
    EXPECT_NEAR(value, 3, 100 * std::numeric_limits<Real>::epsilon());
}
