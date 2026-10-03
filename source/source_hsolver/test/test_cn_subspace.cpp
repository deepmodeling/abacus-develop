#include "source_base/parallel_device.h"
#include "source_hsolver/cn_subspace.h"
#include "source_hsolver/test/linear_test_utils.h"

#include <cmath>
#include <gtest/gtest.h>
#include <type_traits>

namespace
{
template <typename Real>
class CNSubspaceTest : public testing::Test
{
};
using Precisions = testing::Types<float, double>;
TYPED_TEST_SUITE(CNSubspaceTest, Precisions);

TYPED_TEST(CNSubspaceTest, ProjectedSeedAndResidualOrRankDeficientFallback)
{
    using T = std::complex<TypeParam>;
    using Wide = linear_test::Complex;
    const hsolver::diag_comm_info comm = linear_test::world_comm();
    const int n = 5;
    const int bands = 2;
    const int start = n * comm.rank / comm.nproc;
    const int dim = n * (comm.rank + 1) / comm.nproc - start;
    const int ld = dim + 2;
    const double tolerance = std::is_same<TypeParam, float>::value ? 1e-5 : 1e-12;
    const std::vector<T> h = linear_test::hermitian<T>(n);
    std::vector<T> a(n * n);
    std::vector<T> right(n * n);
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            const T identity(i == j ? 1 : 0);
            a[j * n + i] = identity + T(0, 0.15) * h[j * n + i];
            right[j * n + i] = identity - T(0, 0.15) * h[j * n + i];
        }
    }
    hsolver::LinearAlgebra<T, base_device::DEVICE_CPU> algebra(comm);
    hsolver::CNSubspace<T, base_device::DEVICE_CPU> projection;
    for (const bool singular: {false, true})
    {
        SCOPED_TRACE(singular);
        std::vector<T> u(n * bands);
        for (int j = 0; j < bands; ++j)
        {
            for (int i = 0; i < n; ++i)
            {
                const double phase = singular ? 0.0 : 2.0 * std::acos(-1.0) * i * j / n;
                u[j * n + i] = std::polar(TypeParam(1.0 / std::sqrt(n)), TypeParam(phase));
            }
        }
        const std::vector<T> b = linear_test::multiply(right, u, n, bands);
        const std::vector<T> au = linear_test::multiply(a, u, n, bands);
        const std::vector<T> local_u = linear_test::local_columns(u, n, bands, start, dim, ld);
        const std::vector<T> local_b = linear_test::local_columns(b, n, bands, start, dim, ld);
        const bool prepared = projection.prepare(algebra, ld, dim, bands, local_u.data(), local_b.data());
        if (singular)
        {
            EXPECT_FALSE(prepared);
            continue;
        }
        ASSERT_TRUE(prepared);
        std::vector<Wide> projected(bands * bands, Wide(0));
        std::vector<Wide> projected_rhs(bands * bands, Wide(0));
        for (int j = 0; j < bands; ++j)
        {
            for (int i = 0; i < bands; ++i)
            {
                for (int row = 0; row < n; ++row)
                {
                    projected[j * bands + i] += std::conj(Wide(u[i * n + row])) * Wide(au[j * n + row]);
                    projected_rhs[j * bands + i] += std::conj(Wide(u[i * n + row])) * Wide(b[j * n + row]);
                }
            }
        }
        const std::vector<Wide> coefficients = linear_test::lapack_solve(projected, projected_rhs, bands, bands);
        std::vector<Wide> actual(n * bands, Wide(0));
        for (int j = 0; j < bands; ++j)
        {
            for (int i = 0; i < dim; ++i)
            {
                Wide expected(0);
                for (int k = 0; k < bands; ++k)
                {
                    expected += Wide(u[k * n + start + i]) * coefficients[j * bands + k];
                }
                EXPECT_LT(std::abs(Wide(projection.seed()[j * ld + i]) - expected), tolerance);
                EXPECT_LT(std::abs(Wide(projection.image()[j * ld + i]) - Wide(au[j * n + start + i])), tolerance);
                actual[j * n + start + i] = projection.seed()[j * ld + i];
            }
        }
#ifdef __MPI
        const int count = actual.size();
        Parallel_Common::reduce_data(actual.data(), count, comm.comm);
#endif
        const std::vector<Wide> wide_a(a.begin(), a.end());
        const std::vector<Wide> ax = linear_test::multiply(wide_a, actual, n, bands);
        for (int j = 0; j < bands; ++j)
        {
            for (int i = 0; i < dim; ++i)
            {
                const Wide residual = Wide(b[j * n + start + i]) - ax[j * n + start + i];
                EXPECT_LT(std::abs(Wide(projection.residual()[j * ld + i]) - residual), tolerance);
            }
            for (int k = 0; k < bands; ++k)
            {
                Wide overlap(0);
                for (int i = 0; i < n; ++i)
                {
                    overlap += std::conj(Wide(u[k * n + i])) * (Wide(b[j * n + i]) - ax[j * n + i]);
                }
                EXPECT_LT(std::abs(overlap), tolerance);
            }
        }
    }
}
} // namespace
