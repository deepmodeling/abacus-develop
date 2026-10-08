#include "source_hsolver/linear_low_rank.h"
#include "source_hsolver/test/linear_test_utils.h"

#include <gtest/gtest.h>
#include <type_traits>

namespace
{
template <typename Real>
class LinearLowRankTest : public testing::Test
{
};
using Precisions = testing::Types<float, double>;
TYPED_TEST_SUITE(LinearLowRankTest, Precisions);

TYPED_TEST(LinearLowRankTest, CoarseAndHistoryCorrectionsMatchDenseReference)
{
    using T = std::complex<TypeParam>;
    using Wide = linear_test::Complex;
    using Device = base_device::DEVICE_CPU;
    const hsolver::diag_comm_info comm = linear_test::world_comm();
    const int n = 5;
    const int rank = 2;
    const int start = n * comm.rank / comm.nproc;
    const int dim = n * (comm.rank + 1) / comm.nproc - start;
    const int ld = dim + 1;
    const double tolerance = std::is_same<TypeParam, float>::value ? 2e-5 : 1e-12;
    std::vector<T> z(n * rank);
    std::vector<T> w(n * rank);
    std::vector<T> diagonal(n);
    std::vector<T> x(n);
    for (int i = 0; i < n; ++i)
    {
        diagonal[i] = T(0.7 + 0.02 * i, 0.1);
        x[i] = T(0.3 - 0.12 * i, 0.04 * i * i);
        for (int j = 0; j < rank; ++j)
        {
            z[j * n + i] = T(j == 0 ? 1.0 : 0.1 * i, 0.1 * (j + i));
            w[j * n + i] = T(1.0, 0.2 * (i + 1)) * z[j * n + i];
        }
    }
    const std::vector<T> local_z = linear_test::local_columns(z, n, rank, start, dim, ld);
    const std::vector<T> local_w = linear_test::local_columns(w, n, rank, start, dim, ld);
    const std::vector<T> local_d = linear_test::local_columns(diagonal, n, 1, start, dim, ld);
    const std::vector<T> local_x = linear_test::local_columns(x, n, 1, start, dim, ld);
    const std::vector<T> zero(ld * rank, T(0));
    hsolver::LinearAlgebra<T, Device> algebra(comm);
    hsolver::LinearResponse<T, Device> response;
    ct::Tensor workspace;
    response.update(algebra, ld, dim, rank, local_z.data(), zero.data(), local_w.data(), 1e-8, &workspace);
    ASSERT_EQ(response.rank(), rank);
    ct::Tensor correction_workspace;
    for (const bool history: {false, true})
    {
        SCOPED_TRACE(history);
        const std::vector<T>& test = history ? w : z;
        std::vector<Wide> gram(rank * rank, Wide(0));
        std::vector<Wide> projected(rank, Wide(0));
        for (int i = 0; i < rank; ++i)
        {
            for (int row = 0; row < n; ++row)
            {
                projected[i] += std::conj(Wide(test[i * n + row])) * Wide(x[row]);
                for (int j = 0; j < rank; ++j)
                {
                    gram[j * rank + i] += std::conj(Wide(test[i * n + row])) * Wide(w[j * n + row]);
                }
            }
        }
        const std::vector<Wide> coefficients = linear_test::lapack_solve(gram, projected, rank, 1);
        hsolver::LinearSmallLU factor;
        ASSERT_TRUE(factor.factor(gram, rank));
        hsolver::LinearLowRank<T, Device> preconditioner(algebra, local_d.data(), dim);
        if (history)
        {
            preconditioner.prepare_response(ld, response.rank(), response.directions(), response.images(), &correction_workspace);
        }
        else
        {
            preconditioner.prepare_subspace(ld, rank, local_z.data(), local_w.data(), factor, &correction_workspace);
        }
        ASSERT_EQ(preconditioner.rank(), rank);
        std::vector<T> actual(ld * rank, T(0));
        // Both corrections must satisfy P*(A*Z) = Z on the represented space.
        preconditioner.apply(local_w.data(), actual.data(), ld, rank);
        for (int j = 0; j < rank; ++j)
        {
            for (int i = 0; i < dim; ++i)
            {
                EXPECT_LT(std::abs(Wide(actual[j * ld + i]) - Wide(local_z[j * ld + i])), tolerance);
            }
        }
        preconditioner.apply(local_x.data(), actual.data(), ld, 1);
        for (int i = 0; i < dim; ++i)
        {
            const int row = start + i;
            Wide expected = Wide(diagonal[row]) * Wide(x[row]);
            for (int j = 0; j < rank; ++j)
            {
                expected += (Wide(z[j * n + row]) - Wide(diagonal[row]) * Wide(w[j * n + row])) * coefficients[j];
            }
            EXPECT_LT(std::abs(Wide(actual[i]) - expected), tolerance);
        }
        if (history)
        {
            preconditioner.prepare_response(ld, 0, nullptr, nullptr, &correction_workspace);
        }
        else
        {
            preconditioner.prepare_subspace(ld, 0, nullptr, nullptr, factor, &correction_workspace);
        }
        EXPECT_EQ(preconditioner.rank(), 0);
        preconditioner.apply(local_x.data(), actual.data(), ld, 1);
        for (int i = 0; i < dim; ++i)
        {
            EXPECT_LT(std::abs(actual[i] - local_d[i] * local_x[i]), tolerance);
        }
    }
    // A new unresolved response must discard an older usable history.
    response.update(algebra, ld, dim, rank, zero.data(), zero.data(), zero.data(), 1e-8, &workspace);
    EXPECT_EQ(response.rank(), 0);
}
} // namespace
