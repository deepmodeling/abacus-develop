#include "source_hsolver/orthonormal.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

namespace
{
using Wide = std::complex<double>;
using Device = base_device::DEVICE_CPU;
constexpr int dim = 4;
constexpr int ld = 6;
constexpr int bands = 2;

hsolver::diag_comm_info self_comm()
{
#ifdef __MPI
    return hsolver::diag_comm_info(MPI_COMM_SELF, 0, 1);
#else
    return hsolver::diag_comm_info(0, 1);
#endif
}

template <typename T>
std::vector<Wide> overlaps(const std::vector<T>& psi)
{
    std::vector<Wide> gram(bands * bands, Wide(0));
    for (int j = 0; j < bands; ++j)
    {
        for (int i = 0; i < bands; ++i)
        {
            for (int row = 0; row < dim; ++row)
            {
                gram[i + j * bands] += std::conj(Wide(psi[i * ld + row])) * Wide(psi[j * ld + row]);
            }
        }
    }
    return gram;
}

template <typename T>
std::vector<T> orbitals()
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<T> psi(ld * bands, T(nan, nan));
    for (int band = 0; band < bands; ++band)
    {
        std::fill_n(psi.begin() + band * ld, dim, T(0));
    }
    psi[0] = T(1.001);
    psi[ld] = T(0.002, -0.003);
    psi[ld + 1] = T(0.999);
    return psi;
}

template <typename T>
void check_orthogonal(const std::vector<T>& psi)
{
    const std::vector<Wide> gram = overlaps(psi);
    const double tolerance = hsolver::orth_tolerance<T>();
    for (int j = 0; j < bands; ++j)
    {
        for (int i = 0; i < bands; ++i)
        {
            const Wide expected = i == j ? Wide(1) : Wide(0);
            EXPECT_LE(std::abs(gram[i + j * bands] - expected), tolerance);
        }
    }
}

template <typename Real>
class OrthonormalTest : public testing::Test
{
};
using Precisions = testing::Types<double, float>;
TYPED_TEST_SUITE(OrthonormalTest, Precisions);

TYPED_TEST(OrthonormalTest, CorrectsComplexOverlapsWithoutFallback)
{
    using T = std::complex<TypeParam>;
    const hsolver::diag_comm_info comm = self_comm();
    hsolver::LinearAlgebra<T, Device> algebra(comm);
    hsolver::Orthonormal<T, Device> orth(comm, algebra);
    for (const hsolver::OrthMethod method: {hsolver::OrthMethod::cholesky, hsolver::OrthMethod::lowdin, hsolver::OrthMethod::newton_schulz})
    {
        SCOPED_TRACE(hsolver::orth_method_name(method));
        std::vector<T> psi = orbitals<T>();
        const std::vector<T> original(psi);
        const hsolver::OrthResult result = orth.apply(psi.data(), ld, dim, bands, method, false);
        ASSERT_EQ(result.status, hsolver::OrthStatus::accepted);
        EXPECT_EQ(result.actual, method);
        EXPECT_EQ(result.fallbacks, 0);
        EXPECT_EQ(result.failure, hsolver::OrthFailure::none);
        check_orthogonal(psi);
        for (int band = 0; band < bands; ++band)
        {
            const int padding = band * ld + dim;
            EXPECT_EQ(std::memcmp(psi.data() + padding, original.data() + padding, (ld - dim) * sizeof(T)), 0);
        }
    }
}

TYPED_TEST(OrthonormalTest, NewtonSchulzFallsBackOutsideItsRange)
{
    using T = std::complex<TypeParam>;
    const hsolver::diag_comm_info comm = self_comm();
    hsolver::LinearAlgebra<T, Device> algebra(comm);
    hsolver::Orthonormal<T, Device> orth(comm, algebra);
    std::vector<T> psi(ld * bands, T(0));
    psi[0] = T(2);
    psi[ld + 1] = T(2);
    const hsolver::OrthResult result = orth.apply(psi.data(), ld, dim, bands, hsolver::OrthMethod::newton_schulz, false);
    ASSERT_EQ(result.status, hsolver::OrthStatus::accepted);
    EXPECT_EQ(result.actual, hsolver::OrthMethod::cholesky);
    EXPECT_EQ(result.fallbacks, 1);
    EXPECT_EQ(result.failure, hsolver::OrthFailure::none);
    check_orthogonal(psi);
}

TYPED_TEST(OrthonormalTest, RejectsInvalidOrbitalsWithoutOverwritingThem)
{
    using T = std::complex<TypeParam>;
    const hsolver::diag_comm_info comm = self_comm();
    hsolver::LinearAlgebra<T, Device> algebra(comm);
    hsolver::Orthonormal<T, Device> orth(comm, algebra);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for (const hsolver::OrthMethod method: {hsolver::OrthMethod::cholesky, hsolver::OrthMethod::none})
    {
        for (const double value: {0.0, nan, inf})
        {
            SCOPED_TRACE(hsolver::orth_method_name(method));
            SCOPED_TRACE(value);
            std::vector<T> psi(ld * bands, T(0));
            psi[0] = T(value);
            psi[ld + 1] = T(1);
            const std::vector<T> original(psi);
            const hsolver::OrthResult result = orth.apply(psi.data(), ld, dim, bands, method, false);
            EXPECT_EQ(result.status, hsolver::OrthStatus::failed);
            const hsolver::OrthFailure nonfinite
                = method == hsolver::OrthMethod::none ? hsolver::OrthFailure::nonfinite_norm : hsolver::OrthFailure::nonfinite_gram;
            const hsolver::OrthFailure expected = value == 0.0 ? hsolver::OrthFailure::nonpositive_norm : nonfinite;
            EXPECT_EQ(result.failure, expected);
            EXPECT_EQ(std::memcmp(psi.data(), original.data(), psi.size() * sizeof(T)), 0);
        }
    }
    std::vector<T> duplicate(ld * bands, T(0));
    duplicate[0] = T(1);
    duplicate[ld] = T(1);
    const std::vector<T> original(duplicate);
    const hsolver::OrthResult result = orth.apply(duplicate.data(), ld, dim, bands, hsolver::OrthMethod::cholesky, false);
    EXPECT_EQ(result.status, hsolver::OrthStatus::failed);
    EXPECT_EQ(result.failure, hsolver::OrthFailure::factorization_failed);
    EXPECT_EQ(std::memcmp(duplicate.data(), original.data(), duplicate.size() * sizeof(T)), 0);
}

TYPED_TEST(OrthonormalTest, InitialInspectionAndNoneAreReadOnly)
{
    using T = std::complex<TypeParam>;
    const hsolver::diag_comm_info comm = self_comm();
    hsolver::LinearAlgebra<T, Device> algebra(comm);
    hsolver::Orthonormal<T, Device> orth(comm, algebra);
    std::vector<T> psi = orbitals<T>();
    const std::vector<T> original(psi);
    const std::vector<Wide> gram = overlaps(psi);
    for (const bool collect: {false, true})
    {
        const hsolver::OrthResult inspected = orth.inspect(psi.data(), ld, dim, bands, true, collect);
        const hsolver::OrthResult disabled = orth.apply(psi.data(), ld, dim, bands, hsolver::OrthMethod::none, collect);
        EXPECT_EQ(inspected.status, hsolver::OrthStatus::inspected);
        EXPECT_GT(inspected.after, hsolver::orth_tolerance<T>());
        EXPECT_EQ(disabled.status, hsolver::OrthStatus::disabled);
        EXPECT_FALSE(disabled.gram_checked);
        EXPECT_TRUE(std::isnan(disabled.after));
        const std::size_t expected_size = collect ? bands : 0;
        ASSERT_EQ(inspected.norms.size(), expected_size);
        ASSERT_EQ(disabled.norms.size(), expected_size);
        for (int band = 0; collect && band < bands; ++band)
        {
            EXPECT_NEAR(inspected.norms[band], gram[band + band * bands].real(), 1e-14);
            EXPECT_NEAR(disabled.norms[band], gram[band + band * bands].real(), 1e-14);
        }
        EXPECT_EQ(std::memcmp(psi.data(), original.data(), psi.size() * sizeof(T)), 0);
    }
}
} // namespace
