#include "source_hsolver/kernels/cuda/diag_cusolver.cuh"
#include "gtest/gtest.h"
#include <cuda_runtime_api.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace
{
void check_retained_density(const bool subset)
{
    const int sizes[3] = {11, 67, 11};
    Diag_Cusolver_gvd retained;
    std::string error;
    for (int size_index = 0; size_index < 3; ++size_index)
    {
        const int n = sizes[size_index];
        const int bands = n - 3;
        const std::size_t elements = static_cast<std::size_t>(n) * n;
        const std::size_t wave_elements = static_cast<std::size_t>(n) * bands;
        std::vector<double> h(elements, 0.0);
        std::vector<double> s(elements, 0.0);
        for (int col = 0; col < n; ++col)
        {
            for (int row = 0; row <= col; ++row)
            {
                const int index = col * n + row;
                h[index] = row == col ? row + 0.75 : 0.03;
                s[index] = row == col ? 1.0 + 0.001 * row : 0.001;
            }
        }
        std::vector<double> legacy_values(n);
        std::vector<double> retained_values(n);
        std::vector<double> legacy_wave(elements);
        // H can receive the returned eigenvectors once its input has been uploaded.
        double* host_wave = h.data();
        std::vector<double> density(elements);
        std::vector<double> weights(bands);
        Diag_Cusolver_gvd legacy;
        legacy.solve_full(n, n, h.data(), s.data(), legacy_values.data(), legacy_wave.data(), false);
        if (subset)
        {
            const bool solved = retained.solve_subset(n, bands, h.data(), s.data(), retained_values.data(), host_wave, true);
            ASSERT_TRUE(solved);
        }
        else
        {
            retained.solve_full(n, bands, h.data(), s.data(), retained_values.data(), host_wave, true);
        }
        for (int band = 0; band < bands; ++band)
        {
            const double phase = 0.37 * band;
            weights[band] = std::cos(phase);
            EXPECT_NEAR(legacy_values[band], retained_values[band], 1.0e-12);
        }
        weights[0] = 2.0;
        weights[1] = 0.0;
        weights[2] = 0.125;
        // The direct DM must not depend on the host eigenvector materialization.
        std::fill_n(host_wave, wave_elements, 0.0);
        const bool ready = retained.density_from_retained(bands, weights.data(), density.data(), error);
        ASSERT_TRUE(ready);
        EXPECT_TRUE(error.empty());
        for (int col = 0; col < n; ++col)
        {
            for (int row = 0; row < n; ++row)
            {
                double reference = 0.0;
                for (int band = 0; band < bands; ++band)
                    reference += weights[band] * legacy_wave[band * n + row] * legacy_wave[band * n + col];
                const int index = col * n + row;
                EXPECT_NEAR(density[index], reference, 2.0e-12);
            }
        }
        const bool invalid = retained.density_from_retained(0, weights.data(), density.data(), error);
        EXPECT_FALSE(invalid);
        EXPECT_NE(error.find("validate retained eigenvectors"), std::string::npos);
    }
}
}

TEST(CusolverRetainedDm, FractionalSignedWeightsAfterHostMaterialization)
{
    check_retained_density(false);
}

#if CUDART_VERSION >= 10010
TEST(CusolverRetainedDm, SubsetFractionalSignedWeightsAfterHostMaterialization)
{
    check_retained_density(true);
}
#endif

TEST(CusolverRetainedDm, ComplexSolveInvalidatesRealEigenvectors)
{
    const int n = 2;
    double h[4] = {1.0, 0.0, 0.0, 2.0};
    double s[4] = {1.0, 0.0, 0.0, 1.0};
    double values[2] = {};
    double wave[4] = {};
    double density[4] = {};
    const double weights[2] = {2.0, 0.125};
    Diag_Cusolver_gvd solver;
    std::string error;
    const bool empty = solver.density_from_retained(n, weights, density, error);
    EXPECT_FALSE(empty);
    solver.solve_full(n, n, h, s, values, wave, true);
    const bool real = solver.density_from_retained(n, weights, density, error);
    ASSERT_TRUE(real);
    EXPECT_TRUE(error.empty());
    const std::complex<double> complex_h[4] = {1.0, 0.0, 0.0, 2.0};
    const std::complex<double> complex_s[4] = {1.0, 0.0, 0.0, 1.0};
    std::complex<double> complex_wave[4] = {};
    solver.solve_full(n, n, complex_h, complex_s, values, complex_wave);
    const bool stale = solver.density_from_retained(n, weights, density, error);
    EXPECT_FALSE(stale);
}

TEST(CusolverRetainedDm, HostOnlySolveReleasesRetainedEigenvectors)
{
    const int n = 2;
    const double h[4] = {1.0, 0.0, 0.0, 2.0};
    const double s[4] = {1.0, 0.0, 0.0, 1.0};
    const double weights[2] = {2.0, 0.125};
    double values[2] = {};
    double wave[4] = {};
    double density[4] = {};
    std::string error;
    Diag_Cusolver_gvd solver;
    solver.solve_full(n, n, h, s, values, wave, true);
    ASSERT_TRUE(solver.density_from_retained(n, weights, density, error));
    solver.solve_full(n, n, h, s, values, wave, false);
    EXPECT_FALSE(solver.density_from_retained(n, weights, density, error));
#if CUDART_VERSION >= 10010
    ASSERT_TRUE(solver.solve_subset(n, 1, h, s, values, wave, true));
    ASSERT_TRUE(solver.density_from_retained(1, weights, density, error));
    ASSERT_TRUE(solver.solve_subset(n, 1, h, s, values, wave, false));
    EXPECT_FALSE(solver.density_from_retained(1, weights, density, error));
#endif
}
