#include "source_hsolver/kernels/cuda/diag_cusolver.cuh"
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <vector>

TEST(CusolverRetainedDm, FractionalSignedWeightsAfterHostMaterialization)
{
    const int sizes[2] = {11, 67};
    for (int size_index = 0; size_index < 2; ++size_index)
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
        std::vector<double> host_wave(wave_elements);
        std::vector<double> density(elements);
        std::vector<double> weights(bands);
        Diag_Cusolver_gvd legacy;
        legacy.Dngvd_double(n, n, h.data(), s.data(), legacy_values.data(), legacy_wave.data());
        Diag_Cusolver_gvd retained;
        retained.Dngvd_double_retained(n, h.data(), s.data(), retained_values.data(), host_wave.data(), bands);
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
        std::fill(host_wave.begin(), host_wave.end(), 0.0);
        ASSERT_TRUE(retained.density_from_retained(bands, weights.data(), density.data()));
        for (int col = 0; col < n; ++col)
        {
            for (int row = 0; row < n; ++row)
            {
                double reference = 0.0;
                for (int band = 0; band < bands; ++band)
                    reference += weights[band] * legacy_wave[band * n + row] * legacy_wave[band * n + col];
                EXPECT_NEAR(density[col * n + row], reference, 2.0e-12);
            }
        }
        EXPECT_FALSE(retained.density_from_retained(0, weights.data(), density.data()));
    }
}
