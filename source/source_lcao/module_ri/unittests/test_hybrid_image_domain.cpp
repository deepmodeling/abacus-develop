#include "source_lcao/module_ri/hybrid_image_domain.h"
#include "source_lcao/module_ri/hybrid_gauge_phase.h"
#include "source_base/matrix3.h"
#include "source_base/vector3.h"
#include "source_lcao/module_ri/mix_dmk_2d.h"
#include <gtest/gtest.h>

namespace
{
using Cell = RI_2D_Comm::ImageCell;
using Complex = std::complex<double>;
using Density = RI_2D_Comm::HybridGaugeComplexDensity;

Density source_density()
{
    Density density(1);
    for (int R = -2; R < 2; ++R)
    {
        RI::Tensor<Complex> tensor({1, 1});
        tensor.ptr()[0] = Complex(0.4 + 0.02 * R, 0.03 * R);
        const Cell cell{{R, 0, 0}};
        density[0][0][{1, cell}] = tensor;
    }
    return density;
}
}

TEST(HybridImageDomain, CentersOddAndEvenPeriods)
{
    const Cell cell{{3, -3, 1}};
    const Cell physical{{4, 4, 1}};
    const Cell expected{{-1, 1, 0}};
    EXPECT_EQ(RI_2D_Comm::centered_image(cell, physical), expected);
    const Cell images{{7, 7, 3}};
    EXPECT_EQ(RI_2D_Comm::centered_image(cell, images), cell);
    const Cell invalid{{0, 4, 1}};
    EXPECT_THROW(RI_2D_Comm::centered_image(cell, invalid), std::invalid_argument);
}

TEST(HybridImageDomain, BoundsAllFourContractions)
{
    const Cell C{{1, 2, 0}};
    const Cell V{{3, 1, 0}};
    const Cell expected{{11, 11, 1}};
    EXPECT_EQ(RI_2D_Comm::hybrid_image_period(C, V), expected);
}

TEST(HybridImageDomain, RetainsExistingMaxNormThreshold)
{
    Density density = source_density();
    const Cell extent = RI_2D_Comm::retained_image_extent(density[0], 0.39);
    const Cell expected{{1, 0, 0}};
    EXPECT_EQ(extent, expected);
    EXPECT_EQ(density[0][0].size(), 2u);
}

TEST(HybridImageDomain, ExtendsEachSourceBeforeRealMixerAndReset)
{
    const Cell physical{{4, 1, 1}};
    const Cell images{{7, 1, 1}};
    const ModuleBase::Matrix3 lattice(1, 0, 0, 0, 1, 0, 0, 0, 1);
    Mix_DMk_2D<Complex> mixer;
    mixer.set_mixing_plain(0.3);
    RI_2D_Comm::HybridDensityLayout layout;
    std::map<int, Complex> previous;
    for (int iteration = 0; iteration < 5; ++iteration)
    {
        Density density = source_density();
        const Density unmodified = density;
        const bool restart = iteration == 0 || iteration == 3;
        const double A = 0.17 + 0.11 * iteration;
        const ModuleBase::Vector3<double> source_A(A, 0, 0);
        const std::array<double, 3> mesh_origin{{0.0, 0.0, 0.0}};
        const std::array<double, 3> input_phase = RI_2D_Comm::hybrid_input_phase_per_cell(lattice, 1.0, source_A);
        RI_2D_Comm::unfold_hybrid_images(density, physical, images, mesh_origin, input_phase);
        RI_2D_Comm::mix_common_frame_density(restart, mixer, layout, density);
        ASSERT_EQ(density[0][0].size(), 7u);
        for (int R = -3; R <= 3; ++R)
        {
            const Cell cell{{R, 0, 0}};
            const Cell representative = RI_2D_Comm::centered_image(cell, physical);
            const Complex raw = unmodified[0].at(0).at({1, representative}).ptr()[0];
            const Complex converted = raw * std::exp(Complex(0, -A * R));
            const Complex expected = restart ? converted : 0.7 * previous[R] + 0.3 * converted;
            EXPECT_NEAR(std::abs(density[0][0].at({1, cell}).ptr()[0] - expected), 0, 1e-14);
            previous[R] = expected;
        }
        const Complex original = unmodified[0].at(0).at({1, Cell{{1, 0, 0}}}).ptr()[0];
        EXPECT_NEAR(std::abs(original - Complex(0.42, 0.03)), 0, 1e-15);
    }
}

TEST(HybridImageDomain, DecodesOutputUsingInternalPeriod)
{
    Density density = source_density();
    RI_2D_Comm::HybridSpinDensity& blocks = density[0];
    const Cell images{{7, 1, 1}};
    const Cell external{{5, 0, 0}};
    blocks[0][{2, external}] = blocks[0].begin()->second;
    RI_2D_Comm::decode_hybrid_images(blocks, images);
    EXPECT_EQ(blocks[0].count({2, Cell{{-2, 0, 0}}}), 1u);
}

// Independent finite-mesh Fourier sums cover Gamma, even MP, and shifted grids.
TEST(HybridImageDomain, MatchesDirectFourierSumAtActualImages)
{
    const Cell physical{{4, 2, 1}};
    const Cell images{{13, 7, 3}};
    const double pi = std::acos(-1.0);
    const ModuleBase::Matrix3 lattice(1.2, 0.3, 0.1, -0.2, 0.9, 0.4, 0.05, -0.1, 1.7);
    const double lat0 = 2.3;
    const ModuleBase::Vector3<double> source_A(0.37, -0.19, 0.23);
    for (const std::array<double, 3>& origin : {std::array<double, 3>{{0.0, 0.0, 0.0}},
                               std::array<double, 3>{{-0.375, -0.25, 0.0}},
                               std::array<double, 3>{{0.07, 0.13, -0.09}}})
    {
        auto direct = [&](const Cell& cell) {
            Complex value = 0.0;
            for (int ix = 0; ix < physical[0]; ++ix)
            {
                for (int iy = 0; iy < physical[1]; ++iy)
                {
                    const double kx = origin[0] + double(ix) / physical[0];
                    const double ky = origin[1] + double(iy) / physical[1];
                    const double cycles = kx * cell[0] + ky * cell[1] + origin[2] * cell[2];
                    const Complex weight(0.4 + 0.03 * ix, -0.02 * iy);
                    value += weight * std::exp(Complex(0, -2.0 * pi * cycles));
                }
            }
            return value;
        };
        Density density(1);
        for (int x = -2; x < 2; ++x)
        {
            for (int y = -1; y < 1; ++y)
            {
                const Cell cell{{x, y, 0}};
                RI::Tensor<Complex> tensor({1, 1});
                tensor.ptr()[0] = direct(cell);
                density[0][0][{0, cell}] = tensor;
            }
        }
        const Density original = density;
        const std::array<double, 3> input_phase = RI_2D_Comm::hybrid_input_phase_per_cell(lattice, lat0, source_A);
        RI_2D_Comm::unfold_hybrid_images(density, physical, images, origin, input_phase);
        Mix_DMk_2D<Complex> mixer;
        mixer.set_mixing_plain(1.0);
        RI_2D_Comm::HybridDensityLayout layout;
        RI_2D_Comm::mix_common_frame_density(true, mixer, layout, density);
        ASSERT_EQ(density[0][0].size(), 273u);
        for (const RI_2D_Comm::HybridDensityRow::value_type& block : density[0][0])
        {
            const Cell& cell = block.first.second;
            // Independent Cartesian projection checks row convention and lat0.
            const double rx = lat0 * (1.2 * cell[0] - 0.2 * cell[1] + 0.05 * cell[2]);
            const double ry = lat0 * (0.3 * cell[0] + 0.9 * cell[1] - 0.1 * cell[2]);
            const double rz = lat0 * (0.1 * cell[0] + 0.4 * cell[1] + 1.7 * cell[2]);
            const double phase = -(source_A.x * rx + source_A.y * ry + source_A.z * rz);
            const Complex expected = direct(cell) * std::exp(Complex(0, phase));
            EXPECT_NEAR(std::abs(block.second.ptr()[0] - expected), 0.0, 2e-14);
        }
        // Expanded tensors and input phase must not mutate the source/history.
        for (const RI_2D_Comm::HybridDensityRow::value_type& block : original[0].at(0))
        {
            EXPECT_NEAR(std::abs(block.second.ptr()[0] - direct(block.first.second)), 0.0, 1e-15);
        }
    }
}

TEST(HybridImageDomain, ZeroFieldStillPreservesShiftedMeshTwist)
{
    Density density = source_density();
    const Density original = density;
    const Cell physical{{4, 1, 1}};
    const Cell images{{7, 1, 1}};
    const std::array<double, 3> origin{{0.125, 0.0, 0.0}};
    const std::array<double, 3> input_phase{{0.0, 0.0, 0.0}};
    RI_2D_Comm::unfold_hybrid_images(density, physical, images, origin, input_phase);
    const Cell actual{{3, 0, 0}};
    const Cell representative{{-1, 0, 0}};
    const Complex source = original[0].at(0).at({1, representative}).ptr()[0];
    const Complex result = density[0].at(0).at({1, actual}).ptr()[0];
    // Half a mesh step: exp(-i 2*pi*(1/8)*4) = -1.
    EXPECT_NEAR(std::abs(result + source), 0.0, 1e-15);
    EXPECT_EQ(original[0].at(0).at({1, representative}).ptr()[0], source);
}
